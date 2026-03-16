# ==============================================================================

# Name: merge_MCSA.R
# Author: Braden Yang
# Created: 04/04/2025
# Description: Load MCSA CSVs, preprocess and merge

# ==============================================================================

rm(list = ls())

SAVE_CSV <- TRUE

# ===========================
# ===== IMPORT PACKAGES =====
# ===========================

library(tidyverse)
library(optparse)
library(devtools)

# ===========================
# ===== PARSE ARGUMENTS =====
# ===========================

option_list <- list(
    make_option(c("-w", "--wdir"), action="store", default=NULL,
        type="character", help="Path to project directory (if none, uses current working directory)"),
    make_option(c("-o", "--odir"), action="store", default="data/tidy/mcsa",
        type="character", help="Path to output directory"),
    make_option(c("-c", "--csv"), action="store_true", default=FALSE,
        help="output in .csv format instead of .RDS")
)

opt <- parse_args(OptionParser(option_list = option_list))

# set working directory
if (!is.null(opt$wdir)) {
    setwd(opt$wdir)
}

# load byyfunctions
load_all("submodules/byyfunctions")

# make output directory
byyfunctions::make_dir(opt$odir)

# ============================
# ===== DEFINE VARIABLES =====
# ============================

mcsa_dir <- "/ceph/chpc/shared/aristeidis_sotiras_group/aris_data/MCSA"
mcsa_img_dir <- file.path(mcsa_dir, "images/bids")
mcsa_tabular_dir <- file.path(mcsa_dir, "tabular")

# ====================
# ===== LOAD CSV =====
# ====================

data_df <- read_csv(file.path(mcsa_tabular_dir, "MCSA_Data_04Apr2025.csv"))
# data_dict <- read_csv(file.path(mcsa_tabular_dir, "MCSA_Data_Dictionary_04Apr2025.csv"))

# ==============================
# ===== TIDY PET AND MERGE =====
# ==============================

data_df_tidy <- data_df %>%
    mutate(
        site = "MCSA",
        subj = mcsa_id,
        visit = visit_num,
        visit_date = visitdate,
        scan_date = imagingdate,
        age.visit = calc_age_vis,
        age.scan = age.visit + (as.numeric(scan_date - visit_date) / 365.25),
        sex = factor(if_else(male==1, "M", "F")),
        apoe = factor(any_e4),   # NOTE: equals one if any allele is 4, but doesn't exclude 2
        cdr = CDRGLOB,
        amyloid_positive = ABNORMAL_SPM12_PIB == 1,
        summary_cortical_amyloid.centiloid = SPM12_PIB_CENTILOID,
        .keep = "none"
    )

pet_df <- data_df_tidy %>%
    filter(!is.na(amyloid_positive)) %>%
    select(site, subj, scan_date, age.scan, sex, apoe, amyloid_positive, summary_cortical_amyloid.centiloid) %>%
    rename(age = age.scan) %>%
    mutate(tracer = "PIB")
cdr_df <- data_df_tidy %>%
    filter(!is.na(cdr)) %>%
    select(subj, visit_date, age.visit, cdr) %>%
    rename(age = age.visit)

# merge with raw PIB PET scans
pet_raw_df <- read_csv("data/raw_img_paths/mcsa_amyloid.txt", col_names = "filepath.raw_amyloid", col_types = "c") %>%
    mutate(.temp = str_replace(filepath.raw_amyloid, str_c(mcsa_img_dir, "/"), "")) %>%
    separate(.temp, into = c("sub", "ses", "mod", "file"), sep = "/") %>%
    mutate(
        subj = sub %>% str_sub(start = 5),
        scan_date = ses %>% str_sub(start = 5) %>% ymd(),
        filepath.raw_amyloid.json = str_replace(filepath.raw_amyloid, ".nii.gz", ".json")
    ) %>%
    select(subj, scan_date, filepath.raw_amyloid, filepath.raw_amyloid.json)
pet_df <- pet_df %>% left_join(pet_raw_df, by = c("subj", "scan_date"))

# merge with raw T1 scans
mri_raw_df <- read_csv("data/raw_img_paths/mcsa_t1.txt", col_names = "filepath", col_types = "c") %>%
    mutate(.temp = str_replace(filepath, str_c(mcsa_img_dir, "/"), "")) %>%
    separate(.temp, into = c("sub", "ses", "mod", "file"), sep = "/") %>%
    mutate(
        subj = sub %>% str_sub(start = 5),
        scan_date = ses %>% str_sub(start = 5) %>% ymd()
    ) %>%
    select(subj, scan_date, filepath)
pet_df <- pet_df %>% byyfunctions::match_data(
    match_df = mri_raw_df,
    ref_col = "scan_date",
    match_col = "scan_date",
    group_col = "subj",
    select_col = c("filepath"),
    suffix = ".raw_t1",
    diff_col = "age_diff.raw_t1"
) %>%
    mutate(
        age_diff.raw_t1 = as.numeric(age_diff.raw_t1) / 365.25,
        filepath.raw_t1.json = str_replace(filepath.raw_t1, ".nii.gz", ".json"),  # get json file
        age.raw_t1 = age + (as.numeric(scan_date.raw_t1 - scan_date) / 365.25),
    )

# match CDR
pet_df <- pet_df %>% byyfunctions::match_data(
    match_df = cdr_df,
    ref_col = "scan_date",
    match_col = "visit_date",
    group_col = "subj",
    select_col = "cdr",
    suffix = ".cdr",
    diff_col = "age_diff.cdr"
) %>%
    mutate(age_diff.cdr = as.numeric(age_diff.cdr) / 365.25)

# ========================
# ===== SAVE OUTPUTS =====
# ========================

if (SAVE_CSV) {
    byyfunctions::make_dir(opt$odir)
    if (opt$csv) {
        write_csv(pet_df, file.path(opt$odir, "mcsa_amyloid.csv"))
        write_csv(cdr_df, file.path(opt$odir, "mcsa_cdr.csv"))
    } else {
        write_rds(pet_df, file.path(opt$odir, "mcsa_amyloid.RDS"))
        write_rds(cdr_df, file.path(opt$odir, "mcsa_cdr.RDS"))
    }
}
