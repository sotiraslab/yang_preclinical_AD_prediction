# ==============================================================================

# Name: merge_HABS.R
# Author: Braden Yang
# Created: 03/17/2025
# Description: Load HABS CSVs, preprocess and merge

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
    make_option(c("-o", "--odir"), action="store", default="data/tidy/habs",
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
# ===== DEFINE FUNCTIONS =====
# ============================

habs_get_raw_img_list <- function(raw_img_dir, scan_type) {

    if (scan_type == "pet") {
        img_subj_pattern <- "(P_\\w{6})_(\\d{4}-\\d{2}-\\d{2})_(HAB_\\d\\.\\d)_PIB_40-60\\.nii\\.gz"
    } else if (scan_type == "mri") {
        img_subj_pattern <- "(P_\\w{6})_(\\d{4}-\\d{2}-\\d{2})_(HAB_\\d\\.\\d)_ADNI_\\dX-MPRAGE\\.nii\\.gz"
    } else if (scan_type == "tau") {
        img_subj_pattern <- "(P_\\w{6})_(\\d{4}-\\d{2}-\\d{2})_(HAB_\\d\\.\\d)_FTP_\\d{2}-\\d{3}\\.nii\\.gz"
    } else {
        stop("Invalid scan type")
    }
    
    list.files(path = raw_img_dir, pattern = "*.nii.gz", full.names = TRUE) %>%
        as_tibble() %>%
        rename(filepath = value) %>%
        mutate(filename = basename(filepath)) %>%
        mutate(
            matches = str_match(filename, img_subj_pattern),
            subj = matches[,2],
            session = matches[,3],
            StudyArc = matches[,4]
        ) %>%
        select(-matches)
    
}

# ============================
# ===== DEFINE VARIABLES =====
# ============================

habs_dir <- "/ceph/chpc/shared/aristeidis_sotiras_group/aris_data/HABS"
habs_img_dir <- file.path(habs_dir, "images")
habs_tabular_dir <- file.path(habs_dir, "tabular")

# ====================
# ===== LOAD CSV =====
# ====================

# tabular data
csv_files <- list.files(path = habs_tabular_dir, pattern = "*.csv", full.names = TRUE)
csv_list <- csv_files %>%
    set_names(basename) %>%
    map(~read_csv(.x, show_col_types = FALSE))

# raw imaging paths
pet_raw_df <- habs_get_raw_img_list(file.path(habs_img_dir, "pib"), "pet")
mri_raw_df <- habs_get_raw_img_list(file.path(habs_img_dir, "mprage"), "mri")
tau_raw_df <- habs_get_raw_img_list(file.path(habs_img_dir, "ftp"), "tau")

# ==============================
# ===== TIDY PET AND MERGE =====
# ==============================

# NOTE: there are 4 types of PET regional analysis columns: 
# "PIB_FS_DVR", "PIB_FS_SUVR", "PIB_T1_DVR", "PIB_T1_SUVR"
# not sure what "FS" and "T1" mean, but for now we will use "PIB_T1_SUVR"

# get base PET df
pet_df <- csv_list[["PIB_FS6_DVR_HABS_DataRelease_2.0.csv"]] %>%
    rename(
        subj = SubjIDshort,
        age = PIB_Age,
        session = PIB_SessionDate
    ) %>%
    select(subj, age, session, StudyArc, PIB_T1_SUVR_Group) %>%
    mutate(
        session = as.character(session),
        site = "HABS",
        tracer = "PIB",
        amyloid_positive = case_match(PIB_T1_SUVR_Group, "PIB+"~TRUE, "PIB-"~FALSE, .default = NA)
    ) %>%
    full_join(pet_raw_df %>% select(-c(StudyArc, filename)), by = c("subj", "session")) %>%
    rename(
        filepath.raw_amyloid = filepath,
    )

# get sex and APOE
demog_df <- csv_list[["Demographics_HABS_DataRelease_2.0.csv"]] %>%
    mutate(
        sex = factor(BiologicalSex, levels = c("F", "M")),
        apoe = case_when(APOE_haplotype %in% c(34, 44) ~ 1, is.na(APOE_haplotype) ~ NA, .default = 0) %>% as.factor()
    ) %>%
    rename(subj = SubjID) %>%
    select(subj, sex, apoe) %>%
    distinct()
pet_df <- pet_df %>% left_join(demog_df, by = "subj")

# match with raw MRI paths
mri_df <- csv_list[["ADNI_MRI_FS6_XSec_HABS_DataRelease_2.0.csv"]] %>%
    rename(
        subj = SubjIDshort,
        age = MRI_Age,
        session = MRI_SessionDate
    ) %>%
    select(subj, age, session, StudyArc) %>%
    mutate(session = as.character(session)) %>%
    full_join(mri_raw_df %>% select(-StudyArc), by = c("subj", "session"))
pet_df <- pet_df %>% byyfunctions::match_data(
    match_df = mri_df,
    ref_col = "age",
    match_col = "age",
    group_col = "subj",
    select_col = c("filepath", "session", "StudyArc"),
    suffix = ".raw_t1",
    diff_col = "age_diff.raw_t1"
)

# ====================
# ===== TIDY CDR =====
# ====================

# NOTE: there is some disagreement between global CDR and the HABS diagnosis
# (ex. 113 cases of CDR=0.5 have diagnosis of CN). Here we will use the CDR
# score as our measure of cognitive impairment to keep consistent with other
# datasets

cdr_df <- csv_list[["ClinicalMeasures_HABS_DataRelease_2.0.csv"]] %>%
    select(SubjIDshort, NP_SessionDate, CDR_Global) %>%
    rename(
        subj = SubjIDshort,
        session = NP_SessionDate,
        cdr = CDR_Global
    ) %>%
    mutate(session = as.character(session))

# get age at each CDR session
session_age <- csv_list[["Demographics_HABS_DataRelease_2.0.csv"]] %>% 
    select(SubjID, NP_Age, NP_SessionDate) %>%
    rename(subj = SubjID, age = NP_Age, session = NP_SessionDate) %>%
    mutate(session = as.character(session))
cdr_df <- cdr_df %>% left_join(session_age, by = c("subj", "session")) %>%
    drop_na()

# match CDR
pet_df <- pet_df %>% byyfunctions::match_data(
    match_df = cdr_df,
    ref_col = "age",
    match_col = "age",
    group_col = "subj",
    select_col = "cdr",
    suffix = ".cdr",
    diff_col = "age_diff.cdr"
)

# ========================
# ===== SAVE OUTPUTS =====
# ========================

if (SAVE_CSV) {
    byyfunctions::make_dir(opt$odir)
    if (opt$csv) {
        write_csv(pet_df, file.path(opt$odir, "habs_amyloid.csv"))
        write_csv(cdr_df, file.path(opt$odir, "habs_cdr.csv"))
    } else {
        write_rds(pet_df, file.path(opt$odir, "habs_amyloid.RDS"))
        write_rds(cdr_df, file.path(opt$odir, "habs_cdr.RDS"))
    }
}
