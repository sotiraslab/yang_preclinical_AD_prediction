# ==============================================================================

# Name: merge_A4.R
# Author: Braden Yang
# Created: 11/11/2024
# Description: Load A4 CSVs, preprocess and merge

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
    make_option(c("-o", "--odir"), action="store", default="data/tidy/a4",
        type="character", help="Path to output directory"),
    make_option(c("-c", "--csv"), action="store_true", default=FALSE,
        help="output in .csv format instead of .RDS"),
    make_option(c("-t", "--treatment"), action="store_true", default=FALSE,
        help="include treatment arm subjects in output")
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
# ===== SOURCE FUNCTIONS =====
# ============================

source("code/a4_functions.R")

# ============================
# ===== DEFINE FUNCTIONS =====
# ============================

a4_get_raw_img_list <- function(raw_dir) {
    list.files(path = raw_dir, pattern = "*.nii.gz", full.names = TRUE) %>%
        as_tibble() %>%
        rename(filepath = value) %>%
        mutate(filename = basename(filepath)) %>%
        separate_wider_delim(filename, names = c("SUBSTUDY", "MODALITY", "SUBMODALITY", "BID", "VISCODE"), delim = "_") %>%
        mutate(VISCODE = str_replace(VISCODE, ".nii.gz", ""))
}

# ============================
# ===== DEFINE VARIABLES =====
# ============================

a4_data_dir <- "/ceph/chpc/shared/aristeidis_sotiras_group/aris_data/A4/rawdata"
a4_tabular_dir <- file.path(a4_data_dir, "Clinical")

# ====================
# ===== LOAD CSV =====
# ====================

# tabular data
df_paths <- list(
    "subjinfo" = "Derived Data/SUBJINFO.csv",
    "sv" = "Derived Data/SV.csv",
    "amyloid" = "External Data/imaging_SUVR_amyloid.csv",
    "amyloid_VA" = "External Data/imaging_PET_VA.csv",
    "tau" = "External Data/imaging_SUVR_tau.csv",
    "mri" = "External Data/imaging_volumetric_mri.csv",
    "cdr" = "Raw Data/cdr.csv",
    "ptdemog" = "Raw Data/ptdemog.csv"
)
df_list <- purrr::map(
    df_paths,
    ~file.path(a4_tabular_dir, .x) %>%
        readr::read_csv(show_col_types = FALSE)
)

# raw imaging paths
pet_raw_df <- a4_get_raw_img_list(file.path(a4_data_dir, "FBP")) %>%
    a4_get_age() %>% select(BID, VISCODE, age, filepath)
mri_raw_df <- a4_get_raw_img_list(file.path(a4_data_dir, "T1")) %>%
    a4_get_age() %>% select(BID, VISCODE, age, filepath)
tau_raw_df <- a4_get_raw_img_list(file.path(a4_data_dir, "FTP")) %>%
    a4_get_age() %>% select(BID, VISCODE, age, filepath)

# ==============================
# ===== TIDY PET AND MERGE =====
# ==============================

# NOTE: we use the following procedure to compute age:
# - if `scan_date_DAYS_CONSENT` is not NA, then add this to age at consent to compute age at scan
# - if `scan_date_DAYS_CONSENT` is NA, then refer to SV.csv to get # of days since consent of the corresponding visit
# Note that there are discrepancies between `scan_date_DAYS_CONSENT` and SV.csv in some instances (about 800 scans)

pet_with_rescreens <- df_list[["amyloid"]] %>% 
    filter(brain_region != '') %>%
    pivot_wider(names_from=brain_region, values_from=suvr_cer) %>%
    left_join(df_list[["amyloid_VA"]] %>% select(SUBSTUDY, BID, VISCODE, elig_vi_1, elig_vi_2), by=c('SUBSTUDY','BID','VISCODE')) %>%
    mutate(
        amyloid_positive.a4_thresh = Composite_Summary > 1.15,
        amyloid_positive.visual_read = case_when(
            is.na(elig_vi_2) ~ case_match(elig_vi_1, "positive" ~ TRUE, "negative" ~ FALSE, .default = NA),
            elig_vi_1 == "positive" & elig_vi_2 == "positive" ~ TRUE,
            elig_vi_1 == "negative" & elig_vi_2 == "negative" ~ FALSE,
            .default = NA
        )
    ) %>%
    a4_get_age("scan_date_DAYS_CONSENT") %>%
    a4_merge_subjinfo("TX")
if (!opt$treatment) {  # remove subjects in the treatment arm
    pet_with_rescreens <- pet_with_rescreens %>%
        filter(TX != "Solanezumab" | is.na(TX))
}
pet_df <- pet_with_rescreens %>%
    a4_remove_rescreens() %>%  # remove rescreens
    rename_with(~str_c(.x, ".amyloid"), lprecuneus_gm:xlaal_frontal_med_orb) %>%
    select(SUBSTUDY, TX, BID, VISCODE, starts_with("age"), Composite_Summary, starts_with("amyloid_positive"), ends_with(".amyloid"))

# get sex and APOE
# NOTE: in SUBJINFO.csv, 1 corresponds to Female and 2 corresponds to Male; however, it is the opposite in ptdemog.csv. Here we use SUBJINFO.csv
pet_df <- pet_df %>%
    a4_merge_subjinfo(c("SEX", "APOEGN")) %>%
    rename(sex = SEX, apoe = APOEGN) %>%
    mutate(
        sex = case_match(sex, 1~"F", 2~"M", .default=NA) %>% as.factor(),
        apoe = case_when(apoe %in% c("E4/E4","E3/E4") ~ 1, is.na(apoe) ~ NA, .default = 0) %>% as.factor()
    )

# get ethnicity
ethnic_levels <- get_levels("PTETHNIC")
ethnic_labels <- get_labels("PTETHNIC")
ethnic_df <- df_list$ptdemog %>% select(BID, PTETHNIC) %>%
    mutate(ethnicity = factor(PTETHNIC, levels = ethnic_levels, labels = ethnic_labels)) %>%
    group_by(BID) %>% summarise(ethnicity = unique(ethnicity))
pet_df <- pet_df %>% left_join(ethnic_df, by = "BID")

race_levels <- get_levels("PTRACE")
race_labels <- get_labels("PTRACE")
race_df <- df_list$ptdemog %>% select(BID, PTRACE) %>%
    mutate(race = factor(PTRACE, levels = race_levels, labels = race_labels)) %>%
    group_by(BID) %>% summarise(race = unique(race))
pet_df <- pet_df %>% left_join(race_df, by = "BID")

# match with raw FBP paths
pet_df <- pet_df %>%
    left_join(
        pet_raw_df %>% rename(filepath.raw_amyloid = filepath) %>% select(-age),
        by = c("BID", "VISCODE")
    ) %>%
    mutate(filepath.raw_amyloid.json = str_replace(filepath.raw_amyloid, ".nii.gz", ".json"))  # get json file

# merge tabular MRI volume
mri_df <- df_list[["mri"]] %>%
    a4_remove_rescreens() %>%
    a4_get_age("Date_DAYS_CONSENT") %>%
    mutate(
        across(LeftCorticalGrayMatter:ForebrainParenchyma,
               ~ .x / IntraCranialVolume,
               .names = "{.col}.icvnorm"),
        filepath.from_tabular = file.path(a4_data_dir, "T1", str_glue("{SUBSTUDY}_MR_T1_{BID}_{VISCODE}.nii.gz")),
        filepath.from_tabular = if_else(file.exists(filepath.from_tabular), filepath.from_tabular, NA)
    )
mri_select <- c("filepath.from_tabular", mri_df %>% select(ends_with("icvnorm")) %>% colnames)
pet_df <- pet_df %>% byyfunctions::match_data(
    match_df = mri_df,
    ref_col = "age",
    match_col = "age",
    group_col = "BID",
    select_col = mri_select,
    suffix = ".mri",
    diff_col = "age_diff.mri"
)

# match with raw MRI paths
pet_df <- pet_df %>% byyfunctions::match_data(
    match_df = mri_raw_df,
    ref_col = "age",
    match_col = "age",
    group_col = "BID",
    select_col = "filepath",
    suffix = ".raw_t1",
    diff_col = "age_diff.raw_t1"
) %>%
    mutate(filepath.raw_t1.json = str_replace(filepath.raw_t1, ".nii.gz", ".json"))  # get json file

# merge tabular tau
tau_df <- df_list[["tau"]] %>%
    filter(brain_region != '') %>%
    pivot_wider(id_cols=c('SUBSTUDY','BID','VISCODE','scan_date_DAYS_CONSENT'), names_from=brain_region, values_from=suvr_persi) %>%  # persi white matter reference region
    a4_remove_rescreens() %>%
    a4_get_age("scan_date_DAYS_CONSENT") %>%
    mutate(
        filepath.from_tabular = file.path(a4_data_dir, "FTP", str_glue("{SUBSTUDY}_PET_FTP_{BID}_{VISCODE}.nii.gz")),
        filepath.from_tabular = if_else(file.exists(filepath.from_tabular), filepath.from_tabular, NA)
    )
tau_select <- c("filepath.from_tabular", tau_df %>% select(ends_with("_VOI")) %>% colnames)
pet_df <- pet_df %>% byyfunctions::match_data(
    match_df = tau_df,
    ref_col = "age",
    match_col = "age",
    group_col = "BID",
    select_col = tau_select,
    suffix = ".tau",
    diff_col = "age_diff.tau"
)

# match with raw FTP paths
pet_df <- pet_df %>% byyfunctions::match_data(
    match_df = tau_raw_df,
    ref_col = "age",
    match_col = "age",
    group_col = "BID",
    select_col = "filepath",
    suffix = ".raw_tau",
    diff_col = "age_diff.raw_tau"
) %>%
    mutate(filepath.raw_tau.json = str_replace(filepath.raw_tau, ".nii.gz", ".json"))  # get json file

# tidy
pet_df <- pet_df %>%
    rename(
        subj = BID,
        a4_substudy = SUBSTUDY,
        summary_cortical_amyloid = Composite_Summary
    ) %>%
    mutate(site = "A4", tracer = "FBP")

# ====================
# ===== TIDY CDR =====
# ====================

cdr_df <- df_list[["cdr"]] %>%
    a4_get_age("CDADTC_DAYS_CONSENT") %>%
    a4_remove_rescreens()
if (!opt$treatment) {  # include treatment arm subjects
    cdr_df <- cdr_df %>%
        a4_merge_subjinfo("TX") %>%
        filter(TX != "Solanezumab" | is.na(TX))  # remove subjects in the treatment arm
}
cdr_df <- cdr_df %>%
    select(BID, starts_with("age"), CDGLOBAL) %>%
    rename(subj = BID, cdr = CDGLOBAL)

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
        write_csv(pet_df, file.path(opt$odir, "a4_amyloid.csv"))
        write_csv(cdr_df, file.path(opt$odir, "a4_cdr.csv"))
    } else {
        write_rds(pet_df, file.path(opt$odir, "a4_amyloid.RDS"))
        write_rds(cdr_df, file.path(opt$odir, "a4_cdr.RDS"))
    }
}
