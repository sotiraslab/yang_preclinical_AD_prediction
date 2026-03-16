# ==============================================================================

# Name: merge_ADNI.R
# Author: Braden Yang
# Created: 08/23/22
# Description: merge ADNI CSV tables into a single table, tidy data, and 
#   perform additional preprocessing

# Notes:
# - in the CDR table, there are several entries which contain negative values;
#   for these cases, I am replacing them with NA
# - according to the amyloid PET processing document, all regional SUVRs have
#   already been normalized by a whole cerebellum SUVR

# ==============================================================================
# =========================== BEGIN SCRIPT ===================================== 
# ==============================================================================

# clear environment
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
    make_option(c("-o", "--odir"), action="store", default="data/tidy/adni",
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

adni_mutate_df <- function(.data) {

    # if VISCODE2 exists, rename to VISCODE
    if ("VISCODE2" %in% colnames(.data)) {
        .data <- .data %>%
            mutate(VISCODE = if_else(!is.na(VISCODE2), VISCODE2, VISCODE))
    }

    # replace "sc" and "scmri" with "bl"
    if ("VISCODE" %in% colnames(.data)){
        .data <- .data %>%
            mutate(VISCODE = ifelse(VISCODE %in% c("sc", "scmri"), "bl", VISCODE) %>% as.character)
    }

    # convert RID to numeric
    if ("RID" %in% colnames(.data)) {
        .data <- .data %>%
            mutate(RID = as.numeric(RID))
    }
    return(.data)
    
}

adni_compute_age <- function(.data, date_col = "EXAMDATE", group_col = "RID") {

    dob_merge <- dob_df %>%
        rename_with(~ group_col, RID)

    .data %>%
        left_join(dob_merge, by = group_col) %>%
        mutate(age = as.numeric(.data[[date_col]] - .DOB) / 365.25) %>%
        select(-.DOB)
    
}

adni_match_demo <- function(.data, group_col = "RID") {
    .data %>%
        left_join(sex_df %>% rename(!!group_col := RID), by = group_col) %>%  # sex
        left_join(apoe_df %>% rename(!!group_col := RID), by = group_col)  # APOE
}

rename_with_list <- function(.data, l) {
    .data %>% rename_with(~ l[.x] %>% unlist(), .cols = any_of(names(l)))
}

# ============================
# ===== DEFINE VARIABLES =====
# ============================

adni_csv_dir <- "/ceph/chpc/shared/aristeidis_sotiras_group/aris_data/ADNI/studydata"
adni_bids_dir <- "/ceph/chpc/shared/aristeidis_sotiras_group/aris_data/ADNI/bids/rawdata"

# =====================
# ===== LOAD DATA =====
# =====================

# load table of ROI feature names
roi_features_df <- readr::read_csv("data/adni_oasis_roi_features.csv", col_types="ccccc")

# load ADNI CSVs (ADNIMERGE, AV45 PET, CDR, and tau) as dataframes
df_paths <- list(
    "adnimerge" = "ADNIMERGE_29Aug2023.csv",
    "amyloid" = "UCBERKELEY_AMY_6MM_29Aug2023.csv",
    "cdr" = "CDR_29Aug2023.csv",
    "tau" = "UCBERKELEY_TAU_6MM_29Aug2023.csv",
    "tau_pvc" = "UCBERKELEY_TAUPVC_6MM_29Aug2023.csv",
    "registry" = "REGISTRY_29Aug2023.csv",
    "demo" = "PTDEMOG_29Aug2023.csv",
    "apoe" = "APOERES_29Aug2023.csv"
)
df_list <- purrr::map(
    df_paths,
    ~file.path(adni_csv_dir, .x) %>%
        readr::read_csv(show_col_types = FALSE) %>%
        adni_mutate_df()
)

# get DOB, sex and APOE tables
dob_df <- df_list[["demo"]] %>%
    drop_na(contains("PTDOB")) %>%
    select(RID, contains("PTDOB")) %>%
    mutate(.DOB = lubridate::ymd(str_c(PTDOBYY, PTDOBMM, "15", sep = "-"))) %>%
    select(-contains("PTDOB")) %>%
    unique()
sex_df <- df_list[["demo"]] %>%
    select(RID, PTGENDER) %>%
    filter(PTGENDER %in% c(1,2)) %>%
    mutate(sex = factor(PTGENDER, levels = c(1,2), labels = c("M","F"))) %>%
    select(-PTGENDER) %>%
    unique()
apoe_df <- df_list[["apoe"]] %>%
    mutate(apoe = case_when(
        APGEN1 == 2 | APGEN2 == 2 ~ 0,
        APGEN1 == 4 | APGEN2 == 4 ~ 1,
        APGEN1 == 3 & APGEN2 == 3 ~ 0,
        .default = NA,
    ) %>% as.factor()) %>%
    select(RID, apoe) %>%
    unique()

# get list of raw imaging
pet_raw_df <- read_csv("data/raw_img_paths/adni_amyloid.txt", col_names = "filepath.raw_amyloid", col_types = "c") %>%
    mutate(.temp = str_replace(filepath.raw_amyloid, str_c(adni_bids_dir, "/"), "")) %>%
    separate(.temp, into = c("sub", "ses", "mod", "file"), sep = "/") %>%
    mutate(
        RID = sub %>% str_sub(start = 9) %>% as.numeric,
        SCANDATE = ses %>% str_sub(start = 5) %>% ymd(),
        TRACER = case_when(
            str_detect(file, "AV45") ~ "FBP",
            str_detect(file, "FBB") ~ "FBB",
            .default = NA
        )
    ) %>%
    select(RID, SCANDATE, TRACER, filepath.raw_amyloid)
mri_raw_df <- read_csv("data/raw_img_paths/adni_t1.txt", col_names = "filepath", col_types = "c") %>%
    mutate(.temp = str_replace(filepath, str_c(adni_bids_dir, "/"), "")) %>%
    separate(.temp, into = c("sub", "ses", "mod", "file"), sep = "/") %>%
    mutate(
        subj = sub %>% str_sub(start = 9) %>% as.numeric,
        SCANDATE = ses %>% str_sub(start = 5) %>% ymd(),
        .accel_tesla = filepath %>% basename %>% str_split_fixed("_", 5) %>% .[,3],
        accel = case_when(
            str_detect(.accel_tesla, "noaccel") ~ FALSE,
            str_detect(.accel_tesla, "accel") ~ TRUE,
            .default = NA
        ),
        tesla = case_when(
            str_detect(.accel_tesla, "3T") ~ 3,
            str_detect(.accel_tesla, "1.5T") ~ 1.5,
            .default = NA
        )
    ) %>%
    adni_compute_age("SCANDATE", group_col = "subj") %>%
    select(subj, age, SCANDATE, filepath, accel, tesla)

# ====================
# ===== TIDY PET =====
# ====================

# get column remappings
pet_col_rename <- roi_features_df %>% filter(!is.na(adni_suvr)) %>% pull(fs_label) %>% str_c(".amyloid")
names(pet_col_rename) <- roi_features_df %>% filter(!is.na(adni_suvr)) %>% pull(adni_suvr)
vol_col_rename <- roi_features_df %>% filter(!is.na(adni_vol)) %>% pull(fs_label) %>% str_c(".vol")
names(vol_col_rename) <- roi_features_df %>% filter(!is.na(adni_vol)) %>% pull(adni_vol)
col_rename <- c(pet_col_rename, vol_col_rename)

# tidy PET
pet_df <- df_list[["amyloid"]] %>%
    adni_compute_age("SCANDATE") %>%
    adni_match_demo() %>%
    left_join(pet_raw_df, by = c("RID", "SCANDATE", "TRACER")) %>%
    mutate(
        subj = as.character(RID),
        tracer = as.factor(TRACER),
        site = as.factor("ADNI"),
        summary_cortical_amyloid = SUMMARY_SUVR,
        summary_cortical_volume = SUMMARY_VOLUME,
        amyloid_positive = AMYLOID_STATUS == 1,
        amyloid_positive.compref = AMYLOID_STATUS_COMPOSITE_REF == 1
    ) %>%
    select(subj, age, tracer, site, sex, apoe,
           summary_cortical_amyloid, summary_cortical_volume,
           amyloid_positive, amyloid_positive.compref,
           all_of(names(col_rename)),
           filepath.raw_amyloid, SCANDATE) %>%
    rename_with_list(col_rename) %>%
    mutate(filepath.raw_amyloid.json = str_replace(filepath.raw_amyloid, ".nii.gz", ".json"))  # get json file

# ====================
# ===== TIDY CDR =====
# ====================

cdr_df <- df_list[["cdr"]] %>%
    mutate(
        EXAMDATE = if_else(!is.na(EXAMDATE), EXAMDATE, USERDATE),
        cdr = if_else(CDGLOBAL < 0, NA, CDGLOBAL)
    ) %>%
    adni_compute_age("EXAMDATE") %>%
    rename(subj = RID, cdrsob = CDRSB) %>%
    mutate(subj = as.character(subj)) %>%
    select(subj, age, cdr, cdrsob)

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

# =================================
# ===== MERGE WITH OTHER DATA =====
# =================================

# match with raw T1 scans
pet_df <- pet_df %>% byyfunctions::match_data(
    match_df = mri_raw_df,
    ref_col = "SCANDATE",
    match_col = "SCANDATE",
    group_col = "subj",
    select_col = c("filepath", "accel", "tesla", "age"),
    suffix = ".raw_t1",
    diff_col = "age_diff.raw_t1"
) %>%
    mutate(
        age_diff.raw_t1 = as.numeric(age_diff.raw_t1) / 365.25,
        filepath.raw_t1.json = str_replace(filepath.raw_t1, ".nii.gz", ".json")  # get json file
    )

# ========================
# ===== SAVE OUTPUTS =====
# ========================

if (SAVE_CSV) {
    byyfunctions::make_dir(opt$odir)
    if (opt$csv) {
        write_csv(pet_df, file.path(opt$odir, "adni_amyloid.csv"))
        write_csv(cdr_df, file.path(opt$odir, "adni_cdr.csv"))
    } else {
        write_rds(pet_df, file.path(opt$odir, "adni_amyloid.RDS"))
        write_rds(cdr_df, file.path(opt$odir, "adni_cdr.RDS"))
    }
}
