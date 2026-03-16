# ==============================================================================

# Name: merge_OASIS.R
# Author: Braden Yang
# Created: 06/05/2023
# Description: get cross-sectional amyloid PET features, tidy, and match clinical
#   and other non-imaging data

# ==============================================================================

# clear environment
rm(list = ls())

# *** toggle variables ***
INTERACTIVE <- FALSE
SAVE_CSV <- TRUE
# *** toggle variables ***

# ===========================
# ===== IMPORT PACKAGES =====
# ===========================

library(rjson)
library(tidyverse)
library(optparse)
library(devtools)

# ===========================
# ===== PARSE ARGUMENTS =====
# ===========================

option_list <- list(
    make_option(c("-w", "--wdir"), action="store", default=NULL,
        type="character", help="Path to project directory (if none, uses current working directory)"),
    make_option(c("-o", "--odir"), action="store", default="data/tidy/oasis",
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

rename_with_list <- function(.data, l) {
    .data %>% rename_with(~ l[.x] %>% unlist(), .cols = any_of(names(l)))
}

oasis_metadata_from_id <- function(.data, id_col, sep = "_", num_sep = 3) {

    # Return dataframe with new columns "Subject" and "daysFromEntry",
    # infered from the column "PUP_PUPTIMECOURSEDATA.ID"
    # 
    # Parameters
    # ----------
    #   .data:   OASIS-3 dataframe of PUP data (loaded from `oasis3_pup_full.csv`)
    #   id_col:  column containing unique subject/session ID (use tidyverse data
    #            masking when inputing column, i.e. no quotes)
    #   sep:     separating character (default = "_")
    #   num_sep: how many substrings result after splitting by sep?
    # 
    # Returns
    # -------
    #   dataframe with new columns

    # split string into components
    id_split <- .data %>%
        pull({{id_col}}) %>%
        str_split_fixed(pattern = sep, n = num_sep)

    # get subject ID as new column "Subject"
    .data$Subject <- id_split[, 1]

    # get days from entry as new column "daysFromEntry"
    # remove first character of every string: https://stackoverflow.com/questions/35113553/r-remove-first-character-from-string
    .data$daysFromEntry <- as.integer(sub(".", "", id_split[, num_sep]))

    return(.data)

}

oasis_match_features <- function(pet_df, cdr_df, fs_df, demog_df) {

    # ===== 1. Get APOE, age and sex =====

    pet_df <- demog_df %>%
        rename(Subject = OASISID, apoe = APOE, ageAtEntry = AgeatEntry, sex = GENDER) %>%
        mutate(
            sex = factor(sex, levels = c(1, 2), labels = c("M", "F")),
            apoe = as.factor(
                case_when(
                    apoe %in% c("34", "44") ~ 1,
                    apoe %in% c("22", "23", "24", "33") ~ 0,
                    .default = NA
                )
            ),
        ) %>%
        group_by(Subject) %>%
        slice_head(n = 1) %>%
        select(Subject, apoe, ageAtEntry, sex) %>%
        right_join(pet_df, by = "Subject", multiple = "all")

    # ===== 2. Merge FreeSurfer volumetric data =====

    pet_df <- fs_df %>%  # then merge
        select(FSId, all_of(names(vol_col_rename))) %>%
        right_join(pet_df, by = "FSId", multiple = "all")

    # ===== 3. Compute age at visit =====

    pet_df <- pet_df %>%
        mutate(ageAtVisit = ageAtEntry + (daysFromEntry) / (365.25))

    # ===== 4. Match CDR and MMSE =====

    pet_df <- pet_df %>%
        byyfunctions::match_data(
            match_df = cdr_df,
            ref_col = "daysFromEntry",
            match_col = "daysFromEntry",
            select_col = c("cdr", "mmse"),
            group_col = "Subject",
            suffix = ".cdr"
        ) %>%
        rename(
            cdr = cdr.cdr,
            mmse = mmse.cdr
        )

    # ===== 5. Rename columns =====

    pet_df <- pet_df %>%
        rename_with_list(col_rename) %>%
        rename(
            subj = Subject,
            age = ageAtVisit,
            summary_cortical_amyloid = PET_fSUVR_TOT_CORTMEAN
        ) %>%
        select(-c(ageAtEntry, FSId), -starts_with("daysFromEntry"), daysFromEntry)

    return(pet_df)

}

oasis_compute_age <- function(
    .data,
    demog_df,
    days_from_entry_col = "daysFromEntry",
    group_col = "Subject"
) {

    demog_merge <- demog_df %>%
        rename_with(~ group_col, OASISID) %>%
        select(all_of(group_col), AgeatEntry)
    .data <- .data %>%
        left_join(demog_merge, by = group_col) %>%
        mutate(age = AgeatEntry + (.data[[days_from_entry_col]] / 365.25)) %>%
        select(-AgeatEntry)

    return(.data)

}

mark_amyloid_positive <- function(suvr, tracer, site) {

    # for ADNI, use the summary SUVR w/whole cerebellum normalization
    # for OASIS, use mean cortical SUVR

    amyloid_positive <- case_when(
        site == "ADNI" & tracer == "AV45" ~ suvr >= 1.11,
        site == "ADNI" & tracer == "FBB" ~ suvr >= 1.08,
        site == "OASIS" & tracer == "AV45" ~ suvr >= 1.24,
        site == "OASIS" & tracer == "PIB" ~ suvr >= 1.31,
        .default = NA
    )

    return(amyloid_positive)

}

assign_clinical_group <- function(cdr_vec) {

    group_vec <- dplyr::case_when(
        {{cdr_vec}} == 0 ~ "CN",
        {{cdr_vec}} == 0.5 ~ "MCI",
        {{cdr_vec}} > 0.5 ~ "AD",
        .default = NA 
    ) %>%
    factor(levels = c("CN", "MCI", "AD"))

    return(group_vec)

}

extend_clinical_group <- function(clinical_group, amyloid_positive) {

    group_vec <- dplyr::if_else(
        {{clinical_group}} == "CN" & {{amyloid_positive}},
        "preclinical",
        {{clinical_group}}
    ) %>%
    factor(levels = c("CN", "preclinical", "MCI", "AD"))

    return(group_vec)

}

normalize_by_reference <- function(.data, cols_to_norm, ref_cols) {

    # Divide a set of columns (cols_to_norm) by the average of the list of
    # columns defined in ref_cols, return the normalized table

    # summing multiple columns: https://stackoverflow.com/questions/47759347/create-a-new-column-which-is-the-sum-of-specific-columns-selected-by-their-name

    .data <- .data %>%
        mutate(ref = reduce(select(., all_of(ref_cols)), `+`) / length(ref_cols)) %>%
        mutate_at(cols_to_norm, ~ . / ref) %>%
        select(-ref)
    
    return(.data)

}

# =======================
# ===== LOAD PARAMS =====
# =======================

params <- fromJSON(file = "params.json")
oasis_paths <- params$OASIS_paths

# =======================
# ===== LOAD TABLES =====
# =======================

# get regions to select
oasis_df <- purrr::map(oasis_paths, ~read_csv(.x, show_col_types = FALSE))
roi_features_df <- read_csv("data/adni_oasis_roi_features.csv", col_types = "ccccc")
pet_col_rename <- roi_features_df %>% filter(!is.na(oasis_suvr)) %>% pull(fs_label) %>% str_c(".amyloid")
names(pet_col_rename) <- roi_features_df %>% filter(!is.na(oasis_suvr)) %>% pull(oasis_suvr)
vol_col_rename <- roi_features_df %>% filter(!is.na(oasis_vol)) %>% pull(fs_label) %>% str_c(".vol")
names(vol_col_rename) <- roi_features_df %>% filter(!is.na(oasis_vol)) %>% pull(oasis_vol)
col_rename <- c(pet_col_rename, vol_col_rename)

# get list of raw PET imaging
oasis_data_dir <- "/ceph/chpc/rcif_datasets/oasis/OASIS3"
pet_raw_df <- read_csv("data/raw_img_paths/oasis_amyloid.txt", col_names = "filepath", show_col_types = FALSE) %>%
    mutate(.temp = str_replace(filepath, str_c(oasis_data_dir, "/"), "")) %>%
    separate(.temp, into = c("scan_id", ".rm1", ".rm2", ".rm3", "filename"), sep = "/") %>%
    oasis_metadata_from_id("scan_id") %>%
    mutate(
        tracer = case_when(
            str_detect(filename, "AV45") ~ "AV45",
            str_detect(filename, "PIB") ~ "PIB",
            .default = NA
        )
    ) %>%
    rename(subj = Subject) %>%
    select(subj, daysFromEntry, tracer, filepath)
mri_raw_df <- read_csv("data/raw_img_paths/oasis_t1.txt", col_names = "filepath", show_col_types = FALSE) %>%
    filter(!str_detect(filepath, "hippocampus"), !str_detect(filepath, "echo")) %>%  # remove hippocampal acquisitions and separate echo scans
    mutate(.temp = str_replace(filepath, str_c(oasis_data_dir, "/"), "")) %>%
    separate(.temp, into = c("scan_id", ".rm1", ".rm2", ".rm3", "filename"), sep = "/") %>%
    oasis_metadata_from_id("scan_id") %>%
    rename(subj = Subject) %>%
    mutate(
        run = str_extract(filename, "(?<=_run-)\\d{1,2}") %>% as.numeric,
        run = if_else(is.na(run), 1, run)
    ) %>%
    select(subj, daysFromEntry, run, filepath) %>%
    group_by(subj, daysFromEntry) %>%
    slice_max(run) %>%  # select the latest run 
    ungroup() %>%
    oasis_compute_age(demog_df = oasis_df[["demog"]], group_col = "subj")

# ==================================
# ===== EXTRACT OASIS FEATURES =====
# ==================================

# get new columns (Subject, daysFromEntry) in clinical_df and pup_full_df
pet_df <- oasis_df[["pet"]] %>%
    oasis_metadata_from_id(id_col = `PUP_PUPTIMECOURSEDATA ID`, num_sep = 4)
cdr_df <- oasis_df[["cdr"]] %>%
    rename(Subject = OASISID, daysFromEntry = days_to_visit, cdr = CDRTOT, mmse = MMSE)

# get new columns in FS df and rename columns
oasis_fs_df <- oasis_df[["freesurfer"]] %>%
    rename(FSId = `FS_FSDATA ID`) %>%
    mutate(
        daysFromEntry = stringr::str_split_fixed(FSId, pattern = "_", n=3)[,3] %>% sub(".", "", .) %>% as.integer()
    )

# select relevant cols in pet_df
pet_df <- pet_df %>%
    select(
        Subject, daysFromEntry, FSId, tracer,
        all_of(names(pet_col_rename)),
        PET_fSUVR_TOT_CORTMEAN
    )

# match features
pet_df <- pet_df %>%
    oasis_match_features(
        cdr_df = cdr_df,
        fs_df = oasis_fs_df,
        demog_df = oasis_df[["demog"]]
    )

# tidy CDR dataframe
cdr_df <- cdr_df %>%
    rename(subj = Subject, age = `age at visit`, cdrsob = CDRSUM) %>%
    select(subj, age, cdr, cdrsob)

# match with raw images
pet_df <- pet_df %>%
    left_join(
        pet_raw_df %>% rename(filepath.raw_amyloid = filepath),
        by = c("subj", "daysFromEntry", "tracer")
    ) %>%
    mutate(filepath.raw_amyloid.json = 
        str_replace(filepath.raw_amyloid, ".nii.gz", ".json") %>%
        str_replace("/NIFTI/", "/BIDS/")
    )  # get json file
pet_df <- pet_df %>% byyfunctions::match_data(
    match_df = mri_raw_df,
    ref_col = "daysFromEntry",
    match_col = "daysFromEntry",
    group_col = "subj",
    select_col = c("filepath", "age"),
    suffix = ".raw_t1",
    diff_col = "age_diff.raw_t1"
) %>%
    mutate(
        age_diff.raw_t1 = as.numeric(age_diff.raw_t1) / 365.25,
        filepath.raw_t1.json = str_replace(filepath.raw_t1, ".nii.gz", ".json")  # get json file
    )

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

# ====================================
# ===== ADDITIONAL PREPROCESSING =====
# ====================================

pet_df <- pet_df %>%
    mutate(site = "OASIS") %>%
    filter(if_all(c(age, sex, cdr), ~!is.na(.x))) %>%
    mutate(
        amyloid_positive = mark_amyloid_positive(summary_cortical_amyloid, tracer, site),  # mark amyloid-positivity by summary cortical region
        clinical_group = assign_clinical_group(cdr),  # mark clinical groups by CDR
        clinical_group_extended = extend_clinical_group(clinical_group, amyloid_positive)  # mark preclinical subjects
    ) %>%
    normalize_by_reference(  # normalize by reference region 
        unname(pet_col_rename),  # don't normalize summary; this is already normalized for us
        c("left.cerebellum.cortex.amyloid", "right.cerebellum.cortex.amyloid")  # average of left and right cerebellar cortex as ref
    ) %>%
    mutate(tracer = if_else(tracer == "AV45", "FBP", tracer))  # replace "AV45" with "FBP"

# =======================
# ===== SAVE TABLES =====
# =======================

if (SAVE_CSV) {
    byyfunctions::make_dir(opt$odir)
    if (opt$csv) {
        write_csv(pet_df, file.path(opt$odir, "oasis_amyloid.csv"))
        write_csv(cdr_df, file.path(opt$odir, "oasis_cdr.csv"))
    } else {
        write_rds(pet_df, file.path(opt$odir, "oasis_amyloid.RDS"))
        write_rds(cdr_df, file.path(opt$odir, "oasis_cdr.RDS"))
    }
}