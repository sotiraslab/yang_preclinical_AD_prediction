# ==============================================================================

# Name: merge_PAC.R
# Author: Braden Yang
# Created: 02/19/2024
# Description: Load PAC CSVs, preprocess and merge

# ==============================================================================

rm(list = ls())

SAVE_CSV <- TRUE

# ===========================
# ===== IMPORT PACKAGES =====
# ===========================

library(tidyverse)
library(optparse)
library(rjson)

library(devtools)

# ===========================
# ===== PARSE ARGUMENTS =====
# ===========================

option_list <- list(
    make_option(c("-w", "--wdir"), action="store", default=NULL,
        type="character", help="Path to project directory (if none, uses current working directory)"),
    make_option(c("-o", "--odir"), action="store", default="data/tidy/pac",
        type="character", help="Path to output directory"),
    make_option(c("-m", "--muse_odir"), action="store", default="data",
        type="character", help="Path to output MUSE feature names"),
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

pac_amyloid_positive <- function(.data, summary_col, site_col, analysis_col) {

    # mark amyloid positivity for PAC data, using both
    # site-specific thresholds and a pooled threshold

    .data <- .data %>%
        mutate(
            amyloid_positive.site = case_match(
                {{analysis_col}},
                "suvr" ~ case_match(
                        {{site_col}},
                        "AIBL" ~ {{summary_col}} >= 1.220,
                        "BIOCARD" ~ {{summary_col}} >= 1.220,
                        "BLSA" ~ {{summary_col}} >= 1.229,
                        "OASIS" ~ {{summary_col}} >= 1.236,
                        "WRAP" ~ {{summary_col}} >= 1.219
                    ),
                "dvr" ~ case_match(
                        {{site_col}},
                        "BIOCARD" ~ {{summary_col}} >= 1.043,
                        "BLSA" ~ {{summary_col}} >= 1.054,
                        "OASIS" ~ {{summary_col}} >= 1.076,
                        "WRAP" ~ {{summary_col}} >= 1.086
                    )
            ),
            amyloid_positive.pooled = case_match(
                {{analysis_col}},
                "suvr" ~ {{summary_col}} >= 1.237,
                "dvr" ~ {{summary_col}} >= 1.069
            )
        )

    return(.data)

}

# ============================
# ===== DEFINE VARIABLES =====
# ============================

pac_dir <- "/ceph/chpc/shared/aristeidis_sotiras_group/aris_data/PAC/tabular"
clinical_csv_path <- file.path(pac_dir, "ClinicalCognitiveMedical", "PAC_Demo_Clinical_Cognitive.csv")
pet_path <- file.path(pac_dir, "PET", "pac_pet_pib_50-70min_nopvc_cbgmref_suvr_imputed_20220307.csv")
mri_xlsx_path <- file.path(pac_dir, "MRI", "PAC_sMRI-Results_v1.3_Share.xlsx")

# ====================
# ===== MUSE ROI =====
# ====================

# here, we get important cortical and subcortical regions for input into predictive models, categorize them,
# and store them in an index table to use in downstream analyses

roi_dict <-
    readxl::read_excel(mri_xlsx_path, "Dictionary_ROI_Hierarchy") %>%  # read Dictionary_ROI_Hierarchy from MRI xlsx
    rename_with(str_to_lower, everything()) %>%  # lowercase colnames
    select(roi_index, hemisphere, tissue_seg, roi_name) %>%
    filter(tissue_seg == "GM" | is.na(tissue_seg)) %>%  # remove WM ROI
    filter(roi_name != "Brain Stem") %>%  # remove brain stem
    mutate(roi_type = case_when(  # categorize ROI
        str_detect(roi_name, "Cerebellar|Cerebellum") ~ "cerebellum",
        str_detect(roi_name, "Accumbens Area|Caudate|Pallidum|Putamen|Thalamus Proper|Amygdala|Basal Forebrain|Hippocampus") ~ "subcortical",
        str_detect(roi_name, "Ventricle|Lat Vent") ~ "ventricle",
        .default = "cortical"
    ))

# =================================
# ===== DEMOGRAPHIC/COGNITIVE =====
# =================================

clinical_df <- readr::read_csv(
    clinical_csv_path,
    show_col_types = FALSE
) %>%
    filter(ID != ".") %>%  # remove ID = "." from WRAP dataset
    rename(sex = Sex) %>%
    mutate(
        site = if_else(site == "ACS", "OASIS", site),  # replace "ACS" with "OASIS"
        ID = case_match(   # remove "OAS" prefix from IDs and convert to integer
            site,
            "OASIS" ~ str_remove(ID, "OAS"),
            "BIOCARD" ~ str_remove(ID, "JHU"),
            .default = ID
        ) %>% as.integer,
        subj = interaction(site, ID),
        age = AgeAtEntry + (MOFROMBL / 12),
        sex = factor(sex, levels = c(1, 2), labels = c("M", "F")),
        apoe = as.factor(
            case_when(
                apoe %in% c(34, 44) ~ 1,
                apoe %in% c(22, 23, 24, 33) ~ 0,
                .default = NA
            )
        ),
        pac_baseline = age - AgeAtPACBaseline,
        mci_ad = diagnosis %in% c("MCI", "Dementia")
    )

# 03/11/2024: I noticed a problem with the `MOFROMBL` column in certain rows of the clinical dataframe. It seems that in 9 entries from OASIS subjects, the value of `MOFROMBL` is too high to be plausible (in the range of 1382 to 1398). The data here was most likely entered incorrectly, so we will filter out any instances where `MOFROMBL` is greater than 500 (arbitrary threshold; the next highest was 345.8)
clinical_df <- clinical_df %>% filter(MOFROMBL < 500)

# remove any subject with diagnosis "Impaired not MCI"; these are likely people with non-AD related dementia
clinical_df <- clinical_df %>%
    group_by(subj) %>%
    filter(all(diagnosis != "Impaired not MCI" | is.na(diagnosis))) %>%
    ungroup()

# ===============
# ===== PET =====
# ===============

# there are slightly different column names for the same regions; make column names same here
#   - missing period in "posterior limb of internal capsule inc{.} ..."
#       - add additional period to suvr col
#   - `age` in suvr, but missing in dvr and r1
#       - we just repeat `AgeAtScan` as `age`
#   - note that dvr and r1 have exactly the same columns
#   - missing in suvr
#       - Scanner, Reconstruction, AC method, AC scan, Injected radioactivity (MBq), suvr5070.how
#   - missing in dvr/r1
#       - tracer, temporal-parietal-occipital GM, tau metaROI Jack, MTL, Braak I/II, Braak III/IV, Braak V/VI, FDG Landau metaROI
pet_df <- 
    read_csv(pet_path, show_col_types = FALSE) %>%
    mutate(subj = interaction(site, id)) %>%
    rename_with(
        .fn = ~ str_replace(.x, " inc ", " inc. "),
        .cols = starts_with("posterior limb of internal capsule")
    ) %>%
    mutate(
        pet_analysis = "suvr",
        tracer = "PIB",
        site = if_else(site == "OASIS3", "OASIS", site),
        petdate.mdy = mdy(petdate)
    ) %>%
    rename(summary_cortical_amyloid = `mean cortical`) %>%
    pac_amyloid_positive(  # mark amyloid positivity
        summary_cortical_amyloid,
        site,
        pet_analysis
    )

# ===============
# ===== MRI =====
# ===============

sheets <- readxl::excel_sheets(mri_xlsx_path)
mri_df <- map(
    sheets,
    ~ readxl::read_excel(mri_xlsx_path, .x, guess_max = 1048576)
) %>% set_names(sheets)
# wrong col data type: https://stackoverflow.com/questions/31633891/specifying-column-types-when-importing-xlsx-data-to-r-with-package-readxl

# get ROI volume df and tidy
vol_df <- mri_df[["ROI_Volume"]] %>%
    pivot_longer(
        cols = starts_with("H_MUSE_Volume_"),
        names_pattern = "H_MUSE_Volume_(.*)",
        names_to = "roi_idx",
        values_to = "vol"
    )

# map ROI idx to feature names
idx2roi_df <- mri_df[["Dictionary_ROI"]] %>%
    select(ROI_INDEX, ROI_NAME) %>%
    filter(ROI_INDEX != "ICV") %>%
    mutate(across(all_of("ROI_INDEX"), as.integer))
vol_df <- vol_df %>%
    mutate(
        roi_label = plyr::mapvalues(roi_idx, from = idx2roi_df$ROI_INDEX, to = idx2roi_df$ROI_NAME),
        PAC_ID = case_match(
            Study,
            "OASIS" ~ str_remove(PAC_ID, "OAS"),
            "BIOCARD" ~ str_remove(PAC_ID, "JHU"),
            .default = PAC_ID
        ) %>% as.integer,
        subj = interaction(Study, PAC_ID)
    ) %>%
    rename(age = AgeAtScan)

# convert back to wide format, standardize subj ID and column names
vol_wide <- vol_df %>%
    select(-roi_idx) %>%
    filter(roi_label %in% roi_dict$roi_name) %>%
    pivot_wider(
        names_from = roi_label,
        values_from = vol,
        names_glue = "{roi_label}.vol"
    ) %>%
    rename(ICV.vol = ICV)

# ===========================
# ===== MATCH AND MERGE =====
# ===========================

# match all subjects' demographics
subj_demographics <- clinical_df %>%
    select(subj, sex, race, apoe) %>%
    group_by(subj) %>%
    summarise(across(
        everything(),
        unique
    ))
pet_merge <- pet_df %>%
    left_join(subj_demographics, by = "subj")

# match MRI to each PET
# NOTE: PET metadata CSV has some info on the difference in dates between MRI and PET, but there
#   wasn't enough information to completely map each of these rows back to the PET data tables
#   - therefore, to match MRIs, we are just going to match by closest date from the CSVs, rather
#     than try to figure out which MRIs were used to process each PET (which would have been the
#     more ideal way to do it uwu)
pet_merge <- pet_merge %>%
    byyfunctions::match_data(
        match_df = vol_wide %>% filter(ROI_QCFlag != "FAIL"),  # exclude scans which failed QC (their data is NA anyways)
        ref_col = "age",
        match_col = "age",
        group_col = "subj",
        select_col = vol_wide %>% select(-c(subj, age)) %>% colnames,
        suffix = ".mri",
        max_diff = 1,
        diff_col = "age_diff.mri"
    )

# match diagnoses to each PET
pet_merge <- pet_merge %>%
    byyfunctions::match_data(
        match_df = clinical_df %>% filter(!is.na(diagnosis)),
        ref_col = "age",
        match_col = "age",
        group_col = "subj",
        select_col = "diagnosis",
        suffix = ".diagnosis",
        max_diff = 1,
        diff_col = "age_diff.diagnosis"
    )

# match cognitive scores to each PET (need to have all 3 composites + MMSE)
pet_merge <- pet_merge %>%
    byyfunctions::match_data(
        match_df = clinical_df %>%
            filter(
                if_all(starts_with("zf"), ~ !is.na(.x)),
                !is.na(MMSE)
            ),
        ref_col = "age",
        match_col = "age",
        group_col = "subj",
        select_col = c("MMSE", "zfgcp", "zfmem", "zfexf"),
        suffix = ".cognitive",
        max_diff = 1,
        diff_col = "age_diff.cognitive"
    )

# ===============================
# ===== MATCH RAW PET SCANS =====
# ===============================

oasis_age_at_entry <- read_csv(file.path(pac_dir, "ClinicalCognitiveMedical", "PAC_ACS_Demo_Clinical_Cognitive.csv"), show_col_types = FALSE) %>%
    mutate(
        id = OASIS_Subject_ID %>% str_replace("OAS","") %>% as.numeric,
        AgeAtEntry.oasis = AgeAtEntry
    ) %>%
    select(id, AgeAtEntry.oasis) %>%
    unique()

pac_amyloid_dir <- "/ceph/chpc/shared/aristeidis_sotiras_group/aris_data/PAC/images/nifti/pet/unprocessed"
pac_t1_dir <- "/ceph/chpc/shared/aristeidis_sotiras_group/aris_data/PAC/images/nifti/mri/rawdata"

add_age <- function(v) {
    if (!any(is.na(v))) {
        a <- as.numeric(v[[1]]) + as.numeric(v[[2]])/10
        return(a)
    } else {
        return(NA)
    }
}

raw_amyloid_scans <- read_csv("data/raw_img_paths/pac_amyloid.txt", col_names=FALSE, show_col_types = FALSE) %>%
    rename(filepath.raw_amyloid = X1) %>%
    mutate(
        filepath.raw_amyloid.json = filepath.raw_amyloid %>% str_replace(".nii","") %>% str_replace(".gz","") %>% str_c(".json"),
        site = filepath.raw_amyloid %>% str_replace(str_c(pac_amyloid_dir,"/"),"") %>% str_extract("[^/]+") %>% str_to_upper(),
        site = if_else(site == "OASIS3", "OASIS", site),
        filename = filepath.raw_amyloid %>% basename(),
        id = case_when(
            site %in% c("AIBL", "BLSA") ~ str_extract(filename, "^[^_]+_[^_]+") %>% str_split_i("_", 2),
            site %in% c("BIOCARD", "OASIS", "WRAP") ~ str_extract(filename, "^[^_]+")
        ),
        id = case_when(
            site == "BIOCARD" ~ str_replace(id, "JHU", ""),
            site == "OASIS" ~ str_replace(id, "sub-OAS", ""),
            .default = id
        ) %>% as.numeric,
        petdate = str_extract(filename, "\\d{2}-\\d{2}-\\d{4}"),
        petdate.mdy = mdy(petdate),
        age.wrap = str_extract(filename, "\\d+p\\d") %>% str_split("p") %>% map_dbl(add_age),
        ses.oasis = str_extract(filename, "ses-d\\d+") %>% str_replace("ses-", ""),
        ses.oasis.numeric = as.numeric(ses.oasis %>% str_replace("d",""))
    ) %>%
    left_join(oasis_age_at_entry, by = join_by(id == id)) %>%
    mutate(
        age.oasis = if_else(site=="OASIS", AgeAtEntry.oasis + (ses.oasis.numeric/365) %>% round(1), NA),
        age.join = coalesce(age.wrap, age.oasis)
    )

petdate_match <- pet_merge %>%
    filter(site %in% c("AIBL", "BLSA", "BIOCARD")) %>%
    left_join(
        raw_amyloid_scans %>% select(filepath.raw_amyloid, filepath.raw_amyloid.json, site, id, petdate.mdy),
        by = join_by(
            site == site,
            id == id,
            petdate.mdy == petdate.mdy
        )
    )
age_match <- pet_merge %>%
    filter(site %in% c("WRAP")) %>%
    mutate(age.join = round(age, 1)) %>%
    left_join(
        raw_amyloid_scans %>% select(filepath.raw_amyloid, filepath.raw_amyloid.json, site, id, age.join),
        by = join_by(
            site == site,
            id == id,
            age.join == age.join
        )
    )
fuzzy_age_match <- pet_merge %>%
    filter(site %in% c("OASIS")) %>%
    mutate(age.join = round(age, 1)) %>%
    byyfunctions::match_data(
        match_df = raw_amyloid_scans %>% filter(site %in% c("OASIS"), !str_detect(filepath.raw_amyloid, "AV45")),   # for OASIS, exclude AV45 scans from matching
        ref_col = "age.join",
        match_col = "age.join",
        group_col = "id",
        select_col = c("filepath.raw_amyloid", "filepath.raw_amyloid.json"),
        diff_col = "age_diff.raw_amyloid",
        suffix = ".raw_amyloid"
    ) %>%
    rename(
        filepath.raw_amyloid = filepath.raw_amyloid.raw_amyloid,
        filepath.raw_amyloid.json = filepath.raw_amyloid.json.raw_amyloid
    ) %>%
    select(-c(age.join.raw_amyloid, age_diff.raw_amyloid))

pet_merge <- bind_rows(petdate_match, age_match, fuzzy_age_match) %>% select(-age.join)

# ===============================
# ===== MATCH RAW MRI SCANS =====
# ===============================

# NOTE: BLSA MRI scans are named differently from other sites, and as a result we're not
# able to pull either a PAC subject ID or scan date; all matching scans are "NA" for now 

raw_t1_scans <- read_csv("data/raw_img_paths/pac_t1.txt", col_names=FALSE, show_col_types = FALSE) %>%
    rename(filepath.raw_t1 = X1) %>%
    mutate(
        site = filepath.raw_t1 %>% str_replace(str_c(pac_t1_dir,"/"),"") %>% str_extract("[^/]+") %>% str_to_upper(),
        filename = filepath.raw_t1 %>% basename,
        folder = filepath.raw_t1 %>% dirname
    ) %>%
    mutate(
        id = case_when(
            site %in% c("BLSA") ~ str_extract(filename, "^[^_]+_[^_]+") %>% str_split_i("_", 2),
            site %in% c("AIBL", "BIOCARD", "OASIS", "WRAP") ~ str_extract(filename, "^[^_]+")
        ),
        id = case_when(
            site == "BIOCARD" ~ str_replace(id, "JHU", ""),
            site == "OASIS" ~ str_replace(id, "OAS", ""),
            .default = id
        ) %>% as.numeric,
        mridate = str_extract(filename, "\\d{8}"),
        mridate.mdy = ymd(mridate),
        age.wrap = str_extract(filename, "\\d+p\\d") %>% str_split("p") %>% map_dbl(add_age),
        ses.oasis = str_extract(filename, "d\\d+"),
        ses.oasis.numeric = as.numeric(ses.oasis %>% str_replace("d",""))
    ) %>%
    left_join(oasis_age_at_entry, by = join_by(id == id)) %>%
    mutate(
        age.oasis = if_else(site=="OASIS", AgeAtEntry.oasis + (ses.oasis.numeric/365) %>% round(1), NA),
        age.raw_t1 = coalesce(age.wrap, age.oasis)
    )
muse_segmentations <- read_csv("data/raw_img_paths/pac_muse.txt", col_names=FALSE, show_col_types = FALSE) %>%
    rename(filepath.muse_segmentations = X1) %>%
    mutate(folder = filepath.muse_segmentations %>% dirname)
dlicv_masks <- read_csv("data/raw_img_paths/pac_dlicv.txt", col_names=FALSE, show_col_types = FALSE) %>%
    rename(filepath.dlicv_mask = X1) %>%
    mutate(folder = filepath.dlicv_mask %>% dirname)

raw_t1_scans <- raw_t1_scans %>%
    full_join(muse_segmentations, by = "folder") %>%
    full_join(dlicv_masks, by = "folder")

petdate_match <- pet_merge %>%
    filter(site %in% c("AIBL", "BIOCARD")) %>%
    byyfunctions::match_data(
        match_df = raw_t1_scans %>% filter(site %in% c("AIBL", "BIOCARD")),
        ref_col = "petdate.mdy",
        match_col = "mridate.mdy",
        group_col = "id",
        select_col = c("filepath.raw_t1", "filepath.muse_segmentations", "filepath.dlicv_mask"),
        diff_col = "age_diff.raw_t1",
        suffix = ""
    ) %>%
    mutate(
        age_diff.raw_t1 = as.numeric(age_diff.raw_t1) / 365.25,
        age.raw_t1 = age + (as.numeric(mridate.mdy - petdate.mdy) / 365.25)
    )
age_match <- pet_merge %>%
    filter(site %in% c("WRAP", "OASIS")) %>%
    mutate(age.join = round(age, 1)) %>%
    byyfunctions::match_data(
        match_df = raw_t1_scans %>% filter(site %in% c("WRAP", "OASIS")),
        ref_col = "age.join",
        match_col = "age.raw_t1",
        group_col = "id",
        select_col = c("filepath.raw_t1", "filepath.muse_segmentations", "filepath.dlicv_mask"),
        diff_col = "age_diff.raw_t1",
        suffix = ""
    )
blsa_codes <- read_csv("/home/b.y.yang/BLSA_PACIDs.csv")
subj_match <- pet_merge %>%
    filter(site %in% c("BLSA")) %>%
    left_join(blsa_codes, by = c("id" = "PAC_ID")) %>%
    left_join(
        raw_t1_scans %>% select(site, id, contains("filepath")),
        by = c("site" = "site", "BLSA_ID_NoPadding" = "id")
    ) %>%
    mutate(age.raw_t1 = id)    # no age at T1 provided; use ID as placeholder (only one T1 per BLSA subject, so should be unique)

pet_merge <- bind_rows(petdate_match, age_match, subj_match) %>% select(-age.join)

# =================================
# ===== PROCESS AND TIDY DATA =====
# =================================

#' Here is the rough procedure for tidying and processing data
#' 
#' 1. select the non-imaging and imaging columns of interest (use regex)
#' 2. filter to only select SUVR (also repeat with DVR and R1)
#' 3. rename columns

# non-imaging columns to select
nonimg_col <- c("subj", "site", "tracer", "pet_analysis", "PAC_VisitNo",
    "age", "sex", "apoe", "diagnosis.diagnosis", "amyloid_positive.site", "amyloid_positive.pooled")
cog_col <- c("MMSE.cognitive", "zfgcp.cognitive", "zfmem.cognitive", "zfexf.cognitive")

# for imaging features, create regex
ctx_subctx_roi <- roi_dict %>% filter(roi_type %in% c("cortical", "subcortical")) %>% pull(roi_name)
ctx_subctx_regex <- str_c("^", ctx_subctx_roi, "(|\\.vol\\.mri)$")
ctx_subctx_regex <- str_c(ctx_subctx_regex, collapse = "|")

# cerebellum_roi <- roi_dict %>% filter(roi_type == "cerebellum") %>% pull(roi_name)
# cerebellum_regex <- str_c("^", cerebellum_roi, "(|\\.vol\\.mri)$")
# cerebellum_regex <- str_c(cerebellum_regex, collapse = "|")

ventricle_roi <- roi_dict %>% filter(roi_type == "ventricle") %>% pull(roi_name)
ventricle_regex <- str_c("^", ventricle_roi, "\\.vol\\.mri$")
ventricle_regex <- str_c(ventricle_regex, collapse = "|")

# select_regex <- str_c(ctx_subctx_regex, cerebellum_regex, ventricle_regex, sep = "|")
select_regex <- str_c(ctx_subctx_regex, ventricle_regex, sep = "|")

tidy_df <- pet_merge %>%
    select(  # select relevant features
        all_of(nonimg_col),
        all_of(cog_col),
        contains("filepath"),
        contains(".raw"),
        summary_cortical_amyloid,
        ICV.vol.mri,
        matches(select_regex)
    ) %>%
    rename_with(
        ~ str_c(.x, ".amyloid"),
        all_of(ctx_subctx_roi)
    ) %>%
    rename_with(
        ~ str_replace_all(.x, ".mri", ""),
        ends_with(".vol.mri")
    ) %>%
    # rename_with(
    #     ~ str_c(.x, ".ref"),
    #     all_of(cerebellum_roi)
    # ) %>%
    rename_with(~ str_remove(.x, ".cognitive"), ends_with(".cognitive")) %>%
    rename_with(~ str_remove(.x, ".diagnosis"), ends_with(".diagnosis"))

# # normalize SUVR
# pet_col <- tidy_df %>% select(ends_with(".amyloid")) %>% colnames
# tidy_df <- tidy_df %>%
#     mutate(
#         .reference = rowMeans(pick(ends_with(".ref"))),
#         across(
#             all_of(c("summary_cortical_amyloid", pet_col)),
#             ~ if_else(  # just pick out SUVR to normalize, don't normalize DVR or R1
#                 pet_analysis == "suvr",
#                 .x / .reference,
#                 .x
#             )
#         )
#     ) %>%
#     select(-.reference) %>%
#     rename(mmse = MMSE)

# tidy diagnosis df
clinical_tidy <- clinical_df %>%
    select(site, subj, age, diagnosis)

# remove OASIS from tables (since we're using the openly available OASIS)
tidy_df <- tidy_df %>% filter(site != "OASIS")
clinical_tidy <- clinical_tidy %>% filter(site != "OASIS")

# match diagnosis
tidy_df <- tidy_df %>% byyfunctions::match_data(
    match_df = clinical_tidy,
    ref_col = "age",
    match_col = "age",
    group_col = "subj",
    select_col = "diagnosis",
    suffix = ".diagnosis",
    diff_col = "age_diff.diagnosis"
)

# add column to indicate site is in PAC
tidy_df <- tidy_df %>% mutate(site_in_pac = TRUE)

# =====================
# ===== SAVE DATA =====
# =====================

if (SAVE_CSV) {
    if (opt$csv) {
        write_csv(tidy_df, file.path(opt$odir, "pac_amyloid.csv"))
        write_csv(clinical_tidy, file.path(opt$odir, "pac_diagnosis.csv"))
    } else {
        write_rds(tidy_df, file.path(opt$odir, "pac_amyloid.RDS"))
        write_rds(clinical_tidy, file.path(opt$odir, "pac_diagnosis.RDS"))
    }

    # output column names
    amyloid_col <- tidy_df %>% select(ends_with(".amyloid")) %>% colnames %>% str_replace(".amyloid", "") %>% tibble()
    vol_col <- tidy_df %>% select(ends_with(".vol"), -ICV.vol) %>% colnames %>% str_replace(".vol", "") %>% tibble()
    write_csv(amyloid_col, file.path(opt$muse_odir, "MUSE_regions_amyloid.csv"), col_names = FALSE)
    write_csv(vol_col, file.path(opt$muse_odir, "MUSE_regions_vol.csv"), col_names = FALSE)
}
