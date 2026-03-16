# ==============================================================================

# Name: merge_outputs.R
# Author: Braden Yang
# Created: 12/31/2024
# Description: Load MUSE and PET processing outputs, merge with main multisite
#   dataframe

# ==============================================================================

rm(list = ls()) 

SAVE_CSV <- TRUE
INTERACTIVE <- FALSE

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
    make_option(c("-o", "--odir"), action="store", default="data/features",
        type="character", help="Path to output directory"),
    make_option(c("-r", "--reload"), action="store_true", default=FALSE,
        help="if specified, reload all generated SUVR and MUSE CSVs")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (INTERACTIVE) {
    opt$wdir <- "/home/b.y.yang/sotiraslab/BradenADLongitudinalPrediction"
    opt$odir <- file.path(opt$wdir, "data/features")
    opt$reload <- TRUE
}

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

source("code/functions.R")

read_csv_batches <- function(read_csv_func, csv_list, batch_size = 1000) {
    # split read_csv into batches (https://stackoverflow.com/questions/74107723/read-csv-crashes-because-too-many-files-open)    
    batch_no <- (seq_along(csv_list) - 1) %/% batch_size
    split(csv_list, batch_no) %>%
        map(read_csv_func) %>%
        list_rbind()
}

# =====================
# ===== LOAD DATA =====
# =====================

multisite_df <- read_rds("/home/b.y.yang/sotiraslab/BradenADLongitudinalPrediction/data/subj/multisite.RDS")
muse_regions_amyloid <- read_csv("/home/b.y.yang/sotiraslab/BradenADLongitudinalPrediction/data/MUSE_regions_amyloid.csv", col_names=FALSE, show_col_types=FALSE) %>% pull()
muse_regions_vol <- read_csv("/home/b.y.yang/sotiraslab/BradenADLongitudinalPrediction/data/MUSE_regions_vol.csv", col_names=FALSE, show_col_types=FALSE) %>% pull()

# ================================
# ===== PET PIPELINE OUTPUTS =====
# ================================

get_suvr_csv <- function(d) {
    f <- list.files(path=d, pattern=".*ROI_stats_mean\\.csv", full.names=TRUE)
    if (length(f) == 0) {
        return(NA)
    } else {
        return(f[[1]])
    }
}

pet_analysis_wdir <- "/ceph/chpc/shared/aristeidis_sotiras_group/b.y.yang_scratch/pet_analysis"
suvr_csv_path <- file.path(pet_analysis_wdir, "suvr.csv")

if ( opt$reload | !file.exists(suvr_csv_path) ) {

    subj_dir_list <- list.files(path = file.path(pet_analysis_wdir, "nipype"), full.names = TRUE)
    suvr_csv_list <- map(file.path(subj_dir_list, "wf/suvr_csv/suvr_csv"), get_suvr_csv) %>% unlist()

    # failed instances
    # TODO: double check these subjects
    pet_failed_idx <- is.na(suvr_csv_list)
    pet_failed_subj <- subj_dir_list[pet_failed_idx]
    suvr_csv_list_passed <- suvr_csv_list[!pet_failed_idx]

    suvr_df <- read_csv_batches(
        ~read_csv(.x, skip=1, id="suvr_csv_path", show_col_types=FALSE),
        suvr_csv_list_passed
    )

    write_csv(suvr_df, suvr_csv_path)

} else {

    suvr_df <- read_csv(suvr_csv_path)

}

# filter and tidy SUVR outputs from PET pipeline
suvr_df_select <- suvr_df %>%
    select(id, `mean cortical`, all_of(muse_regions_amyloid)) %>%
    rename(idvi = id, summary_cortical_amyloid = `mean cortical`) %>%
    mutate(.proc_pass_amyloid = if_all(all_of(muse_regions_amyloid), ~ !is.na(.))) %>%
    rename_with(~str_c(.x, ".amyloid"), all_of(muse_regions_amyloid))

# ========================
# ===== MUSE OUTPUTS =====
# ========================

muse_wdir <- "/ceph/chpc/shared/aristeidis_sotiras_group/b.y.yang_scratch/T1_MUSE"
muse_seg_dir <- file.path(muse_wdir, "Protocols/MUSE")
mri_xlsx_path <- "/ceph/chpc/shared/aristeidis_sotiras_group/aris_data/PAC/tabular/MRI/PAC_sMRI-Results_v1.3_Share.xlsx"

# get MUSE idx to ROI mappings
sheets <- readxl::excel_sheets(mri_xlsx_path)
mri_df <- map(
    sheets,
    ~ readxl::read_excel(mri_xlsx_path, .x, guess_max = 1048576)
) %>% set_names(sheets)
idx2roi_df <- mri_df[["Dictionary_ROI"]] %>%
    select(ROI_INDEX, ROI_NAME) %>%
    filter(ROI_INDEX != "ICV") %>%
    mutate(across(all_of("ROI_INDEX"), as.integer)) %>%
    add_row(ROI_INDEX = 702, ROI_NAME = "ICV")

muse_csv_path <- file.path(muse_wdir, "muse_vol.csv")
if ( opt$reload | !file.exists(muse_csv_path) ) {

    # load all MUSE output files and tidy
    muse_subj_list <- list.files(path = muse_seg_dir, full.names = TRUE)
    muse_subj_list <- muse_subj_list[!(str_detect(muse_subj_list, "_reprocess"))]
    muse_list <- map(muse_subj_list, ~list.files(.x, pattern = "*.csv", recursive = TRUE, full.names = TRUE)) %>% unlist()

    muse_df <- read_csv_batches(
        ~read_csv(.x, id="muse_csv_path", show_col_types=FALSE),
        muse_list
    )

    muse_df_newnames <- muse_df %>%
        mutate(ID = muse_csv_path %>% str_replace(str_c(muse_seg_dir,"/"), "") %>% dirname) %>%
        pivot_longer(-c(ID, muse_csv_path)) %>%
        mutate(name = plyr::mapvalues(name, from = idx2roi_df$ROI_INDEX, to = idx2roi_df$ROI_NAME)) %>%  # rename columns from idx to ROI name
        pivot_wider(names_from = name, values_from = value)

    write_csv(muse_df_newnames, muse_csv_path)

} else {

    muse_df_newnames <- read_csv(muse_csv_path)

}

muse_df_select <- muse_df_newnames %>%
    rename(site.subj.age = ID) %>%
    select(site.subj.age, ICV, all_of(muse_regions_vol)) %>%
    mutate(.proc_pass_t1 = if_all(all_of(muse_regions_vol), ~ !is.na(.))) %>%
    rename_with(~str_c(.x,".vol"), -c(site.subj.age, .proc_pass_t1))

# ===========================
# ===== QUALITY CONTROL =====
# ===========================

qc_df <- read_csv("code/2_Preprocessing/qc/qc.csv", show_col_types = FALSE) %>%
    rename(idvi = id)

# =================
# ===== MERGE =====
# =================

# merge PET pipeline SUVR outputs
multisite_merge <- left_join(
    multisite_df,
    suvr_df_select %>% rename_with(
        ~str_c(.x,".pipeline"),
        -c(idvi, .proc_pass_amyloid)
    ),
    by="idvi"
)

# merge MUSE volume outputs
multisite_merge <- left_join(
    multisite_merge,
    muse_df_select %>% rename_with(
        ~str_c(.x,".pipeline"),
        -c(site.subj.age, .proc_pass_t1)
    ),
    by=c("site.subj.age.raw_t1" = "site.subj.age")  # merge on T1's site.subj.age code
)

# merge QC
multisite_merge <- left_join(
    multisite_merge, qc_df,
    by="idvi"
)

# =====================================
# ===== TIDY FOR 3_TRAINTESTMODEL =====
# =====================================

# get final data for python scripts
features_df <- multisite_merge %>%
    select(site, site_in_pac, idvi, subj, age, age_diff.raw_t1, tracer, starts_with("qc"), starts_with("."),
           sex, apoe,
           summary_cortical_amyloid.pipeline, ICV.vol.pipeline,
           contains("amyloid_positive"),
           all_of(str_c(muse_regions_amyloid,".amyloid.pipeline")),
           all_of(str_c(muse_regions_vol,".vol.pipeline"))) %>%
    rename_with(  # remove ".pipeline" suffix
        ~str_replace(.x, "\\.pipeline$", ""),
        ends_with(".pipeline")
    ) %>%
    rename(ICV = ICV.vol)  # rename ICV column

# filtering steps

# 1. drop rows with missing values in any of the PET SUVR, MUSE volume columns, or [age, sex, apoe]
features_df_filter <- features_df %>%
    drop_na(all_of(str_c(muse_regions_amyloid,".amyloid"))) %>%
    drop_na(all_of(str_c(muse_regions_vol,".vol"))) %>%
    drop_na(age, sex, apoe)

# 2. select scans that passed QC (excluding PAC sites)
features_df_filter <- features_df_filter %>%
    filter(qc > 0)  # include qc==1 and qc==0.5

# 3. remove scans where matched T1 is greater than 1yr apart from PET
features_df_filter <- features_df_filter %>%
    filter(age_diff.raw_t1 <= 1 | site_in_pac)

# 4. remove BIOCARD and WRAP (too few scans, just exclude these sites)
features_df_filter <- features_df_filter %>%
    filter(!(site %in% c("BIOCARD", "WRAP")))

# additional processing
features_df_filter <- features_df_filter %>%
    mutate(across(ends_with(".vol"), ~ . / ICV))  # normalize volumetric features by ICV

# get cross-sectional data
features_cross_sectional <- features_df_filter %>%
    group_by(subj) %>%
    slice_min(age) %>%
    ungroup()

# =====================
# ===== SAVE DATA =====
# =====================

if (SAVE_CSV) {
    # tidied features: contains SUVR, volume and non-imaging features; also removed QC=0 scans
    write_csv(features_df_filter, file.path(opt$odir, "features.csv"))
    write_csv(features_cross_sectional, file.path(opt$odir, "features_cross_sectional.csv"))

    byyfunctions::make_dir("figures/flowchart")
    write_rds(multisite_merge, "figures/flowchart/for_flowchart.RDS")
}
