# ==============================================================================

# Name: get_proc_tables.R
# Author: Braden Yang
# Created: 02/23/2024
# Description: Create CSVs for MUSE and PET processing

# ==============================================================================

rm(list = ls())

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
    make_option(c("-o", "--odir"), action="store", default="data/subj/proc",
        type="character", help="Path to output directory")
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

# ======================
# ===== CREATE CSV =====
# ======================

# load data tables
multisite_df <- read_rds("data/subj/multisite.RDS")
multisite_progressor_stable <- read_rds("data/subj/multisite_progressor_stable.RDS")

# get all scans that need to be processed
idvi_to_process <- multisite_progressor_stable %>% pull(idvi) %>% unique()
multisite_df_to_proc <- multisite_df %>%
    filter(
        idvi %in% idvi_to_process,     # only scans that need to be processed
        !is.na(filepath.raw_amyloid),  # remove subjects with missing PET
        !is.na(filepath.raw_t1),       # remove subjects with missing T1
        (age_diff.raw_t1 <= 1 | site_in_pac)           # matching T1 within 1 year of PET (only for non-PAC sites)
    )

# get MUSE processing table
muse_proc_table <- multisite_df_to_proc %>%
    select(site.subj.age.raw_t1, filepath.raw_t1) %>%
    unique()

# get PET-MRI processing table
muse_dir <- "/ceph/chpc/shared/aristeidis_sotiras_group/b.y.yang_scratch/T1_MUSE"
pet_proc_table <- multisite_df_to_proc %>%
    select(site.subj.age, idvi, site, subj, tracer, age.round, site.subj.age.raw_t1, age.raw_t1.round, filepath.raw_amyloid, filepath.raw_amyloid.json) %>%
    mutate(
        musemaskpath = file.path(muse_dir, "Protocols/Skull-Stripped", site.subj.age.raw_t1, str_glue("{site.subj.age.raw_t1}_T1_LPS_N4_brainmask_muse-ss.nii.gz")),
        musemriskullpath = file.path(muse_dir, "Protocols/Skull-Stripped", site.subj.age.raw_t1, str_glue("{site.subj.age.raw_t1}_T1_LPS_N4_brain_muse-ss.nii.gz")),
        muselabelpath = file.path(muse_dir, "Protocols/MUSE", site.subj.age.raw_t1, str_glue("{site.subj.age.raw_t1}_T1_LPS_N4_brain_muse-ss_fastbc_muse.nii.gz"))
    )

# save outputs
write_delim(muse_proc_table, file.path(opt$odir, "muse.txt"), col_names = FALSE, delim = " ")
write_csv(pet_proc_table, file.path(opt$odir, "petanalysis.csv"))

if (FALSE) {
    # TEMP: process PAC images separately
    muse_proc_table_pac <- multisite_df_to_proc %>%
        filter(site_in_pac) %>%
        select(site.subj.age.raw_t1, filepath.raw_t1) %>%
        unique()
    write_delim(muse_proc_table_pac, file.path(opt$odir, "muse_pac.txt"), col_names = FALSE, delim = " ")
}
