rm(list = ls())

library(tidyverse)

# Limitations of using Centiloid:
# - since we are using the conversion equations provided by each dataset and not calibrating
#   our own, there will be slight differences between the SUVRs of the dataset and SUVRs obtained
#   from our pipeline due to differences in the pipeline operations (e.g. smoothing, acquisition
#   window, reference region, etc.)

# ==============
# ===== A4 =====
# ==============

# FBP equation: 183.07 * SUVR - 177.26

a4_csv_dir <- "/ceph/chpc/shared/aristeidis_sotiras_group/aris_data/A4/rawdata/Clinical"
a4_data_dict <- read_csv(file.path(a4_csv_dir, "Documents/Data Dictionaries/derived_datadic.csv"))
a4_data_dict %>%
    filter(
        FILE_NAME == "ADQS.csv",
        FIELD_NAME == "AMYLCENT"
    ) %>%
    pull(FIELD_DESC)

# ================
# ===== ADNI =====
# ================

# FBP equation: 188.22 * SUVR - 189.16

# reference region: whole cerebellum

# see UC Berkeley Amyloid PET Processing Methods document (June 20 2025 version) pg. 5 for equations

# ================
# ===== HABS =====
# ================

# Centiloid equations not readily available for HABS

# ================
# ===== MCSA =====
# ================

# PIB equation: 88.72 * SUVR - 109.75

mcsa_df <- read_csv("/ceph/chpc/shared/aristeidis_sotiras_group/aris_data/MCSA/tabular/MCSA_Data_04Apr2025.csv")

df <- mcsa_df %>%
    select(SPM12_PIB_RATIO, SPM12_PIB_CENTILOID) %>%
    drop_na()
suvr_to_cl <- lm(SPM12_PIB_CENTILOID ~ SPM12_PIB_RATIO, data = df)
suvr_to_cl$coefficients

# =================
# ===== OASIS =====
# =================

# FBP equation: 163.6 * SUVR - 181.0
# PIB equation: 111.8 * SUVR - 119.3

# FBP acquisition window: 50-70 minutes
# PIB acquisition window: 30-60 minutes

# reference region: cerebellar cortex

# see OASIS-3 data dictionary pg. 27 (v2.3) for the equations

# ===============
# ===== PAC =====
# ===============

# Centiloid equations not readily available for PAC datasets (AIBL, BLSA)

# ================================
# ===== CONVERT TO CENTILOID =====
# ================================

features_df <- read_csv("data/features/features_cross_sectional.csv")
features_cl <- features_df %>%
    filter(site %in% c("A4", "ADNI", "MCSA", "OASIS")) %>%
    mutate(across(
        ends_with(".amyloid"),
        ~ case_match(site,
            "A4" ~ 183.07 * .x - 177.26,
            "ADNI" ~ 188.22 * .x - 189.16,
            "MCSA" ~ 88.72 * .x - 109.75,
            "OASIS" ~ if_else(
                tracer == "FBP",
                163.6 * .x - 181.0,
                111.8 * .x - 119.3
            )
        )
    ))
write_csv(features_cl, "data/features/features_cross_sectional_CL.csv")
