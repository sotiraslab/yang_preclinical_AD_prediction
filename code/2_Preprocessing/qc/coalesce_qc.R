library(tidyverse)

setwd("/home/b.y.yang/sotiraslab/BradenADLongitudinalPrediction/code/2_Preprocessing/qc")

pet_analysis_dir <- "/ceph/chpc/shared/aristeidis_sotiras_group/b.y.yang_scratch/pet_analysis"

# get complete QC csv
qc_failed_df <- read_csv("qc_failed.csv", col_types=c("cnc"))
all_scans <- list.files(file.path(pet_analysis_dir, "nipype")) %>%
    as_tibble() %>%
    rename(id = value)
qc_df <- all_scans %>%
    left_join(qc_failed_df, by = "id") %>%
    replace_na(list(qc = 1))
qc_df %>% write_csv("qc.csv")
