# ==============================================================================

# Name: flowchart_and_demographics.R
# Author: Braden Yang
# Created: 04/25/2024
# Description: Create inclusion/exclusion flowchart and demographics table

# ==============================================================================

rm(list = ls())

# ===========================
# ===== IMPORT PACKAGES =====
# ===========================

library(tidyverse)
library(optparse)
library(flowchart)
library(gt)
library(gtsummary)
library(rstatix)

# ===========================
# ===== PARSE ARGUMENTS =====
# ===========================

option_list <- list(
    make_option(c("-w", "--wdir"), action="store", default=NULL,
        type="character", help="Path to project directory (if none, uses current working directory)")
)

opt <- parse_args(OptionParser(option_list = option_list))

# set working directory
if (!is.null(opt$wdir)) {
    setwd(opt$wdir)
}

# ============================
# ===== SOURCE FUNCTIONS =====
# ============================

source("code/functions.R")

# =====================
# ===== LOAD DATA =====
# =====================

# data for flowchart
flowchart_df <- read_rds("figures/archive/flowchart/for_flowchart.RDS")

# data for demographics table
features_cross <- read_csv("data/features/features_cross_sectional.csv")

# get list of scans to include in flowcharts & demographic table
progressor_idvi <- map(1:5, ~filter_progressor(flowchart_df, "amyloid_positive", .x) %>% pull(idvi) %>% unique) %>% unlist %>% unique
stable_idvi <- filter_stable(flowchart_df, "amyloid_positive", 5) %>% pull(idvi) %>% unique

# ==============================
# ===== DEMOGRAPHICS TABLE =====
# ==============================

demo_df <-
    bind_rows(features_cross %>% mutate(site = "pooled"), features_cross) %>%
    filter(idvi %in% c(progressor_idvi, stable_idvi)) %>%
    mutate(
        .progressor.fct = factor(.progressor, levels = c("stable", "progressor")),
        .time_stable = if_else(.progressor == "progressor", NA, .time_stable),
        sex = case_match(sex, "M"~"male", "F"~"female", .default=sex),
        sex.F.bool = sex == "female",
        apoe.bool = apoe == 1,
        n_subj = 1,
        time_to_progression_category = case_when(
            .time_to_progression <= 1 ~ "0-1 year",
            .time_to_progression <= 2 ~ "1-2 years",
            .time_to_progression <= 3 ~ "2-3 years",
            .time_to_progression <= 4 ~ "3-4 years",
            .time_to_progression <= 5 ~ "4-5 years",
            .default = NA
        )
    )

# run statistical tests
age_ttest <- demo_df %>%
    group_by(site) %>%
    t_test(age ~ .progressor.fct) %>%
    add_significance("p") %>%
    select(site, .y., group1, group2, p)
sex_fisher <- demo_df %>%
    group_by(site) %>%
    reframe(
        fisher_test(table(pick(sex.F.bool, .progressor.fct)))
    ) %>%
    mutate(.y. = "sex.F.bool")
apoe_fisher <- demo_df %>%
    group_by(site) %>%
    reframe(
        fisher_test(table(pick(apoe.bool, .progressor.fct)))
    ) %>%
    mutate(.y. = "apoe.bool")
test_results_all <- bind_rows(age_ttest, sex_fisher, apoe_fisher) %>%
    mutate(p.adj = p.adjust(p, method = "bonferroni")) %>%
    add_significance("p.adj")
test_results_all %>% write_csv("tables/demographics_tests.csv")

demo_tbl <- demo_df %>%
    tbl_strata(
        strata = site,
        .tbl_fun = ~ tbl_summary(
            data = .x,
            by = .progressor.fct,
            include = c(n_subj, age, sex.F.bool, apoe.bool, tracer, .time_to_progression, .time_stable, time_to_progression_category),
            label = list(
                n_subj = "number of subjects",
                age = "age, mean (SD)",
                sex.F.bool = "female, count (%)",
                apoe.bool = "APOE4 carriership, count (%)",
                tracer = "tracer, count (%)",
                .time_to_progression = "years to MCI/AD progression, mean (SD)",
                .time_stable = "years of CN-stable, mean (SD)",
                time_to_progression_category = "years to MCI/AD progression, count (%)"
            ),
            statistic = list(
                n_subj ~ "{N}",
                all_continuous() ~ "{mean} ({sd})",
                c(sex.F.bool, apoe.bool, tracer, time_to_progression_category) ~ "{n} ({p})"
            ),
            type = list(c(where(is.numeric),-n_subj) ~ "continuous"),
            digits = list(c(where(is.numeric),-n_subj) ~ 2)
        ) %>%
            remove_row_type(variables=c(.time_to_progression, .time_stable, time_to_progression_category), type="missing") %>%
            modify_header(all_stat_cols() ~ "**{level}**", label = "") %>%
            modify_footnote(everything() ~ NA_character_)
    ) %>%
    gtsummary::modify_table_body(  # replace NA entires with "-" (https://stackoverflow.com/questions/66351409/gtsummary-table-replace-empty-cell-information-with)
        ~.x %>% mutate(
            across(all_stat_cols(), ~gsub("^NA.*", "-", .)),
            across(all_stat_cols(), ~gsub(".*(NA).*", "-", .))
        )
    ) %>%
    as_gt()

# save tables
odir_tbl <- "tables"
if (!dir.exists(odir_tbl)) {dir.create(odir_tbl)}
gtsave(demo_tbl, "demographics.html", odir_tbl)
gtsave(demo_tbl, "demographics.rtf", odir_tbl)

# ==========================
# ===== SUBJECT COUNTS =====
# ==========================

# # get subject counts at each time window
# progressor_count <- map(1:5, ~filter_progressor(features_cross, "amyloid_positive", .x) %>% count(tracer) %>% rename(!!paste0("progressor_", .x, "yr") := n)) %>%
#     reduce(full_join, by = "tracer") %>%
#     mutate(across(starts_with("progressor_"), ~replace_na(.x, 0)))
# stable_count <- filter_stable(features_cross, "amyloid_positive", 5) %>% count(tracer) %>%
#     rename(stable_5yr = n)
# subj_count <- full_join(stable_count, progressor_count, by = "tracer") %>%
#     bind_rows(
#         summarise(
#             .,
#             across(where(is.numeric), sum),
#             across(where(is.character), ~'Total')
#         )
#     )

# ==================================
# ===== INFO FOR MANUAL PRISMA =====
# ==================================

# baseline_scans <- flowchart_df %>%
#     mutate(subj_in_study = idvi %in% features_cross$idvi) %>%
#     filter(tracer != "FBB")
# baseline_scans1 <- baseline_scans %>%
#     filter(subj_in_study)
# baseline_scans2 <- baseline_scans %>%
#     filter(!(subj %in% features_cross$subj)) %>%
#     group_by(site, subj) %>%
#     slice_min(age) %>%
#     ungroup() %>%
#     filter(!(site %in% c("BIOCARD", "WRAP")))
# baseline_scans <- bind_rows(baseline_scans1, baseline_scans2)

# reverter <- reverter %>% mutate(subj.site = str_c(site, subj, sep = "."))

# baseline_scans %>%
#     mutate(subj.site = str_c(site, subj, sep = ".")) %>%
#     filter(
#         amyloid_positive == TRUE,
#         .is_valid_baseline,
#         ((!is.na(filepath.raw_t1) & age_diff.raw_t1 <= 1) | (site_in_pac)) &
#         !is.na(age) & !is.na(sex) & !is.na(apoe),
#         subj.site %in% reverter$subj.site,
#         # qc > 0 & .proc_pass_amyloid & .proc_pass_t1,
#     ) %>%
#     count(site)

# baseline_scans %>%
#     filter(.progressor == "progressor") %>%
#     select(.time_stable, .time_to_progression)
#     # mutate(.revert = .time_stable > .time_to_progression) %>%
#     # count(.revert)

# baseline_scans %>% 
#     filter(is.na(.progressor))

# # number of progressors that revert
# foo %>%
#     filter(.progressor == "progressor") %>%
#     mutate(.revert = .time_stable > .time_to_progression) %>%
#     count(.revert)

# fc2 <- baseline_scans %>%
#     as_fc(label = "baseline amyloid PET scans") %>%
#     fc_split(site) %>%
#     fc_filter(
#         amyloid_positive == TRUE,
#         label="amyloid-positive",
#         show_exc=TRUE,
#         label_exc="not amyloid-positive"
#     ) %>%
#     fc_filter(
#         .is_valid_baseline,
#         label="baseline CDR = 0",
#         show_exc=TRUE,
#         label_exc="baseline CDR > 0 or is missing"
#     ) %>%
#     fc_filter(
#         ((!is.na(filepath.raw_t1) & age_diff.raw_t1 <= 1) | (site_in_pac)) &
#         !is.na(age) & !is.na(sex) & !is.na(apoe),
#         label="matching T1 & demographics",
#         show_exc=TRUE,
#         label_exc="missing T1, demographics and APOE4"
#     ) %>%
#     fc_filter(
#         idvi %in% c(progressor_idvi, stable_idvi),
#         label="stable or progressor",
#         show_exc=TRUE,
#         label_exc="baseline CDR > 0, didn't progress within 5 years"
#     ) %>%
#     fc_filter(
#         qc > 0 & .proc_pass_amyloid & .proc_pass_t1,
#         label="PET-MRI processing passed",
#         show_exc=TRUE,
#         label_exc="PET-MRI processing failed"
#     ) %>%
#     fc_split(.progressor) %>%
#     fc_split(tracer) %>%
#     fc_draw()

# # save figure
# odir_fc <- file.path("figures/flowchart")
# if (!dir.exists(odir_fc)) {dir.create(odir_fc)}
# fc_export(fc2, file.path(odir_fc, "prisma2.png"), width = 20, height = 10, res = 300, units = "in")
