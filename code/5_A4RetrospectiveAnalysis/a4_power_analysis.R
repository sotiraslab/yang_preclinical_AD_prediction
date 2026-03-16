# ==============================================================================

# Name: a4_power_analysis.R
# Author: Braden Yang
# Created: 07/07/2025
# Description: Perform power analysis to detect significant differences in FBP
#   PET SUVR using ANCOVA with or without stratification by ML model predictions

# ==============================================================================

rm(list = ls())

INTERACTIVE <- TRUE

library(tidyverse)
library(optparse)
library(devtools)
library(emmeans)

# ===========================
# ===== PARSE ARGUMENTS =====
# ===========================

option_list <- list(
    make_option(c("-w", "--wdir"), action="store", default=NULL,
        type="character", help="Path to project directory (if none, uses current working directory)"),
    make_option(c("-o", "--odir_fig"), action="store", default="figures",
        type="character", help="Path to output directory for figures"),
    make_option(c("-t", "--odir_table"), action="store", default="tables",
        type="character", help="Path to output directory for tables")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (INTERACTIVE) {
    opt$wdir <- "/home/b.y.yang/sotiraslab/BradenADLongitudinalPrediction"
    opt$odir_fig <- "/home/b.y.yang/sotiraslab/BradenADLongitudinalPrediction/figures"
    opt$odir_table <- "/home/b.y.yang/sotiraslab/BradenADLongitudinalPrediction/tables"
}

# set working directory
if (!is.null(opt$wdir)) {
    setwd(opt$wdir)
}

# load byyfunctions
load_all("/home/b.y.yang/sotiraslab/BradenADLongitudinalPrediction/submodules/byyfunctions")

# make output directories
byyfunctions::make_dir(opt$odir_fig)
byyfunctions::make_dir(opt$odir_table)

# ============================
# ===== DEFINE FUNCTIONS =====
# ============================

# procedure:
# - set a fixed sample size per group
# - select a random sample from the cohort, either unenriched or enriched
# - compute the proportion of significant differences out of the total number of random samples
# - repeat for a list of sample sizes

source("/home/b.y.yang/sotiraslab/BradenADLongitudinalPrediction/code/a4_functions.R")

stratify_df <- function(.data, t, stratifications, class = "progressor") {
    .data %>%
        filter(BID %in% (stratifications %>% filter(time_window == t, predicted_class == class) %>% pull(subj)))
}

fit_amyloid_ancova <- function(.data, time_window = NA) {

    # contrasts in emmeans: https://cran.r-project.org/web/packages/emmeans/vignettes/comparisons.html

    n_df <- .data %>% count(TX)
    n_total = nrow(.data)
    model <- aov(
        Composite_Summary.diff ~ age.start + apoe + Composite_Summary.start + TX,
        data = .data
    )
    emm <- emmeans(model, ~ TX)
    emm.contrast <- pairs(emm, reverse = TRUE)
    emm.eff_size <- eff_size(emm.contrast, sigma = sigma(model), edf = df.residual(model), method = "identity")

    emm.df <- emm %>%
        as_tibble() %>%
        mutate(stratification = if_else(is.na(time_window), "orig", as.character(time_window))) %>%
        left_join(n_df, by = "TX")
    emm.contrast.df <- emm.contrast %>%
        as_tibble() %>%
        mutate(
            stratification = if_else(is.na(time_window), "orig", as.character(time_window)),
            n = n_total,
            lower.CL = estimate - 1.96*SE,
            upper.CL = estimate + 1.96*SE
        ) %>%
        bind_cols(emm.eff_size %>% as_tibble() %>% select(effect.size, SE) %>% rename(SE.effect.size = SE))
    
    return(list(emm.df, emm.contrast.df))

}

fit_amyloid_ancova_plus_stratification <- function(.data, time_window = NA, class = "progressor") {
    if (!is.na(time_window)) {
        .data <- .data %>% stratify_df(time_window, a4_stratifications, class)
    }
    fit_amyloid_ancova(.data, time_window)
}

random_sample_df <- function(.data, n_sample) {
    .data %>%
        group_by(TX) %>%
        slice_sample(n = n_sample) %>%
        ungroup()
}

single_bootstrap_run <- function(.data, n_sample) {
    amyloid_ancova_results <- .data %>%
        random_sample_df(n_sample) %>%
        fit_amyloid_ancova()
    return(amyloid_ancova_results[[2]]$p.value)
}

compute_power <- function(.data, n_sample, time_window = NA, n_iter = 1000, alpha = 0.05) {

    print("++ working on sample size = {n_sample} ++" %>% str_glue())

    if (!is.na(time_window)) {
        .data <- .data %>% stratify_df(time_window, a4_stratifications)
    }
    pval <- map(1:n_iter, ~single_bootstrap_run(.data, n_sample))
    sig <- unlist(pval) < alpha
    sum(sig) / n_iter

}

# ====================================
# ===== LOAD AND PREPROCESS DATA =====
# ====================================

a4_predict <- read_csv(file.path(opt$wdir, "features_A4_treatment_with_prediction_ANN.csv")) %>%
    mutate(age.round = round(age, 2) %>% format(nsmall = 2, trim = TRUE) %>% as.numeric) %>%
    filter(!is.na(TX))  # remove missing TX rows

# get stratifications by ML model
a4_stratifications <- a4_predict %>%
    select(subj, starts_with("y_score.svm")) %>%
    mutate(across(
        starts_with("y_score"),
        ~if_else(.x > 0, "progressor", "stable")
    )) %>%
    rename_with(~str_replace(.x, "y_score.svm", "predicted_class"), starts_with("y_score")) %>%

    pivot_longer(
        cols = -subj,
        names_pattern = "predicted_class\\.(\\d)",
        names_to = "time_window",
        values_to = "predicted_class"
    )

subj_with_stratifications <- a4_stratifications %>% pull(subj) %>% unique()

amyloid_pet_df <- read_csv("/ceph/chpc/shared/aristeidis_sotiras_group/aris_data/A4/rawdata/Clinical/External Data/imaging_SUVR_amyloid.csv")

amyloid_composite <- amyloid_pet_df %>%
    filter(BID %in% subj_with_stratifications) %>%
    a4_get_age() %>%
    a4_merge_subjinfo("TX") %>%
    filter(!is.na(TX)) %>%
    left_join(
        subjinfo %>% select(BID, APOEGN),
        by = "BID"
    ) %>%
    mutate(apoe = case_when(APOEGN %in% c("E4/E4","E3/E4") ~ 1, is.na(APOEGN) ~ NA, .default = 0) %>% as.factor()) %>%
    filter(brain_region == "Composite_Summary") %>%
    pivot_wider(
        id_cols = c(SUBSTUDY, TX, BID, VISCODE, age, apoe),
        names_from = brain_region,
        values_from = suvr_cer
    )

amyloid_df_for_ancova <- amyloid_composite %>%
    filter(VISCODE %in% c("002", "066")) %>%
    mutate(trial_period = case_match(
        VISCODE,
        "002" ~ "start",
        "066" ~ "end"
    )) %>%
    select(-VISCODE) %>%
    pivot_wider(
        names_from = trial_period,
        names_sep = ".",
        values_from = c(age, Composite_Summary)
    ) %>%
    filter(if_all(c(age.start, starts_with("Composite_Summary")), ~!is.na(.x))) %>%
    mutate(Composite_Summary.diff = Composite_Summary.end - Composite_Summary.start)

# ========================================
# ===== BOOTSTRAPPING POWER ANALYSIS =====
# ========================================

# define variables
alpha_for_significance <- 0.05
n_iter <- 1000

# set random seed
set.seed(42)

get_power_df <- function(n_sample_list, time_window, label) {
    power <- map(
        n_sample_list,
        ~compute_power(
            amyloid_df_for_ancova,
            .x,
            time_window,
            n_iter,
            alpha_for_significance
        )
    )
    tibble(
        sample_size = n_sample_list,
        power = unlist(power),
        stratification = label
    )
}

# unenriched power
power_unenriched <- get_power_df(seq(10, 80, 10), NA, "original cohort")

# enriched power
power_enriched_1 <- get_power_df(seq(10, 30, 10), 1, "1-year progressor")
power_enriched_2 <- get_power_df(seq(10, 70, 10), 2, "2-year progressor")
power_enriched_3 <- get_power_df(seq(10, 80, 10), 3, "3-year progressor")
power_enriched_4 <- get_power_df(seq(10, 80, 10), 4, "4-year progressor")
power_enriched_5 <- get_power_df(seq(10, 80, 10), 5, "5-year progressor")

# combine df
power_df <- bind_rows(
    power_unenriched,
    power_enriched_1,
    power_enriched_2,
    power_enriched_3,
    power_enriched_4,
    power_enriched_5
)

# # load previous power analysis results
# power_df <- read_csv("tables/amyloid_power_analysis.csv")

# plot
colors <- RColorBrewer::brewer.pal(5, "Set1")
colors_list <- c(
    "original cohort" = "black",
    "1-year progressor" = colors[1],
    "2-year progressor" = colors[2],
    "3-year progressor" = colors[3],
    "4-year progressor" = colors[4],
    "5-year progressor" = colors[5]
)
shape_list <- c(
    "original cohort" = 16,
    "1-year progressor" = 17,
    "2-year progressor" = 15,
    "3-year progressor" = 18,
    "4-year progressor" = 3,
    "5-year progressor" = 4
)
power_plot <- power_df %>%
    mutate(stratification = factor(stratification, levels = names(colors_list))) %>%
    ggplot(aes(x = sample_size, y = power, color = stratification)) +
        geom_line() +
        geom_point(aes(shape = stratification)) +
        labs(
            x = "Sample size per trial group",
            y = "Power"
        ) +
        scale_x_continuous(breaks = seq(10,80,10)) +
        scale_y_continuous(labels = scales::percent) +
        scale_color_manual(values = colors_list, name = "Cohort") +
        scale_shape_manual(values = shape_list, name = "Cohort") +
        theme_minimal() +
        theme(
            legend.position = "inside",
            legend.position.inside = c(0.85, 0.2),
            legend.key.size = unit(0.75, "lines"),
            panel.background = element_rect(fill = "white", color = NA),
            plot.background  = element_rect(fill = "white", color = NA),
            legend.background = element_rect(fill = alpha("grey", 0.15), color = NA)
        )

power_diff_plot <- power_df %>%
    pivot_wider(names_from = stratification, values_from = power) %>%
    mutate(across(ends_with("progressor"), ~ .x - `original cohort`)) %>%
    pivot_longer(cols = ends_with("progressor"), names_to = "stratification", values_to = "power_diff") %>%
    mutate(stratification = factor(stratification, levels = names(colors_list))) %>%
    filter(!is.na(power_diff)) %>%
    ggplot(aes(x = sample_size, y = power_diff, color = stratification)) +
        geom_abline(slope = 0, intercept = 0, color = "black") +
        geom_line() +
        geom_point(aes(shape = stratification)) +
        labs(
            x = "Sample size per trial group",
            y = "Difference in power"
        ) +
        scale_x_continuous(breaks = seq(10,80,10)) +
        scale_y_continuous(labels = scales::percent) +
        scale_color_manual(values = colors_list, name = "Cohort") +
        scale_shape_manual(values = shape_list, name = "Cohort") +
        theme_minimal() +
        theme(
            legend.position = "inside",
            legend.position.inside = c(0.15, 0.8),
            legend.key.size = unit(0.75, "lines"),
            panel.background = element_rect(fill = "white", color = NA),
            plot.background  = element_rect(fill = "white", color = NA),
            legend.background = element_rect(fill = alpha("grey", 0.15), color = NA)
        )

# save results
ggsave(file.path(opt$odir_fig, "amyloid_power_analysis.png"), power_plot, width = 7, height = 4, units = "in", dpi = 300)
ggsave(file.path(opt$odir_fig, "amyloid_power_analysis_diff.png"), power_diff_plot, width = 7, height = 4, units = "in", dpi = 300)
write_csv(power_df, file.path(opt$odir_table, "amyloid_power_analysis.csv"))
