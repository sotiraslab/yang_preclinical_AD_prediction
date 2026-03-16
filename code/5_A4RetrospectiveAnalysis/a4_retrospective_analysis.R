# ==============================================================================

# Name: a4_retrospective_analysis.R
# Author: Braden Yang
# Created: 07/07/2025
# Description: Perform a retrospective analysis on the A4 dataset using predictions
#   from ML classifiers to stratify participants, retroactively select only
#   progressors, and rerun statistical analyses performed by the original A4 study
#   paper

# ==============================================================================

rm(list = ls())

INTERACTIVE <- TRUE

library(tidyverse)
library(optparse)
library(gt)
library(gtsummary)
library(devtools)
library(splines)
library(nlme)
library(emmeans)
library(clubSandwich)
library(lme4)
library(broom.mixed)
library(knitr)
library(kableExtra)
library(patchwork)
library(rstatix)

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

source("/home/b.y.yang/sotiraslab/BradenADLongitudinalPrediction/code/a4_functions.R")

formatp <- function(x) case_when(
  x < 0.001 ~ "p<0.001",
  x > 0.01 ~ Hmisc::format.pval(x, digits=2, eps=0.01, nsmall=2),
  TRUE ~ Hmisc::format.pval(x, digits=3, eps=0.001, nsmall=3))

stratify_df <- function(.data, t, stratifications, class = "progressor") {
    .data %>%
        filter(BID %in% (stratifications %>% filter(time_window == t, predicted_class == class) %>% pull(subj)))
}

sig_level <- function(p) {
    ifelse(
        p <= 0.0001, '****', ifelse(
        p <= 0.001, '***', ifelse(
        p <= 0.01, '**', ifelse(
        p <= 0.05, '*', ''))))
}

# =====================
# ===== LOAD DATA =====
# =====================

# data dictionary
derived_data_dict <- read_csv("/ceph/chpc/shared/aristeidis_sotiras_group/aris_data/A4/rawdata/Clinical/Documents/Data Dictionaries/derived_datadic.csv")
external_data_dict <- read_csv("/ceph/chpc/shared/aristeidis_sotiras_group/aris_data/A4/rawdata/Clinical/Documents/Data Dictionaries/external_datadic.csv")
raw_data_dict <- read_csv("/ceph/chpc/shared/aristeidis_sotiras_group/aris_data/A4/rawdata/Clinical/Documents/Data Dictionaries/clinical_datadic.csv")

# load stable/progressor predictions by ML models
a4_predict <- read_csv("data/a4_retospective_analysis/features_A4_treatment_with_prediction.csv") %>%
    mutate(age.round = round(age, 2) %>% format(nsmall = 2, trim = TRUE) %>% as.numeric) %>%
    filter(!is.na(TX))  # remove missing TX rows

# subject list
subj_list <- a4_predict %>% pull(subj) %>% unique

# load PACC observations and preprocess (code is obtained from `Intro-to-A4-data.Rmd` provided by the A4 dataset)
ADQS_raw <- read_csv("/ceph/chpc/shared/aristeidis_sotiras_group/aris_data/A4/rawdata/Clinical/Derived Data/ADQS.csv")
ADQS_PACC <- ADQS_raw %>%
    filter(MITTFL== 1) %>%
    filter(EPOCH == "BLINDED TREATMENT") %>%
    filter(QSTESTCD == "PACC") %>%
    rename(PACC = QSSTRESN) %>%
    select(BID, ASEQNCS, TX, ADURW, TX, AGEYR, 
    AAPOEGNPRSNFLG, EDCCNTU, SUVRCER, QSVERSION, PACC) %>%
    mutate(TX = factor(TX, levels = c("Placebo", "Solanezumab"))) %>%
    na.omit() %>%
    filter(BID %in% subj_list)  # filter to only include subjects that were processed by PET-MRI pipeline

# imaging biomarkers
amyloid_pet_df <- read_csv("/ceph/chpc/shared/aristeidis_sotiras_group/aris_data/A4/rawdata/Clinical/External Data/imaging_SUVR_amyloid.csv")

# ===========================
# ===== PREPROCESS DATA =====
# ===========================

# compute PACC rate-of-change using linear regression
pacc_slope <- ADQS_PACC %>%
    group_by(BID) %>%
    summarise(
        pacc.bl = PACC[ADURW == min(ADURW)],
        pacc.slope = lm(PACC ~ ADURW)$coefficients[2] * 52,  # convert 1/weeks to 1/yrs
        pacc.n = n()
    )

# merge baseline PACC and change in PACC
a4_predict_with_pacc <- a4_predict %>%
    mutate(num_subj = 1) %>%
    left_join(pacc_slope, by = c("subj" = "BID"))

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

a4_strata_list <- map(1:5, ~stratify_df(ADQS_PACC, .x, a4_stratifications))
a4_strata_df <- tibble(
    stratify = c("orig", str_c("stratify_", 1:5, "yr_progressor")),
    data = c(list(ADQS_PACC), a4_strata_list)
)

# preprocess amyloid PET
subj_with_stratifications <- a4_stratifications %>% pull(subj) %>% unique()

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

# ==================================================
# ===== FIT NATURAL CUBIC SPLINE MODEL ON PACC =====
# ==================================================

# fit a natural cubic spline model to longitudinal PACC scores during the placebo-controlled
# period of the clinical trial, then estimate the model-adjusted mean and 95% CI of PACC
# at the end of the period
# NOTE: we're only using the 5-year progressor stratification here, since the placebo-controlled period
#   of the A4 study was 240 weeks, or about 5 years

a4_strata_df_filter <- a4_strata_df %>%
    mutate(stratify = factor(
        stratify,
        levels = c(
            "orig",
            "stratify_1yr_progressor",
            "stratify_2yr_progressor",
            "stratify_3yr_progressor",
            "stratify_4yr_progressor",
            "stratify_5yr_progressor"
        ),
        labels = c(
            "original cohort",
            "1-year progressors",
            "2-year progressors",
            "3-year progressors",
            "4-year progressors",
            "5-year progressors"
        )
    ))
emm_df_list <- list()
for (i in 1:6) {

    foo <- a4_strata_df_filter %>% slice(i)
    s <- foo %>% pull(stratify)
    df <- foo %>% pull(data) %>% .[[1]]

    print("working on {s}" %>% str_glue)

    ns21 <- function(t) {
        as.numeric(predict(splines::ns(df$ADURW, df=2,
            Boundary.knots = c(0, max(df$ADURW))), t)[,1])
    }
    ns22 <- function(t) {
        as.numeric(predict(splines::ns(df$ADURW, df=2,
            Boundary.knots = c(0, max(df$ADURW))), t)[,2])
    }

    # fit GLS model
    # - following the implementation given in the `Intro-to-A4-data.Rmd` document
    model <- nlme::gls(PACC ~ 
        I(ns21(ADURW)) + I(ns22(ADURW)) +
        (I(ns21(ADURW)) + I(ns22(ADURW))):TX + 
        AGEYR + AAPOEGNPRSNFLG + EDCCNTU + SUVRCER + QSVERSION,
        data = df,
        weights = varIdent(form = ~ 1 | ASEQNCS),
        correlation = corARMA(form = ~ ASEQNCS | BID, p = 10)
    )

    # estimate PACC scores at 240 weeks and contrast between Solanezumab and Placebo using `emmeans` package
    # - following the implementation given in the `Intro-to-A4-data.Rmd` document
    n_df <- df %>% select(BID, TX) %>% unique() %>% count(TX)
    emm <- ref_grid(
        model, 
        at = list(ADURW = 240, TX = levels(df$TX)),
        vcov. = clubSandwich::vcovCR(model, type = "CR2") %>% as.matrix(), 
        data = df, 
        mode = "satterthwaite"
    ) %>%
        emmeans(specs = "TX", by = "ADURW")
    emm.contrast <- emm %>% pairs(reverse = TRUE, adjust = "none")

    # estimate effect size
    model.df = model$dims$N - length(model$coefficients)  # df = num observations minus num of fixed effects
    model.sigma = sigma(model)  # SD of model residuals
    emm.eff_size <- eff_size(emm.contrast, sigma = model.sigma, edf = model.df, method = "identity")

    emm.df <- emm %>%
        as_tibble() %>%
        left_join(n_df, by = "TX")
    emm.contrast.df <- emm.contrast %>%
        as_tibble() %>%
        mutate(
            n = n_df %>% pull(n) %>% sum(),
            lower.CL = estimate - 1.96*SE,
            upper.CL = estimate + 1.96*SE
        ) %>%
        bind_cols(emm.eff_size %>% as_tibble() %>% select(effect.size, SE) %>% rename(SE.effect.size = SE))

    df_for_gt <- emm.df %>%
        select(-c(df, ADURW)) %>%
        rename(estimate = emmean) %>%
        bind_rows(emm.contrast.df %>%
                    select(-c(df, t.ratio, ADURW)) %>%
                    rename(TX = contrast)) %>%
        mutate(
            p.value.sig = sig_level(p.value),
            p.value.str = if_else(!is.na(p.value),
                str_c(
                    formatC(p.value, format = "e", digits = 2),
                    "<sup>{p.value.sig}</sup>" %>% str_glue()
                ),
                NA
            ),
            across(c(estimate, SE), ~if_else(is.na(.x), NA, formatC(.x, digits = 3, format = "fg"))),
            estimate.SE = str_glue("{estimate} \u00B1 {SE}")
        ) %>%
        mutate(
            across(p.value.str, ~ map(.x, html)),
            p.value.str = if_else(is.na(p.value), NA, p.value.str)
        )
    emm_df_list[[as.character(s)]] <- df_for_gt

}

# combine tables and save
df_for_gt <- bind_rows(emm_df_list, .id = "stratification") %>%
    mutate(p.value = p.adjust(p.value, method = "bonferroni")) %>%
    add_significance("p.value") %>%
    mutate(
        p.value.str = if_else(is.na(p.value), NA, str_c(p.value, p.value.sig)),
        across(c(effect.size, SE.effect.size), ~if_else(is.na(.x), NA, formatC(.x, digits = 3, format = "fg"))),
        effect.size.SE = if_else(is.na(effect.size), NA, str_glue("{effect.size} \u00B1 {SE.effect.size}"))
    )

pacc_gt <- df_for_gt %>%
    select(stratification, TX, n, estimate.SE, lower.CL, upper.CL, effect.size.SE, p.value.str) %>%
    gt(
        rowname_col = "TX",
        groupname_col = "stratification"
    ) %>%
    cols_label(
        TX = "trial arm",
        estimate.SE = "model-adjusted mean PACC (\u00B1 SE)",
        lower.CL = "lower 95% CI",
        upper.CL = "upper 95% CI",
        effect.size.SE = "effect size (\u00B1 SE)",
        p.value.str = "p-value"
    ) %>%
    fmt_number(
        columns = all_of(c("estimate.SE", "lower.CL", "upper.CL")),
        n_sigfig = 4
    ) %>%
    sub_missing(
        columns = everything(),
        missing_text = ""
    )
gtsave(pacc_gt, file.path(opt$odir_table, "pacc_ncs_tx.html"))
gtsave(pacc_gt, file.path(opt$odir_table, "pacc_ncs_tx.rtf"))

# ===================================================================
# ===== REFIT NCS MODEL ON PACC TO COMPARE STABLE VS PROGRESSOR =====
# ===================================================================

# refit NCS models, but instead compare stables vs. progressors

stable_progressor_emm_df_list <- list()
stable_progressor_ncs_plot_list <- list()
for (time_window_i in seq(1,5)) {
    print("+ working on time window {time_window_i} +" %>% str_glue())
    emm_df_list <- list(); plot_list <- list()
    for (TX_group in ADQS_PACC$TX %>% unique()) {

        print("working on {TX_group}" %>% str_glue)

        ns21 <- function(t) {
            as.numeric(predict(splines::ns(df$ADURW, df=2,
                Boundary.knots = c(0, max(df$ADURW))), t)[,1])
        }
        ns22 <- function(t) {
            as.numeric(predict(splines::ns(df$ADURW, df=2,
                Boundary.knots = c(0, max(df$ADURW))), t)[,2])
        }

        df <- ADQS_PACC %>%
            left_join(
                a4_stratifications %>% filter(time_window == time_window_i),
                by = c("BID" = "subj")
            ) %>%
            mutate(across(predicted_class, as.factor)) %>%
            filter(TX == TX_group)

        # fit GLS model
        # - following the implementation given in the `Intro-to-A4-data.Rmd` document
        model <- nlme::gls(PACC ~ 
            I(ns21(ADURW)) + I(ns22(ADURW)) +
            (I(ns21(ADURW)) + I(ns22(ADURW))):predicted_class + 
            AGEYR + AAPOEGNPRSNFLG + EDCCNTU + SUVRCER + QSVERSION,
            data = df,
            weights = varIdent(form = ~ 1 | ASEQNCS),
            correlation = corARMA(form = ~ ASEQNCS | BID, p = 10)
        )

        # estimate PACC scores at 240 weeks and contrast between Solanezumab and Placebo using `emmeans` package
        # - following the implementation given in the `Intro-to-A4-data.Rmd` document
        n_df <- df %>% select(BID, predicted_class) %>% unique() %>% count(predicted_class)
        emm <- ref_grid(
            model, 
            at = list(ADURW = 240, predicted_class = levels(df$predicted_class)),
            vcov. = clubSandwich::vcovCR(model, type = "CR2") %>% as.matrix(), 
            data = df, 
            mode = "satterthwaite"
        ) %>%
            emmeans(specs = "predicted_class", by = "ADURW")
        emm.contrast <- emm %>% pairs(reverse = TRUE, adjust = "none")

        # estimate effect size
        model.df = model$dims$N - length(model$coefficients)  # df = num observations minus num of fixed effects
        model.sigma = sigma(model)  # SD of model residuals
        emm.eff_size <- eff_size(emm.contrast, sigma = model.sigma, edf = model.df, method = "identity")

        emm.df <- emm %>%
            as_tibble() %>%
            left_join(n_df, by = "predicted_class")
        emm.contrast.df <- emm.contrast %>%
            as_tibble() %>%
            mutate(
                n = n_df %>% pull(n) %>% sum(),
                lower.CL = estimate - 1.96*SE,
                upper.CL = estimate + 1.96*SE
            ) %>%
            bind_cols(emm.eff_size %>% as_tibble() %>% select(effect.size, SE) %>% rename(SE.effect.size = SE))

        # estimate overall emmean
        n_overall <- df %>% select(BID) %>% unique() %>% nrow()
        emm_overall.df <- ref_grid(
            model, 
            at = list(ADURW = 240),
            vcov. = clubSandwich::vcovCR(model, type = "CR2") %>% as.matrix(), 
            data = df, 
            mode = "satterthwaite"
        ) %>%
            emmeans(spec = ~1, by = "ADURW") %>%
            as_tibble() %>%
            mutate(n = n_overall) %>% rename(predicted_class = `1`)

        # store data to plot model-adjusted PACC
        data_for_g <- ref_grid(
            model,
            at = list(ADURW = seq(0, 312, by=12), predicted_class = levels(df$predicted_class)),
            vcov. = clubSandwich::vcovCR(model, type = "CR2") %>% as.matrix(), 
            data = df, 
            mode = "satterthwaite"
        ) %>%
            emmeans(specs = "predicted_class", by = "ADURW") %>%
            as_tibble()
        data_for_g_overall <- ref_grid(
            model,
            at = list(ADURW = seq(0, 312, by=12)),
            vcov. = clubSandwich::vcovCR(model, type = "CR2") %>% as.matrix(), 
            data = df, 
            mode = "satterthwaite"
        ) %>%
            emmeans(specs = ~1, by = "ADURW") %>%
            as_tibble() %>%
            rename(predicted_class = `1`)
        data_for_g <- data_for_g %>%
            bind_rows(data_for_g_overall) %>%
            mutate(
                time_window = time_window_i,
                TX = TX_group
            )
        plot_list[[TX_group]] <- data_for_g        
        
        # create gt table of statistical test results
        df_for_gt <- bind_rows(emm_overall.df, emm.df) %>%
            select(-c(df, ADURW)) %>%
            rename(estimate = emmean) %>%
            bind_rows(emm.contrast.df %>%
                        select(-c(df, t.ratio, ADURW)) %>%
                        rename(predicted_class = contrast)) %>%
            mutate(
                p.value.sig = sig_level(p.value),
                p.value.str = if_else(!is.na(p.value),
                    str_c(
                        formatC(p.value, format = "e", digits = 2),
                        "<sup>{p.value.sig}</sup>" %>% str_glue()
                    ),
                    NA
                ),
                across(c(estimate, SE), ~if_else(is.na(.x), NA, formatC(.x, digits = 3, format = "fg"))),
                estimate.SE = str_glue("{estimate} \u00B1 {SE}")
            ) %>%
            mutate(
                across(p.value.str, ~ map(.x, html)),
                p.value.str = if_else(is.na(p.value), NA, p.value.str),
                time_window = time_window_i,
                TX = TX_group
            )
        emm_df_list[[TX_group]] <- df_for_gt
        
    }

    stable_progressor_emm_df_list[[time_window_i]] <- emm_df_list
    stable_progressor_ncs_plot_list[[time_window_i]] <- plot_list

}

# combine tables and save
df_for_gt <- stable_progressor_emm_df_list %>% purrr::flatten() %>%
    bind_rows() %>%
    mutate(p.value = p.adjust(p.value, method = "bonferroni")) %>%
    add_significance("p.value") %>%
    mutate(
        p.value.signif = if_else(p.value.signif=="ns", "", p.value.signif),
        p.value.str = if_else(!is.na(p.value),
            str_c(
                formatC(p.value, format = "e", digits = 2),
                p.value.signif
            ),
            NA
        ),
        across(c(effect.size, SE.effect.size), ~if_else(is.na(.x), NA, formatC(.x, digits = 3, format = "fg"))),
        effect.size.SE = if_else(is.na(effect.size), NA, str_glue("{effect.size} \u00B1 {SE.effect.size}"))
    )
stable_progressor_pacc_gt <- df_for_gt %>%
    select(TX, time_window, predicted_class, n, estimate.SE, lower.CL, upper.CL, effect.size.SE, p.value.str) %>%
    gt(
        rowname_col = "time_window",
        groupname_col = "TX"
    ) %>%
    tab_stubhead(label = "progression time-horizon (yrs)") %>%
    cols_label(
        predicted_class = "predicted class",
        TX = "trial arm",
        estimate.SE = "model-adjusted mean PACC (\u00B1 SE)",
        lower.CL = "lower 95% CI",
        upper.CL = "upper 95% CI",
        effect.size.SE = "effect size (\u00B1 SE)",
        p.value.str = "p-value"
    ) %>%
    fmt_number(
        columns = all_of(c("estimate.SE", "lower.CL", "upper.CL")),
        n_sigfig = 4
    ) %>%
    sub_missing(
        columns = everything(),
        missing_text = ""
    )
gtsave(stable_progressor_pacc_gt, file.path(opt$odir_table, "pacc_ncs_stable_progressor_overall.html"))
gtsave(stable_progressor_pacc_gt, file.path(opt$odir_table, "pacc_ncs_stable_progressor_overall.rtf"))

# combine plots and save
text_df <- df_for_gt %>%
    filter(predicted_class == "stable - progressor") %>%
    mutate(
        p.value.round = case_when(
            p.value < 0.001 ~ "< 0.001",
            p.value < 0.01 ~ str_c("= ", sprintf("%.3f", p.value)),
            TRUE ~ str_c("= ", sprintf("%.2f", p.value))
        ),
        p.value.signif = if_else(p.value.signif=="****", "***", p.value.signif),
        p.value.str = str_c(p.value.round, p.value.signif)
    ) %>%
    select(TX, time_window, effect.size.SE, p.value.str) %>%
    mutate(
        time_window = str_glue("{time_window}-year progressors"),
        text = str_glue("effect size = {effect.size.SE}\np {p.value.str}")
    )
df_for_plots <- stable_progressor_ncs_plot_list %>% purrr::flatten() %>% bind_rows()
colors <- RColorBrewer::brewer.pal(3, "Set1")
colors_list <- c(
    "overall" = "gray30",
    "progressor" = colors[1],
    "stable" = colors[2]
)
stable_progressor_ncs_plot_all <- df_for_plots %>%
    mutate(time_window = str_glue("{time_window}-year progressors")) %>%
    ggplot(aes(x=ADURW, y=emmean)) +
        geom_line(aes(color=predicted_class)) +
        geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = predicted_class), alpha = 0.2) +
        geom_vline(xintercept = 240, linetype = "dotted") +
        scale_x_continuous(breaks = seq(0, 288, by=48)) +
        scale_color_manual(values = colors_list, name = "predicted class") +
        scale_fill_manual(values = colors_list, name = "predicted class") +
        facet_grid(rows = vars(time_window), cols = vars(TX)) +
        labs(
            x = "Time from randomization (weeks)",
            y = "Model-adjusted mean PACC"
        ) +
        theme_minimal() +
        theme(
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12),
            strip.text.y = element_text(size = 12),
            strip.text.x = element_text(size = 12),
            legend.position = "bottom",
            panel.background = element_rect(fill = "white", color = "gray"),
            plot.background  = element_rect(fill = "white", color = NA)
        ) +
        geom_label(
            data = text_df,
            mapping = aes(label = text),
            x = 12, y = -7.5,
            fill = "white", alpha = 0.5,
            size = 3, hjust = 0, vjust = 0.5, family = "Arial",
            label.size = NA
        )

ggsave(file.path(opt$odir_fig, "pacc_ncs_stable_progressor_overall_all.png"), stable_progressor_ncs_plot_all, width = 8, height = 10, dpi = 300)

# =====================================
# ===== STRATA DEMOGRAPHICS TABLE =====
# =====================================

# get demographics of each group separated by ML strata
a4_predict_with_pacc_strata <- a4_predict_with_pacc %>%
    full_join(a4_stratifications, by = "subj") %>%
    filter(predicted_class == "progressor") %>%
    mutate(time_window = str_c(time_window, "-year progressor"))
a4_predict_with_pacc_strata <- bind_rows(
    a4_predict_with_pacc %>% mutate(time_window = "original cohort"),
    a4_predict_with_pacc_strata
) %>%
    mutate(
        time_window = factor(time_window, levels = c("original cohort", str_c(seq(1,5), "-year progressor"))),
        sex.F.bool = sex == "F"
    ) %>%
    left_join(amyloid_composite %>% filter(VISCODE == "002") %>% select(BID, Composite_Summary), by = c("subj" = "BID")) %>%
    mutate(
        has_amyloid_pet = subj %in% (amyloid_composite %>%
            group_by(BID) %>%
            filter(any(VISCODE == "002"), any(VISCODE == "066")) %>%
            pull(BID) %>%
            unique())
    )

strata_gt <- a4_predict_with_pacc_strata %>%
    tbl_strata(
        strata = TX,
        .tbl_fun = ~ tbl_summary(
            data = .x,
            by = time_window,
            include = c(num_subj, has_amyloid_pet, age, sex.F.bool, apoe, Composite_Summary, pacc.bl, pacc.slope),
            label = list(
                num_subj = "number of subjects",
                has_amyloid_pet = "subjects with FBP PET, count (%)",
                age = "baseline age, mean (SD)",
                sex.F.bool = "female, count (%)",
                apoe = "APOE4 carriership, count (%)",
                Composite_Summary = "baseline cortical FBP SUVR, mean (SD)",
                pacc.bl = "baseline PACC, mean (SD)",
                pacc.slope = "annualized PACC rate-of-change, mean (SD)"
            ),
            statistic = list(
                all_categorical() ~ "{n} ({p})",
                all_continuous() ~ "{mean} ({sd})",
                num_subj ~ "{N}"
            )
        ) %>%
            modify_header(all_stat_cols() ~ "**{level}**", label = "") %>%
            modify_footnote(everything() ~ NA_character_)
    ) %>%
    as_gt()
gtsave(strata_gt, file.path(opt$odir_table, "a4_enriched_demographics.html"))
gtsave(strata_gt, file.path(opt$odir_table, "a4_enriched_demographics.rtf"))

# run statistical tests
tests_df <- a4_predict_with_pacc_strata %>% group_by(TX)
age_ttest <- tests_df %>% t_test(age ~ time_window, ref.group = "original cohort", p.adjust.method = "none") %>% select(TX, .y., group1, group2, p, p.adj.signif)
suvr_ttest <- tests_df %>% t_test(Composite_Summary ~ time_window, ref.group = "original cohort", p.adjust.method = "none") %>% select(TX, .y., group1, group2, p, p.adj.signif)
bl_pacc_ttest <- tests_df %>% t_test(pacc.bl ~ time_window, ref.group = "original cohort", p.adjust.method = "none") %>% select(TX, .y., group1, group2, p, p.adj.signif)
slope_pacc_ttest <- tests_df %>% t_test(pacc.slope ~ time_window, ref.group = "original cohort", p.adjust.method = "none") %>% select(TX, .y., group1, group2, p, p.adj.signif)
sex_fisher <- tests_df %>%
    reframe(
        pairwise_fisher_test(table(pick(sex.F.bool, time_window)), p.adjust.method = "none")
    ) %>%
    filter(group1 == "original cohort") %>%
    select(TX, group1, group2, p, p.adj.signif) %>%
    mutate(.y. = "sex.F.bool")
apoe_fisher <- tests_df %>%
    reframe(
        pairwise_fisher_test(table(pick(apoe, time_window)), p.adjust.method = "none")
    ) %>%
    filter(group1 == "original cohort") %>%
    select(TX, group1, group2, p, p.adj.signif) %>%
    mutate(.y. = "apoe")

# combine test results, apply Bonferroni correction
test_results_all <- bind_rows(age_ttest, suvr_ttest, bl_pacc_ttest, slope_pacc_ttest, sex_fisher, apoe_fisher) %>%
    rename(p.signif = p.adj.signif) %>%
    mutate(p.adj = p.adjust(p, method = "bonferroni")) %>%
    add_significance("p.adj")
test_results_all %>% write_csv(file.path(opt$odir_table, "a4_enriched_demographics_tests.csv"))

# ================================================
# ===== PERFORM ANCOVA ON IMAGING BIOMARKERS =====
# ================================================

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

amyloid_ancova_emm <- map(c(NA, 1:5), ~ fit_amyloid_ancova_plus_stratification(amyloid_df_for_ancova, time_window = .x))
amyloid_ancova_emm_df <- map(amyloid_ancova_emm, ~.x[[1]]) %>% bind_rows()
amyloid_ancova_emm_contrast_df <- map(amyloid_ancova_emm, ~.x[[2]]) %>%
    bind_rows() %>%
    mutate(p.value = p.adjust(p.value, method = "bonferroni"))

# # plot model-adjusted change in SUVR
# dodge <- position_dodge(width=0.9)
# amyloid_ancova_plot <- amyloid_ancova_emm_df %>%
#     ggplot(aes(x = stratification, y = emmean, group = TX, color = TX)) +
#         geom_point(position = dodge) +
#         geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), position = dodge) +
#         labs(
#             title = "model-adjusted change in FBP composite SUVR from 0-weeks to 240-weeks - FBP amyloid PET",
#             y = "change in SUVR"
#         )
# ggsave(file.path(opt$odir_fig, "amyloid_ancova.png"), amyloid_ancova_plot, dpi = 300, width = 8, height = 4, units = "in")

# # plot contrasts in change in SUVR between placebo and Solanezumab
# amyloid_ancova_contrast_plot <- amyloid_ancova_emm_contrast_df %>%
#     ggplot(aes(x = stratification, y = estimate)) +
#         geom_point() +
#         geom_errorbar(aes(ymin = estimate - 1.96*SE, ymax = estimate + 1.96*SE)) +
#         labs(
#             title = str_wrap("difference in model-adjusted change in SUVR between placebo and Solanezumab groups - FBP amyloid PET", width = 60),
#             y = str_wrap("difference in change in SUVR (Placebo - Solanezumab)", width = 35)
#         )
# ggsave(file.path(opt$odir_fig, "amyloid_ancova_contrast.png"), amyloid_ancova_contrast_plot, dpi = 300, width = 8, height = 4, units = "in")

# ================================================
# ===== SAVE IMAGING BIOMARKERS STATS TABLES =====
# ================================================

generate_ancova_table <- function(emm_df, emm_contrast_df) {
    df_for_gt <- emm_df %>%
        select(-df) %>%
        rename(estimate = emmean) %>%
        bind_rows(emm_contrast_df %>%
                    select(-c(df, t.ratio)) %>%
                    rename(TX = contrast)) %>%
        mutate(
            stratification = case_match(stratification,
                "orig"~"original cohort",
                "1"~"1-year progressors",
                "2"~"2-year progressors",
                "3"~"3-year progressors",
                "4"~"4-year progressors",
                "5"~"5-year progressors"
            ),
            p.value.sig = sig_level(p.value),
            p.value.str = if_else(!is.na(p.value),
                str_c(
                    formatC(p.value, format = "e", digits = 2),
                    "<sup>{p.value.sig}</sup>" %>% str_glue()
                ),
                NA
            ),
            across(c(estimate, SE, effect.size, SE.effect.size), ~if_else(is.na(.x), NA, formatC(.x, digits = 3, format = "fg"))),
            estimate.SE = str_glue("{estimate} \u00B1 {SE}"),
            effect.size.SE = if_else(is.na(effect.size), NA, str_glue("{effect.size} \u00B1 {SE.effect.size}"))
        ) %>%
        mutate(
            across(p.value.str, ~ map(.x, html)),
            p.value.str = if_else(is.na(p.value), NA, p.value.str)
        )
    df_for_gt %>%
        select(TX, stratification, n, estimate.SE, lower.CL, upper.CL, effect.size.SE, p.value.str) %>%
        gt(
            rowname_col = "TX",
            groupname_col = "stratification"
        ) %>%
        cols_label(
            TX = "group",
            estimate.SE = "model-adjusted mean (\u00B1 SE)",
            lower.CL = "lower 95% CI",
            upper.CL = "upper 95% CI",
            effect.size.SE = "effect size (\u00B1 SE)",
            p.value.str = "p-value"
        ) %>%
        fmt_number(
            columns = all_of(c("estimate.SE", "effect.size.SE", "lower.CL", "upper.CL")),
            n_sigfig = 4
        ) %>%
        sub_missing(
            columns = everything(),
            missing_text = ""
        )
}

amyloid_ancova_gt <- generate_ancova_table(amyloid_ancova_emm_df, amyloid_ancova_emm_contrast_df)
gtsave(amyloid_ancova_gt, file.path(opt$odir_table, "amyloid_ancova.html"))
gtsave(amyloid_ancova_gt, file.path(opt$odir_table, "amyloid_ancova.rtf"))
