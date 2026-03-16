# ==============================================================================

# Name: centiloid_eval_tables.R
# Author: Braden Yang
# Created: 07/27/2025
# Description: Create gt tables for comparing Centiloid vs. no Centiloid model
#   performance

# ==============================================================================

rm(list = ls())

library(tidyverse)
library(gt)
library(gtExtras)
library(pROC)
library(broom)
library(rstatix)

wdir <- "/home/b.y.yang/sotiraslab/BradenADLongitudinalPrediction"
odir <- file.path(wdir, "tables/centiloid")

acc_df <- read_csv(file.path(wdir, "code/3_TrainTestModel/eval/acc_metrics_CL.csv"))

# run Delong test
y_score_df <- read_csv("code/3_TrainTestModel/eval/y_score_centiloid_experiment.csv")

run_delong <- function(.data) {

    if (length(unique(.data$y_true)) < 2) {
        return(NA)
    }

    df1 <- .data %>% filter(experiment == "base_reduced_sites")
    df2 <- .data %>% filter(experiment == "centiloid")

    roc1 <- roc(df1$y_true, df1$y_score, direction = "<")
    roc2 <- roc(df2$y_true, df2$y_score, direction = "<")

    delong <- roc.test(roc1, roc2, method = "delong", paired = TRUE) %>%
        tidy()
    
    return(delong)

}

delong_df <- y_score_df %>%
    group_by(time_window, leave_out_group) %>%
    reframe(
        delong = run_delong(pick(experiment, y_true, y_score))
    ) %>%
    unnest(delong)
delong_sig <- delong_df %>%
    mutate(p.value.adj = p.adjust(p.value, method = "bonferroni")) %>%
    add_significance("p.value.adj") %>%
    mutate(
        delong.p.value.adj.signif = if_else(p.value.adj.signif == "ns", "", p.value.adj.signif)
    ) %>%
    select(time_window, leave_out_group, delong.p.value.adj.signif)

# compute difference in accuracy metrics
acc_wide <- acc_df %>%
    filter(
        features == "avn",
        eval_dataset == "test_metrics"
    ) %>%
    select(-c(features, eval_dataset, model)) %>%
    left_join(delong_sig, by = c("time_window", "leave_out_group")) %>%
    pivot_wider(
        names_from = experiment,
        values_from = value
    ) %>%
    mutate(
        base.str = as.character(formatC(base, format = "f", digits = 4)),
        centiloid.diff = centiloid - base,
        centiloid.diff.str = as.character(formatC(centiloid.diff, format = "f", digits = 4)),
        centiloid.diff.str = if_else(centiloid.diff > 0, str_c("+", centiloid.diff.str), centiloid.diff.str),
        centiloid.diff.sig = if_else(
            metric == "roc_auc",
            str_c(centiloid.diff.str, delong.p.value.adj.signif),
            centiloid.diff.str
        )
    ) %>%
    select(-c(base, centiloid, centiloid.diff, centiloid.diff.str, delong.p.value.adj.signif)) %>%
    rename(base = base.str, centiloid = centiloid.diff.sig) %>%
    pivot_longer(
        cols = c(base, centiloid),
        names_to = "experiment",
        values_to = "value"
    ) %>%
    pivot_wider(
        names_from = leave_out_group,
        values_from = value
    ) %>%
    mutate(metric = factor(
        metric,
        levels = c("roc_auc", "accuracy", "balanced_accuracy", "f1", "sensitivity", "specificity", "positive_predictive_value", "negative_predictive_value"),
        labels = c("ROC AUC", "Accuracy", "Balanced Accuracy", "F1 Score", "Sensitivity", "Specificity", "PPV", "NPV")
    )) %>%
    group_by(metric) %>%
    mutate(metric_str = if_else(row_number() == 1, metric, "")) %>%
    ungroup()

roc_gt <- acc_wide %>%
    filter(metric == "ROC AUC") %>%
    select(-c(metric, metric_str, starts_with("aggregated"))) %>%
    gt(rowname_col = "time_window") %>%
    tab_stubhead(label = "progression time-horizon (yrs)") %>%
    fmt_number(
        columns = c(where(is.numeric), -time_window),
        n_sigfig = 4
    ) %>%
    sub_missing(
        columns = everything(),
        missing_text = "-"
    ) %>%
    gt_add_divider(column = "OASIS", include_labels = FALSE)

acc_gt <- acc_wide %>%
    filter(metric != "ROC AUC") %>%
    select(-metric) %>%
    gt(rowname_col = "metric_str") %>%
    cols_label(
        time_window = "progression time-horizon (yrs)",
        aggregated_site = "aggregated (site)",
        aggregated_tracer = "aggregated (tracer)"
    ) %>%
    fmt_number(
        columns = c(where(is.numeric), -time_window),
        n_sigfig = 4
    ) %>%
    gt_add_divider(column = "aggregated_site", include_labels = FALSE)

gtsave(acc_gt, file.path(odir, "acc_metrics_CL.html"))
gtsave(acc_gt, file.path(odir, "acc_metrics_CL.rtf"))
gtsave(roc_gt, file.path(odir, "roc_auc_metrics_CL.html"))
gtsave(roc_gt, file.path(odir, "roc_auc_metrics_CL.rtf"))
