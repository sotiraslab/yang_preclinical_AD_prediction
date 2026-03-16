# ==============================================================================

# Name: eval_tables.R
# Author: Braden Yang
# Created: 02/23/2024
# Description: Create gt tables for model evaluation metrics

# ==============================================================================

rm(list = ls())

library(tidyverse)
library(gt)
library(gtExtras)

wdir <- "/home/b.y.yang/sotiraslab/BradenADLongitudinalPrediction"
odir <- file.path(wdir, "tables")

acc_df <- read_csv(file.path(wdir, "code/3_TrainTestModel/eval/acc_metrics.csv"))
acc_table <- acc_df %>%
    mutate(leave_out_group = factor(
        leave_out_group,
        levels = c("A4", "ADNI", "AIBL", "BLSA", "HABS", "MCSA", "OASIS", "aggregated_site", "FBP", "PIB", "aggregated_tracer")
    )) %>%
    arrange(leave_out_group) %>%
    filter(
        features == "avn",
        eval_dataset == "test_metrics",
        experiment == "base"
    ) %>%
    select(-c(features, eval_dataset, experiment, model)) %>%
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
acc_gt <- acc_table %>%
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
roc_gt <- acc_table %>%
    filter(metric == "ROC AUC") %>%
    select(-c(metric, metric_str, aggregated_site, aggregated_tracer)) %>%
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

gtsave(acc_gt, file.path(odir, "acc_metrics.html"))
gtsave(acc_gt, file.path(odir, "acc_metrics.rtf"))
gtsave(roc_gt, file.path(odir, "roc_auc_metrics.html"))
gtsave(roc_gt, file.path(odir, "roc_auc_metrics.rtf"))
