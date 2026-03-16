# ==============================================================================

# Name: nested_model_delong.R
# Author: Braden Yang
# Created: 07/27/2025
# Description: Create gt tables to report ROC-AUC differences and Delong test
#   results of nested model analysis

# ==============================================================================

rm(list = ls())

library(tidyverse)
library(pROC)
library(broom)
library(rstatix)
library(gt)
# library(gtExtras)

y_score_df <- read_csv("code/3_TrainTestModel/eval/revisions/y_score_nested.csv")

run_delong <- function(.data) {

    if (length(unique(.data$y_true)) < 2) {
        return(NA)
    }

    delong_list <- list()
    for (f in c("an", "vn", "av", "a", "v", "n")) {
    # for (f in c("an", "vn", "n")) {
        df1 <- .data %>% filter(features == "avn")
        df2 <- .data %>% filter(features == f)

        roc1 <- roc(df1$y_true, df1$y_score, direction = "<")
        roc2 <- roc(df2$y_true, df2$y_score, direction = "<")

        delong <- roc.test(roc1, roc2, method = "delong", paired = TRUE) %>%
            tidy() %>%
            mutate(nested_model = f)
        delong_list[[length(delong_list) + 1]] <- delong
    }
    delong_df <- bind_rows(delong_list)
    
    return(delong_df)

}

# filter
y_score_df <- y_score_df %>%
    filter(
        leave_out_group %in% c("A4", "ADNI", "AIBL", "BLSA", "HABS", "MCSA", "OASIS")
    )

delong_df <- y_score_df %>%
    group_by(time_window, leave_out_group) %>%
    reframe(
        delong = run_delong(pick(features, y_true, y_score))
    ) %>%
    unnest(delong)
delong_df <- delong_df %>%
    mutate(p.value.adj = p.adjust(p.value, method = "bonferroni")) %>%
    add_significance("p.value.adj") %>%
    add_significance("p.value") %>%
    mutate(
        p.value.signif = if_else(p.value.signif == "ns", "", p.value.signif),
        p.value.adj.signif = if_else(p.value.adj.signif == "ns", "", p.value.adj.signif),
        estimate.diff = estimate2 - estimate1,
        across(c(estimate1, estimate2, estimate.diff), ~ if_else(is.na(.x), NA, as.character(formatC(.x, format = "f", digits = 4))), .names = "{.col}.str"),
        estimate_p.value = case_when(
            is.na(estimate2) ~ NA,
            p.value.adj.signif == "" ~ estimate2.str,
            .default = str_glue("{estimate2.str}{p.value.adj.signif}")
        ),
        estimate.diff_p.value.adj = case_when(
            is.na(estimate.diff) ~ NA,
            .default = if_else(estimate.diff > 0,
                str_glue("+{estimate.diff.str}{p.value.adj.signif}"),
                str_glue("{estimate.diff.str}{p.value.adj.signif}"))
        ),
        estimate.diff_p.value = case_when(
            is.na(estimate.diff) ~ NA,
            .default = if_else(estimate.diff > 0,
                str_glue("+{estimate.diff.str}{p.value.signif}"),
                str_glue("{estimate.diff.str}{p.value.signif}"))
        )
    )
avn_df <- delong_df %>% select(time_window, leave_out_group, estimate1.str) %>%
    unique() %>%
    mutate(nested_model = "avn") %>%
    rename(estimate.diff_p.value.adj = estimate1.str)
delong_table <- bind_rows(
    delong_df %>% select(time_window, leave_out_group, nested_model, estimate.diff_p.value.adj),
    avn_df
) %>%
    mutate(nested_model = factor(
        nested_model,
        levels = c("avn", "an", "vn", "av", "a", "v", "n"),
        labels = c("amyloid + volume + non-imaging", "amyloid + non-imaging", "volume + non-imaging", "amyloid + volume", "amyloid", "volume", "non-imaging")
    )) %>%
    arrange(time_window, leave_out_group, nested_model) %>%
    pivot_wider(
        names_from = leave_out_group,
        values_from = estimate.diff_p.value.adj
    ) %>%
    drop_na(nested_model)
    # %>%
    # relocate(FBP, .after = "OASIS")

delong_gt <- delong_table %>%
    gt(rowname_col = "time_window") %>%
    tab_stubhead(label = "progression time-horizon (yrs)") %>%
    cols_label(
        nested_model = "feature set",
    ) %>%
    tab_spanner(
        label = "leave-out group",
        columns = delong_table %>% select(A4:OASIS) %>% colnames
    ) %>%
    # gtExtras::gt_add_divider(column = "OASIS", include_labels = FALSE) %>%
    sub_missing(missing_text = "-")
    
gtsave(delong_gt, "tables/nested_model_delong_kleen.html")
gtsave(delong_gt, "tables/nested_model_delong_kleen.rtf")
