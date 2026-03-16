rm(list = ls())

library(tidyverse)
library(patchwork)
library(tidytext)

# plot Haufe coefficients
haufe <- read_csv("code/4_FeatureImportance/haufe.csv") %>%
    mutate(
        feature_name_raw = if_else(
            feature_type == "nonimg",
            str_replace(feature_name, "nonimg__", ""),
            str_match(feature_name, ".*__(.*)\\..*")[,2]
        ),
        feature_type_fct = factor(feature_type, levels = c("amyloid", "volume", "nonimg"), labels = c("amyloid PET SUVR", "volumetric", "non-imaging"))
    )
haufe_boxplots <- haufe %>%   # all time windows aggregated
    ggplot(aes(x = tidytext::reorder_within(feature_name_raw, haufe, feature_type_fct, median), y = haufe, color = feature_type_fct)) +
        geom_hline(yintercept = 0) +
        geom_boxplot(alpha = 0.75) +
        labs(x = "features", y = "Haufe-corrected weights") +
        facet_grid(cols = vars(feature_type_fct), scales = "free_x", space = "free_x") + 
        scale_color_brewer(palette = "Set1", name = "feature type") +
        scale_x_reordered() + 
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.background = element_rect(fill = "white", color = NA),
            plot.background  = element_rect(fill = "white", color = NA),
            legend.position = "top",
            strip.text = element_blank(), strip.background = element_blank()
        )
ggsave(
    "figures/haufe/haufe_boxplot.png",
    haufe_boxplots,
    width = 25, height = 10, units = "in",
    dpi = 500
)

# plot slope
haufe_slope <- read_csv("code/4_FeatureImportance/haufe_slope.csv") %>%
    mutate(
        feature_name_raw = if_else(
            feature_type == "nonimg",
            str_replace(feature_name, "nonimg__", ""),
            str_match(feature_name, ".*__(.*)\\..*")[,2]
        ),
        feature_type_fct = factor(feature_type, levels = c("amyloid", "volume", "nonimg"), labels = c("amyloid PET SUVR", "volumetric", "non-imaging"))
    )

haufe_slope_boxplots <- haufe_slope %>%
    ggplot(aes(x = tidytext::reorder_within(feature_name_raw, slope, feature_type_fct, median), y = slope, color = feature_type_fct)) +
        geom_hline(yintercept = 0) +
        geom_boxplot(alpha = 0.75) +
        labs(x = "features", y = "rate-of-change of absolute value Haufe-corrected weights") +
        facet_grid(cols = vars(feature_type_fct), scales = "free_x", space = "free_x") + 
        scale_color_brewer(palette = "Set1", name = "feature type") +
        scale_x_reordered() + 
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.background = element_rect(fill = "white", color = NA),
            plot.background  = element_rect(fill = "white", color = NA),
            legend.position = "top",
            strip.text = element_blank(), strip.background = element_blank()
        )
ggsave(
    "figures/haufe/haufe_slope_boxplot.png",
    haufe_slope_boxplots,
    width = 25, height = 10, units = "in",
    dpi = 500
)
