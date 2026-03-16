rm(list = ls())

library(tidyverse)
library(ggnewscale)

# ======================
# ===== ROC FIGURE =====
# ======================

roc_df <- read_csv("code/3_TrainTestModel/eval/roc_svm_for_figure.csv")

roc_df <- roc_df %>%
    mutate(
        time_window.str = str_c(time_window, "-year progressors")
    )

roc_plot <- ggplot() +
    geom_line(
        data = roc_df %>% filter(leave_out_group_type == "site"),
        mapping = aes(x = fpr, y = tpr, color = leave_out_group),
        linetype = "solid",
        size = 1, alpha = 0.75
    ) +
    scale_color_brewer(name = "leave-out-site", palette = "Set1") +
    guides(color = guide_legend(nrow = 1)) +
    new_scale_color() +
    geom_line(
        data = roc_df %>% filter(leave_out_group_type == "tracer"),
        mapping = aes(x = fpr, y = tpr, color = leave_out_group),
        linetype = "dotted",
        size = 0.75, alpha = 0.75
    ) +
    scale_color_brewer(name = "leave-out-tracer", palette = "Set1") +
    facet_wrap(vars(time_window.str), nrow = 1) + 
    labs(x = "False positive rate", y = "True positive rate") +
    coord_equal() +
    theme_minimal() +
    theme(
        text = element_text(family = "Arial"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        legend.position = "bottom",
        legend.box.spacing = unit(0.01, "lines")
    )
ggsave("figures/roc_plot.png", roc_plot, width = 15, height = 4, units = "in", dpi = 300)
