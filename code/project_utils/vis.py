import numpy as np
import pandas as pd
from plotnine import *


def spaghetti(
    cdr_df,
    pet_df,
    group_col = "subj",
    time_col = "age",
    event_col = "cdr",
    tracer_col = "tracer",
    event_levels = [0, 0.5, 1, 2, 3],
    facet_by_site = False,
    figure_size = (15, 5)
):

    """
    Plot spaghetti plot of subjects' longitudinal CDR and their baseline
    amyloid PET scan(s)
    """
    
    cdr_df = cdr_df.copy(); pet_df = pet_df.copy()
    
    # filter subjects
    common_subj = np.intersect1d(cdr_df[group_col], pet_df[group_col])
    cdr_df = cdr_df[cdr_df[group_col].isin(common_subj)]
    pet_df = pet_df[pet_df[group_col].isin(common_subj)]
    
    # get ordered columns
    time_min_df = pet_df.groupby(group_col)[time_col].min()
    progressor_df = pet_df.groupby(group_col)["y"].unique().apply(lambda x: x[0])
    subj_order = pd.concat([time_min_df, progressor_df], axis = 1).sort_values(["y", "age"]).index
    
    cdr_df[".event"] = pd.Categorical(cdr_df[event_col], categories = event_levels, ordered = True)
    cdr_df[".group"] = pd.Categorical(cdr_df[group_col]).set_categories(subj_order, ordered=True)
    pet_df[".group"] = pd.Categorical(pet_df[group_col]).set_categories(subj_order, ordered=True)
    
    # remove NA columns
    cdr_df.dropna(subset=[".event", ".group"], inplace=True)
    
    # merge site and converter status
    pet_merge = pet_df.loc[:, ["site", group_col, "y"]]
    cdr_df = pd.merge(cdr_df, pet_merge, how = "left", on = group_col)

    # plot
    g = ggplot(data = cdr_df, mapping = aes(x = ".group", y = time_col, group = "subj")) \
        + geom_point(mapping = aes(color = ".event")) \
        + geom_line(mapping = aes(color = ".event")) \
        + geom_point(data = pet_df, mapping = aes(shape = "tracer"), color = "red", fill = None, size = 0.5, alpha = 0.5) \
        + scale_y_reverse() \
        + theme(figure_size = figure_size, axis_text_x = element_text(angle = 90))

    if facet_by_site:
        g = g + facet_wrap("site", nrow = 1, scales = "free_x")

    return g