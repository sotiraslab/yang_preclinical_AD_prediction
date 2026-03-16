#' @title Generate roi_features.RData for byyfunctions package
#' 
#' @description Internal function used to generate the dataset `roi_features.RData`,
#' which contains a tibble of FreeSurfer cortical/subcortical ROI from the DK atlas,
#' along with their corresponding rows in the ADNI/OASIS merged tables and ggseg-
#' compatible labels
#' 
#' @param roi_features_csv path to ROI features CSV file
#' @param outfile path of output RData file containing tibble
#' 
#' @noRd
.create_roi_features <- function(

    roi_features_csv = "/home/b.y.yang/BradenADLongitudinalPrediction/data/adni_oasis_roi_features.csv",
    summary_cort_roi_csv = "/home/b.y.yang/BradenADLongitudinalPrediction/data/summary_cortical_roi.csv",
    outfile = "/home/b.y.yang/BradenADLongitudinalPrediction/code/byyfunctions/data/roi_features.RData"
    ) {

    # load table of ROI feature names
    roi_features_raw <- readr::read_csv(roi_features_csv, col_types="ccccc") %>%
        dplyr::mutate(
            adni_suvr = byyfunctions::paste0_na(adni_suvr, ".av45"),
            oasis_suvr = byyfunctions::paste0_na(oasis_suvr, ".pup"),
            adni_vol = byyfunctions::paste0_na(adni_vol, ".av45"),
            oasis_vol = byyfunctions::paste0_na(oasis_vol, ".fs")
        ) # append suffix to each column label

    # tidy data
    roi_features <- roi_features_raw %>%
        tidyr::pivot_longer(
            cols = -"fs_label",
            names_pattern = "(.*)_(.*)",
            names_to = c("study", "feature_type"),
            values_to = "col"
        ) %>%
        tidyr::drop_na() %>%  # remove rows with NA
        byyfunctions::add_ggseg_label("fs_label")  # map ggseg labels to each ROI feature name

    # categorize each ROI (ADNI summary cortical, other cortical, subcortical)
    summary_roi_df <- readr::read_csv(summary_cort_roi_csv) %>%
        dplyr::mutate(fs_label = roi)
    roi_type_labels <- c("Summary cortical ROI", "Other cortical ROI", "Subcortical ROI")
    roi_features <- roi_features %>%
        dplyr::mutate(
            roi_type = ifelse(
                fs_label %in% summary_roi_df$fs_label,
                1,
                ifelse(startsWith(fs_label, "ctx."), 2, 3))
        ) %>%
        dplyr::mutate(roi_type = factor(roi_type, levels = c(1, 2, 3), labels = roi_type_labels))
    
    # if summary cortical, add lobe
    roi_features <- roi_features %>%
        dplyr::left_join(summary_roi_df %>% dplyr::select(fs_label, lobe), by = "fs_label")

    # save in "code/byyfunctions/data"
    save(roi_features, file = outfile)

}

#' @title Match data by date of entry (slow version)
#' 
#' @description Match a recorded data point to another recorded data point that is closest
#' in time to it
#' 
#' @details For example, match the CDR value closest in time to each subject's
#' PET scan. Each row containing non-NA data in the reference column will be
#' matched with a row containins non-NA data in the matching column which is closest
#' in value to the reference value. User may also specify additional columns to pull
#' data from the matched rows.
#' 
#' @param .data dataframe containing the reference data and matching data
#' @param ref_col name of column to serve as reference value to match to
#' @param match_col name of column to get values to match
#' @param select_col names of column(s) to pull additional data from matched rows
#' @param group_col name of group column, usually a subject identifier column
#' @param max_diff (optional) the maximum allowed time difference for a valid match;
#'     if no entires fall within this allowed time difference, NA will be matched
#' @param match_df (optional) a separate dataframe or tibble from which matched data
#'     will be pulled from. If specified, then `match_col` and `select_col` will be
#'     pulled from `match_df`. If not specified, then the aforementioned columns
#'     are pulled from `.data`
#' @param idx_col (optional) string of name of column to store index from matched
#'     dataframe
#' @param diff_col (optional) string of name of column to store absolute difference
#' @param suffix (optional) suffix appended to matched columns. If not specified,
#'     defaults to ".match"
#' 
#' @return .data with additional columns containing matched data
#' 
#' @section Notes:
#' When using dplyr::group_by, all subsequent functions on a dataframe
#' are performed per group. For our purposes, since we want to find the row index
#' corresponding to the matched entry for each reference value, the returned
#' index is of the group dataframe rather than the entire dataframe, which is not
#' the desired result. To combat this, we make a new column which indicates the
#' global index; that way, our matching function can reference this column and
#' return the global index rather than the group-wise index
#' 
#' @section WARNING:
#' This is an alternate implementation of `match_data`, but now using tidyverse
#' functions. The code itself is much more elegant, but in practice this function
#' takes more than 10x the amount of time to execute compared to the original.
#' Therefore I would recommend to just use the original, as it does exactly the
#' same thing as this implementation.
.match_data_slow <- function(
    ref_df,
    ref_col,
    match_col,
    select_col,
    group_col,
    max_diff = Inf,
    match_df = NULL,
    idx_col = NULL,
    diff_col = NULL,
    suffix = ".match"
) {

    match_row2row <- function(v, df) {

        # helper function to match a single value v to the closest row in df, then
        # return that row (only the specified columns to select)

        # if v is NA, return NA
        if (is.na(v)) {
            return(NA)
        # else return matched row in df
        } else {
            
            result <- df %>%
                mutate(abs_diff = abs(v - .data[[match_col]])) %>%  # compute absolute difference
                filter(abs_diff <= max_diff) %>%  # filter rows where abs diff is less than max allowed difference
                slice_min(abs_diff)  # get row with minimum difference

            # select relevant columns
            result <- result %>%
                select(
                    match_col,
                    all_of(select_col),
                    if (!is.null(idx_col)) {all_of(idx_col)},
                    if (!is.null(diff_col)) {all_of("abs_diff")}
                )
            
            # rename absolute diff column if necessary
            if (!is.null(diff_col)) {
                result <- result %>% rename({{diff_col}} := abs_diff)
            }

            # rename with suffix
            result <- result %>% rename_with(~str_c(.x, suffix), everything())

            return(result)

        }

    }

    # if match_df specified, use match_df
    if (!is.null(match_df)) {

        if (!is.null(idx_col)) {
            match_df <- match_df %>% mutate({{idx_col}} := row_number())
        }

        ref_df <- ref_df %>% nest(.by = all_of(group_col))
        match_df <- match_df %>% nest(.by = all_of(group_col))
        ref_df <- ref_df %>% left_join(match_df, by = group_col)

    # else repeat column from ref_df
    } else {

        if (!is.null(idx_col)) {
            ref_df <- ref_df %>% mutate({{idx_col}} := row_number())
        }

        ref_df <- ref_df %>%
            nest(.by = all_of(group_col)) %>%
            mutate(data.y = data) %>%
            rename(data.x = data)

    }

    # match data
    match_df <- ref_df %>%
        mutate(
            data.x = pmap(
                list(x = data.x, y = data.y),
                function(x, y) {
                    x %>% mutate(
                        match = map(.data[[ref_col]], ~ match_row2row(.x, y))
                    )
                }
            )
        ) %>%
        select(-data.y) %>%
        unnest(data.x) %>%
        unnest(match, keep_empty = TRUE)

    return(match_df)

}