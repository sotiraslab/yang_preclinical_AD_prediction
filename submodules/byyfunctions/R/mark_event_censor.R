#' @title Create new column indicating events/right censor
#' 
#' @description This function takes in 3 input columns: time, event and grouping.
#' Within each group, the function marks whether each entry is a non-event (0) or
#' event (1). In the case that no event exists for a particular group, the last
#' entry in time is marked as the right censor (2)
#' 
#' @param .data dataframe with event data
#' @param time_col double column containing times of each entry (dplyr data-mask)
#' @param event_col logical column containing events, where TRUE indicates the
#' presence of an event (dplyr data-mask)
#' @param group_col name of grouping column, usually a subject ID column (string)
#' @param consecutive_events (optional) number of consecutive events that must occur
#' in order for it to be marked as a true event; if not specified, defaults to 1
#' 
#' @return .data with additional columns containing event-censor markings
#' 
#' @export
mark_event_censor <- function(.data, time_col, event_col, group_col, consecutive_events = 1) {

    # TODO: write unit tests for function
    # TODO: make function robust to missing values, without having to explicitly remove NAs
    # TODO: tidyselect interface for group_col
    # cases:
    # - converter
    # - censor
    # - NA data in age and/or diagnosis
    # - not enough conversion instances

    arrange_by_group <- function(.data, arrange_col, group_col) {

        .data <- .data %>%
            dplyr::group_by(dplyr::pick(dplyr::all_of(group_col))) %>%
            dplyr::arrange({{arrange_col}}, .by_group = TRUE) %>%
            dplyr::ungroup()

        return(.data)

    }

    run_length_df <- .data %>%
        arrange_by_group({{time_col}}, group_col = group_col) %>%
        dplyr::reframe(
            run_length = rle({{event_col}}) %>% unclass %>% data.frame,
            .by = group_col
        ) %>%
        tidyr::unnest(run_length)

    event_df <- run_length_df %>%
        dplyr::group_by(dplyr::pick(dplyr::all_of(group_col))) %>%
        dplyr::mutate(
            cumulative_length = cumsum(lengths),
            start_idx = cumulative_length - lengths + consecutive_events
        ) %>%
        dplyr::filter(values == TRUE, lengths >= consecutive_events) %>%
        dplyr::slice_min(start_idx)
    
    .data <- .data %>%
        dplyr::left_join(event_df %>% dplyr::select(dplyr::all_of(group_col), start_idx), by = group_col) %>%
        dplyr::group_by(dplyr::pick(dplyr::all_of(group_col))) %>%
        dplyr::mutate(
            event = dplyr::row_number() == start_idx,
            censor = dplyr::if_else(
                all(is.na(event)) & age == max(age),
                TRUE,
                FALSE
            ),
            event_censor = dplyr::case_when(
                event ~ 1,
                censor ~ 2,
                .default = NA
            )
        ) %>%
        dplyr::select(-start_idx) %>%
        dplyr::ungroup()
    
    return(.data)

}