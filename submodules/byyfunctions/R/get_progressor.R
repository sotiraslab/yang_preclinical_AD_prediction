#' @title Mark subjects as stable or progressors based on chronological ordering of binary events
#' 
#' @description This function creates a new logical column named `progressor` which contains
#' subject-wise indications of progressor status. Subjects can fall into one of 3 categories:
#' - progressor: subjects who start at event = FALSE and record an event = TRUE at a later time;
#'   marked as TRUE
#' - stable: subjects who start at event = FALSE and never record an event = TRUE; marked as FALSE
#' - NA: subjects who are uncategorizable (see details)
#' 
#' @details - Uncategorizable subjects can be due to one of three cases: (1) if a stable subject hasn't
#' remained stable for at least `stable_min_time`; (2) if a progressor subject doesn't start as
#' stable at baseline. These cases are marked as NA in the `progressor` column; (3) if a subject
#' progresses to an event but then at a later time reverts to the non-event (this is done only if
#' `exclude_revert = TRUE`)
#' - If a subject has more than one baseline specified in the baseline column, then the function defaults
#' to the first baseline in time
#' - extra columns that can be kept via the `keep_extra_col` argument:
#'   - .event = logical column indicating locations of events (useful when `event_col` is a character col)
#'   - .revert = logical column indicating subject-wise reverter status
#'   - .bl = logical column indicating the baseline row
#'   - .time_to_progression = numeric column indicating time from baseline to first event
#' 
#' @param .data dataframe or tibble containing subject data to categorize as stable/progressor
#' @param event_col name of column containing event markings. This can either be a boolean column
#' where events are marked as TRUE, or a character column with multiple event/non-event labels. In
#' the latter case, specify argument `event_label` to indicate which of the labels correspond to a
#' true event
#' @param time_col name of column containing time information
#' @param group_col name of grouping column, usually a subject identifier column
#' @param event_label (optional) character vector containing the labels that correspond to a positive
#' event; specify this if `event_col` is not a boolean column
#' @param baseline_col (optional) name of column containing locations of the baseline to consider; if
#' not specified, defaults to baseline = subject's first entry in time
#' @param time_window (optional) length of time window from baseline to search for an event to check
#' for progressor status (default = Inf, i.e. look indefinitely into future)
#' @param stable_min_time (optional) length of minimum time that a stable subject must remain stable
#' to be labeled as such (default = -Inf, i.e. no minimum)
#' @param progressor_nonevent_start (optional) whether to impose constraint that progressors must start
#' as stable at baseline (default = TRUE)
#' @param exclude_revert (optional) whether to exclude reverters and categorize them as NA (default = TRUE)
#' @param exclude_stable_event (optional) whether to categorize a stable as NA if they have a positive
#' event at any time, even past the stable_min_time window
#' @param keep_extra_col (optional) whether to keep extra cols `.event`, which is the boolean event col
#' (useful when `event_col` is a character col), `.revert`, which indicates which subjects are reverters,
#' `.bl`, which indicates the baseline row, and `.time_to_progression`, which is the time from baseline to
#' first event
#' @param progressor_colname (optional) name of new column containing progressor labels (default =
#' "progressor")
#' 
#' @export

# TODO: make robust to NA values in any of the input columns (especially age and event cols)
# TODO: make tests for NA values

get_progressor <- function(
    .data,
    event_col,
    time_col,
    group_col,
    event_label = NULL,
    baseline_col = NULL,
    time_window = Inf,
    stable_min_time = -Inf,
    progressor_nonevent_start = TRUE,
    exclude_revert = TRUE,
    exclude_stable_event = TRUE,
    keep_extra_col = FALSE,
    progressor_colname = "progressor"
) {

    max_na <- function(x, na_val = -Inf) {
        # wrapper of `max` which checks if all values of vector
        # x are NA, and if so return na_val
        # - defaults to -Inf if NA so that the max cannot be
        #   greater than any other number
        if (all(is.na(x))) {
            return(na_val)
        } else {
            return(max(x))
        }
    }

    min_na <- function(x, na_val = Inf) {
        # wrapper of `min` which checks if all values of vector
        # x are NA, and if so return na_val
        # - defaults to Inf if NA so that the min cannot be
        #   less than any other number
        if (all(is.na(x))) {
            return(na_val)
        } else {
            return(min(x))
        }
    }

    which_first <- function(x) {
        # wrapper of `which` which returns only the first element
        # when multiple matches are made
        return(which(x)[1])
    }

    df <- .data %>%
        dplyr::group_by(dplyr::pick(dplyr::all_of(group_col))) %>%  # group by
        dplyr::arrange(.data[[time_col]], .by_group = TRUE)  # arrange by time column

    # create event col
    if (typeof(.data[[event_col]]) == "logical") {
        df <- df %>% dplyr::mutate(.event = .data[[event_col]])
    } else {
        df <- df %>% dplyr::mutate(.event = .data[[event_col]] %in% event_label)
    }

    # mark subjects that revert to normal after event occurs
    df <- df %>% dplyr::mutate(
        .revert = min_na(which(.event)) < max_na(which(!.event))
    )

    # create baseline column
    if (is.null(baseline_col)) {
        # baseline = first row
        df <- df %>% dplyr::mutate(.bl = .data[[time_col]] == min(.data[[time_col]], na.rm = TRUE))
    } else {
        # baseline = baseline_col
        df <- df %>% dplyr::mutate(.bl = .data[[baseline_col]])
    }

    # select progressors
    df <- df %>%
        dplyr::mutate(
            .progressor = dplyr::if_else(  # only select times within the time window
                any(.event[
                    .data[[time_col]] >= .data[[time_col]][which_first(.bl)] &
                    .data[[time_col]] < .data[[time_col]][which_first(.bl)] + time_window
                ]),
                TRUE, FALSE
            ),
            .progressor = dplyr::if_else(  # stable min time criterion
                !.progressor & (max(.data[[time_col]]) - .data[[time_col]][which_first(.bl)]) < stable_min_time,
                NA, .progressor
            ),
            .progressor = dplyr::if_else(  # progressor starts at .event == 0
                .progressor & progressor_nonevent_start & .event[which_first(.bl)],
                NA, .progressor
            ),
            .progressor = dplyr::if_else(  # progressor reverts
                exclude_revert & .revert,
                NA, .progressor
            ),
            .progressor = dplyr::if_else(  # stable experiences event
                !.progressor & exclude_stable_event & any(.event),
                NA, .progressor
            ),
            .time_to_progression = .data[[time_col]][which_first(.event)] - .data[[time_col]][which_first(.bl)]
        )

    # rename progressor column
    df <- df %>% rename_with(~ progressor_colname, .progressor)

    # remove extra columns
    if (!keep_extra_col) {
        df <- df %>% dplyr::select(-c(.bl, .event, .revert, .time_to_progression))
    }

    # ungroup
    df <- df %>% dplyr::ungroup()

    return(df)

}
