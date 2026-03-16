
a <- dplyr::tribble(
    ~group, ~time, ~event,
    "a", 1, FALSE,
    "a", 2, FALSE,
    "a", 3.5, TRUE,
    "b", 1, FALSE,
    "b", 2, FALSE,
    "b", 4, FALSE
) %>% sample_frac(1)

test_that("correct labeling of stable/progressor", {
    expect_equal(
        get_progressor(
            a,
            event_col = "event",
            time_col = "time",
            group_col = "group"
        ),
        dplyr::tribble(
            ~group, ~time, ~event, ~progressor,
            "a", 1, FALSE, TRUE,
            "a", 2, FALSE, TRUE,
            "a", 3.5, TRUE, TRUE,
            "b", 1, FALSE, FALSE,
            "b", 2, FALSE, FALSE,
            "b", 4, FALSE, FALSE,
        )
    )
})

test_that("time window test", {
    expect_equal(
        get_progressor(
            a,
            event_col = "event",
            time_col = "time",
            group_col = "group",
            time_window = 2,
            exclude_stable_event = FALSE
        ),
        dplyr::tribble(
            ~group, ~time, ~event, ~progressor,
            "a", 1, FALSE, FALSE,
            "a", 2, FALSE, FALSE,
            "a", 3.5, TRUE, FALSE,
            "b", 1, FALSE, FALSE,
            "b", 2, FALSE, FALSE,
            "b", 4, FALSE, FALSE,
        )
    )
})

test_that("stable minimum time test", {
    expect_equal(
        get_progressor(
            a,
            event_col = "event",
            time_col = "time",
            group_col = "group",
            stable_min_time = 5
        ),
        dplyr::tribble(
            ~group, ~time, ~event, ~progressor,
            "a", 1, FALSE, TRUE,
            "a", 2, FALSE, TRUE,
            "a", 3.5, TRUE, TRUE,
            "b", 1, FALSE, NA,
            "b", 2, FALSE, NA,
            "b", 4, FALSE, NA,
        )
    )
})

test_that("extra columns kept", {
    expect_equal(
        get_progressor(
            a,
            event_col = "event",
            time_col = "time",
            group_col = "group",
            keep_extra_col = TRUE
        ),
        dplyr::tribble(
            ~group, ~time, ~event, ~.event, ~.revert, ~.bl, ~progressor, ~.time_to_progression,
            "a", 1, FALSE, FALSE, FALSE, TRUE, TRUE, 2.5,
            "a", 2, FALSE, FALSE, FALSE, FALSE, TRUE, 2.5,
            "a", 3.5, TRUE, TRUE, FALSE, FALSE, TRUE, 2.5,
            "b", 1, FALSE, FALSE, FALSE, TRUE, FALSE, NA,
            "b", 2, FALSE, FALSE, FALSE, FALSE, FALSE, NA,
            "b", 4, FALSE, FALSE, FALSE, FALSE, FALSE, NA,
        )
    )
})

test_that("new column name", {
    expect_equal(
        get_progressor(
            a,
            event_col = "event",
            time_col = "time",
            group_col = "group",
            progressor_colname = "new_col"
        ),
        dplyr::tribble(
            ~group, ~time, ~event, ~new_col,
            "a", 1, FALSE, TRUE,
            "a", 2, FALSE, TRUE,
            "a", 3.5, TRUE, TRUE,
            "b", 1, FALSE, FALSE,
            "b", 2, FALSE, FALSE,
            "b", 4, FALSE, FALSE,
        )
    )
})


b <- dplyr::tribble(
    ~group, ~time, ~event,
    "a", 1, TRUE,
    "a", 2, TRUE,
    "a", 3.5, TRUE,
    "b", 1, FALSE,
    "b", 2, FALSE,
    "b", 4, FALSE
) %>% sample_frac(1)

test_that("progressor starts at non-event", {
    expect_equal(
        get_progressor(
            b,
            event_col = "event",
            time_col = "time",
            group_col = "group"
        ),
        dplyr::tribble(
            ~group, ~time, ~event, ~progressor,
            "a", 1, TRUE, NA,
            "a", 2, TRUE, NA,
            "a", 3.5, TRUE, NA,
            "b", 1, FALSE, FALSE,
            "b", 2, FALSE, FALSE,
            "b", 4, FALSE, FALSE,
        )
    )
})

test_that("progressor starts at non-event (allow)", {
    expect_equal(
        get_progressor(
            b,
            event_col = "event",
            time_col = "time",
            group_col = "group",
            progressor_nonevent_start = FALSE
        ),
        dplyr::tribble(
            ~group, ~time, ~event, ~progressor,
            "a", 1, TRUE, TRUE,
            "a", 2, TRUE, TRUE,
            "a", 3.5, TRUE, TRUE,
            "b", 1, FALSE, FALSE,
            "b", 2, FALSE, FALSE,
            "b", 4, FALSE, FALSE,
        )
    )
})


c <- dplyr::tribble(
    ~group, ~time, ~event, ~bl,
    "a", 1, FALSE, FALSE,
    "a", 2, FALSE, TRUE,
    "a", 3.5, TRUE, FALSE,
    "b", 1, FALSE, FALSE,
    "b", 2, FALSE, FALSE,
    "b", 4, FALSE, TRUE
) %>% sample_frac(1)

test_that("new baseline", {
    expect_equal(
        get_progressor(
            c,
            event_col = "event",
            time_col = "time",
            group_col = "group",
            baseline_col = "bl",
            time_window = 2
        ),
        dplyr::tribble(
            ~group, ~time, ~event, ~bl, ~progressor,
            "a", 1, FALSE, FALSE, TRUE,
            "a", 2, FALSE, TRUE,  TRUE,
            "a", 3.5, TRUE, FALSE,  TRUE,
            "b", 1, FALSE, FALSE, FALSE,
            "b", 2, FALSE, FALSE, FALSE,
            "b", 4, FALSE, TRUE, FALSE,
        )
    )
    expect_equal(
        get_progressor(
            c,
            event_col = "event",
            time_col = "time",
            group_col = "group",
            baseline_col = "bl",
            time_window = 2,
            stable_min_time = 1
        ),
        dplyr::tribble(
            ~group, ~time, ~event, ~bl, ~progressor,
            "a", 1, FALSE, FALSE, TRUE,
            "a", 2, FALSE, TRUE,  TRUE,
            "a", 3.5, TRUE, FALSE,  TRUE,
            "b", 1, FALSE, FALSE, NA,
            "b", 2, FALSE, FALSE, NA,
            "b", 4, FALSE, TRUE, NA,
        )
    )
})


d <- dplyr::tribble(
    ~group, ~time, ~event, ~bl,
    "a", 0.5, FALSE, FALSE,
    "a", 1, TRUE, FALSE,
    "a", 2, FALSE, TRUE,
    "a", 3.5, TRUE, FALSE,
    "b", 1, FALSE, FALSE,
    "b", 2, TRUE, FALSE,
    "b", 4, TRUE, TRUE,
) %>% sample_frac(1)

test_that("reverter test", {
    # reverter is excluded
    expect_equal(
        get_progressor(
            d,
            event_col = "event",
            time_col = "time",
            group_col = "group"
        ),
        dplyr::tribble(
            ~group, ~time, ~event, ~bl, ~progressor,
            "a", 0.5, FALSE, FALSE, NA,
            "a", 1, TRUE, FALSE, NA,
            "a", 2, FALSE, TRUE, NA,
            "a", 3.5, TRUE, FALSE, NA,
            "b", 1, FALSE, FALSE, TRUE,
            "b", 2, TRUE, FALSE, TRUE,
            "b", 4, TRUE, TRUE, TRUE,
        )
    )
    # reverter is allowed
    expect_equal(
        get_progressor(
            d,
            event_col = "event",
            time_col = "time",
            group_col = "group",
            exclude_revert = FALSE
        ),
        dplyr::tribble(
            ~group, ~time, ~event, ~bl, ~progressor,
            "a", 0.5, FALSE, FALSE, TRUE,
            "a", 1, TRUE, FALSE, TRUE,
            "a", 2, FALSE, TRUE, TRUE,
            "a", 3.5, TRUE, FALSE, TRUE,
            "b", 1, FALSE, FALSE, TRUE,
            "b", 2, TRUE, FALSE, TRUE,
            "b", 4, TRUE, TRUE, TRUE,
        )
    )
    # reverter found even with different baseline
    expect_equal(
        get_progressor(
            d,
            event_col = "event",
            time_col = "time",
            group_col = "group",
            baseline_col = "bl"
        ),
        dplyr::tribble(
            ~group, ~time, ~event, ~bl, ~progressor,
            "a", 0.5, FALSE, FALSE, NA,
            "a", 1, TRUE, FALSE, NA,
            "a", 2, FALSE, TRUE, NA,
            "a", 3.5, TRUE, FALSE, NA,
            "b", 1, FALSE, FALSE, NA,
            "b", 2, TRUE, FALSE, NA,
            "b", 4, TRUE, TRUE, NA,
        )
    )
    expect_equal(
        get_progressor(
            d,
            event_col = "event",
            time_col = "time",
            group_col = "group",
            baseline_col = "bl",
            progressor_nonevent_start = FALSE
        ),
        dplyr::tribble(
            ~group, ~time, ~event, ~bl, ~progressor,
            "a", 0.5, FALSE, FALSE, NA,
            "a", 1, TRUE, FALSE, NA,
            "a", 2, FALSE, TRUE, NA,
            "a", 3.5, TRUE, FALSE, NA,
            "b", 1, FALSE, FALSE, TRUE,
            "b", 2, TRUE, FALSE, TRUE,
            "b", 4, TRUE, TRUE, TRUE,
        )
    )
})


e <- dplyr::tribble(
    ~group, ~time, ~label,
    "a", 1, "label1",
    "a", 2, "label3",
    "a", 3.5, "label2",
    "b", 1, "label3",
    "b", 2, "label2",
    "b", 4, "label1",
) %>% sample_frac(1)

test_that("event labels", {
    expect_equal(
        get_progressor(
            e,
            event_col = "label",
            time_col = "time",
            group_col = "group",
            event_label = c("label2", "label3")
        ),
        dplyr::tribble(
            ~group, ~time, ~label, ~progressor,
            "a", 1, "label1", TRUE,
            "a", 2, "label3", TRUE,
            "a", 3.5, "label2", TRUE,
            "b", 1, "label3", NA,
            "b", 2, "label2", NA,
            "b", 4, "label1", NA,
        )
    )
})


f <- dplyr::tribble(
    ~group, ~time, ~event, ~bl,
    "a", 1, FALSE, FALSE,
    "a", 2, FALSE, TRUE,
    "a", 3.5, TRUE, FALSE,
    "a", 4, TRUE, TRUE, 
    "b", 1, FALSE, FALSE,
    "b", 2, FALSE, FALSE,
    "b", 4, FALSE, TRUE
) %>% sample_frac(1)

test_that("multiple baselines", {
    expect_equal(
        get_progressor(
            f,
            event_col = "event",
            time_col = "time",
            group_col = "group",
            baseline_col = "bl"
        ),
        dplyr::tribble(
            ~group, ~time, ~event, ~bl, ~progressor,
            "a", 1, FALSE, FALSE, TRUE,
            "a", 2, FALSE, TRUE, TRUE,
            "a", 3.5, TRUE, FALSE, TRUE,
            "a", 4, TRUE, TRUE, TRUE,
            "b", 1, FALSE, FALSE, FALSE,
            "b", 2, FALSE, FALSE, FALSE,
            "b", 4, FALSE, TRUE, FALSE,
        )
    )
})

g <- dplyr::tribble(
    ~group, ~time, ~event,
    "a", 1, FALSE,
    "a", 2, FALSE,
    "a", 3.5, FALSE,
    "a", 4, FALSE,
    "a", 5.5, TRUE,
    "a", 6, TRUE,
    "b", 1, FALSE,
    "b", 2, TRUE,
    "b", 4, TRUE,
) %>% sample_frac(1)

test_that("exclude stable if event", {
    expect_equal(
        get_progressor(
            g,
            event_col = "event",
            time_col = "time",
            group_col = "group",
            time_window = 2
        ),
        dplyr::tribble(
            ~group, ~time, ~event, ~progressor,
            "a", 1, FALSE, NA,
            "a", 2, FALSE, NA,
            "a", 3.5, FALSE, NA,
            "a", 4, FALSE, NA,
            "a", 5.5, TRUE, NA,
            "a", 6, TRUE, NA,
            "b", 1, FALSE, TRUE,
            "b", 2, TRUE, TRUE,
            "b", 4, TRUE, TRUE,
        )
    )
    expect_equal(
        get_progressor(
            g,
            event_col = "event",
            time_col = "time",
            group_col = "group",
            time_window = 2,
            exclude_stable_event = FALSE,
        ),
        dplyr::tribble(
            ~group, ~time, ~event, ~progressor,
            "a", 1, FALSE, FALSE,
            "a", 2, FALSE, FALSE,
            "a", 3.5, FALSE, FALSE,
            "a", 4, FALSE, FALSE,
            "a", 5.5, TRUE, FALSE,
            "a", 6, TRUE, FALSE,
            "b", 1, FALSE, TRUE,
            "b", 2, TRUE, TRUE,
            "b", 4, TRUE, TRUE,
        )
    )
})


h <- dplyr::tribble(
    ~group, ~time, ~event, ~bl,
    "a", 1, FALSE, FALSE,
    "a", 2, FALSE, TRUE,
    "a", 3.5, TRUE, FALSE,
    "b", 1, TRUE, FALSE,
    "b", 2, TRUE, FALSE,
    "b", 4, TRUE, TRUE
) %>% sample_frac(1)

test_that("time to conversion", {
    expect_equal(
        get_progressor(
            h,
            event_col = "event",
            time_col = "time",
            group_col = "group",
            keep_extra_col = TRUE,
            time_window = 5
        ) %>%
        select(-c(.event, .revert, .bl)),
        dplyr::tribble(
            ~group, ~time, ~event, ~bl, ~progressor, ~.time_to_progression,
            "a", 1, FALSE, FALSE, TRUE, 2.5,
            "a", 2, FALSE, TRUE,  TRUE, 2.5,
            "a", 3.5, TRUE, FALSE,  TRUE, 2.5,
            "b", 1, TRUE, FALSE, NA, 0.0,
            "b", 2, TRUE, FALSE, NA, 0.0,
            "b", 4, TRUE, TRUE, NA, 0.0,
        )
    )
    expect_equal(
        get_progressor(
            h,
            event_col = "event",
            time_col = "time",
            group_col = "group",
            baseline_col = "bl",
            keep_extra_col = TRUE,
            time_window = 5
        ) %>%
        select(-c(.event, .revert, .bl)),
        dplyr::tribble(
            ~group, ~time, ~event, ~bl, ~progressor, ~.time_to_progression,
            "a", 1, FALSE, FALSE, TRUE, 1.5,
            "a", 2, FALSE, TRUE,  TRUE, 1.5,
            "a", 3.5, TRUE, FALSE,  TRUE, 1.5,
            "b", 1, TRUE, FALSE, NA, -3.0,
            "b", 2, TRUE, FALSE, NA, -3.0,
            "b", 4, TRUE, TRUE, NA, -3.0,
        )
    )
})