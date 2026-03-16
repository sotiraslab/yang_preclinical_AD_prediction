# TODO: make unit test to test that function works in the presense of NA values

# =======================
# ===== DEFINE DATA =====
# =======================

# define reference tibble (combined)
a <- dplyr::tribble(
    ~group, ~date1, ~x1, ~date2, ~x2, ~y2,
    "a", "2000-01-01", 1, "2000-01-15", 2, 3,
    "a", "2000-02-01", 4, "2000-03-05", 5, 6,
    "b", "2000-03-01", 7, "2000-04-02", 8, 9,
    "b", "2000-04-01", 10, "2000-03-02", 11, 12,
    "c", "2000-05-01", 13, "2001-10-18", 14, 15,
) %>%
    dplyr::mutate(
        date1 = as.Date(date1),
        date2 = as.Date(date2)
    )

# define reference tibble (separated)
a_ref <- dplyr::tribble(
    ~group, ~date1, ~x1,
    "a", "2000-01-01", 1,
    "a", "2000-02-01", 4,
    "b", "2000-03-01", 7,
    "b", "2000-04-01", 10,
    "c", "2000-05-01", 13,
) %>% dplyr::mutate(date1 = as.Date(date1))
a_match <- dplyr::tribble(
    ~group, ~date2, ~x2, ~y2,
    "a", "2000-01-15", 2, 3,
    "a", "2000-03-05", 5, 6,
    "b", "2000-04-02", 8, 9,
    "b", "2000-03-02", 11, 12,
    "c", "2001-10-18", 14, 15,
) %>% dplyr::mutate(date2 = as.Date(date2))
a_match2 <- dplyr::tribble(
    ~group, ~date2, ~x2, ~y2,
    "a", "2000-01-15", 2, 3,
    "a", "2000-03-05", 5, 6,
    "b", "2000-04-02", 8, 9,
    "b", "2000-03-02", 11, 12,
    "b", "2000-05-02", 100, 200,
    "c", "2001-10-18", 14, 15,
    "c", "2002-01-18", 300, 400,
) %>% dplyr::mutate(date2 = as.Date(date2))
a_match3 <- dplyr::tribble(
    ~group, ~date2, ~x2, ~y2,
    "a", "2003-01-15", 2, 3,
    "a", "2003-03-05", 5, 6,
    "b", "2003-04-02", 8, 9,
    "b", "2003-03-02", 11, 12,
    "b", "2003-05-02", 100, 200,
    "c", "2001-10-18", 14, 15,
    "c", "2002-01-18", 300, 400,
) %>% dplyr::mutate(date2 = as.Date(date2))

# ======================
# ===== UNIT TESTS =====
# ======================


b_merge <- dplyr::tribble(
    ~group, ~date1, ~x1, ~date2, ~x2, ~y2, ~x2.match, ~y2.match, ~date2.match,
    "a", "2000-01-01", 1, "2000-01-15", 2, 3, 2, 3, "2000-01-15",
    "a", "2000-02-01", 4, "2000-03-05", 5, 6, 2, 3, "2000-01-15",
    "b", "2000-03-01", 7, "2000-04-02", 8, 9, 11, 12, "2000-03-02",
    "b", "2000-04-01", 10, "2000-03-02", 11, 12, 8, 9, "2000-04-02",
    "c", "2000-05-01", 13, "2001-10-18", 14, 15, NA, NA, NA,
) %>%
    dplyr::mutate(
        date1 = as.Date(date1),
        date2 = as.Date(date2),
        date2.match = as.Date(date2.match)
    )
test_that("match with merged df is correct", {
    expect_equal(
        match_data(
            a,
            ref_col = "date1",
            match_col = "date2",
            select_col = c("x2", "y2"),
            group_col = "group",
            max_diff = 365.25
        ),
        b_merge
    )
})


b_sep <- dplyr::tribble(
    ~group, ~date1, ~x1, ~x2.match, ~y2.match, ~date2.match,
    "a", "2000-01-01", 1, 2, 3, "2000-01-15",
    "a", "2000-02-01", 4, 2, 3, "2000-01-15",
    "b", "2000-03-01", 7, 11, 12, "2000-03-02",
    "b", "2000-04-01", 10, 8, 9, "2000-04-02",
    "c", "2000-05-01", 13, NA, NA, NA,
) %>%
    dplyr::mutate(
        date1 = as.Date(date1),
        date2.match = as.Date(date2.match)
    )
test_that("match with separate dfs is correct", {
    expect_equal(
        match_data(
            a_ref,
            match_df = a_match,
            ref_col = "date1",
            match_col = "date2",
            select_col = c("x2", "y2"),
            group_col = "group",
            max_diff = 365.25
        ),
        b_sep
    )
    expect_equal(
        match_data(
            a_ref,
            match_df = a_match2,
            ref_col = "date1",
            match_col = "date2",
            select_col = c("x2", "y2"),
            group_col = "group",
            max_diff = 365.25
        ),
        b_sep
    )
})


b_merge_with_idx <- dplyr::tribble(
    ~group, ~date1, ~x1, ~date2, ~x2, ~y2, ~x2.match, ~y2.match, ~date2.match, ~idx.match,
    "a", "2000-01-01", 1, "2000-01-15", 2, 3, 2, 3, "2000-01-15", 1,
    "a", "2000-02-01", 4, "2000-03-05", 5, 6, 2, 3, "2000-01-15", 1,
    "b", "2000-03-01", 7, "2000-04-02", 8, 9, 11, 12, "2000-03-02", 4,
    "b", "2000-04-01", 10, "2000-03-02", 11, 12, 8, 9, "2000-04-02", 3,
    "c", "2000-05-01", 13, "2001-10-18", 14, 15, NA, NA, NA, NA
) %>%
    dplyr::mutate(
        date1 = as.Date(date1),
        date2 = as.Date(date2),
        date2.match = as.Date(date2.match)
    )
test_that("index column appears", {
    expect_equal(
        match_data(
            a,
            ref_col = "date1",
            match_col = "date2",
            select_col = c("x2", "y2"),
            group_col = "group",
            max_diff = 365.25,
            idx_col = "idx.match"
        ),
        b_merge_with_idx
    )
})


b_merge_with_idx_diff <- dplyr::tribble(
    ~group, ~date1, ~x1, ~date2, ~x2, ~y2, ~x2.match, ~y2.match, ~date2.match, ~idx.match,
    "a", "2000-01-01", 1, "2000-01-15", 2, 3, 2, 3, "2000-01-15", 1,
    "a", "2000-02-01", 4, "2000-03-05", 5, 6, 2, 3, "2000-01-15", 1,
    "b", "2000-03-01", 7, "2000-04-02", 8, 9, 11, 12, "2000-03-02", 4,
    "b", "2000-04-01", 10, "2000-03-02", 11, 12, 8, 9, "2000-04-02", 3,
    "c", "2000-05-01", 13, "2001-10-18", 14, 15, NA, NA, NA, NA
) %>%
    dplyr::mutate(
        date1 = as.Date(date1),
        date2 = as.Date(date2),
        date2.match = as.Date(date2.match),
        diff.match = abs(difftime(date1, date2.match, units = "days"))
    )
test_that("diff column appears", {
    expect_equal(
        match_data(
            a,
            ref_col = "date1",
            match_col = "date2",
            select_col = c("x2", "y2"),
            group_col = "group",
            max_diff = 365.25,
            idx_col = "idx.match",
            diff_col = "diff.match"
        ),
        b_merge_with_idx_diff
    )
})


b_no_max_diff <- dplyr::tribble(
    ~group, ~date1, ~x1, ~x2.match, ~y2.match, ~date2.match,
    "a", "2000-01-01", 1, 2, 3, "2000-01-15",
    "a", "2000-02-01", 4, 2, 3, "2000-01-15",
    "b", "2000-03-01", 7, 11, 12, "2000-03-02",
    "b", "2000-04-01", 10, 8, 9, "2000-04-02",
    "c", "2000-05-01", 13, 14, 15, "2001-10-18",
) %>%
    dplyr::mutate(
        date1 = as.Date(date1),
        date2.match = as.Date(date2.match)
    )
test_that("no max diff constraint", {
    expect_equal(
        match_data(
            a_ref,
            match_df = a_match,
            ref_col = "date1",
            match_col = "date2",
            select_col = c("x2", "y2"),
            group_col = "group"
        ),
        b_no_max_diff
    )
})


b_new_suffix <- dplyr::tribble(
    ~group, ~date1, ~x1, ~date2, ~x2, ~y2, ~x2.new, ~y2.new, ~date2.new, ~idx.match,
    "a", "2000-01-01", 1, "2000-01-15", 2, 3, 2, 3, "2000-01-15", 1,
    "a", "2000-02-01", 4, "2000-03-05", 5, 6, 2, 3, "2000-01-15", 1,
    "b", "2000-03-01", 7, "2000-04-02", 8, 9, 11, 12, "2000-03-02", 4,
    "b", "2000-04-01", 10, "2000-03-02", 11, 12, 8, 9, "2000-04-02", 3,
    "c", "2000-05-01", 13, "2001-10-18", 14, 15, NA, NA, NA, NA
) %>%
    dplyr::mutate(
        date1 = as.Date(date1),
        date2 = as.Date(date2),
        date2.new = as.Date(date2.new),
        idx.match = as.integer(idx.match)
    )
test_that("suffix is correct", {
    expect_equal(
        match_data(
            a,
            ref_col = "date1",
            match_col = "date2",
            select_col = c("x2", "y2"),
            group_col = "group",
            max_diff = 365.25,
            idx_col = "idx.match",
            suffix = ".new"
        ),
        b_new_suffix
    )
})


b_no_match <- dplyr::tribble(
    ~group, ~date1, ~x1, ~x2.match, ~y2.match, ~date2.match, ~idx.match,
    "a", "2000-01-01", 1, NA, NA, NA, NA,
    "a", "2000-02-01", 4, NA, NA, NA, NA,
    "b", "2000-03-01", 7, NA, NA, NA, NA,
    "b", "2000-04-01", 10, NA, NA, NA, NA,
    "c", "2000-05-01", 13, NA, NA, NA, NA,
) %>%
    dplyr::mutate(
        date1 = as.Date(date1),
        date2.match = as.Date(date2.match),
        diff = difftime(date1, date2.match, units = "days"),
        across(c(x2.match, y2.match), as.numeric),
        idx.match = as.integer(idx.match)
    )
test_that("no match", {
    expect_equal(
        match_data(
            a_ref,
            match_df = a_match3,
            ref_col = "date1",
            match_col = "date2",
            select_col = c("x2", "y2"),
            group_col = "group",
            max_diff = 365.25,
            idx_col = "idx.match",
            diff_col = "diff"
        ),
        b_no_match
    )
})
