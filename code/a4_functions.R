
# functions from "a4_learn_data_primer.Rmd"
get_levels <- function(x){
  as.numeric(unlist(lapply(strsplit(unlist(strsplit(subset(
    datadic, FIELD_NAME==x)$FIELD_CODE, ';')), '='), function(y) y[1])))
}
get_labels <- function(x){
  unlist(lapply(strsplit(unlist(strsplit(subset(
    datadic, FIELD_NAME==x)$FIELD_CODE, ';')), '='), function(y) y[2]))
}
a4_convert_to_factor <- function(x, name) {
    factor(x, levels = get_levels(name), labels = get_labels(name))
}

a4_remove_rescreens <- function(.data, subjinfo_df = NULL) {

    # remove BIDs which were rescreened at a later date; this
    # removes the first BID and only keeps the second BID

    if (is.null(subjinfo_df)) {
        subjinfo_df <- subjinfo
    }

    # get BID of rescreens (pulled from `a4_learn_data_primer.Rmd`)
    rescreens <- subjinfo_df %>%
        filter(!is.na(PREVBID)) %>%
        rename(BID1 = PREVBID, BID2 = BID) %>%
        select(BID1, BID2)

    # remove any BID that appears in BID1 (only if BID2 also exists in the frame)
    subj_list <- .data %>% pull(BID) %>% unique
    rescreens_remove <- rescreens %>%
        filter(BID2 %in% subj_list) %>%
        pull(BID1)
    .data %>%
        filter(!(BID %in% rescreens_remove))

}

a4_merge_subjinfo <- function(.data, cols = NULL, subjinfo_df = NULL) {

    if (is.null(subjinfo_df)) {
        subjinfo_df <- subjinfo
    }

    # merge everything
    if (is.null(cols)) {
        .data %>%
            left_join(subjinfo_df, by = c("SUBSTUDY", "BID"))
    } else { 
        .data %>% left_join(
            subjinfo_df %>% select(SUBSTUDY, BID, all_of(cols)),
            by = c("SUBSTUDY", "BID")
        )
    }

}

a4_get_age <- function(.data, days_from_consent = NULL, sv_df = NULL, subjinfo_df = NULL, from_col_only = FALSE) {

    # TODO: document function
    # NOTE: if BID/VISCODE has non-NA age from both methods, the age from the days-from-consent method is always chosen over the SV.csv age

    if (is.null(subjinfo_df)) {
        subjinfo_df <- subjinfo
    }
    if (is.null(sv_df)) {
        sv_df <- sv
    }

    from_col <- function(.data, days_from_consent) {
        # compute age using a days-from-consent column and SUBJINFO.csv
        .data %>%
            left_join(subjinfo_df %>% select(BID, AGEYR) %>% rename(.age_at_consent = AGEYR), by = 'BID') %>%
            mutate(
                age.col = .age_at_consent + (.data[[days_from_consent]]) / 365.25
            ) %>%
            select(-.age_at_consent)
    }
    from_sv <- function(.data) {
        # use SV.csv to compute age, then merge with .data
        age_df <- sv_df %>%
            select(BID, VISITCD, SVSTDTC_DAYS_CONSENT) %>%
            left_join(subjinfo_df %>% select(BID, AGEYR), by = "BID") %>%
            mutate(age.sv = AGEYR + SVSTDTC_DAYS_CONSENT/365.25) %>%
            select(BID, VISITCD, age.sv)
        .data %>%
            left_join(
                age_df,
                by = c("BID" = "BID", "VISCODE" = "VISITCD")
            )
    }

    
    if (from_col_only) {
    # use days-from-consent method only

        .data %>%
            from_col(days_from_consent) %>%
            rename(age = age.col)
        
    } else if (is.null(days_from_consent)) {
    # use SV.csv method only

        .data %>%
            from_sv() %>%
            rename(age = age.sv)

    } else {
    # compute age using both methods, then combine results

        .data <- .data %>%
            from_col(days_from_consent) %>% from_sv()
        
        # check if both ages are equal (only from non-na group)
        ages_are_equal <- .data %>% mutate(.equal = age.col == age.sv) %>% pull(.equal) %>% all(., na.rm = TRUE)
        if (!ages_are_equal) {
            warning("WARNING: both methods of computing age were performed, but different ages were computed across the two methods for certain BID/VISCODE's")
        }

        .data %>%
            mutate(age = case_when(
                !is.na(age.col) ~ age.col,
                is.na(age.col) & !is.na(age.sv) ~ age.sv,
                .default = NA
            ))

    }
    
}

subjinfo <- read_csv("/ceph/chpc/shared/aristeidis_sotiras_group/aris_data/A4/rawdata/Clinical/Derived Data/SUBJINFO.csv", show_col_types = FALSE)
sv <- read_csv("/ceph/chpc/shared/aristeidis_sotiras_group/aris_data/A4/rawdata/Clinical/Derived Data/SV.csv", show_col_types = FALSE)
datadic <- read_csv("/ceph/chpc/shared/aristeidis_sotiras_group/aris_data/A4/rawdata/Clinical/Documents/Data Dictionaries/clinical_datadic.csv")
