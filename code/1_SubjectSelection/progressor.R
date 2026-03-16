# ==============================================================================

# Name: progressor.R
# Author: Braden Yang
# Created: 02/27/2024
# Description: Get stables and progressors

# ==============================================================================

rm(list = ls())

# *** toggle variables ***
INTERACTIVE <- FALSE
SAVE_CSV <- TRUE
SAVE_FIG <- FALSE
CONTROL_SUBJ <- FALSE
# *** toggle variables ***

# ===========================
# ===== IMPORT PACKAGES =====
# ===========================

library(tidyverse)
library(optparse)
library(viridisLite)
library(ggnewscale)
library(patchwork)

library(devtools)

# ===========================
# ===== PARSE ARGUMENTS =====
# ===========================

option_list <- list(
    make_option(c("-w", "--wdir"), action="store", default=NULL,
        type="character", help="Path to project directory (if none, uses current working directory)"),
    make_option(c("-o", "--odir"), action="store", default="data/subj",
        type="character", help="Path to output directory"),
    make_option(c("-f", "--odir_fig"), action="store", default="figures/subj",
        type="character", help="Output directory for subject list figures")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (INTERACTIVE) {
    opt$wdir <- "/home/b.y.yang/sotiraslab/BradenADLongitudinalPrediction"
    opt$odir <- "/home/b.y.yang/sotiraslab/BradenADLongitudinalPrediction/data/subj"
    opt$odir_fig <- "/home/b.y.yang/sotiraslab/BradenADLongitudinalPrediction/figures/subj"
}

# set working directory
if (!is.null(opt$wdir)) {
    setwd(opt$wdir)
}

# load byyfunctions
load_all("submodules/byyfunctions")

# make output directory
byyfunctions::make_dir(opt$odir)

# ============================
# ===== DEFINE FUNCTIONS =====
# ============================

# source A4 functions
source("code/a4_functions.R")
source("code/functions.R")

get_all_baselines <- function(
    .data,
    clinical_df,
    time_col = "age",
    event_col = "diagnosis",
    group_col = "subj",
    event_label = c("MCI", "Dementia"),
    normal_event_label = c("Normal", "Impaired not MCI"),
    baseline_diagnosis_method = "closest",
    baseline_max_time_diff_match_cdr = 1,
    return_full = FALSE,
    exclude_revert_bool = TRUE,
    recover_reverter = FALSE
) {

    # TODO: indicate warnings (e..g subj == 2229, missing args to min, returning Inf)

    # *function was pulled from progressor.R in BradenPAC
    # use amyloid-positive baseline to find MCI/AD progressors;
    # this function will first select stable amyloid+ subjects,
    # then use each of their amyloid-positive scans to check
    # whether they are valid baselines (i.e. starts at normal
    # clinical diagnosis); finally, the time-to-progression is
    # computed


    match_single_group <- function(df, time, time_col = "age", method = "closest") {
        
        # compute diff
        diff_df <- df %>% mutate(.diff = .data[[time_col]] - time)
        
        # extract row relative to time (depending on method)
        if (method == "previous") {
            r <- diff_df %>% filter(.diff < 0) %>% slice_max(.diff)
        } else if (method == "next") {
            r <- diff_df %>% filter(.diff > 0) %>% slice_min(.diff)
        } else if (method == "both") {
            r1 <- diff_df %>% filter(.diff < 0) %>% slice_max(.diff)
            r2 <- diff_df %>% filter(.diff > 0) %>% slice_min(.diff)
            r <- bind_rows(r1, r2)
        } else if (method == "closest") {
            r <- diff_df %>% slice_min(abs(.diff))
        } else {
            stop("invalid method")
        }
        
        return(r)
        
    }

    is_valid_baseline <- function(
        df,
        time,
        time_col,
        event_col,
        normal_event_label,
        method = "closest",
        max_time_diff = 1
    ) {
        
        # check if a scan is a valid baseline (i.e. starts as normal clinical diagnosis)
        # - note that one can choose what counts as a "cognitively normal scan" from 3 options:
        #   1. the closest entry = normal
        #   2. only the previous entry = normal
        #   3. both previous and next entry = normal
        # if closest clinical diagnosis is past the max_time_diff, then return False
        
        if (is.null(df)) {
            return(FALSE)
        }
        
        r <- match_single_group(df, time, time_col, method)
        return(all(r[event_col] %in% normal_event_label) & all(abs(r[".diff"]) <= max_time_diff))
        
    }

    get_time_to_conversion <- function(
        df,
        time,
        time_col,
        event_col,
        event_label
    ) {
        
        if (is.null(df)) {
            return(NA)
        }
        
        foo <- df %>%
            filter(.data[[event_col]] %in% event_label) %>%
            slice_min(.data[[time_col]], with_ties = FALSE) %>%
            mutate(.diff = .data[[time_col]] - time)
        
        if (length(foo[[".diff"]]) == 0) {
            return(NA)
        } else {
            return(foo[[".diff"]])
        }
        
    }

    get_max_stable <- function(
        df,
        time,
        time_col,
        event_col,
        normal_event_label
    ) {
        
        if (is.null(df)) {
            return(NA)
        }
        
        foo <- df %>%
            filter(.data[[event_col]] %in% normal_event_label) %>%
            slice_max(.data[[time_col]], with_ties = FALSE) %>%
            mutate(.diff = .data[[time_col]] - time)
        
        if (length(foo[[".diff"]]) == 0) {
            return(NA)
        } else {
            return(foo[[".diff"]])
        }
        
    }

    # first select all valid progressors & stables (+ excluded subjects)
    valid_progressor_df <- clinical_df %>%
        byyfunctions::get_progressor(
            event_col,
            time_col,
            group_col,
            event_label,
            time_window = Inf,
            progressor_nonevent_start = FALSE,
            exclude_revert = exclude_revert_bool,
            keep_extra_col = TRUE
        )
    progressor_subj <- valid_progressor_df %>%
        filter(progressor) %>%
        pull(all_of(group_col)) %>%
        unique()
    stable_subj <- valid_progressor_df %>%
        filter(!progressor) %>%
        pull(all_of(group_col)) %>%
        unique()
    exclude_subj <- valid_progressor_df %>%
        filter(is.na(progressor)) %>%
        pull(all_of(group_col)) %>%
        unique()

    # select all scans that are valid baselines and compute
    # time-to-conversion/max-stable-time
    df <- .data %>%
        left_join(
            clinical_df %>% nest(.by = group_col),
            by = group_col
        ) %>%
        mutate(
            .is_valid_baseline = map2_lgl(  # indicates whether scan starts at cognitively-normal (if so, is valid baseline)
                data, .data[[time_col]],
                ~ is_valid_baseline(.x, .y, time_col, event_col, normal_event_label, baseline_diagnosis_method, baseline_max_time_diff_match_cdr)
            ),
            .time_to_progression = map2_dbl(  # time to conversion to MCI or AD, in years
                data, .data[[time_col]],
                ~ get_time_to_conversion(.x, .y, time_col, event_col, event_label)
            ),
            .time_stable = map2_dbl(  # time remained stable
                data, .data[[time_col]],
                ~ get_max_stable(.x, .y, time_col, event_col, normal_event_label)
            ),
            .progressor = case_when(
                .data[[group_col]] %in% progressor_subj ~ "progressor",
                .data[[group_col]] %in% stable_subj ~ "stable",
                .data[[group_col]] %in% exclude_subj ~ "exclude",
                .default = NA
            )
        )

    # recover some progressors who may have reverted to CDR = 0
    if (recover_reverter) {
        recover_excluded <- function(a, x) {
            all(x[x[time_col] <= a, event_col] == 0)
        }
        df <- df %>%
            mutate(
                .recover_excluded = map2_lgl(.data[[time_col]], data, recover_excluded),
                .progressor = if_else(
                    .progressor == "exclude" & .recover_excluded,
                    "progressor_recovered",
                    .progressor
                )
            )
    }

    return(df)

}

spaghetti <- function(
    cdr_df,
    pet_df,
    group_col = "subj",
    time_col = "age",
    event_col = "cdr",
    tracer_col = "tracer",
    event_levels = c(0, 0.5, 1, 2, 3)
) {

    viridis_cdr <- viridis(5)
    names(viridis_cdr) <- event_levels

    time_min_df <- pet_df %>%
        group_by(pick(all_of(group_col))) %>%
        summarise(.time_min = min(.data[[time_col]]))

    p_list <- list()
    for (t in pet_df %>% pull(.data[[tracer_col]] %>% unique())) {

        pet_tracer_df <- pet_df %>% filter(.data[[tracer_col]] == t)
        n_group <- length(unique(pet_tracer_df %>% pull(.data[[group_col]])))

        p <- cdr_df %>%
            filter(
                .data[[group_col]] %in% (pet_tracer_df %>% pull(.data[[group_col]]) %>% unique())
            ) %>%
            left_join(time_min_df, by = group_col) %>%
            mutate(
                .event = factor(.data[[event_col]], levels = event_levels, ordered = TRUE),
                .group = fct_reorder(.data[[group_col]], .time_min),
                .time = .data[[time_col]]
            ) %>%
            ggplot(aes(x = .group, y = .time, group = .group)) + 
                geom_point(aes(color = .event)) +
                geom_line(aes(color = .event)) +
                scale_color_manual(values = viridis_cdr) +
                new_scale_color() +
                geom_point(
                    data = pet_tracer_df %>%
                        mutate(.group = .data[[group_col]], .time = .data[[time_col]], Tracer = factor(.data[[tracer_col]], levels = pet_df %>% pull(.data[[tracer_col]]) %>% unique())),
                    mapping = aes(color = Tracer),
                    size = 4, shape = 1
                ) +
                labs(x = "subject", y = "age", title = str_glue("{t} (n={n_group})")) +
                theme(axis.text.x = element_text(angle = 90)) +
                scale_y_reverse()
        p_list[[t]] <- p

    }

    wrap_plots(p_list)

}

get_control_subj <- function(pet_df, cdr_df, amypos_col = "amyloid_positive", dx_col = "cdr", cn_val = 0) {

    # criteria for control subjects:
    # - amyloid-negative in *all* PET scans
    # - CDR=0 (or equivalent diagnosis of CN) in *all* entries

    amyneg_subj <- pet_df %>%
        group_by(subj) %>%
        drop_na(all_of(amypos_col)) %>%
        filter(all(!.data[[amypos_col]])) %>%
        pull(subj) %>% unique()

    cdr_df %>%
        filter(subj %in% amyneg_subj) %>%
        drop_na(all_of(dx_col)) %>%
        group_by(subj) %>%
        filter(all(.data[[dx_col]] == cn_val)) %>%
        select(subj) %>% unique
    
}

# =====================
# ===== LOAD DATA =====
# =====================

adni_amyloid <- read_rds("data/tidy/adni/adni_amyloid.RDS")
adni_cdr <- read_rds("data/tidy/adni/adni_cdr.RDS")

oasis_amyloid <- read_rds("data/tidy/oasis/oasis_amyloid.RDS")
oasis_cdr <- read_rds("data/tidy/oasis/oasis_cdr.RDS")

pac_amyloid <- read_rds("data/tidy/pac/pac_amyloid.RDS")
pac_diagnosis <- read_rds("data/tidy/pac/pac_diagnosis.RDS")

a4_amyloid <- read_rds("data/tidy/a4/a4_amyloid.RDS")
a4_cdr <- read_rds("data/tidy/a4/a4_cdr.RDS")

habs_amyloid <- read_rds("data/tidy/habs/habs_amyloid.RDS")
habs_cdr <- read_rds("data/tidy/habs/habs_cdr.RDS")

# habshd_amyloid <- read_rds("data/tidy/habshd/habshd_amyloid.RDS")
# habshd_cdr <- read_rds("data/tidy/habshd/habshd_cdr.RDS")

mcsa_amyloid <- read_rds("data/tidy/mcsa/mcsa_amyloid.RDS")
mcsa_cdr <- read_rds("data/tidy/mcsa/mcsa_cdr.RDS")

# =======================================
# ===== GET STABLES AND PROGRESSORS =====
# =======================================

if (file.exists(file.path(opt$odir, "adni.RDS"))) {
    adni_progressor <- read_rds(file.path(opt$odir, "adni.RDS"))
} else {
    print("+++ Finding stables and progressors for ADNI +++")
    adni_progressor <- get_all_baselines(
        adni_amyloid,
        adni_cdr,
        time_col = "age",
        event_col = "cdr",
        group_col = "subj",
        event_label = c(0.5, 1, 2, 3),
        normal_event_label = 0,
        baseline_diagnosis_method = "closest",
        return_full = TRUE,
        exclude_revert_bool = TRUE,
        recover_reverter = TRUE
    )
    if (SAVE_CSV) {
        write_rds(adni_progressor, file.path(opt$odir, "adni.RDS"))
        write_csv(adni_progressor, file.path(opt$odir, "adni.csv"))
    }
}

if (file.exists(file.path(opt$odir, "oasis.RDS"))) {
    oasis_progressor <- read_rds(file.path(opt$odir, "oasis.RDS"))
} else {
    print("+++ Finding stables and progressors for OASIS +++")
    oasis_progressor <- get_all_baselines(
        oasis_amyloid,
        oasis_cdr,
        time_col = "age",
        event_col = "cdr",
        group_col = "subj",
        event_label = c(0.5, 1, 2, 3),
        normal_event_label = 0,
        baseline_diagnosis_method = "closest",
        return_full = TRUE,
        exclude_revert_bool = TRUE,
        recover_reverter = TRUE
    )
    if (SAVE_CSV) {
        write_rds(oasis_progressor, file.path(opt$odir, "oasis.RDS"))
        write_csv(oasis_progressor, file.path(opt$odir, "oasis.csv"))
    }
}

if (file.exists(file.path(opt$odir, "pac.RDS"))) {
    pac_progressor <- read_rds(file.path(opt$odir, "pac.RDS"))
} else {
    print("+++ Finding stables and progressors for PAC +++")
    pac_progressor <- get_all_baselines(
        pac_amyloid,
        pac_diagnosis %>% drop_na(diagnosis),
        time_col = "age",
        event_col = "diagnosis",
        group_col = "subj",
        event_label = c("MCI", "Dementia"),
        normal_event_label = "Normal",
        baseline_diagnosis_method = "closest",
        return_full = TRUE,
        exclude_revert_bool = TRUE,
        recover_reverter = TRUE
    )
    if (SAVE_CSV) {
        write_rds(pac_progressor, file.path(opt$odir, "pac.RDS"))
        write_csv(pac_progressor, file.path(opt$odir, "pac.csv"))
    }
}

if (file.exists(file.path(opt$odir, "a4.RDS"))) {
    a4_progressor <- read_rds(file.path(opt$odir, "a4.RDS"))
} else {
    print("+++ Finding stables and progressors for A4 +++")
    a4_progressor <- get_all_baselines(
        a4_amyloid,
        a4_cdr %>% drop_na(cdr),
        time_col = "age",
        event_col = "cdr",
        group_col = "subj",
        event_label = c(0.5, 1, 2, 3),
        normal_event_label = 0,
        baseline_diagnosis_method = "closest",
        return_full = TRUE,
        exclude_revert_bool = TRUE,
        recover_reverter = TRUE
    )
    if (SAVE_CSV) {
        write_rds(a4_progressor, file.path(opt$odir, "a4.RDS"))
        write_csv(a4_progressor, file.path(opt$odir, "a4.csv"))
    }
}

if (file.exists(file.path(opt$odir, "habs.RDS"))) {
    habs_progressor <- read_rds(file.path(opt$odir, "habs.RDS"))
} else {
    print("+++ Finding stables and progressors for HABS +++")
    habs_progressor <- get_all_baselines(
        habs_amyloid,
        habs_cdr %>% drop_na(cdr),
        time_col = "age",
        event_col = "cdr",
        group_col = "subj",
        event_label = c(0.5, 1, 2, 3),
        normal_event_label = 0,
        baseline_diagnosis_method = "closest",
        return_full = TRUE,
        exclude_revert_bool = TRUE,
        recover_reverter = TRUE
    )
    if (SAVE_CSV) {
        write_rds(habs_progressor, file.path(opt$odir, "habs.RDS"))
        write_csv(habs_progressor, file.path(opt$odir, "habs.csv"))
    }
}

# if (file.exists(file.path(opt$odir, "habshd.RDS"))) {
#     habshd_progressor <- read_rds(file.path(opt$odir, "habshd.RDS"))
# } else {
#     print("+++ Finding stables and progressors for HABS-HD +++")
#     habshd_progressor <- get_all_baselines(
#         habshd_amyloid,
#         habshd_cdr %>% drop_na(cdr),
#         time_col = "age",
#         event_col = "cdr",
#         group_col = "subj",
#         event_label = c(0.5, 1, 2, 3),
#         normal_event_label = 0,
#         baseline_diagnosis_method = "closest",
#         return_full = TRUE,
#         exclude_revert_bool = TRUE,
#         recover_reverter = TRUE
#     )
#     if (SAVE_CSV) {
#         write_rds(habshd_progressor, file.path(opt$odir, "habshd.RDS"))
#         write_csv(habshd_progressor, file.path(opt$odir, "habshd.csv"))
#     }
# }

if (file.exists(file.path(opt$odir, "mcsa.RDS"))) {
    mcsa_progressor <- read_rds(file.path(opt$odir, "mcsa.RDS"))
} else {
    print("+++ Finding stables and progressors for MCSA +++")
    mcsa_progressor <- get_all_baselines(
        mcsa_amyloid,
        mcsa_cdr,
        time_col = "age",
        event_col = "cdr",
        group_col = "subj",
        event_label = c(0.5, 1, 2, 3),
        normal_event_label = 0,
        baseline_diagnosis_method = "closest",
        return_full = TRUE,
        exclude_revert_bool = TRUE,
        recover_reverter = TRUE
    )
    if (SAVE_CSV) {
        write_rds(mcsa_progressor, file.path(opt$odir, "mcsa.RDS"))
        write_csv(mcsa_progressor, file.path(opt$odir, "mcsa.csv"))
    }
}

# ==============================
# ===== GET CONTROL COHORT =====
# ==============================

if (CONTROL_SUBJ) {

    print("+++ Getting control subjects +++")

    control_subj <- list()
    control_subj[["adni"]] <- get_control_subj(adni_amyloid, adni_cdr)
    control_subj[["oasis"]] <- get_control_subj(oasis_amyloid, oasis_cdr)
    control_subj[["a4"]] <- get_control_subj(a4_amyloid, a4_cdr, amypos_col = "amyloid_positive.a4_thresh")
    control_subj[["pac"]] <- get_control_subj(pac_amyloid, pac_diagnosis, amypos_col = "amyloid_positive.site", dx_col = "diagnosis", cn_val = "Normal")
    control_subj[["habs"]] <- get_control_subj(habs_amyloid, habs_cdr)
    # control_subj[["habshd"]] <- get_control_subj(habshd_amyloid, habshd_cdr, amypos_col = "amyloid_positive.fbb")
    control_subj[["mcsa"]] <- get_control_subj(mcsa_amyloid, mcsa_cdr)

    if (SAVE_CSV) {
        odir_control <- file.path(opt$odir, "control")
        byyfunctions::make_dir(odir_control)
        for (n in names(control_subj)) {
            write_csv(control_subj[[n]], file.path(odir_control, str_glue("{n}.csv")), col_names = FALSE)
        }
    }

}

# ===================================
# ===== GET COMBINED DATAFRAMES =====
# ===================================

print("+++ Merging into multisite dataframe +++")

multisite_df <- bind_rows(adni_progressor, oasis_progressor, pac_progressor, a4_progressor, habs_progressor, mcsa_progressor) %>%
    mutate(
        site_in_pac = if_else(is.na(site_in_pac), FALSE, site_in_pac),
        amyloid_positive = case_when(  # combine amyloid-positive column to a single column
            site == "A4" ~ amyloid_positive.a4_thresh,  # A4 -> use threshold set by A4
            site_in_pac ~ amyloid_positive.site,  # PAC -> use site-specific thresholds
            site == "HABS-HD" ~ amyloid_positive.fbb,  # HABS-HD -> use PET quantitative positivity
            .default = amyloid_positive
        ),
        age.round = round(age, 2) %>% format(nsmall = 2, trim = TRUE),
        age.raw_t1.round = round(age.raw_t1, 2) %>% format(nsmall = 2, trim = TRUE),
        site.subj.age = str_c(site, subj, age.round, sep = "__"),
        site.subj.age.raw_t1 = str_c(site, subj, age.raw_t1.round, sep = "__"),
        idvi = str_c(site, subj, age.round, tracer, sep = "_"),
        idvi.raw_t1 = str_c(site, subj, age.raw_t1.round, "T1", sep = "_"),
    )

# split into general stables and progressors (no time window); this is for deep learning
multisite_progressor_dl <- multisite_df %>% filter_progressor("amyloid_positive", max_time_to_progression = Inf)
multisite_stable_dl <- multisite_df %>% filter_stable("amyloid_positive", min_time_stable = 0)
multisite_progressor_stable <- bind_rows(multisite_progressor_dl, multisite_stable_dl)

if (SAVE_CSV) {

    write_rds(multisite_df, file.path(opt$odir, "multisite.RDS"))
    write_csv(multisite_df, file.path(opt$odir, "multisite.csv"))

    write_rds(multisite_progressor_stable, file.path(opt$odir, "multisite_progressor_stable.RDS"))
    write_csv(multisite_progressor_stable, file.path(opt$odir, "multisite_progressor_stable.csv"))

}

# =====================
# ===== VISUALIZE =====
# =====================

visualize_stable_progressor <- function(progressor_df, cdr_df, odir, prefix, event_col = "cdr", event_levels = c(0, 0.5, 1, 2, 3), recovered_bool = FALSE) {

    p <- progressor_df %>% filter_progressor("amyloid_positive", recovered = recovered_bool, max_time_to_progression = Inf)
    s <- progressor_df %>% filter_stable("amyloid_positive", min_time_stable = 0)

    gp <- spaghetti(cdr_df, p, "subj", "age", event_col, "tracer", event_levels)
    gs <- spaghetti(cdr_df, s, "subj", "age", event_col, "tracer", event_levels)

    byyfunctions::make_dir(odir)
    ggsave(file.path(odir, str_glue("{prefix}_progressor.png")), gp, width = 15, height = 5)
    ggsave(file.path(odir, str_glue("{prefix}_stable.png")), gs, width = 15, height = 5)

}

if (SAVE_FIG) {

    print("+++ Generating spaghetti plots +++")

    visualize_stable_progressor(adni_progressor, adni_cdr, opt$odir_fig, "adni")
    visualize_stable_progressor(oasis_progressor, oasis_cdr, opt$odir_fig, "oasis")
    visualize_stable_progressor(pac_progressor %>% mutate(amyloid_positive = amyloid_positive.site, tracer = "PIB"), pac_diagnosis, opt$odir_fig, "pac", "diagnosis", c("Normal", "MCI", "Dementia"))
    visualize_stable_progressor(a4_progressor %>% mutate(amyloid_positive = amyloid_positive.a4_thresh), a4_cdr, opt$odir_fig, "a4")
    visualize_stable_progressor(habs_progressor, habs_cdr, opt$odir_fig, "habs")
    # visualize_stable_progressor(habshd_progressor %>% mutate(amyloid_positive = amyloid_positive.fbb), habshd_cdr, opt$odir_fig, "habshd")
    visualize_stable_progressor(mcsa_progressor, mcsa_cdr, opt$odir_fig, "mcsa")

    odir_recovered <- file.path(opt$odir_fig, "recovered")
    visualize_stable_progressor(adni_progressor, adni_cdr, odir_recovered, "adni", recovered_bool = TRUE)
    visualize_stable_progressor(oasis_progressor, oasis_cdr, odir_recovered, "oasis", recovered_bool = TRUE)
    visualize_stable_progressor(pac_progressor %>% mutate(amyloid_positive = amyloid_positive.site, tracer = "PIB"), pac_diagnosis, odir_recovered, "pac", "diagnosis", c("Normal", "MCI", "Dementia"), recovered_bool = TRUE)
    visualize_stable_progressor(a4_progressor %>% mutate(amyloid_positive = amyloid_positive.a4_thresh), a4_cdr, odir_recovered, "a4", recovered_bool = TRUE)
    visualize_stable_progressor(habs_progressor, habs_cdr, odir_recovered, "habs", recovered_bool = TRUE)
    # visualize_stable_progressor(habshd_progressor %>% mutate(amyloid_positive = amyloid_positive.fbb), habshd_cdr, odir_recovered, "habshd", recovered_bool = TRUE)
    visualize_stable_progressor(mcsa_progressor, mcsa_cdr, odir_recovered, "mcsa", recovered_bool = TRUE)
}
