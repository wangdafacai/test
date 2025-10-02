#!/usr/bin/env Rscript
set.seed(123)

ensure_packages <- function(pkgs) {
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_install) > 0) {
    install.packages(to_install, repos = "https://cloud.r-project.org", dependencies = TRUE)
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

ensure_packages(c(
  "readr", "dplyr", "tidyr", "stringr", "purrr", "tibble", "readxl", "splines", "glue", "janitor"
))

`%||%` <- function(x, y) {
  if (is.null(x) || (is.atomic(x) && length(x) == 0)) y else x
}

build_age_grid <- function() {
  start <- c(seq(10, 75, by = 5), 80)
  end <- c(seq(14, 79, by = 5), Inf)
  tibble::tibble(
    age_start = start,
    age_end = end,
    age_group = factor(paste0(start, "-", ifelse(is.finite(end), end, "+")),
                       levels = paste0(start, "-", ifelse(is.finite(end), end, "+"))),
    age_mid = ifelse(is.finite(end), (start + end) / 2, start + 5)
  )
}

parse_age_weights <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    return(NULL)
  }
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    df <- readr::read_csv(path, show_col_types = FALSE) %>% janitor::clean_names()
  } else {
    sheets <- readxl::excel_sheets(path)
    sheet <- dplyr::case_when(
      "AgeWeights_Long" %in% sheets ~ "AgeWeights_Long",
      "AgeWeights" %in% sheets ~ "AgeWeights",
      TRUE ~ sheets[1]
    )
    df <- readxl::read_excel(path, sheet = sheet) %>% janitor::clean_names()
  }
  if (all(c("age_start", "age_end", "weight") %in% names(df))) {
    df <- df %>%
      dplyr::mutate(
        profile_id = dplyr::coalesce(
          df$profile_id %||% df$weight_profile_id %||% df$pattern_id %||% df$study_id %||% df$nid %||% "default"
        )
      )
  } else {
    age_cols <- grep("^age_", names(df), value = TRUE)
    if (length(age_cols) == 0) {
      warning(glue::glue("No age columns recognised in {path}."))
      return(NULL)
    }
    id_col <- setdiff(names(df), age_cols)
    if (length(id_col) == 0) {
      id_col <- "profile_id"
      df[[id_col]] <- "default"
    } else {
      id_col <- id_col[1]
    }
    df <- df %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(age_cols),
        names_to = "age_label",
        values_to = "weight"
      ) %>%
      dplyr::mutate(
        age_label = stringr::str_remove(age_label, "^age_"),
        age_label = stringr::str_replace_all(age_label, "\\.", "-"),
        age_label = stringr::str_replace(age_label, "+", "80+"),
        age_start = suppressWarnings(as.numeric(stringr::str_extract(age_label, "^\\d+"))),
        age_end = suppressWarnings(as.numeric(stringr::str_extract(age_label, "(?<=-)\\d+"))),
        age_end = dplyr::coalesce(age_end, Inf),
        profile_id = as.character(.data[[id_col]])
      ) %>%
      dplyr::select(profile_id, age_start, age_end, weight)
  }
  df %>%
    dplyr::mutate(
      age_start = as.numeric(age_start),
      age_end = as.numeric(age_end),
      weight = as.numeric(weight)
    ) %>%
    dplyr::filter(!is.na(age_start) & !is.na(weight))
}

normalise_profile_id <- function(df) {
  profile_cols <- intersect(names(df), c("weight_profile_id", "profile_id", "age_pattern_id", "study_id_clean", "study_id", "nid"))
  if (length(profile_cols) == 0) {
    rep("default", nrow(df))
  } else {
    first_col <- profile_cols[1]
    out <- df[[first_col]]
    out[is.na(out)] <- "default"
    as.character(out)
  }
}

expand_weights_to_grid <- function(weights, grid) {
  weights %>%
    dplyr::mutate(age_end = dplyr::coalesce(age_end, age_start + 4)) %>%
    dplyr::group_by(profile_id, age_start) %>%
    dplyr::summarise(weight = mean(weight, na.rm = TRUE), .groups = "drop") %>%
    tidyr::complete(profile_id, age_start = grid$age_start, fill = list(weight = 0)) %>%
    dplyr::arrange(profile_id, age_start)
}

compute_uniform_weights <- function(age_start, age_end, grid) {
  if (is.na(age_start)) age_start <- min(grid$age_start)
  if (is.na(age_end)) age_end <- max(grid$age_end[is.finite(grid$age_end)], na.rm = TRUE)
  if (!is.finite(age_end)) {
    age_end <- max(grid$age_end[is.finite(grid$age_end)], na.rm = TRUE)
  }
  overlaps <- purrr::map2_dbl(grid$age_start, grid$age_end, function(gs, ge) {
    ge <- ifelse(is.infinite(ge), max(grid$age_end[is.finite(grid$age_end)], na.rm = TRUE), ge)
    overlap_start <- max(age_start, gs)
    overlap_end <- min(age_end, ge)
    length <- overlap_end - overlap_start
    if (length < 0) 0 else length + 1
  })
  if (all(overlaps == 0)) {
    overlaps <- rep(1, length(overlaps))
  }
  overlaps / sum(overlaps)
}

weight_matrix_for_observations <- function(data, weights, grid, log_con) {
  profiles <- normalise_profile_id(data)
  mat <- matrix(0, nrow = nrow(data), ncol = nrow(grid))
  used_profiles <- character(0)
  if (is.null(weights)) {
    for (i in seq_len(nrow(data))) {
      mat[i, ] <- compute_uniform_weights(
        suppressWarnings(as.numeric(data$age_start[i])),
        suppressWarnings(as.numeric(data$age_end[i])),
        grid
      )
    }
    return(list(weights = mat, profile_lookup = profiles, used_profiles = used_profiles))
  }
  for (pid in unique(weights$profile_id)) {
    total <- sum(weights$weight[weights$profile_id == pid], na.rm = TRUE)
    if (total == 0 || is.na(total)) {
      weights$weight[weights$profile_id == pid] <- 0
    } else {
      weights$weight[weights$profile_id == pid] <- weights$weight[weights$profile_id == pid] / total
    }
  }
  for (i in seq_len(nrow(data))) {
    pid <- profiles[i]
    profile_weights <- weights %>% dplyr::filter(profile_id == pid)
    if (nrow(profile_weights) == 0) {
      if (!is.null(log_con)) writeLines(glue::glue("Row {i}: profile {pid} not found; using uniform weights."), log_con)
      mat[i, ] <- compute_uniform_weights(
        suppressWarnings(as.numeric(data$age_start[i])),
        suppressWarnings(as.numeric(data$age_end[i])),
        grid
      )
    } else {
      used_profiles <- union(used_profiles, pid)
      for (j in seq_along(grid$age_start)) {
        age_val <- grid$age_start[j]
        match_weight <- profile_weights$weight[profile_weights$age_start == age_val]
        mat[i, j] <- ifelse(length(match_weight) == 0, 0, match_weight)
      }
      if (sum(mat[i, ]) <= 0) {
        if (!is.null(log_con)) writeLines(glue::glue("Row {i}: profile {pid} weights sum to zero; using uniform."), log_con)
        mat[i, ] <- compute_uniform_weights(
          suppressWarnings(as.numeric(data$age_start[i])),
          suppressWarnings(as.numeric(data$age_end[i])),
          grid
        )
      }
    }
  }
  list(weights = mat, profile_lookup = profiles, used_profiles = used_profiles)
}

build_bs_basis <- function(grid, df) {
  ages <- grid$age_mid
  basis <- splines::bs(ages, df = df, intercept = TRUE)
  colnames(basis) <- paste0("bs_", seq_len(ncol(basis)))
  basis
}

save_outputs <- function(output_dir, obs_weights, weighted_basis, basis_grid, grid, profiles, log_messages) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  base_dir <- file.path(output_dir, "age_basis")
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(obs_weights, file.path(base_dir, "observation_age_weights.rds"))
  saveRDS(weighted_basis, file.path(base_dir, "observation_bspline_basis.rds"))
  saveRDS(basis_grid, file.path(base_dir, "bspline_basis_grid.rds"))
  readr::write_csv(grid, file.path(base_dir, "age_grid.csv"))
  readr::write_csv(tibble::tibble(row_id = seq_along(profiles), profile_id = profiles),
                   file.path(base_dir, "profile_lookup.csv"))
  if (length(log_messages) > 0) {
    writeLines(log_messages, con = file.path(base_dir, "age_basis_warnings.txt"))
  }
  list(
    weight_matrix = file.path(base_dir, "observation_age_weights.rds"),
    weighted_basis = file.path(base_dir, "observation_bspline_basis.rds"),
    basis_matrix = file.path(base_dir, "bspline_basis_grid.rds"),
    age_grid = file.path(base_dir, "age_grid.csv"),
    profile_lookup = file.path(base_dir, "profile_lookup.csv"),
    log = if (length(log_messages) > 0) file.path(base_dir, "age_basis_warnings.txt") else NA_character_
  )
}

build_age_basis <- function(clean_file = file.path("clean", "ed_clean_long.csv"),
                            age_weight_path = NULL,
                            output_dir = "clean",
                            spline_df = 6,
                            log_path = file.path(output_dir, "diagnostics", "age_basis_log.txt")) {
  if (!file.exists(clean_file)) {
    stop(glue::glue("Clean file {clean_file} not found. Run 01_cleaning.R first."))
  }
  data <- readr::read_csv(clean_file, show_col_types = FALSE)
  grid <- build_age_grid()
  log_messages <- character(0)
  log_con <- NULL
  if (!is.null(log_path)) {
    dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
    log_con <- file(log_path, open = "w")
  }
  on.exit({
    if (!is.null(log_con)) close(log_con)
  }, add = TRUE)
  weights <- parse_age_weights(age_weight_path)
  if (is.null(weights)) {
    log_messages <- c(log_messages, "Age weight file missing or unparseable; using uniform weights.")
    if (!is.null(log_con)) writeLines("Age weight file missing or unparseable; using uniform weights.", log_con)
  } else {
    weights <- expand_weights_to_grid(weights, grid)
  }
  weight_info <- weight_matrix_for_observations(data, weights, grid, log_con)
  obs_weights <- weight_info$weights
  zero_rows <- which(rowSums(obs_weights) == 0)
  if (length(zero_rows) > 0) {
    for (idx in zero_rows) {
      obs_weights[idx, ] <- compute_uniform_weights(
        suppressWarnings(as.numeric(data$age_start[idx])),
        suppressWarnings(as.numeric(data$age_end[idx])),
        grid
      )
    }
    msg <- glue::glue("{length(zero_rows)} observations defaulted to uniform weights.")
    log_messages <- c(log_messages, msg)
    if (!is.null(log_con)) writeLines(msg, log_con)
  }
  basis_grid <- build_bs_basis(grid, spline_df)
  weighted_basis <- obs_weights %*% basis_grid
  output_paths <- save_outputs(output_dir, obs_weights, weighted_basis, basis_grid, grid, weight_info$profile_lookup, log_messages)
  if (!is.null(log_path)) {
    output_paths$log_file <- log_path
  }
  output_paths
}

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  clean_file <- if (length(args) >= 1) args[1] else file.path("clean", "ed_clean_long.csv")
  age_weight <- if (length(args) >= 2) args[2] else NULL
  build_age_basis(clean_file = clean_file, age_weight_path = age_weight)
}
