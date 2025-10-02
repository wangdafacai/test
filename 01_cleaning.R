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
  "readxl", "dplyr", "tidyr", "stringr", "readr", "purrr", "lubridate",
  "janitor", "glue", "tibble", "rlang"
))

`%||%` <- function(x, y) {
  if (is.null(x) || (is.atomic(x) && length(x) == 0)) y else x
}

choose_sheet <- function(path, candidates) {
  sheets <- readxl::excel_sheets(path)
  match <- candidates[candidates %in% sheets][1]
  if (is.na(match)) {
    stop(glue::glue(
      "None of the candidate sheets ({paste(candidates, collapse = ', ')}) were found in {path}."
    ))
  }
  match
}

normalise_province <- function(x) {
  x <- stringr::str_trim(tidyr::replace_na(x, ""))
  x <- stringr::str_replace_all(x, "\\s+", " ")
  x <- stringr::str_to_title(x)
  dplyr::if_else(x == "", NA_character_, x)
}

coalesce_columns <- function(df, columns) {
  if (length(columns) == 0) {
    return(rep(NA_character_, nrow(df)))
  }
  out <- df[[columns[1]]]
  for (col in columns[-1]) {
    out <- dplyr::coalesce(out, df[[col]])
  }
  out
}

compute_year_mid <- function(df) {
  if (!"year_mid" %in% names(df)) {
    df$year_mid <- NA_real_
  }
  df <- df %>%
    dplyr::mutate(
      year_start = suppressWarnings(as.numeric(.data$year_start)),
      year_end = suppressWarnings(as.numeric(dplyr::coalesce(.data$year_end, .data$year_start))),
      year = suppressWarnings(as.numeric(dplyr::coalesce(.data$year, .data$year_start, .data$year_end, .data$year_mid))),
      year_mid = dplyr::coalesce(.data$year_mid, (.data$year_start + .data$year_end) / 2, .data$year)
    )
  df$year_mid
}

harmonise_prevalence_units <- function(df) {
  unit_cols <- intersect(names(df), c("measure_unit", "unit", "measure_units", "prevalence_unit"))
  unit_col <- if (length(unit_cols) > 0) unit_cols[1] else NA_character_
  lower_cols <- intersect(names(df), c("lower", "lower_ci", "lcl"))
  upper_cols <- intersect(names(df), c("upper", "upper_ci", "ucl"))
  lower_col <- if (length(lower_cols) > 0) lower_cols[1] else NA_character_
  upper_col <- if (length(upper_cols) > 0) upper_cols[1] else NA_character_
  mean_cols <- intersect(names(df), c("mean", "prevalence", "value"))
  mean_col <- if (length(mean_cols) > 0) mean_cols[1] else NA_character_
  if (is.na(unit_col) || is.na(mean_col)) {
    df$prevalence_unit <- "proportion"
    return(df)
  }
  unit_values <- tolower(dplyr::pull(df, unit_col))
  convert <- function(values, units) {
    purrr::map2_dbl(values, units, function(value, unit) {
      if (is.na(value)) {
        return(NA_real_)
      }
      if (is.na(unit)) {
        return(as.numeric(value))
      }
      unit <- stringr::str_replace_all(unit, "\\s", "")
      if (stringr::str_detect(unit, "per100000") || stringr::str_detect(unit, "per1e5")) {
        as.numeric(value) / 100000
      } else if (stringr::str_detect(unit, "per10000") || stringr::str_detect(unit, "per1e4")) {
        as.numeric(value) / 10000
      } else if (stringr::str_detect(unit, "per1000") || stringr::str_detect(unit, "per1e3")) {
        as.numeric(value) / 1000
      } else if (stringr::str_detect(unit, "per100") || stringr::str_detect(unit, "percent") ||
                 stringr::str_detect(unit, "%")) {
        as.numeric(value) / 100
      } else if (stringr::str_detect(unit, "per1000000") || stringr::str_detect(unit, "per1e6")) {
        as.numeric(value) / 1e6
      } else {
        as.numeric(value)
      }
    })
  }
  df[[mean_col]] <- convert(df[[mean_col]], unit_values)
  if (!is.na(lower_col)) {
    df[[lower_col]] <- convert(df[[lower_col]], unit_values)
  }
  if (!is.na(upper_col)) {
    df[[upper_col]] <- convert(df[[upper_col]], unit_values)
  }
  df$prevalence_unit <- "proportion"
  df
}

back_calculate_cases <- function(df) {
  mean_cols <- intersect(names(df), c("mean", "prevalence", "value"))
  mean_col <- if (length(mean_cols) > 0) mean_cols[1] else NA_character_
  sample_cols <- intersect(names(df), c("sample_size", "n", "denominator"))
  sample_col <- if (length(sample_cols) > 0) sample_cols[1] else NA_character_
  if (!"cases" %in% names(df)) {
    df$cases <- NA_real_
  }
  if (!is.na(mean_col) && !is.na(sample_col)) {
    missing_cases <- which(is.na(df$cases) & !is.na(df[[mean_col]]) & !is.na(df[[sample_col]]))
    df$cases[missing_cases] <- round(df[[mean_col]][missing_cases] * df[[sample_col]][missing_cases])
    missing_mean <- which(is.na(df[[mean_col]]) & !is.na(df$cases) & !is.na(df[[sample_col]]))
    df[[mean_col]][missing_mean] <- df$cases[missing_mean] / df[[sample_col]][missing_mean]
  }
  df
}

apply_validity_checks <- function(df) {
  mean_cols <- intersect(names(df), c("mean", "prevalence", "value"))
  mean_col <- if (length(mean_cols) > 0) mean_cols[1] else NA_character_
  sample_cols <- intersect(names(df), c("sample_size", "n", "denominator"))
  sample_col <- if (length(sample_cols) > 0) sample_cols[1] else NA_character_
  invalid <- tibble::tibble()
  if (!is.na(mean_col)) {
    outside <- df %>% dplyr::filter(!is.na(.data[[mean_col]]) & (.data[[mean_col]] < 0 | .data[[mean_col]] > 1))
    if (nrow(outside) > 0) {
      outside$issue <- "Prevalence outside [0,1]"
      invalid <- dplyr::bind_rows(invalid, outside)
    }
  }
  if (!is.na(sample_col)) {
    over_ss <- df %>% dplyr::filter(!is.na(.data$cases) & !is.na(.data[[sample_col]]) & .data$cases > .data[[sample_col]])
    if (nrow(over_ss) > 0) {
      over_ss$issue <- "Cases exceed sample size"
      invalid <- dplyr::bind_rows(invalid, over_ss)
    }
  }
  invalid
}

identify_mixed_duplicates <- function(df) {
  sex_col <- intersect(names(df), c("sex", "sex_id"))
  if (length(sex_col) == 0) {
    return(list(data = df, dropped = dplyr::tibble()))
  }
  sex_col <- sex_col[1]
  df <- df %>% dplyr::mutate(sex_clean = stringr::str_to_lower(as.character(.data[[sex_col]])))
  mixed_mask <- stringr::str_detect(df$sex_clean %||% "", "both|all|mixed")
  if (!any(mixed_mask)) {
    df <- dplyr::select(df, -"sex_clean")
    return(list(data = df, dropped = dplyr::tibble()))
  }
  key_cols <- intersect(names(df), c(
    "study_id", "study", "nid", "source", "source_id", "publication_id",
    "province", "admin_1", "admin1", "location", "year_mid", "year_start", "year_end",
    "age_start", "age_end"
  ))
  if (length(key_cols) == 0) {
    df <- dplyr::select(df, -"sex_clean")
    return(list(data = df, dropped = dplyr::tibble()))
  }
  dupes <- df %>%
    dplyr::filter(mixed_mask) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(key_cols))) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::ungroup()
  if (nrow(dupes) == 0) {
    df <- dplyr::select(df, -"sex_clean")
    return(list(data = df, dropped = dplyr::tibble()))
  }
  to_drop <- dupes %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(key_cols))) %>%
    dplyr::slice(-1) %>%
    dplyr::ungroup()
  df <- df %>% dplyr::anti_join(to_drop %>% dplyr::select(dplyr::all_of(key_cols)), by = key_cols)
  df <- dplyr::select(df, -"sex_clean")
  list(data = df, dropped = dplyr::select(to_drop, -"sex_clean"))
}

assign_indices <- function(df) {
  study_cols <- intersect(names(df), c("study_id", "nid", "source", "publication_id"))
  if (length(study_cols) == 0) {
    df$study_id_clean <- paste0("study_", seq_len(nrow(df)))
  } else {
    first_col <- study_cols[1]
    df$study_id_clean <- dplyr::if_else(is.na(df[[first_col]]), paste0("study_", seq_len(nrow(df))), as.character(df[[first_col]]))
  }
  df$study_idx <- as.integer(factor(df$study_id_clean))
  df$row_id <- seq_len(nrow(df))
  df
}

export_diagnostics <- function(df, invalid, dropped, diagnostics_dir) {
  dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
  mean_cols <- intersect(names(df), c("mean", "prevalence", "value"))
  mean_col <- if (length(mean_cols) > 0) mean_cols[1] else NA_character_
  if (!is.na(mean_col)) {
    high <- df %>% dplyr::arrange(dplyr::desc(.data[[mean_col]])) %>% dplyr::slice_head(n = 20)
    low <- df %>% dplyr::arrange(.data[[mean_col]]) %>% dplyr::slice_head(n = 20)
    readr::write_csv(high, file.path(diagnostics_dir, "extreme_prevalence_high.csv"))
    readr::write_csv(low, file.path(diagnostics_dir, "extreme_prevalence_low.csv"))
  }
  if (nrow(invalid) > 0) {
    readr::write_csv(invalid, file.path(diagnostics_dir, "invalid_rows.csv"))
  }
  if (nrow(dropped) > 0) {
    readr::write_csv(dropped, file.path(diagnostics_dir, "dropped_mixed_duplicates.csv"))
  }
  overlap_cols <- intersect(names(df), c("study_idx", "province", "year_mid", "age_start", "age_end", "sex"))
  if (length(overlap_cols) >= 4) {
    overlaps <- df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(overlap_cols))) %>%
      dplyr::filter(dplyr::n() > 1) %>%
      dplyr::ungroup()
    if (nrow(overlaps) > 0) {
      readr::write_csv(overlaps, file.path(diagnostics_dir, "potential_overlaps.csv"))
    }
  }
}

run_cleaning <- function(input_path = "ED_meta_template_v5.xlsx",
                         output_dir = "clean",
                         diagnostics_dir = file.path(output_dir, "diagnostics")) {
  if (!file.exists(input_path)) {
    stop(glue::glue("Input file {input_path} not found."))
  }
  ext <- tolower(tools::file_ext(input_path))
  if (ext == "csv") {
    df_long <- readr::read_csv(input_path, show_col_types = FALSE) %>% janitor::clean_names()
  } else {
    sheet_long <- tryCatch(choose_sheet(input_path, c("Data_Long", "Long", "Data")), error = function(e) NA_character_)
    if (is.na(sheet_long)) {
      stop("No long-format sheet found in the input workbook.")
    }
    df_long <- readxl::read_excel(input_path, sheet = sheet_long) %>%
      janitor::clean_names()
  }
  province_cols <- intersect(names(df_long), c("province", "admin_1", "admin1", "location", "area"))
  df_long$province <- normalise_province(coalesce_columns(df_long, province_cols))
  df_long$year_mid <- compute_year_mid(df_long)
  df_long <- harmonise_prevalence_units(df_long)
  df_long <- back_calculate_cases(df_long)
  invalid <- apply_validity_checks(df_long)
  mixed <- identify_mixed_duplicates(df_long)
  df_long <- mixed$data
  dropped <- mixed$dropped
  df_long <- assign_indices(df_long)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(output_dir, "ed_clean_long.csv")
  readr::write_csv(df_long, out_file)
  export_diagnostics(df_long, invalid, dropped, diagnostics_dir)
  list(clean_file = out_file, diagnostics_dir = diagnostics_dir,
       invalid_file = if (nrow(invalid) > 0) file.path(diagnostics_dir, "invalid_rows.csv") else NA_character_,
       dropped_file = if (nrow(dropped) > 0) file.path(diagnostics_dir, "dropped_mixed_duplicates.csv") else NA_character_)
}

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  args <- commandArgs(trailingOnly = TRUE)
  input <- if (length(args) >= 1) args[1] else "ED_meta_template_v5.xlsx"
  output_dir <- if (length(args) >= 2) args[2] else "clean"
  run_cleaning(input_path = input, output_dir = output_dir)
}
