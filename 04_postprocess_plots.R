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
  "readr", "dplyr", "tidyr", "stringr", "purrr", "tibble", "ggplot2", "glue", "scales"
))

logit_inv <- function(x) 1 / (1 + exp(-x))

extract_draws <- function(samples, pattern) {
  cols <- grep(pattern, colnames(samples), value = TRUE)
  if (length(cols) == 0) {
    return(NULL)
  }
  samples[, cols, drop = FALSE]
}

summarise_matrix_draws <- function(mat) {
  tibble::tibble(
    mean = colMeans(mat),
    lower = apply(mat, 2, stats::quantile, probs = 0.025),
    upper = apply(mat, 2, stats::quantile, probs = 0.975)
  )
}

compute_observation_draws <- function(draws, metadata) {
  S <- nrow(draws)
  N <- length(metadata$sex_idx)
  alpha <- draws[, "alpha"]
  spline_draws <- extract_draws(draws, "^spline_coefs\\[")
  method_draws <- extract_draws(draws, "^method_beta\\[")
  bias_draws <- extract_draws(draws, "^bias_beta\\[")
  province_draws <- extract_draws(draws, "^province_effect\\[")
  study_draws <- extract_draws(draws, "^study_effect\\[")
  sex_draws <- cbind(draws[, "sex_beta[1]"], draws[, "sex_beta[2]"])
  eta_base <- matrix(alpha, nrow = S, ncol = N)
  if (!is.null(spline_draws)) {
    eta_base <- eta_base + t(metadata$basis %*% t(spline_draws))
  }
  if (!is.null(method_draws)) {
    eta_base <- eta_base + t(metadata$method_matrix %*% t(method_draws))
  }
  if (!is.null(bias_draws)) {
    eta_base <- eta_base + t(metadata$bias_matrix %*% t(bias_draws))
  }
  if (!is.null(province_draws)) {
    eta_base <- eta_base + province_draws[, metadata$province_idx, drop = FALSE]
  }
  if (!is.null(study_draws)) {
    eta_base <- eta_base + study_draws[, metadata$study_idx, drop = FALSE]
  }
  eta_female <- eta_base + matrix(sex_draws[, 1], nrow = S, ncol = N)
  eta_male <- eta_base + matrix(sex_draws[, 2], nrow = S, ncol = N)
  prob <- matrix(0, nrow = S, ncol = N)
  female_cols <- which(metadata$sex_idx == 1)
  male_cols <- which(metadata$sex_idx == 2)
  mix_cols <- which(metadata$mix_indicator == 1)
  if (length(female_cols) > 0) {
    prob[, female_cols] <- logit_inv(eta_female[, female_cols, drop = FALSE])
  }
  if (length(male_cols) > 0) {
    prob[, male_cols] <- logit_inv(eta_male[, male_cols, drop = FALSE])
  }
  if (length(mix_cols) > 0) {
    prob[, mix_cols] <- logit_inv(eta_female[, mix_cols, drop = FALSE]) * rep(metadata$female_share[mix_cols], each = S) +
      logit_inv(eta_male[, mix_cols, drop = FALSE]) * rep(1 - metadata$female_share[mix_cols], each = S)
  }
  prob
}

create_observation_summary <- function(prob, metadata) {
  obs <- tibble::tibble(
    row_id = metadata$row_id,
    province = metadata$provinces,
    study_idx = metadata$study_idx,
    sex_idx = metadata$sex_idx,
    year_mid = metadata$year_mid,
    age_start = metadata$age_start,
    age_end = metadata$age_end,
    sample_size = metadata$sample_size,
    cases = metadata$cases
  )
  stats <- summarise_matrix_draws(prob)
  dplyr::bind_cols(obs, stats) %>%
    dplyr::mutate(predicted_cases = .data$mean * .data$sample_size)
}

summarise_group <- function(df, group_cols) {
  if (length(group_cols) == 0) {
    total_ss <- sum(df$sample_size, na.rm = TRUE)
    if (total_ss == 0) {
      return(tibble::tibble(sample_size = 0, weighted_mean = NA_real_, weighted_lower = NA_real_, weighted_upper = NA_real_))
    }
    wm <- sum(df$mean * df$sample_size, na.rm = TRUE) / total_ss
    wl <- sum(df$lower * df$sample_size, na.rm = TRUE) / total_ss
    wu <- sum(df$upper * df$sample_size, na.rm = TRUE) / total_ss
    return(tibble::tibble(sample_size = total_ss, weighted_mean = wm, weighted_lower = wl, weighted_upper = wu))
  }
  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      sample_size = sum(.data$sample_size, na.rm = TRUE),
      weighted_mean = sum(.data$mean * .data$sample_size, na.rm = TRUE) / sum(.data$sample_size, na.rm = TRUE),
      weighted_lower = sum(.data$lower * .data$sample_size, na.rm = TRUE) / sum(.data$sample_size, na.rm = TRUE),
      weighted_upper = sum(.data$upper * .data$sample_size, na.rm = TRUE) / sum(.data$sample_size, na.rm = TRUE),
      .groups = "drop"
    )
}

compute_age_draws <- function(draws, basis_grid, sex_draws) {
  alpha <- draws[, "alpha"]
  spline_draws <- extract_draws(draws, "^spline_coefs\\[")
  eta <- matrix(alpha, nrow = nrow(draws), ncol = nrow(basis_grid))
  if (!is.null(spline_draws)) {
    eta <- eta + t(basis_grid %*% t(spline_draws))
  }
  eta <- eta + matrix(sex_draws, nrow = nrow(draws), ncol = nrow(basis_grid))
  logit_inv(eta)
}

summarise_age_draws <- function(prob_matrix, ages, label) {
  stats <- summarise_matrix_draws(prob_matrix)
  stats$age_group <- ages
  stats$sex <- label
  stats
}

parse_diag_file <- function(path) {
  if (!file.exists(path)) return(tibble::tibble())
  lines <- readr::read_lines(path)
  start_rhat <- which(stringr::str_detect(lines, "^R-hat"))
  start_ess <- which(stringr::str_detect(lines, "^Effective sample size"))
  if (length(start_rhat) == 0 || length(start_ess) == 0) return(tibble::tibble())
  rhat_lines <- lines[(start_rhat + 1):(start_ess - 1)]
  ess_lines <- lines[(start_ess + 1):length(lines)]
  tibble::tibble(
    parameter = stringr::str_extract(rhat_lines, "^[^:]+"),
    rhat = as.numeric(stringr::str_extract(rhat_lines, "[0-9.]+$")),
    ess = as.numeric(stringr::str_extract(ess_lines, "[0-9.]+$"))
  )
}

collect_sensitivity_results <- function(base_summary, output_dir) {
  configs <- list(
    age_cap = file.path(output_dir, "sensitivity_age_cap", "posterior_samples.rds"),
    sex_equal = file.path(output_dir, "sensitivity_sex_equal", "posterior_samples.rds"),
    bias_exclude = file.path(output_dir, "sensitivity_bias_exclude", "posterior_samples.rds")
  )
  out <- list()
  for (name in names(configs)) {
    file <- configs[[name]]
    if (file.exists(file)) {
      samples <- readRDS(file)
      draws <- as.matrix(samples)
      meta <- readRDS(file.path(dirname(file), "model_metadata.rds"))
      prob <- compute_observation_draws(draws, meta)
      summary <- create_observation_summary(prob, meta)
      national <- summarise_group(summary, character(0)) %>% dplyr::mutate(config = name)
      out[[name]] <- national
    }
  }
  base_national <- summarise_group(base_summary, character(0)) %>% dplyr::mutate(config = "baseline")
  if (length(out) == 0) {
    return(base_national)
  }
  dplyr::bind_rows(base_national, dplyr::bind_rows(out))
}

run_postprocess <- function(posterior_file = file.path("ed_bayes_outputs", "model", "posterior_samples.rds"),
                            metadata_file = file.path("ed_bayes_outputs", "model", "model_metadata.rds"),
                            clean_file = file.path("clean", "ed_clean_long.csv"),
                            age_basis_dir = file.path("clean", "age_basis"),
                            output_dir = "ed_bayes_outputs") {
  if (!file.exists(posterior_file)) {
    stop(glue::glue("Posterior samples not found at {posterior_file}."))
  }
  samples <- readRDS(posterior_file)
  draws <- as.matrix(samples)
  metadata <- readRDS(metadata_file)
  basis_grid <- readRDS(file.path(age_basis_dir, "bspline_basis_grid.rds"))
  age_grid <- readr::read_csv(file.path(age_basis_dir, "age_grid.csv"), show_col_types = FALSE)
  prob <- compute_observation_draws(draws, metadata)
  observation_summary <- create_observation_summary(prob, metadata)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  summaries_dir <- file.path(output_dir, "summaries")
  plots_dir <- file.path(output_dir, "plots")
  dir.create(summaries_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(observation_summary, file.path(summaries_dir, "observation_summary.csv"))
  province_summary <- summarise_group(observation_summary, "province")
  readr::write_csv(province_summary, file.path(summaries_dir, "province_summary.csv"))
  national_summary <- summarise_group(observation_summary, character(0)) %>% dplyr::mutate(province = "National")
  readr::write_csv(national_summary, file.path(summaries_dir, "national_summary.csv"))
  time_summary <- summarise_group(observation_summary, "year_mid") %>% dplyr::filter(!is.na(year_mid))
  readr::write_csv(time_summary, file.path(summaries_dir, "time_summary.csv"))
  sex_levels <- c("Female", "Male", "Mixed")
  observation_summary$sex <- factor(sex_levels[observation_summary$sex_idx], levels = sex_levels)
  female_draws <- compute_age_draws(draws, basis_grid, draws[, "sex_beta[1]"])
  male_draws <- compute_age_draws(draws, basis_grid, draws[, "sex_beta[2]"])
  age_curves_f <- summarise_age_draws(female_draws, age_grid$age_group, "Female")
  age_curves_m <- summarise_age_draws(male_draws, age_grid$age_group, "Male")
  female_fraction_vec <- dplyr::case_when(
    metadata$sex_idx == 1 ~ 1,
    metadata$sex_idx == 2 ~ 0,
    TRUE ~ dplyr::coalesce(metadata$female_share, 0.5)
  )
  total_sample <- sum(metadata$sample_size, na.rm = TRUE)
  female_total <- sum(metadata$sample_size * female_fraction_vec, na.rm = TRUE)
  female_weight <- if (is.finite(total_sample) && total_sample > 0) female_total / total_sample else 0.5
  if (!is.finite(female_weight) || female_weight < 0 || female_weight > 1) {
    female_weight <- 0.5
  }
  combined_draws <- female_weight * female_draws + (1 - female_weight) * male_draws
  age_curves_overall <- summarise_age_draws(combined_draws, age_grid$age_group, "Overall")
  age_curve_summary <- dplyr::bind_rows(age_curves_f, age_curves_m, age_curves_overall)
  readr::write_csv(age_curve_summary, file.path(summaries_dir, "age_group_summary.csv"))
  readr::write_csv(dplyr::filter(age_curve_summary, sex %in% c("Female", "Male")),
                   file.path(summaries_dir, "age_curves.csv"))
  sensitivity <- collect_sensitivity_results(observation_summary, output_dir)
  readr::write_csv(sensitivity, file.path(summaries_dir, "sensitivity_summary.csv"))
  # Plots
  province_plot <- ggplot2::ggplot(province_summary, ggplot2::aes(x = reorder(province, weighted_mean), y = weighted_mean)) +
    ggplot2::geom_col(fill = "#3c6e71") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = weighted_lower, ymax = weighted_upper), width = 0.2) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    ggplot2::labs(title = "Province-level prevalence", x = NULL, y = "Prevalence") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(plots_dir, "province_prevalence.png"), province_plot, width = 8, height = 6, dpi = 120)
  national_plot <- ggplot2::ggplot(time_summary, ggplot2::aes(x = year_mid, y = weighted_mean)) +
    ggplot2::geom_line(color = "#284b63") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = weighted_lower, ymax = weighted_upper), fill = "#96c5f7", alpha = 0.4) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    ggplot2::labs(title = "Time trend (2000-2025)", x = "Year", y = "Prevalence") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(plots_dir, "time_trend.png"), national_plot, width = 8, height = 5, dpi = 120)
  age_curves_plot <- dplyr::filter(age_curve_summary, sex %in% c("Female", "Male"))
  age_plot <- ggplot2::ggplot(age_curves_plot, ggplot2::aes(x = age_group, y = mean, color = sex, group = sex)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = sex), alpha = 0.2, colour = NA) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    ggplot2::labs(title = "Age curves", x = "Age group", y = "Prevalence", color = "Sex", fill = "Sex") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(plots_dir, "age_curves.png"), age_plot, width = 9, height = 5, dpi = 120)
  sex_plot <- observation_summary %>%
    dplyr::filter(!is.na(year_mid)) %>%
    dplyr::group_by(year_mid, sex) %>%
    dplyr::summarise(weighted_mean = sum(mean * sample_size, na.rm = TRUE) / sum(sample_size, na.rm = TRUE), .groups = "drop") %>%
    ggplot2::ggplot(ggplot2::aes(x = year_mid, y = weighted_mean, color = sex)) +
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    ggplot2::labs(title = "Sex comparisons over time", x = "Year", y = "Prevalence", color = "Sex") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(plots_dir, "sex_comparison.png"), sex_plot, width = 8, height = 5, dpi = 120)
  diag <- parse_diag_file(file.path(dirname(posterior_file), "mcmc_diag.txt"))
  if (nrow(diag) > 0) {
    diag_plot <- ggplot2::ggplot(diag, ggplot2::aes(x = reorder(parameter, rhat), y = rhat)) +
      ggplot2::geom_col(fill = "#8d99ae") +
      ggplot2::geom_hline(yintercept = 1.1, colour = "red", linetype = "dashed") +
      ggplot2::coord_flip() +
      ggplot2::labs(title = "R-hat diagnostics", x = NULL, y = "R-hat") +
      ggplot2::theme_minimal()
    ggplot2::ggsave(file.path(plots_dir, "diagnostics_rhat.png"), diag_plot, width = 8, height = 6, dpi = 120)
  }
  sensitivity_plot <- ggplot2::ggplot(sensitivity, ggplot2::aes(x = config, y = weighted_mean)) +
    ggplot2::geom_col(fill = "#2a9d8f") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = weighted_lower, ymax = weighted_upper), width = 0.2) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    ggplot2::labs(title = "Sensitivity comparisons", x = "Configuration", y = "Prevalence") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(plots_dir, "sensitivity.png"), sensitivity_plot, width = 7, height = 5, dpi = 120)
  list(
    observation_summary = file.path(summaries_dir, "observation_summary.csv"),
    province_summary = file.path(summaries_dir, "province_summary.csv"),
    national_summary = file.path(summaries_dir, "national_summary.csv"),
    time_summary = file.path(summaries_dir, "time_summary.csv"),
    age_curves = file.path(summaries_dir, "age_curves.csv"),
    sensitivity = file.path(summaries_dir, "sensitivity_summary.csv"),
    plots_dir = plots_dir
  )
}

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  run_postprocess()
}
