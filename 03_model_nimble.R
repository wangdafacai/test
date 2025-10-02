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
  "readr", "dplyr", "tibble", "stringr", "purrr", "nimble", "coda", "bayesplot", "ggplot2", "glue"
))

`%||%` <- function(x, y) {
  if (is.null(x) || (is.atomic(x) && length(x) == 0)) y else x
}

prepare_model_data <- function(clean_file, basis_file, weight_file = NULL, sensitivity = list()) {
  data <- readr::read_csv(clean_file, show_col_types = FALSE)
  if (!is.null(sensitivity$exclude_bias)) {
    bias_cols <- intersect(names(data), c("bias", "bias_flag", "exclude_bias"))
    if (length(bias_cols) > 0 && isTRUE(sensitivity$exclude_bias)) {
      data <- data %>% dplyr::filter(is.na(.data[[bias_cols[1]]]) | .data[[bias_cols[1]]] == 0)
    }
  }
  if (!is.null(sensitivity$age_cap) && isTRUE(sensitivity$age_cap)) {
    if ("age_end" %in% names(data)) {
      data$age_end <- pmin(suppressWarnings(as.numeric(data$age_end)), 80)
    }
  }
  data <- data %>% dplyr::filter(!is.na(.data$cases), !is.na(.data$sample_size))
  data <- data %>% dplyr::mutate(
    cases = as.numeric(.data$cases),
    sample_size = as.numeric(.data$sample_size)
  )
  basis <- readRDS(basis_file)
  if (nrow(basis) != nrow(data)) {
    stop("Basis matrix and data have mismatched rows.")
  }
  obs_weights <- if (!is.null(weight_file) && file.exists(weight_file)) readRDS(weight_file) else NULL
  sex_levels <- c("female", "male", "both")
  sex_col <- intersect(names(data), c("sex", "sex_id"))
  if (length(sex_col) == 0) {
    data$sex <- "female"
  } else {
    sex_col <- sex_col[1]
    data$sex <- stringr::str_to_lower(as.character(data[[sex_col]]))
  }
  data$sex[data$sex %in% c("f", "women", "girl")] <- "female"
  data$sex[data$sex %in% c("m", "men", "boy")] <- "male"
  data$sex[!data$sex %in% sex_levels] <- "both"
  sex_idx <- match(data$sex, sex_levels)
  mix_indicator <- ifelse(sex_idx == 3, 1, 0)
  female_share_cols <- intersect(names(data), c("female_weight", "female_share", "female_ratio", "prop_female"))
  if (length(female_share_cols) == 0) {
    female_share <- rep(0.5, nrow(data))
  } else {
    female_share <- suppressWarnings(as.numeric(data[[female_share_cols[1]]]))
    female_share[is.na(female_share)] <- 0.5
    female_share <- pmin(pmax(female_share, 0.05), 0.95)
  }
  if (!is.null(sensitivity$sex_equal) && isTRUE(sensitivity$sex_equal)) {
    female_share <- rep(0.5, length(female_share))
  }
  provinces <- data$province %||% data$admin_1 %||% data$location
  if (is.null(provinces)) {
    provinces <- rep("Unknown", nrow(data))
  }
  province_idx <- as.integer(factor(provinces))
  studies <- data$study_idx %||% data$study_id_clean %||% data$study_id %||% data$nid
  if (is.null(studies)) {
    studies <- seq_len(nrow(data))
  }
  study_idx <- as.integer(factor(studies))
  method_cols <- intersect(names(data), c("method", "method_group", "instrument", "collection_method"))
  bias_cols <- intersect(names(data), c("bias_category", "bias_type", "risk_of_bias"))
  method_matrix <- matrix(0, nrow = nrow(data), ncol = 0)
  method_names <- character(0)
  if (length(method_cols) > 0) {
    method_factor <- factor(data[[method_cols[1]]])
    method_matrix <- model.matrix(~ method_factor - 1)
    colnames(method_matrix) <- stringr::str_replace(colnames(method_matrix), "method_factor", "method")
    method_names <- colnames(method_matrix)
  }
  bias_matrix <- matrix(0, nrow = nrow(data), ncol = 0)
  bias_names <- character(0)
  if (length(bias_cols) > 0) {
    bias_factor <- factor(data[[bias_cols[1]]])
    bias_matrix <- model.matrix(~ bias_factor - 1)
    colnames(bias_matrix) <- stringr::str_replace(colnames(bias_matrix), "bias_factor", "bias")
    bias_names <- colnames(bias_matrix)
  }
  has_method <- ncol(method_matrix) > 0
  has_bias <- ncol(bias_matrix) > 0
  if (!has_method) {
    method_matrix <- matrix(0, nrow = nrow(data), ncol = 1)
    method_names <- "method_none"
  }
  if (!has_bias) {
    bias_matrix <- matrix(0, nrow = nrow(data), ncol = 1)
    bias_names <- "bias_none"
  }
  list(
    data = data,
    basis = basis,
    obs_weights = obs_weights,
    sex_idx = sex_idx,
    mix_indicator = mix_indicator,
    female_share = female_share,
    province_idx = province_idx,
    study_idx = study_idx,
    method_matrix = method_matrix,
    bias_matrix = bias_matrix,
    method_names = method_names,
    bias_names = bias_names,
    has_method = has_method,
    has_bias = has_bias,
    provinces = provinces,
    studies = studies,
    method_cols = method_cols,
    bias_cols = bias_cols
  )
}

build_model_code <- function(constants) {
  nimble::nimbleCode({
    for (i in 1:N) {
      eta_base[i] <- alpha + province_effect[province_idx[i]] + study_effect[study_idx[i]] + inprod(basis_matrix[i, 1:K], spline_coefs[1:K])
      eta_method[i] <- inprod(method_matrix[i, 1:n_method], method_beta[1:n_method])
      eta_bias[i] <- inprod(bias_matrix[i, 1:n_bias], bias_beta[1:n_bias])
      eta_female[i] <- eta_base[i] + eta_method[i] + eta_bias[i] + sex_beta[1]
      eta_male[i] <- eta_base[i] + eta_method[i] + eta_bias[i] + sex_beta[2]
      p_female[i] <- ilogit(eta_female[i])
      p_male[i] <- ilogit(eta_male[i])
      p_single[i] <- p_female[i] * equals(sex_idx[i], 1) + p_male[i] * equals(sex_idx[i], 2)
      p_mix[i] <- female_share[i] * p_female[i] + (1 - female_share[i]) * p_male[i]
      p[i] <- (1 - mix_indicator[i]) * p_single[i] + mix_indicator[i] * p_mix[i]
      cases[i] ~ dbin(p[i], sample_size[i])
    }
    alpha ~ dnorm(0, sd = 5)
    for (s in 1:2) {
      sex_beta[s] ~ dnorm(0, sd = 2)
    }
    for (k in 1:n_method) {
      method_beta[k] ~ dnorm(0, sd = method_prior_sd[k])
    }
    for (k in 1:n_bias) {
      bias_beta[k] ~ dnorm(0, sd = bias_prior_sd[k])
    }
    for (k in 1:K) {
      spline_coefs[k] ~ dnorm(0, sd = 2)
    }
    for (p in 1:n_province) {
      province_effect[p] ~ dnorm(0, tau = tau_province)
    }
    tau_province ~ dgamma(2, 2)
    for (j in 1:n_study) {
      study_effect[j] ~ dnorm(0, tau = tau_study)
    }
    tau_study ~ dgamma(2, 2)
  })
}

initialise_chains <- function(constants, chains) {
  inits <- vector("list", chains)
  for (i in seq_len(chains)) {
    inits[[i]] <- list(
      alpha = rnorm(1, 0, 0.1),
      sex_beta = rnorm(2, 0, 0.1),
      method_beta = rep(0, constants$n_method),
      bias_beta = rep(0, constants$n_bias),
      spline_coefs = rnorm(constants$K, 0, 0.1),
      province_effect = rnorm(constants$n_province, 0, 0.1),
      study_effect = rnorm(constants$n_study, 0, 0.1),
      tau_province = rgamma(1, 2, 2),
      tau_study = rgamma(1, 2, 2)
    )
  }
  inits
}

run_mcmc_with_diagnostics <- function(model, mcmc_conf, constants, chains = 4, iter = 6000, burnin = 2000, thin = 5,
                                       max_iter = 20000, ess_threshold = 200, rhat_threshold = 1.1, monitors) {
  mcmc_builder <- nimble::buildMCMC(mcmc_conf)
  cmodel <- nimble::compileNimble(model)
  cmcmc <- nimble::compileNimble(mcmc_builder, project = model)
  current_iter <- iter
  samples_list <- NULL
  inits <- initialise_chains(constants, chains)
  repeat {
    samples_list <- nimble::runMCMC(cmcmc, niter = current_iter, nchains = chains, nburnin = burnin,
                                   thin = thin, samplesAsCodaMCMC = TRUE, monitors = monitors, inits = inits)
    combined <- samples_list
    rhat <- tryCatch(coda::gelman.diag(combined, autoburnin = FALSE)$psrf[, "Point est."], error = function(e) rep(NA_real_, length(monitors)))
    ess <- tryCatch(coda::effectiveSize(combined), error = function(e) rep(NA_real_, length(monitors)))
    if (all(!is.na(rhat) & rhat <= rhat_threshold) && all(!is.na(ess) & ess >= ess_threshold)) {
      break
    }
    if (current_iter >= max_iter) {
      warning("Reached maximum iterations without satisfying convergence thresholds.")
      break
    }
    current_iter <- min(current_iter * 2, max_iter)
  }
  list(samples = combined, rhat = rhat, ess = ess, iter = current_iter, burnin = burnin, thin = thin)
}

save_mcmc_outputs <- function(samples, output_dir, meta, monitors) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  posterior_path <- file.path(output_dir, "posterior_samples.rds")
  saveRDS(samples, posterior_path)
  diag_path <- file.path(output_dir, "mcmc_diag.txt")
  diag_lines <- c(
    glue::glue("Total iterations: {meta$iter}"),
    glue::glue("Burn-in: {meta$burnin}"),
    glue::glue("Thin: {meta$thin}"),
    "R-hat:",
    paste(names(meta$rhat), sprintf("%.3f", meta$rhat), sep = ": "),
    "Effective sample size:",
    paste(names(meta$ess), sprintf("%.1f", meta$ess), sep = ": ")
  )
  writeLines(diag_lines, diag_path)
  diag_dir <- file.path(output_dir, "diagnostics")
  dir.create(diag_dir, showWarnings = FALSE)
  trace_plot <- bayesplot::mcmc_trace(samples, pars = monitors)
  ggplot2::ggsave(file.path(diag_dir, "trace_plot.png"), trace_plot, width = 12, height = 8, dpi = 120)
  density_plot <- bayesplot::mcmc_dens_overlay(samples, pars = monitors)
  ggplot2::ggsave(file.path(diag_dir, "density_plot.png"), density_plot, width = 12, height = 8, dpi = 120)
  list(posterior = posterior_path, diag = diag_path, trace = file.path(diag_dir, "trace_plot.png"), density = file.path(diag_dir, "density_plot.png"))
}

run_nimble_model <- function(clean_file = file.path("clean", "ed_clean_long.csv"),
                             basis_file = file.path("clean", "age_basis", "observation_bspline_basis.rds"),
                             weight_file = file.path("clean", "age_basis", "observation_age_weights.rds"),
                             output_dir = file.path("ed_bayes_outputs", "model"),
                             sensitivity = list(),
                             mcmc = list(chains = 4, iter = 6000, burnin = 2000, thin = 5, max_iter = 20000,
                                         ess_threshold = 200, rhat_threshold = 1.1)) {
  prep <- prepare_model_data(clean_file, basis_file, weight_file, sensitivity)
  constants <- list(
    N = nrow(prep$data),
    K = ncol(prep$basis),
    n_method = ncol(prep$method_matrix),
    n_bias = ncol(prep$bias_matrix),
    n_province = max(prep$province_idx),
    n_study = max(prep$study_idx)
  )
  method_prior_sd <- if (prep$has_method) rep(1, constants$n_method) else rep(1e-6, constants$n_method)
  bias_prior_sd <- if (prep$has_bias) rep(1, constants$n_bias) else rep(1e-6, constants$n_bias)
  constants$method_prior_sd <- method_prior_sd
  constants$bias_prior_sd <- bias_prior_sd
  data_list <- list(
    cases = prep$data$cases,
    sample_size = prep$data$sample_size,
    sex_idx = prep$sex_idx,
    mix_indicator = prep$mix_indicator,
    female_share = prep$female_share,
    province_idx = prep$province_idx,
    study_idx = prep$study_idx,
    basis_matrix = prep$basis,
    method_matrix = prep$method_matrix,
    bias_matrix = prep$bias_matrix
  )
  monitors <- c("alpha", "sex_beta", "method_beta", "bias_beta", "spline_coefs", "province_effect", "study_effect", "tau_province", "tau_study")
  model_code <- build_model_code(constants)
  model <- nimble::nimbleModel(code = model_code, data = data_list, constants = constants)
  mcmc_conf <- nimble::configureMCMC(model, monitors = monitors)
  meta <- run_mcmc_with_diagnostics(model, mcmc_conf, constants,
                                    chains = mcmc$chains, iter = mcmc$iter, burnin = mcmc$burnin, thin = mcmc$thin,
                                    max_iter = mcmc$max_iter, ess_threshold = mcmc$ess_threshold, rhat_threshold = mcmc$rhat_threshold,
                                    monitors = monitors)
  metadata <- list(
    method_names = prep$method_names,
    bias_names = prep$bias_names,
    has_method = prep$has_method,
    has_bias = prep$has_bias,
    provinces = prep$provinces,
    province_idx = prep$province_idx,
    province_levels = levels(factor(prep$provinces)),
    studies = prep$studies,
    study_idx = prep$study_idx,
    study_levels = levels(factor(prep$studies)),
    sex_idx = prep$sex_idx,
    female_share = prep$female_share,
    mix_indicator = prep$mix_indicator,
    method_matrix_colnames = colnames(prep$method_matrix),
    bias_matrix_colnames = colnames(prep$bias_matrix),
    method_matrix = prep$method_matrix,
    bias_matrix = prep$bias_matrix,
    basis = prep$basis,
    row_id = prep$data$row_id %||% seq_len(nrow(prep$data)),
    age_start = prep$data$age_start %||% NA,
    age_end = prep$data$age_end %||% NA,
    year_mid = prep$data$year_mid %||% NA,
    sample_size = prep$data$sample_size,
    cases = prep$data$cases
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(metadata, file.path(output_dir, "model_metadata.rds"))
  output_paths <- save_mcmc_outputs(meta$samples, output_dir, meta, monitors)
  output_paths
}

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  run_nimble_model()
}
