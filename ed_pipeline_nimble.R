#!/usr/bin/env Rscript
set.seed(123)

update_config_value <- function(config, key, value) {
  parts <- strsplit(key, "\\.")[[1]]
  key <- parts[1]
  if (length(parts) == 1) {
    config[[key]] <- value
    return(config)
  }
  head <- key
  tail <- paste(parts[-1], collapse = ".")
  if (is.null(config[[head]]) || !is.list(config[[head]])) {
    config[[head]] <- list()
  }
  config[[head]] <- update_config_value(config[[head]], tail, value)
  config
}

coerce_value <- function(value) {
  if (value %in% c("TRUE", "FALSE")) {
    return(as.logical(value))
  }
  num <- suppressWarnings(as.numeric(value))
  if (!is.na(num)) {
    return(num)
  }
  value
}

parse_args <- function(default_config) {
  args <- commandArgs(trailingOnly = TRUE)
  config <- default_config
  for (arg in args) {
    if (!stringr::str_detect(arg, "=")) next
    key <- sub("=.*$", "", arg)
    value <- sub("^[^=]+=", "", arg)
    config <- update_config_value(config, key, coerce_value(value))
  }
  config
}

run_pipeline <- function(config) {
  source("01_cleaning.R")
  source("02_build_age_basis.R")
  source("03_model_nimble.R")
  source("04_postprocess_plots.R")
  message("Starting cleaning...")
  clean_info <- run_cleaning(input_path = config$input_path, output_dir = config$clean_dir)
  message("Building age basis...")
  age_info <- build_age_basis(clean_file = clean_info$clean_file,
                              age_weight_path = config$age_weight_path,
                              output_dir = config$clean_dir,
                              spline_df = config$spline_df)
  run_model <- function(sens, output_path) {
    run_nimble_model(clean_file = clean_info$clean_file,
                     basis_file = file.path(config$clean_dir, "age_basis", "observation_bspline_basis.rds"),
                     weight_file = file.path(config$clean_dir, "age_basis", "observation_age_weights.rds"),
                     output_dir = output_path,
                     sensitivity = sens,
                     mcmc = config$mcmc)
  }
  message("Running baseline model...")
  model_dir <- file.path(config$output_dir, "model")
  run_model(config$sensitivity, model_dir)
  if (!is.null(config$sensitivity_runs) && length(config$sensitivity_runs) > 0) {
    for (name in names(config$sensitivity_runs)) {
      if (!isTRUE(config$sensitivity_runs[[name]])) next
      message(glue::glue("Running sensitivity: {name}"))
      sens_list <- config$sensitivity
      sens_list[[name]] <- TRUE
      sens_dir <- file.path(config$output_dir, glue::glue("sensitivity_{name}"))
      run_model(sens_list, sens_dir)
    }
  }
  message("Post-processing results...")
  run_postprocess(posterior_file = file.path(model_dir, "posterior_samples.rds"),
                  metadata_file = file.path(model_dir, "model_metadata.rds"),
                  clean_file = clean_info$clean_file,
                  age_basis_dir = file.path(config$clean_dir, "age_basis"),
                  output_dir = config$output_dir)
  if (!is.null(config$session_info)) {
    capture.output(sessionInfo(), file = config$session_info)
  }
  message("Pipeline complete.")
}

default_config <- list(
  input_path = "ED_meta_template_v5.xlsx",
  age_weight_path = NULL,
  clean_dir = "clean",
  output_dir = "ed_bayes_outputs",
  spline_df = 6,
  mcmc = list(chains = 4, iter = 6000, burnin = 2000, thin = 5, max_iter = 20000, ess_threshold = 200, rhat_threshold = 1.1),
  sensitivity = list(age_cap = FALSE, sex_equal = FALSE, exclude_bias = FALSE),
  sensitivity_runs = list(age_cap = FALSE, sex_equal = FALSE, exclude_bias = FALSE),
  session_info = "sessionInfo.txt"
)

if (identical(environment(), globalenv()) && !length(sys.calls())) {
  ensure_packages <- function(pkgs) {
    to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
    if (length(to_install) > 0) {
      install.packages(to_install, repos = "https://cloud.r-project.org", dependencies = TRUE)
    }
    invisible(lapply(pkgs, library, character.only = TRUE))
  }
  ensure_packages(c("stringr", "glue"))
  config <- parse_args(default_config)
  run_pipeline(config)
}
