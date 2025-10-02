# Eating Disorder Bayesian Pipeline

This repository contains a modular R pipeline for cleaning eating disorder surveillance data, constructing age bases, fitting a hierarchical NIMBLE model, and generating summaries and graphics.

## Repository structure

```
.
├── 01_cleaning.R                # ingest + clean the source workbook
├── 02_build_age_basis.R         # transform optional age weights + construct B-spline basis
├── 03_model_nimble.R            # specify and run the hierarchical NIMBLE model
├── 04_postprocess_plots.R       # summarise posterior draws + create visualisations
├── ed_pipeline_nimble.R         # orchestration script with configurable parameters
├── ED_meta_template_v5.xlsx     # template workbook expected by the cleaner
├── examples/
│   ├── mock_data_long.csv       # small long-format example for quick testing
│   └── mock_ageweights.csv      # optional example age-weight template
└── README.md
```

Running the full pipeline will produce two main output directories:

* `clean/` – cleaned long-format data, diagnostics, and age-basis artefacts.
* `ed_bayes_outputs/` – model outputs, posterior draws, derived summaries, and plots.

These directories are ignored in version control so that large derived files do not bloat the repository.

## Software requirements

* R (≥ 4.1 recommended).
* Internet access for first-run package installation (packages are auto-installed by each script).

No additional system dependencies are required.

## Core scripts

### 01_cleaning.R
* Reads the workbook (default `ED_meta_template_v5.xlsx`) and searches for a long-format sheet (prefers `Data_Long`).
* Standardises province names, derives mid-year (`year_mid`), harmonises prevalence units to proportions, back-calculates missing case counts, and enforces basic validity checks.
* Drops duplicate mixed-sex rows, assigns `row_id` and `study_idx`, and exports `clean/ed_clean_long.csv` with diagnostics (extreme values, invalid rows, duplicates).

Run standalone:

```bash
Rscript 01_cleaning.R path/to/workbook.xlsx optional/output/dir
```

### 02_build_age_basis.R
* Consumes the cleaned CSV and an optional age-weight workbook or CSV (default `ED_ageweights_template.xlsx`).
* Aligns supplied weights to a fixed 5-year grid (10–14 … 80+), falls back to uniform weights with logged warnings, and constructs B-spline bases.
* Saves observation-specific weight matrices, spline bases, and supporting metadata under `clean/age_basis/`.

Run standalone:

```bash
Rscript 02_build_age_basis.R clean/ed_clean_long.csv optional/age_weight_file.xlsx
```

### 03_model_nimble.R
* Builds a hierarchical logistic model with sex-specific logits, method/bias covariates, province and study random effects, and spline effects.
* Handles mixed-sex observations using sample-size weighting, runs multi-chain MCMC with automatic extension until convergence targets (R̂ ≤ 1.1, ESS ≥ 200) are met or the maximum iteration cap is reached.
* Writes posterior draws (`posterior_samples.rds`), convergence diagnostics (`mcmc_diag.txt`), and trace/density plots under `ed_bayes_outputs/model/`. A metadata RDS captures design matrices used during fitting.

Run standalone (after cleaning + age basis):

```bash
Rscript 03_model_nimble.R
```

### 04_postprocess_plots.R
* Combines posterior draws and metadata to generate observation-level posterior summaries, province/national aggregates, time trends (2000–2025), age curves (10–80), sex comparisons, and sensitivity summaries.
* Exports CSV summaries to `ed_bayes_outputs/summaries/` and ggplot visualisations (province bars, national trends, age curves, sex contrasts, sensitivity comparisons, convergence diagnostics) to `ed_bayes_outputs/plots/`.

Run standalone (after modelling):

```bash
Rscript 04_postprocess_plots.R
```

## Orchestrated pipeline

`ed_pipeline_nimble.R` ties the modules together, exposing command-line overrides for key settings:

* `input_path` – workbook to clean (`ED_meta_template_v5.xlsx` by default).
* `age_weight_path` – optional age-weight file (Excel or CSV interpreted by `02_build_age_basis.R`).
* `spline_df` – degrees of freedom for the B-spline basis (default 6).
* `mcmc.*` – NIMBLE sampler settings (`chains`, `iter`, `burnin`, `thin`, `max_iter`, `ess_threshold`, `rhat_threshold`).
* `sensitivity.*` – baseline toggles passed to the model (`age_cap`, `sex_equal`, `exclude_bias`).
* `sensitivity_runs.*` – set to `TRUE` to run additional sensitivity re-fits saved under `ed_bayes_outputs/sensitivity_<name>/`.

Example:

```bash
Rscript ed_pipeline_nimble.R input_path=ED_meta_template_v5.xlsx spline_df=8 mcmc.iter=8000 sensitivity_runs.age_cap=TRUE
```

The orchestrator always captures `sessionInfo()` in `sessionInfo.txt` once the workflow finishes.

## Minimal example

For a lightweight smoke test, use the files under `examples/`:

1. Copy `examples/mock_data_long.csv` to a workbook containing a sheet named `Data_Long` (the cleaner expects Excel input). Any spreadsheet editor can import the CSV and save as `.xlsx`.
2. Optionally supply `examples/mock_ageweights.csv` to `02_build_age_basis.R` via `age_weight_path`.
3. Run the pipeline with:

```bash
Rscript ed_pipeline_nimble.R input_path=examples/mock_data_long.csv age_weight_path=examples/mock_ageweights.csv
```

The cleaning script can ingest CSV for the mock example; production runs should continue to use the official template workbook.

## Logging & diagnostics

* Cleaning diagnostics live in `clean/diagnostics/` (`invalid_rows.csv`, `extreme_prevalence_*`, etc.).
* Age-basis warnings are written to `clean/diagnostics/age_basis_log.txt` and `clean/age_basis/age_basis_warnings.txt`.
* Model diagnostics include `mcmc_diag.txt` plus trace/density plots under `ed_bayes_outputs/model/diagnostics/`.
* Post-processing writes aggregated summaries and sensitivity comparisons for downstream reporting.

## Re-running components

All scripts are idempotent – rerunning them will overwrite previous outputs in their respective folders. This makes it straightforward to tweak parameters and regenerate artefacts without manual cleanup.

## Questions

If you encounter missing dependencies or want to add new sensitivity scenarios, extend the `sensitivity` list within `ed_pipeline_nimble.R` and provide corresponding flags in the orchestrator call.
