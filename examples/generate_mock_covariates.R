#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (!requireNamespace("readr", quietly = TRUE)) {
    install.packages("readr", repos = "https://cloud.r-project.org", dependencies = TRUE)
  }
})

set.seed(20240524)

provinces <- c(
  "北京", "天津", "河北", "山西", "内蒙古", "辽宁", "吉林", "黑龙江",
  "上海", "江苏", "浙江", "安徽", "福建", "江西", "山东", "河南",
  "湖北", "湖南", "广东", "广西", "海南", "重庆", "四川", "贵州",
  "云南", "西藏", "陕西", "甘肃", "青海", "宁夏", "新疆"
)

n_prov <- length(provinces)

make_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
}

cov_dir <- file.path("examples", "mock_covariates")
make_dir(cov_dir)

# Simulate SDI between 0.45 and 0.85 with mild spatial trend proxy via latitude rank
rank_factor <- seq_len(n_prov)
base_sdi <- 0.45 + (rank_factor - min(rank_factor)) / (max(rank_factor) - min(rank_factor)) * 0.35
noise <- rnorm(n_prov, sd = 0.02)
sdi <- pmin(pmax(base_sdi + noise, 0.35), 0.95)

sdi_tbl <- data.frame(
  province_std = provinces,
  sdi = round(sdi, 3)
)

# Simulate total population (millions) with east-west gradient
base_pop <- exp(rnorm(n_prov, mean = 2, sd = 0.6)) * 1e6
# Upscale coastal provinces by factor
coastal <- provinces %in% c("上海", "江苏", "浙江", "福建", "广东", "山东", "辽宁", "天津", "海南")
base_pop[coastal] <- base_pop[coastal] * runif(sum(coastal), 1.2, 1.6)
# Ensure Tibet, Qinghai smaller
small <- provinces %in% c("西藏", "青海", "宁夏")
base_pop[small] <- base_pop[small] * runif(sum(small), 0.3, 0.5)

population <- round(base_pop)
total_population <- sum(population)
weight <- population / total_population

# Age structure shares for 5-year groups (10-14 ... 80+) using normalized gamma draws
groups <- c("10_14", "15_19", "20_24", "25_29", "30_34", "35_39",
            "40_44", "45_49", "50_54", "55_59", "60_64", "65_69", "70_74", "75_79", "80_plus")

age_matrix <- replicate(length(groups), rgamma(n_prov, shape = 2, rate = 1))
age_shares <- age_matrix / rowSums(age_matrix)
age_df <- as.data.frame(age_shares)
colnames(age_df) <- paste0("share_", groups)
age_df[] <- lapply(age_df, function(col) round(col, 4))

pop_tbl <- data.frame(
  province_std = provinces,
  population = population,
  weight = round(weight, 5)
)
pop_age_tbl <- cbind(pop_tbl["province_std"], age_df)

# Household urban share + GDP per capita as extra optional covariates
urban_share <- pmin(pmax(rbeta(n_prov, 5, 3), 0.2), 0.9)
gdp_pc <- round(rlnorm(n_prov, meanlog = log(60000), sdlog = 0.35))
extra_tbl <- data.frame(
  province_std = provinces,
  urban_share = round(urban_share, 3),
  gdp_per_capita = gdp_pc
)

readr::write_csv(sdi_tbl, file.path(cov_dir, "province_sdi.csv"))
readr::write_csv(pop_tbl, file.path(cov_dir, "province_population_2020.csv"))
readr::write_csv(pop_age_tbl, file.path(cov_dir, "province_age_structure.csv"))
readr::write_csv(extra_tbl, file.path(cov_dir, "province_extra_covariates.csv"))

cat("Mock covariate files written to", cov_dir, "\n")
