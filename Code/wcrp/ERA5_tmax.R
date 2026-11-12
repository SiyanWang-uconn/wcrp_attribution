library(dplyr)

# ============================================================================
# 1. Paths and settings
# ============================================================================
project_root <- "C:/Users/siyan/WCRP/wcrp_common_attribution_assessment"
path_path <- file.path(project_root,"Code", "config", "paths.R")
output_dir <- file.path(project_root,"Results","ERA5")


source(path_path)
source(file.path(project_root,"Code", "helpers", "gev.R"))
source(file.path(project_root,"Code", "helpers", "workflow.R"))


era5_file <- era5_tmax_file

# output_dir <- cpc_tmax_output_dir
# 
# dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

file.exists(era5_file)
# file.exists(era5_file)

#file.exists(gmst_file)

# ============================================================================
# 2.Check variables
# ============================================================================

library(ncdf4)
# source(here::here("Code", "helpers", "gev.R"))
# source(here::here("Code", "helpers", "workflow.R"))
nc <- nc_open(era5_file)

nc
var_names<-names(nc$var)
names(nc$dim)

era5_dat <- as.data.frame(setNames(
  lapply(var_names, function(v) as.vector(ncvar_get(nc, varid = v))),
  var_names
))

date_units <- ncatt_get(nc, "tmax_date", "units")$value 


nc_close(nc) ## CLose ncdf file

## Processing

era5_dat <- era5_dat %>%
  mutate(tmax_date = as.Date(
    tmax_date, origin = sub("^days since ", "", date_units)
  ),
  year = as.integer(format(tmax_date, "%Y")) )


era5_dat <- era5_dat %>%
  mutate(complete=expected_days==valid_days)%>%
  select(-c(expected_days,valid_days))

# ============================================================================
# 3. Read the shared precomputed smoothed GMST
# ============================================================================

stopifnot(file.exists(smoothed_gmst_file))
gmst_complete <- read.csv(smoothed_gmst_file)
gmst_2021 <- gmst_complete$gmst[gmst_complete$year == 2021]

stopifnot(length(gmst_2021) == 1, is.finite(gmst_2021))
gmst_2021

gmst_1951 <- gmst_complete$gmst[gmst_complete$year == 1951]
gmst_1951
#tail(gmst_complete)

# ============================================================================
# 4. Join Tx3x and GMST, then select fitting years
# ============================================================================

analysis_data <- merge(
  era5_dat,
  gmst_complete[, c("year", "gmst")],
  by = "year",
  all.x = TRUE
)

analysis_data$used_in_fit <- analysis_data$complete &
  analysis_data$year <= 2025 &
  is.finite(analysis_data$gmst)

annual_data <- analysis_data[ analysis_data$used_in_fit,]



range(annual_data$year)
nrow(annual_data)


write.csv(
  analysis_data,
  file.path(output_dir, "annual_extremes_and_gmst.csv"),
  row.names = FALSE
)

# ============================================================================
# 5. Fit the nonstationary GEV models
# ============================================================================

fit_annual <- fit_gev(annual_data$tmax, annual_data$gmst)
fit_annual



# ============================================================================
# 6. Calculate Risk Ratio and intensity point estimates
# ============================================================================

# Table 2 uses 25 years for Annual Tx3x and 100 years for June Tx3x.
annual_result <- attribution(fit_annual, gmst_2021,gmst_1951, annual_data)
annual_result$season <- "Annual"
annual_result$n_years <- nrow(annual_data)

# june_result <- attribution(fit_june, gmst_2026, 100)
# june_result$index <- "Tx3x"
# june_result$season <- "June"
# june_result$n_years <- nrow(june_data)

annual_result



write.csv(
  annual_result,
  file.path(output_dir, "era5_tx_risk_ratios_PROVISIONAL.csv"),
  row.names = FALSE
)
# saveRDS(
#   list(fits = list(Annual = fit_annual, June = fit_june),
#        bootstrap = NULL, results = cpc_results),
#   file.path(output_dir, "fits_and_bootstrap.rds")
# )

# ============================================================================
# 7. Confidence intervals
# ============================================================================

# This is the slow section. Inspect the point results first. While testing,
# change bootstrap_B near the top to a small number such as 20.


### 7.1 Optional bootstrap 
bootstrap_B <- 1000L
annual_bootstrap <- bootstrap_attribution(
  fit_annual,
  annual_data,
  gmst_2021,
  gmst_1951,
  bootstrap_B
)

annual_bootstrap_success <- sum(colSums(is.na(annual_bootstrap)) == 0)


annual_bootstrap_success


era5_results <- add_bootstrap_intervals(
  annual_result, annual_bootstrap
)
era5_results

write.csv(
  era5_results,
  file.path(output_dir, "era5_tx_risk_ratios_bootstrap.csv"),
  row.names = FALSE
)
### 7.2 Delta method

delta_ci <- rr_ci_delta(
  fit_annual, annual_data, gmst_2021, gmst_1951,
  point_results = annual_result
)

### LRT method
lrt_ci <- rr_ci_lrt(
  fit_annual, annual_data, gmst_2021, gmst_1951
)
