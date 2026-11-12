# CPC Tmin / Tn3x attribution -- interactive RStudio version

# 1. Load paths and helpers
source(here::here("Code", "config", "paths.R"))
library(ncdf4)
source(here::here("Code", "helpers", "gev.R"))
source(here::here("Code", "helpers", "workflow.R"))

# 2. Dataset settings
dataset_file <- cpc_tmin_file
variable_name <- "tmin"
output_dir <- cpc_tmin_output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
file.exists(dataset_file)

# 3. Open and inspect the original NetCDF file
dataset <- open_temperature_dataset(dataset_file, variable_name)
nc <- dataset$nc
metadata <- dataset$metadata
lon <- dataset$lon
lat <- dataset$lat
date <- dataset$date
sample_tmin <- dataset$sample_grid
nc
metadata
range(date)
dim(sample_tmin)
range(sample_tmin, na.rm = TRUE)
sample_tmin[1:8, 1:8]

# 4. Determine the stable land footprint
footprint_result <- find_persistent_footprint(
  nc, variable_name, date, length(lon), length(lat)
)
land_footprint <- footprint_result$included
valid_fraction <- footprint_result$valid_fraction
table(land_footprint)
summary(valid_fraction)
footprint_table <- data.frame(
  expand.grid(lon = lon, lat = lat), valid_fraction = valid_fraction,
  included = land_footprint
)
write.csv(footprint_table, file.path(output_dir, "grid_footprint.csv"), row.names = FALSE)

# 5. Calculate the daily regional Tmin
daily_result <- calculate_daily_regional_temperature(
  nc, variable_name, date, lon, lat, land_footprint
)
nc_close(nc)
head(daily_result)
summary(daily_result$value)

# 6. Calculate three-day means and Annual/June Tn3x
extreme_result <- calculate_three_day_extremes(daily_result)
daily_cpc_tmin <- extreme_result$daily
cpc_tn_extremes <- extreme_result$combined
names(daily_cpc_tmin)[names(daily_cpc_tmin) == "value"] <- "tmin"
names(daily_cpc_tmin)[names(daily_cpc_tmin) == "value3"] <- "tn3"
head(daily_cpc_tmin)
head(cpc_tn_extremes)
write.csv(daily_cpc_tmin, file.path(output_dir, "daily_regional_tmin.csv"), row.names = FALSE)

# 7. Read the shared precomputed smoothed GMST
gmst_complete <- read.csv(smoothed_gmst_file)
gmst_2026 <- gmst_complete$gmst[gmst_complete$year == 2026]
gmst_2026
tail(gmst_complete)

# 8. Join Tn3x and GMST
joined_result <- join_extremes_and_gmst(cpc_tn_extremes, gmst_complete)
analysis_data <- joined_result$all
annual_data <- joined_result$annual
june_data <- joined_result$june
nrow(annual_data)
nrow(june_data)
write.csv(analysis_data, file.path(output_dir, "annual_extremes_and_gmst.csv"), row.names = FALSE)

# 9. Fit separate nonstationary GEV models
fit_annual <- fit_gev(annual_data$X, annual_data$gmst)
fit_june <- fit_gev(june_data$X, june_data$gmst)
fit_annual
fit_june

# 10. Point estimates: Annual Tn3x uses 20 years; June uses 100 years
annual_result <- attribution(fit_annual, gmst_2026, 20)
june_result <- attribution(fit_june, gmst_2026, 100)
annual_result$index <- "Tn3x"
annual_result$season <- "Annual"
annual_result$n_years <- nrow(annual_data)
june_result$index <- "Tn3x"
june_result$season <- "June"
june_result$n_years <- nrow(june_data)
cpc_tn_results <- rbind(annual_result, june_result)
cpc_tn_results
write.csv(cpc_tn_results, file.path(output_dir, "cpc_tn_risk_ratios_PROVISIONAL.csv"), row.names = FALSE)
saveRDS(list(fits = list(Annual = fit_annual, June = fit_june),
             bootstrap = NULL, results = cpc_tn_results),
        file.path(output_dir, "fits_and_bootstrap.rds"))

# 11. Optional bootstrap; run only after checking the point estimates
annual_bootstrap <- bootstrap_attribution(
  fit_annual, annual_data, gmst_2026, 20, bootstrap_B, random_seed
)
june_bootstrap <- bootstrap_attribution(
  fit_june, june_data, gmst_2026, 100, bootstrap_B, random_seed + 1
)
cpc_tn_results <- add_bootstrap_intervals(
  cpc_tn_results, annual_bootstrap, june_bootstrap
)
cpc_tn_results

write.csv(
  cpc_tn_results,
  file.path(output_dir, "cpc_tn_risk_ratios_PROVISIONAL.csv"),
  row.names = FALSE
)

bootstrap_samples <- write_bootstrap_samples(
  annual_bootstrap,
  june_bootstrap,
  file.path(output_dir, "bootstrap_samples.csv")
)

saveRDS(list(fits = list(Annual = fit_annual, June = fit_june),
             bootstrap = list(Annual = annual_bootstrap, June = june_bootstrap),
             results = cpc_tn_results),
        file.path(output_dir, "fits_and_bootstrap.rds"))
