# ERA5 Tx3x analysis -- interactive RStudio version
# Run one section at a time. This script deliberately keeps intermediate
# objects so you can inspect them before continuing.

# ============================================================================
# 1. Paths and settings
# ============================================================================

source(here::here("Code", "config", "paths.R"))

era5_file <- era5_tmax_file
land_mask_file <- era5_land_mask_file
output_dir <- era5_tx_output_dir

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
file.exists(era5_file)
file.exists(land_mask_file)
file.exists(gmst_file)

# ============================================================================
# 2. Packages and GEV functions
# ============================================================================

library(ncdf4)
source(here::here("Code", "helpers", "gev.R"))
source(here::here("Code", "helpers", "workflow.R"))

# ============================================================================
# 3. Open and inspect the original ERA5 NetCDF file
# ============================================================================

nc <- nc_open(era5_file)
#nc
names(nc$var)
names(nc$dim)

tmax_metadata <- nc$var$tmax
#tmax_metadata

dimension_names <- vapply(tmax_metadata$dim, function(x) x$name, character(1))
temperature_units <- tmax_metadata$units
#dimension_names [1] "lon"  "lat"  "time"
#temperature_units [1] "K"

stopifnot(identical(dimension_names, c("lon", "lat", "time")))
stopifnot(temperature_units %in% c("K", "kelvin", "Kelvin"))
capture.output(nc, file = file.path(output_dir, "netcdf_metadata.txt"))

# ============================================================================
# 4. Read coordinates and dates
# ============================================================================

lon_original <- nc$dim$lon$vals
lon <- (lon_original + 180) %% 360 - 180
lat <- nc$dim$lat$vals

time_values <- nc$dim$time$vals
time_units <- nc$dim$time$units
time_origin <- as.Date(substr(sub("^days since ", "", time_units), 1, 10))
date <- time_origin + time_values

# length(lon)
# range(lon)
# length(lat)
# range(lat)
# length(date) [1] 28002
# range(date) [1] "1950-01-01" "2026-08-31"

stopifnot(!anyDuplicated(date))
stopifnot(all(diff(date) > 0))

# ============================================================================
# 5. Read one raw day before doing any averaging
# ============================================================================

sample_day_number <- 1L
sample_date <- date[sample_day_number]

sample_tmax_K <- ncvar_get(
  nc,
  "tmax",
  start = c(1, 1, sample_day_number),
  count = c(-1, -1, 1),
  collapse_degen = FALSE
)

sample_tmax_K <- sample_tmax_K[, , 1]
sample_tmax_C <- sample_tmax_K - 273.15

# sample_date
# dim(sample_tmax_C)
# range(sample_tmax_C, na.rm = TRUE)
# sample_tmax_C[1:8, 1:8]

# ============================================================================
# 6. Read the approximate land mask
# ============================================================================

# ERA5 has valid values over both land and ocean. We use the previously saved
# E-Obs land footprint as an approximate European land mask.

mask_table <- read.csv(land_mask_file)
names(mask_table)
head(mask_table)
table(mask_table$included, useNA = "ifany")

mask_lon <- unique(mask_table$lon)
mask_lat <- unique(mask_table$lat)

mask_grid <- matrix(
  as.logical(mask_table$included),
  nrow = length(mask_lon),
  ncol = length(mask_lat)
)

nearest_mask_lon <- vapply(
  lon,
  function(x) which.min(abs(mask_lon - x)),
  integer(1)
)

nearest_mask_lat <- vapply(
  lat,
  function(x) which.min(abs(mask_lat - x)),
  integer(1)
)

land_mask <- as.vector(
  mask_grid[nearest_mask_lon, nearest_mask_lat, drop = FALSE]
)

table(land_mask, useNA = "ifany")

# ============================================================================
# 7. Define grid-cell area weights
# ============================================================================

latitude_weights <- cos(lat * pi / 180)
grid_weights <- rep(latitude_weights, each = length(lon))
grid_weights[!land_mask] <- 0

length(grid_weights)
sum(grid_weights > 0)
range(grid_weights[grid_weights > 0])

# ============================================================================
# 8. Read ERA5 in one-year chunks and calculate daily regional Tmax
# ============================================================================

# The file is too large to read comfortably in one operation. Each iteration
# reads at most 366 days. You can stop after a completed chunk and inspect the
# values already stored in daily_tmax.

number_of_days <- length(date)
chunk_starts <- seq(1L, number_of_days, by = 366L)
daily_tmax <- rep(NA_real_, number_of_days)
spatial_coverage <- rep(NA_real_, number_of_days)

for (chunk_start in chunk_starts) {
  chunk_end <- min(chunk_start + 365L, number_of_days)
  chunk_days <- chunk_start:chunk_end

  temperature_K <- ncvar_get(
    nc,
    "tmax",
    start = c(1, 1, chunk_start),
    count = c(-1, -1, length(chunk_days)),
    collapse_degen = FALSE
  )

  temperature_C <- temperature_K - 273.15
  temperature_matrix <- matrix(temperature_C, nrow = length(grid_weights))
  valid_value <- is.finite(temperature_matrix)
  weighted_denominator <- colSums(valid_value * grid_weights)

  temperature_matrix[!valid_value] <- 0

  daily_tmax[chunk_days] <- colSums(
    temperature_matrix * grid_weights
  ) / weighted_denominator

  spatial_coverage[chunk_days] <- weighted_denominator / sum(grid_weights)

  message("Finished ", date[chunk_start], " to ", date[chunk_end])
}

nc_close(nc)
daily_tmax[spatial_coverage < 0.95] <- NA_real_

range(daily_tmax, na.rm = TRUE)
#summary(spatial_coverage)

# ============================================================================
# 9. Inspect and save the daily regional series
# ============================================================================

daily_era5 <- data.frame(
  date = date,
  tmax = daily_tmax,
  spatial_coverage = spatial_coverage
)

daily_era5$year <- as.integer(format(daily_era5$date, "%Y"))
daily_era5$month <- as.integer(format(daily_era5$date, "%m"))

head(daily_era5)
tail(daily_era5)
summary(daily_era5$tmax)

plot(
  daily_era5$date,
  daily_era5$tmax,
  type = "l",
  xlab = "Date",
  ylab = "Regional daily Tmax (degC)"
)

# write.csv(
#   daily_era5,
#   file.path(output_dir, "daily_regional_tmax.csv"),
#   row.names = FALSE
# )

saveRDS(daily_era5, file.path(output_dir, "daily_regional_tmax.rds"))

# ============================================================================
# 10. Calculate the trailing three-day mean
# ============================================================================

daily_era5$tx3 <- as.numeric(
  stats::filter(daily_era5$tmax, rep(1 / 3, 3), sides = 1)
)

head(daily_era5, 10)
summary(daily_era5$tx3)

# ============================================================================
# 11. Extract Annual Tx3x and June Tx3x
# ============================================================================

years <- sort(unique(daily_era5$year))
annual_extremes <- data.frame()
june_extremes <- data.frame()

for (selected_year in years) {
  annual_start <- as.Date(sprintf("%d-01-01", selected_year))
  annual_end <- as.Date(sprintf("%d-12-31", selected_year))

  annual_block <- daily_era5[
    daily_era5$date >= annual_start & daily_era5$date <= annual_end,
  ]

  expected_days <- as.integer(annual_end - annual_start) + 1L
  annual_complete <- nrow(annual_block) == expected_days &&
    mean(is.finite(annual_block$tmax)) >= 0.95

  annual_values <- annual_block$tx3[annual_block$date >= annual_start + 2]
  annual_tx3x <- if (annual_complete) max(annual_values, na.rm = TRUE) else NA_real_

  annual_extremes <- rbind(
    annual_extremes,
    data.frame(
      year = selected_year,
      season = "Annual",
      complete = annual_complete,
      X = annual_tx3x
    )
  )

  june_start <- as.Date(sprintf("%d-06-01", selected_year))
  june_end <- as.Date(sprintf("%d-06-30", selected_year))

  june_block <- daily_era5[
    daily_era5$date >= june_start & daily_era5$date <= june_end,
  ]

  june_complete <- nrow(june_block) == 30L &&
    mean(is.finite(june_block$tmax)) >= 0.95

  # Start on 3 June so all three days used by the rolling mean lie in June.
  june_values <- june_block$tx3[june_block$date >= june_start + 2]
  june_tx3x <- if (june_complete) max(june_values, na.rm = TRUE) else NA_real_

  june_extremes <- rbind(
    june_extremes,
    data.frame(
      year = selected_year,
      season = "June",
      complete = june_complete,
      X = june_tx3x
    )
  )
}

era5_extremes <- rbind(annual_extremes, june_extremes)
head(era5_extremes)
tail(era5_extremes)
table(era5_extremes$season, era5_extremes$complete)

# ============================================================================
# 12. Read the shared precomputed smoothed GMST
# ============================================================================

stopifnot(file.exists(smoothed_gmst_file))
gmst_complete <- read.csv(smoothed_gmst_file)
gmst_2026 <- gmst_complete$gmst[gmst_complete$year == 2026]

stopifnot(length(gmst_2026) == 1, is.finite(gmst_2026))
gmst_2026
tail(gmst_complete)

# ============================================================================
# 13. Join Tx3x and GMST, then select fitting years
# ============================================================================

analysis_data <- merge(
  era5_extremes,
  gmst_complete[, c("year", "gmst")],
  by = "year",
  all.x = TRUE
)

analysis_data$used_in_fit <- analysis_data$complete &
  analysis_data$year <= 2025 &
  is.finite(analysis_data$gmst)

annual_data <- analysis_data[
  analysis_data$season == "Annual" & analysis_data$used_in_fit,
]

june_data <- analysis_data[
  analysis_data$season == "June" & analysis_data$used_in_fit,
]

range(annual_data$year)
nrow(annual_data)
range(june_data$year)
nrow(june_data)

write.csv(
  analysis_data,
  file.path(output_dir, "annual_extremes_and_gmst.csv"),
  row.names = FALSE
)

# ============================================================================
# 14. Fit the nonstationary GEV models
# ============================================================================

fit_annual <- fit_gev(annual_data$X, annual_data$gmst)
fit_annual

fit_june <- fit_gev(june_data$X, june_data$gmst)
fit_june

# ============================================================================
# 15. Calculate point estimates
# ============================================================================

# Table 2: Annual Tx3x uses 25 years; June Tx3x uses 100 years.
annual_result <- attribution(fit_annual, gmst_2026, 25)
annual_result$index <- "Tx3x"
annual_result$season <- "Annual"
annual_result$n_years <- nrow(annual_data)

june_result <- attribution(fit_june, gmst_2026, 100)
june_result$index <- "Tx3x"
june_result$season <- "June"
june_result$n_years <- nrow(june_data)

annual_result
june_result

era5_results <- rbind(annual_result, june_result)
era5_results


write.csv(
  era5_results,
  file.path(output_dir, "era5_tx_risk_ratios_PROVISIONAL.csv"),
  row.names = FALSE
)
saveRDS(
  list(fits = list(Annual = fit_annual, June = fit_june),
       bootstrap = NULL, results = era5_results),
  file.path(output_dir, "fits_and_bootstrap.rds")
)

# ============================================================================
# 16. Optional bootstrap confidence intervals
# ============================================================================

# This is the slow section. Inspect the point estimates above before running it.
# Reduce bootstrap_B at the top while testing, for example bootstrap_B <- 20L.

annual_bootstrap <- bootstrap_attribution(
  fit_annual,
  annual_data,
  gmst_2026,
  25,
  bootstrap_B,
  random_seed
)

june_bootstrap <- bootstrap_attribution(
  fit_june,
  june_data,
  gmst_2026,
  100,
  bootstrap_B,
  random_seed + 1
)

annual_bootstrap_success <- sum(colSums(is.na(annual_bootstrap)) == 0)
june_bootstrap_success <- sum(colSums(is.na(june_bootstrap)) == 0)

annual_bootstrap_success
june_bootstrap_success

era5_results <- add_bootstrap_intervals(
  era5_results, annual_bootstrap, june_bootstrap
)
era5_results

write.csv(
  era5_results,
  file.path(output_dir, "era5_tx_risk_ratios_PROVISIONAL.csv"),
  row.names = FALSE
)

bootstrap_samples <- write_bootstrap_samples(
  annual_bootstrap,
  june_bootstrap,
  file.path(output_dir, "bootstrap_samples.csv")
)

saveRDS(
  list(
    fit_annual = fit_annual,
    fit_june = fit_june,
    annual_bootstrap = annual_bootstrap,
    june_bootstrap = june_bootstrap,
    point_results = era5_results
  ),
  file.path(output_dir, "era5_fits_and_bootstrap.rds")
)
saveRDS(
  list(fits = list(Annual = fit_annual, June = fit_june),
       bootstrap = list(Annual = annual_bootstrap, June = june_bootstrap),
       results = era5_results),
  file.path(output_dir, "fits_and_bootstrap.rds")
)
