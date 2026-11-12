# CPC Tmax / Tx3x attribution -- interactive RStudio version
#
# Run one numbered section at a time. Intermediate objects stay in the Global
# Environment so that you can inspect, plot, or save them before continuing.
# There are no command-line arguments, environment variables, or wrappers.

# ============================================================================
# 1. Paths and settings
# ============================================================================


source(here::here("Code", "config", "paths.R"))

cpc_file <- cpc_tmax_file
output_dir <- cpc_tmax_output_dir

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

file.exists(cpc_file)
file.exists(gmst_file)

# ============================================================================
# 2. Packages and GEV functions
# ============================================================================

library(ncdf4)
source(here::here("Code", "helpers", "gev.R"))
source(here::here("Code", "helpers", "workflow.R"))

# ============================================================================
# 3. Open and inspect the original CPC NetCDF file
# ============================================================================

nc <- nc_open(cpc_file)

nc
names(nc$var)
names(nc$dim)

tmax_metadata <- nc$var$tmax
tmax_metadata

dimension_names <- vapply(tmax_metadata$dim, function(x) x$name, character(1))
temperature_units <- tmax_metadata$units

dimension_names
temperature_units

stopifnot(identical(dimension_names, c("lon", "lat", "time")))
stopifnot(temperature_units %in% c("degC", "C", "degree C", "Celsius"))

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

length(lon)
range(lon)
length(lat)
range(lat)
length(date)
range(date)

stopifnot(!anyDuplicated(date))
stopifnot(all(diff(date) > 0))

# ============================================================================
# 5. Read and inspect one raw day
# ============================================================================

sample_day_number <- 1L
sample_date <- date[sample_day_number]

sample_tmax <- ncvar_get(
  nc,
  "tmax",
  start = c(1, 1, sample_day_number),
  count = c(-1, -1, 1),
  collapse_degen = FALSE
)

sample_tmax <- sample_tmax[, , 1]

sample_date
dim(sample_tmax)
range(sample_tmax, na.rm = TRUE)
sample_tmax[1:8, 1:8]
sum(is.finite(sample_tmax))
sum(!is.finite(sample_tmax))

# Convert the sample grid to a normal table for inspection.
sample_grid <- data.frame(
  date = sample_date,
  expand.grid(lon = lon, lat = lat),
  tmax = as.vector(sample_tmax)
)

head(sample_grid)
head(sample_grid[is.finite(sample_grid$tmax), ])

# ============================================================================
# 6. Identify the persistent CPC land footprint
# ============================================================================

# CPC ocean cells and unavailable cells are stored as missing values. We retain
# cells that contain a finite value on at least 95% of all dates. This creates a
# stable spatial footprint rather than changing the region from day to day.

number_of_days <- length(date)
number_of_grid_cells <- length(lon) * length(lat)
chunk_starts <- seq(1L, number_of_days, by = 366L)

valid_day_count <- rep(0, number_of_grid_cells)

for (chunk_start in chunk_starts) {
  chunk_end <- min(chunk_start + 365L, number_of_days)
  chunk_days <- chunk_start:chunk_end

  temperature_chunk <- ncvar_get(
    nc,
    "tmax",
    start = c(1, 1, chunk_start),
    count = c(-1, -1, length(chunk_days)),
    collapse_degen = FALSE
  )

  temperature_matrix <- matrix(
    temperature_chunk,
    nrow = number_of_grid_cells
  )

  valid_day_count <- valid_day_count + rowSums(is.finite(temperature_matrix))

  message("Footprint: finished ", date[chunk_start], " to ", date[chunk_end])
}

valid_fraction <- valid_day_count / number_of_days
land_footprint <- valid_fraction >= 0.95

summary(valid_fraction)
table(land_footprint)

footprint_table <- data.frame(
  expand.grid(lon = lon, lat = lat),
  valid_fraction = valid_fraction,
  included = land_footprint
)

head(footprint_table)

write.csv(
  footprint_table,
  file.path(output_dir, "grid_footprint.csv"),
  row.names = FALSE
)

# ============================================================================
# 7. Define grid-cell area weights
# ============================================================================

# A regular longitude-latitude cell represents less area at higher latitudes.
# The approximate area weight is cos(latitude).

latitude_weights <- cos(lat * pi / 180)
grid_weights <- rep(latitude_weights, each = length(lon))
grid_weights[!land_footprint] <- 0

length(grid_weights)
sum(grid_weights > 0)
range(grid_weights[grid_weights > 0])

# ============================================================================
# 8. Calculate the daily area-weighted regional Tmax
# ============================================================================

daily_tmax <- rep(NA_real_, number_of_days)
spatial_coverage <- rep(NA_real_, number_of_days)

for (chunk_start in chunk_starts) {
  chunk_end <- min(chunk_start + 365L, number_of_days)
  chunk_days <- chunk_start:chunk_end

  temperature_chunk <- ncvar_get(
    nc,
    "tmax",
    start = c(1, 1, chunk_start),
    count = c(-1, -1, length(chunk_days)),
    collapse_degen = FALSE
  )

  temperature_matrix <- matrix(
    temperature_chunk,
    nrow = number_of_grid_cells
  )

  valid_value <- is.finite(temperature_matrix)
  weighted_denominator <- colSums(valid_value * grid_weights)

  temperature_matrix[!valid_value] <- 0

  daily_tmax[chunk_days] <- colSums(
    temperature_matrix * grid_weights
  ) / weighted_denominator

  spatial_coverage[chunk_days] <- weighted_denominator / sum(grid_weights)

  message("Regional mean: finished ", date[chunk_start], " to ", date[chunk_end])
}

nc_close(nc)

daily_tmax[spatial_coverage < 0.95] <- NA_real_

range(daily_tmax, na.rm = TRUE)
summary(spatial_coverage)

# ============================================================================
# 9. Inspect and save the daily regional series
# ============================================================================

daily_cpc <- data.frame(
  date = date,
  tmax = daily_tmax,
  spatial_coverage = spatial_coverage
)

daily_cpc$year <- as.integer(format(daily_cpc$date, "%Y"))
daily_cpc$month <- as.integer(format(daily_cpc$date, "%m"))

head(daily_cpc)
tail(daily_cpc)
summary(daily_cpc$tmax)

plot(
  daily_cpc$date,
  daily_cpc$tmax,
  type = "l",
  xlab = "Date",
  ylab = "Regional daily Tmax (degC)"
)

write.csv(
  daily_cpc,
  file.path(output_dir, "daily_regional_tmax.csv"),
  row.names = FALSE
)

saveRDS(daily_cpc, file.path(output_dir, "daily_regional_tmax.rds"))

# ============================================================================
# 10. Calculate the trailing three-day mean
# ============================================================================

daily_cpc$tx3 <- as.numeric(
  stats::filter(daily_cpc$tmax, rep(1 / 3, 3), sides = 1)
)

head(daily_cpc, 10)
summary(daily_cpc$tx3)

# ============================================================================
# 11. Extract Annual Tx3x and June Tx3x
# ============================================================================

years <- sort(unique(daily_cpc$year))
annual_extremes <- data.frame()
june_extremes <- data.frame()

for (selected_year in years) {
  annual_start <- as.Date(sprintf("%d-01-01", selected_year))
  annual_end <- as.Date(sprintf("%d-12-31", selected_year))

  annual_block <- daily_cpc[
    daily_cpc$date >= annual_start & daily_cpc$date <= annual_end,
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

  june_block <- daily_cpc[
    daily_cpc$date >= june_start & daily_cpc$date <= june_end,
  ]

  june_complete <- nrow(june_block) == 30L &&
    mean(is.finite(june_block$tmax)) >= 0.95

  # Start on 3 June so every three-day window lies fully within June.
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

cpc_extremes <- rbind(annual_extremes, june_extremes)

head(cpc_extremes)
tail(cpc_extremes)
table(cpc_extremes$season, cpc_extremes$complete)

plot(
  annual_extremes$year,
  annual_extremes$X,
  type = "b",
  xlab = "Year",
  ylab = "Annual Tx3x (degC)"
)

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
  cpc_extremes,
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
# 15. Calculate Risk Ratio and intensity point estimates
# ============================================================================

# Table 2 uses 25 years for Annual Tx3x and 100 years for June Tx3x.
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

cpc_results <- rbind(annual_result, june_result)
cpc_results


write.csv(
  cpc_results,
  file.path(output_dir, "cpc_tx_risk_ratios_PROVISIONAL.csv"),
  row.names = FALSE
)
saveRDS(
  list(fits = list(Annual = fit_annual, June = fit_june),
       bootstrap = NULL, results = cpc_results),
  file.path(output_dir, "fits_and_bootstrap.rds")
)

# ============================================================================
# 16. Optional bootstrap confidence intervals
# ============================================================================

# This is the slow section. Inspect the point results first. While testing,
# change bootstrap_B near the top to a small number such as 20.

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

cpc_results <- add_bootstrap_intervals(
  cpc_results, annual_bootstrap, june_bootstrap
)
cpc_results

write.csv(
  cpc_results,
  file.path(output_dir, "cpc_tx_risk_ratios_PROVISIONAL.csv"),
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
    point_results = cpc_results
  ),
  file.path(output_dir, "cpc_fits_and_bootstrap.rds")
)
saveRDS(
  list(fits = list(Annual = fit_annual, June = fit_june),
       bootstrap = list(Annual = annual_bootstrap, June = june_bootstrap),
       results = cpc_results),
  file.path(output_dir, "fits_and_bootstrap.rds")
)
