# Smith Section 0.6.1 GEV MLE reproduction: Table 0.1 and endpoint checks.
# Run sections in RStudio. No bootstrap, MCMC, or figures are executed here.

# ---- 1. Paths and configuration ----
project_root <- "C:/Users/siyan/WCRP/wcrp_common_attribution_assessment"
smith_dir <- file.path(project_root, "Code/smith_reproduce")
source(file.path(smith_dir, "helpers/smith_pnw_pot.R"))
source(file.path(smith_dir, "helpers/smith_pnw_gev.R"))
daily_file <- file.path(project_root, "outputs/pnw_heatwave/ERA5_reference_daily_tmax.nc")
gmst_file <- file.path(project_root, "data/processed/gmst_smoothed.csv")
output_dir <- file.path(project_root, "Results/ERA5_Smith_PNW_GEV_MLE")
first_year <- 1950L
last_historical_year <- 2020L
event_year <- 2021L
counterfactual_cooling_C <- 1.2
initial_shapes <- c(-0.45, -0.2, 0.05)
optimizer_maxit <- 50000L
paper_event_temperature_C <- 39.6
paper_previous_record_C <- 33.98

stopifnot(requireNamespace("ncdf4", quietly = TRUE), requireNamespace("ismev", quietly = TRUE),
          file.exists(daily_file), file.exists(gmst_file), event_year == last_historical_year + 1L)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 2. Full-calendar annual TXx, not warm-season maxima or cluster peaks ----
daily <- smith_read_daily(daily_file)
gmst <- read.csv(gmst_file)
annual_data <- smith_gev_annual(daily, gmst, first_year, event_year)
annual_data
write.csv(annual_data, file.path(output_dir, "annual_txx_and_coverage.csv"), row.names = FALSE)
fit_excluding <- annual_data[annual_data$used_in_fit & annual_data$year <= last_historical_year, ]
fit_including <- annual_data[annual_data$used_in_fit & annual_data$year <= event_year, ]
event <- annual_data[annual_data$year == event_year, ]
stopifnot(nrow(event) == 1L, event$used_in_fit)
g_factual <- gmst$gmst[gmst$year == event_year]
g_counterfactual <- g_factual - counterfactual_cooling_C
stopifnot(length(g_factual) == 1L, is.finite(g_factual))

# ---- 3. Fit both Table 0.1 models using ismev::gev.fit ----
# mu(g)=beta0+betaG*g, constant sigma and xi. Center GMST for numerical stability;
# undo centering before comparing beta0 with the published table.
fits <- list(
  exclude_2021 = smith_gev_fit(fit_excluding, initial_shapes, optimizer_maxit),
  include_2021 = smith_gev_fit(fit_including, initial_shapes, optimizer_maxit))
parameter_table <- rbind(
  smith_gev_table_parameters(fits$exclude_2021, "exclude_2021"),
  smith_gev_table_parameters(fits$include_2021, "include_2021"))
table_0_1_comparison <- merge(parameter_table, smith_gev_paper_table(),
                              by = c("case", "parameter"), sort = FALSE)
table_0_1_comparison
optimizer_audit <- rbind(
  data.frame(case = "exclude_2021", fits$exclude_2021$audit),
  data.frame(case = "include_2021", fits$include_2021$audit))

# ---- 4. MLE endpoints, delta standard errors and fixed-climate event risks ----
endpoints <- rbind(
  smith_gev_endpoints(fits$exclude_2021, g_factual, g_counterfactual, "exclude_2021"),
  smith_gev_endpoints(fits$include_2021, g_factual, g_counterfactual, "include_2021"))
# Paper excluding-2021 endpoints: 35.05 and 32.94 C (SE 0.61 and 0.40).
# Existing ERA5/GMST differ; these are references, not imposed constraints.
endpoints
thresholds <- sort(unique(c(paper_previous_record_C, paper_event_temperature_C, event$txx)))
attribution <- rbind(
  smith_gev_attribution(fits$exclude_2021, g_factual, g_counterfactual, thresholds, "exclude_2021"),
  smith_gev_attribution(fits$include_2021, g_factual, g_counterfactual, thresholds, "include_2021"))
attribution

# ---- 5. Figure 1 QQ coordinates and goodness-of-fit statistics, no figures ----
gof_excluding <- smith_gev_gof(fits$exclude_2021, "exclude_2021")
gof_including <- smith_gev_gof(fits$include_2021, "include_2021")
gof_statistics <- rbind(gof_excluding$statistics, gof_including$statistics)
qq_coordinates <- rbind(gof_excluding$qq, gof_including$qq)
gof_statistics
# Simulated GOF p-values are not recreated: their source calibration is not fully
# specified. Do not use ordinary fixed-parameter p-values for these fitted tests.

# ---- 6. Save a dedicated snapshot for the GEV Bayesian script ----
settings <- list(daily_file = daily_file, gmst_file = gmst_file,
  first_year = first_year, last_historical_year = last_historical_year, event_year = event_year,
  g_factual = g_factual, g_counterfactual = g_counterfactual,
  counterfactual_cooling_C = counterfactual_cooling_C,
  paper_event_temperature_C = paper_event_temperature_C,
  paper_previous_record_C = paper_previous_record_C,
  observed_event_temperature_C = event$txx, initial_shapes = initial_shapes,
  optimizer_maxit = optimizer_maxit, event_metric = "Full-calendar annual TXx")
write.csv(table_0_1_comparison, file.path(output_dir, "table_0_1_comparison.csv"), row.names = FALSE)
write.csv(optimizer_audit, file.path(output_dir, "optimizer_audit.csv"), row.names = FALSE)
write.csv(endpoints, file.path(output_dir, "mle_endpoints.csv"), row.names = FALSE)
write.csv(attribution, file.path(output_dir, "mle_attribution.csv"), row.names = FALSE)
write.csv(gof_statistics, file.path(output_dir, "gof_statistics.csv"), row.names = FALSE)
write.csv(qq_coordinates, file.path(output_dir, "figure_1_qq_coordinates.csv"), row.names = FALSE)
saveRDS(list(settings = settings, annual_data = annual_data, fits = fits,
             table_0_1_comparison = table_0_1_comparison, endpoints = endpoints,
             attribution = attribution, gof_statistics = gof_statistics,
             qq_coordinates = qq_coordinates),
        file.path(output_dir, "smith_pnw_gev_mle.rds"))
