# Central paths and run settings for the attribution scripts.
# Edit paths here when files move; dataset scripts should not contain paths.

project_dir <- "C:/Users/siyan/WCRP/wcrp_common_attribution_assessment"

# Input temperature datasets
cpc_tmax_file <- file.path(
  project_dir, "outputs/pnw_heatwave",
  "CPC_reference_daily_tmax.nc"
)

# cpc_tmin_file <- file.path(
#   project_dir, "data/raw",
#   "tmin_cpc_daily_-10.5-20E_40-60N_firstyear-lastyear.nc"
# )




berkeley_tmax_file <- file.path(
  project_dir, "outputs/pnw_heatwave",
  "Berkeley_reference_tmax.nc"
)

era5_tmax_file <- file.path(
  project_dir, "outputs/pnw_heatwave",
  "ERA5_reference_mx2t.nc"
)

# GMST input
gmst_file <- file.path(project_dir, "data", "gistemp_global.csv")
smoothed_gmst_file <- file.path(
  project_dir,
  "data",
  "processed",
  "gmst_smoothed.csv"
)

# Output folders
cpc_tmax_output_dir <- file.path(project_dir, "results_cpc_tmax")
#cpc_tmin_output_dir <- file.path(project_dir, "results_cpc_tmin")
berkeley_tx_output_dir <- file.path(project_dir, "results_berkeley_tx")
era5_tx_output_dir <- file.path(project_dir, "results_era5_tx")
figure7_output_dir <- file.path(project_dir, "figure")

# ERA5 needs a land-only mask because its temperature field also covers oceans.
# This is a fixed input copied into data/masks, so ERA5 can run independently
# without first running the E-Obs analysis script.
era5_land_mask_file <- file.path(
  project_dir,
  "data",
  "masks",
  "eobs_v33e_land_footprint.csv"
)

# Settings used by optional bootstrap sections.
bootstrap_B <- 200L
figure7_bootstrap_B <- 300L
random_seed <- 20260915

# Quick path check. Run this line interactively after sourcing this file.
input_file_check <- c(
  cpc_tmax = file.exists(cpc_tmax_file),
  #cpc_tmin = file.exists(cpc_tmin_file),
  berkeley_tmax = file.exists(berkeley_tmax_file),
  era5_tmax = file.exists(era5_tmax_file),
  era5_land_mask = file.exists(era5_land_mask_file),
  gmst = file.exists(gmst_file),
  smoothed_gmst = file.exists(smoothed_gmst_file)
)
