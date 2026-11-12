# Prepare the shared smoothed GMST covariate once.
# Run this script only when the source GISTEMP file or smoothing assumptions
# change. Dataset scripts read the saved CSV directly.

# 1. Load paths and the GMST helper
source(here::here("Code", "config", "paths.R"))
source(here::here("Code", "helpers", "workflow.R"))

# 2. Read the original NASA GISTEMP monthly anomaly table
file.exists(gmst_file)
gmst_result <- prepare_smoothed_gmst(gmst_file)

# 3. Inspect intermediate and final values
gmst_observed <- gmst_result$observed
endpoint_trend <- gmst_result$endpoint_trend
gmst_complete <- gmst_result$complete
gmst_2026 <- gmst_result$gmst_2026

head(gmst_observed)
tail(gmst_observed)
summary(endpoint_trend)
tail(gmst_complete, 10)
gmst_2026

# 4. Save the shared result
dir.create(dirname(smoothed_gmst_file), recursive = TRUE, showWarnings = FALSE)
write.csv(gmst_complete, smoothed_gmst_file, row.names = FALSE)

smoothed_gmst_file
file.exists(smoothed_gmst_file)
