# =============================================================================
# 02_lc_priors.R
#
# PURPOSE
# -------
# Derives moderately informative prior distributions for the reference carbon
# (beta_init_*) and maximum SOC (beta_max_*) coefficients used in
# 03_bayes_analysis.R.  Priors are obtained by fitting linear regressions of
# log SOC density against environmental covariates on three globally gridded
# datasets:
#
#   1. No-land-use scenario  (Sanderman et al. 2017 PNAS) — informs beta_max_*
#      priors, representing the potential maximum SOC under no human disturbance
#   2. Current cropland areas (SoilGrids, masked by Copernicus cropland cover)
#      — informs beta_init_* priors for the cropland [1] land use class
#   3. Current pasture areas  (SoilGrids, masked by SEDAC pasture cover)
#      — informs beta_init_* priors for the pasture [2] land use class
#
# All covariates are centred and scaled using the parameters in sc_df.csv
# (produced by 01_soc_data_preparation.R) so that the resulting regression
# coefficients are directly comparable to the Bayesian model parameters they
# inform.
#
# SCRIPT DEPENDENCIES (must be run first)
# ----------------------------------------
#   01_soc_data_preparation.R  — produces ./data_inputs/sc_df.csv
#
# INPUTS
# ------
#   ./data_inputs/sc_df.csv
#   ./data_inputs/raster/Sanderman2017/OCD_*cm_year_NoLU_10km.tif  (4 depths)
#   ./data_inputs/raster/via_gee/amt.tif          — mean annual temperature
#   ./data_inputs/raster/via_gee/tap.tif          — total annual precipitation
#   ./data_inputs/raster/via_gee/clay_*_mean.tif  (6 depth layers)
#   ./data_inputs/raster/via_gee/nitrogen_*_mean.tif (6 depth layers)
#   ./data_inputs/raster/via_gee/ocd_*_mean.tif   (6 depth layers, SoilGrids)
#   ./data_inputs/raster/via_gee/cfvo_*_mean.tif  (6 depth layers, coarse
#                                                   fragment volume fraction)
#   ./data_inputs/raster/via_gee/cropland_copernicus.tif
#   ./data_inputs/raster/pasture_sedac.tif
#
# OUTPUTS
# -------
#   ./data_output/no_lu_coefs.RData      — summary(lm) for the no-land-use model
#   ./data_output/cropland_coefs.RData   — summary(lm) for the cropland model
#   ./data_output/pasture_coefs.RData    — summary(lm) for the pasture model
#
# NOTE ON PRIOR PRECISION
# -----------------------
# The regression coefficient standard errors are divided by 1e7 when passed to
# JAGS as prior precision values (tau = 1 / SE^2 / 1e7).  This makes the
# priors moderately informative — they regularise the posteriors toward
# ecologically plausible values while allowing the data to dominate.
# =============================================================================


# ── Packages ──────────────────────────────────────────────────────────────────
# ipak() installs and loads packages in one call
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(c("tidyverse", "terra"))


# =============================================================================
# USER FLAG — REFIT MODELS
# =============================================================================
# Set refit_models <- TRUE  to force all three raster regressions (no-land-use,
# cropland, pasture) to be re-run even if their saved .RData files already
# exist on disk.
#
# Set refit_models <- FALSE (the default) to load cached .RData files whenever
# all three are found, skipping the expensive raster loading and regression
# steps.  If any one of the three output files is missing, all three sections
# will be re-run to ensure the outputs are mutually consistent.

refit_models <- FALSE


# =============================================================================
# SECTION 1 — LOAD SCALING PARAMETERS
# =============================================================================
# sc_df provides the centering and scaling constants used throughout this
# script.  It must reflect the same dataset that will be passed to
# 03_bayes_analysis.R, so re-run 01_soc_data_preparation.R first whenever
# the analysis dataset changes.

sc_df <- read.csv('./data_inputs/sc_df.csv', row.names = 1)

# Helper: apply the centering and scaling stored in sc_df to a raw covariate
# value.  If natural_val is NULL, sd_val (a pre-scaled value) is returned
# directly — used for depth, which is set to a fixed reference value rather
# than extracted from a raster.
sc_fun <- function(var, natural_val, sd_val) {
  if (!is.null(natural_val)) {
    natural_val <- ifelse(sc_df[var, "log"], log(natural_val), natural_val)
    out_val     <- (natural_val - sc_df[var, "center"]) / sc_df[var, "scale"]
  } else {
    out_val <- sd_val
  }
  return(out_val)
}


# =============================================================================
# CACHE CHECK — load outputs if all three exist and refit_models is FALSE
# =============================================================================
# All three output files must be present for any one to be trusted as current,
# since they share raster inputs and scaling parameters.  If any is missing,
# all three sections are re-run.

no_lu_outfile   <- "./data_output/no_lu_coefs.RData"
cropland_outfile <- "./data_output/cropland_coefs.RData"
pasture_outfile  <- "./data_output/pasture_coefs.RData"

all_cached <- file.exists(no_lu_outfile) &&
  file.exists(cropland_outfile) &&
  file.exists(pasture_outfile)

if (!refit_models && all_cached) {
  message("All three prior coefficient files found — loading from cache.")
  message("  ", no_lu_outfile)
  message("  ", cropland_outfile)
  message("  ", pasture_outfile)
  load(no_lu_outfile)
  load(cropland_outfile)
  load(pasture_outfile)
} else {
  if (refit_models) {
    message("refit_models = TRUE — re-running all three raster regressions.")
  } else {
    missing_files <- c(no_lu_outfile, cropland_outfile, pasture_outfile)[
      !c(file.exists(no_lu_outfile),
         file.exists(cropland_outfile),
         file.exists(pasture_outfile))
    ]
    message("One or more prior coefficient files not found — re-running all sections.")
    message("  Missing: ", paste(missing_files, collapse = ", "))
  }
  
  # =============================================================================
  # SECTION 2 — NO LAND USE MODEL (informs beta_max_* priors)
  # =============================================================================
  # The Sanderman et al. 2017 PNAS "no land use" scenario provides gridded
  # estimates of SOC density (g C / m²) at four standard depths under the
  # counterfactual of no human land disturbance.  Regressing log SOC density
  # against environmental covariates at these depths — using the same covariate
  # scaling as the Bayesian model — yields coefficient estimates that serve as
  # moderately informative priors for beta_max_*, which represents the potential
  # maximum SOC that abandoned land could recover toward.
  
  # ── 2a. Load Sanderman no-land-use rasters ────────────────────────────────────
  no_lu_0   <- rast("./data_inputs/raster/Sanderman2017/OCD_0cm_year_NoLU_10km.tif")
  no_lu_30  <- rast("./data_inputs/raster/Sanderman2017/OCD_30cm_year_NoLU_10km.tif")
  no_lu_100 <- rast("./data_inputs/raster/Sanderman2017/OCD_100cm_year_NoLU_10km.tif")
  no_lu_200 <- rast("./data_inputs/raster/Sanderman2017/OCD_200cm_year_NoLU_10km.tif")
  
  # ── 2b. Load and reproject climate and soil covariates ────────────────────────
  # All covariates are reprojected to the CRS and resolution of the Sanderman
  # rasters so that values can be extracted at the same pixel locations.
  proj_crs <- crs(no_lu_0)
  
  mat <- rast("./data_inputs/raster/via_gee/amt.tif") %>%
    project(proj_crs) %>%
    resample(no_lu_0, method = "bilinear")
  
  map <- rast("./data_inputs/raster/via_gee/tap.tif") %>%
    project(proj_crs) %>%
    resample(no_lu_0, method = "bilinear")
  
  clay_0   <- rast("./data_inputs/raster/via_gee/clay_0_5cm_mean.tif")   %>% project(proj_crs) %>% resample(no_lu_0, method = "bilinear")
  clay_30  <- rast("./data_inputs/raster/via_gee/clay_15_30cm_mean.tif") %>% project(proj_crs) %>% resample(no_lu_0, method = "bilinear")
  clay_100 <- rast("./data_inputs/raster/via_gee/clay_60_100cm_mean.tif")  %>% project(proj_crs) %>% resample(no_lu_0, method = "bilinear")
  clay_200 <- rast("./data_inputs/raster/via_gee/clay_100_200cm_mean.tif") %>% project(proj_crs) %>% resample(no_lu_0, method = "bilinear")
  
  nitrogen_0   <- rast("./data_inputs/raster/via_gee/nitrogen_0_5cm_mean.tif")   %>% project(proj_crs) %>% resample(no_lu_0, method = "bilinear")
  nitrogen_30  <- rast("./data_inputs/raster/via_gee/nitrogen_15_30cm_mean.tif") %>% project(proj_crs) %>% resample(no_lu_0, method = "bilinear")
  nitrogen_100 <- rast("./data_inputs/raster/via_gee/nitrogen_60_100cm_mean.tif")  %>% project(proj_crs) %>% resample(no_lu_0, method = "bilinear")
  nitrogen_200 <- rast("./data_inputs/raster/via_gee/nitrogen_100_200cm_mean.tif") %>% project(proj_crs) %>% resample(no_lu_0, method = "bilinear")
  
  # ── 2c. Assemble regression data frame ───────────────────────────────────────
  # Each depth layer becomes a block of rows.  The response is log SOC density
  # (g C / m² → Mg C / ha conversion: divide by 100).  All predictors are
  # centred and scaled using sc_df.  Pixels with very low precipitation
  # (map_scaled < -6) are excluded as they represent hyper-arid areas outside
  # the range of the cropland/pasture dataset.
  raster_df <- tibble(
    c_dens   = log(as.vector(no_lu_0) / 100),
    mat      = (as.vector(mat / 10)   - sc_df["mat",      "center"]) / sc_df["mat",      "scale"],
    map      = (log(as.vector(map))   - sc_df["map",      "center"]) / sc_df["map",      "scale"],
    clay     = (as.vector(clay_0)     - sc_df["clay",     "center"]) / sc_df["clay",     "scale"],
    nitrogen = (log(as.vector(nitrogen_0)) - sc_df["nitrogen", "center"]) / sc_df["nitrogen", "scale"],
    depth    = (log(1)                - sc_df["depth",    "center"]) / sc_df["depth",    "scale"]
  ) %>%
    bind_rows(tibble(
      c_dens   = log(as.vector(no_lu_30) / 100),
      mat      = (as.vector(mat / 10)   - sc_df["mat",      "center"]) / sc_df["mat",      "scale"],
      map      = (log(as.vector(map))   - sc_df["map",      "center"]) / sc_df["map",      "scale"],
      clay     = (as.vector(clay_30)    - sc_df["clay",     "center"]) / sc_df["clay",     "scale"],
      nitrogen = (log(as.vector(nitrogen_30)) - sc_df["nitrogen", "center"]) / sc_df["nitrogen", "scale"],
      depth    = (log(30)               - sc_df["depth",    "center"]) / sc_df["depth",    "scale"]
    )) %>%
    bind_rows(tibble(
      c_dens   = log(as.vector(no_lu_100) / 100),
      mat      = (as.vector(mat / 10)    - sc_df["mat",      "center"]) / sc_df["mat",      "scale"],
      map      = (log(as.vector(map))    - sc_df["map",      "center"]) / sc_df["map",      "scale"],
      clay     = (as.vector(clay_100)    - sc_df["clay",     "center"]) / sc_df["clay",     "scale"],
      nitrogen = (log(as.vector(nitrogen_100)) - sc_df["nitrogen", "center"]) / sc_df["nitrogen", "scale"],
      depth    = (log(100)               - sc_df["depth",    "center"]) / sc_df["depth",    "scale"]
    )) %>%
    bind_rows(tibble(
      c_dens   = log(as.vector(no_lu_200) / 100),
      mat      = (as.vector(mat / 10)    - sc_df["mat",      "center"]) / sc_df["mat",      "scale"],
      map      = (log(as.vector(map))    - sc_df["map",      "center"]) / sc_df["map",      "scale"],
      clay     = (as.vector(clay_200)    - sc_df["clay",     "center"]) / sc_df["clay",     "scale"],
      nitrogen = (log(as.vector(nitrogen_200)) - sc_df["nitrogen", "center"]) / sc_df["nitrogen", "scale"],
      depth    = (log(200)               - sc_df["depth",    "center"]) / sc_df["depth",    "scale"]
    )) %>%
    drop_na() %>%
    filter(is.finite(c_dens)) %>%
    filter(map > -6)   # exclude hyper-arid pixels outside cropland/pasture range
  
  # ── 2d. Fit regression and save ───────────────────────────────────────────────
  no_lu_mod <- lm(c_dens ~ depth + map * mat + mat + I(mat^2) +
                    map + I(map^2) + clay + nitrogen,
                  data = raster_df)
  
  no_lu_coefs <- summary(no_lu_mod)
  save(no_lu_coefs, file = "./data_output/no_lu_coefs.RData")
  
  
  # =============================================================================
  # SECTION 3 — CROPLAND MODEL (informs beta_init_* priors, cropland [1])
  # =============================================================================
  # SoilGrids organic carbon density rasters, masked to pixels with > 70%
  # cropland cover (Copernicus land cover), provide the training data for the
  # cropland reference carbon prior.  SoilGrids OCD values are first corrected
  # for coarse fragment volume (cfvo) before log-transforming.
  
  # ── 3a. Load SoilGrids OCD and coarse fragment rasters ───────────────────────
  ocd_0   <- rast("./data_inputs/raster/via_gee/ocd_0_5cm_mean.tif")
  ocd_5   <- rast("./data_inputs/raster/via_gee/ocd_5_15cm_mean.tif")
  ocd_15  <- rast("./data_inputs/raster/via_gee/ocd_15_30cm_mean.tif")
  ocd_30  <- rast("./data_inputs/raster/via_gee/ocd_30_60cm_mean.tif")
  ocd_60  <- rast("./data_inputs/raster/via_gee/ocd_60_100cm_mean.tif")
  ocd_100 <- rast("./data_inputs/raster/via_gee/ocd_100_200cm_mean.tif")
  
  # Coarse fragment volume fraction (cm³/dm³ — divide by 1000 to get
  # dimensionless fraction).  OCD is corrected by multiplying by (1 - cfvo)
  # so that stocks reflect the fine-earth fraction only.
  cfvo_0   <- rast("./data_inputs/raster/via_gee/cfvo_0_5cm_mean.tif")
  cfvo_5   <- rast("./data_inputs/raster/via_gee/cfvo_5_15cm_mean.tif")
  cfvo_15  <- rast("./data_inputs/raster/via_gee/cfvo_15_30cm_mean.tif")
  cfvo_30  <- rast("./data_inputs/raster/via_gee/cfvo_30_60cm_mean.tif")
  cfvo_60  <- rast("./data_inputs/raster/via_gee/cfvo_60_100cm_mean.tif")
  cfvo_100 <- rast("./data_inputs/raster/via_gee/cfvo_100_200cm_mean.tif")
  
  ocd_0   <- ocd_0   * (1 - cfvo_0   / 1000)
  ocd_5   <- ocd_5   * (1 - cfvo_5   / 1000)
  ocd_15  <- ocd_15  * (1 - cfvo_15  / 1000)
  ocd_30  <- ocd_30  * (1 - cfvo_30  / 1000)
  ocd_60  <- ocd_60  * (1 - cfvo_60  / 1000)
  ocd_100 <- ocd_100 * (1 - cfvo_100 / 1000)
  
  # ── 3b. Load climate and soil covariates (native CRS — no reprojection needed)
  proj_crs <- crs(ocd_0)
  
  mat <- rast("./data_inputs/raster/via_gee/amt.tif")
  map <- rast("./data_inputs/raster/via_gee/tap.tif")
  
  clay_0   <- rast("./data_inputs/raster/via_gee/clay_0_5cm_mean.tif")
  clay_5   <- rast("./data_inputs/raster/via_gee/clay_5_15cm_mean.tif")
  clay_15  <- rast("./data_inputs/raster/via_gee/clay_15_30cm_mean.tif")
  clay_30  <- rast("./data_inputs/raster/via_gee/clay_30_60cm_mean.tif")
  clay_60  <- rast("./data_inputs/raster/via_gee/clay_60_100cm_mean.tif")
  clay_100 <- rast("./data_inputs/raster/via_gee/clay_100_200cm_mean.tif")
  
  nitrogen_0   <- rast("./data_inputs/raster/via_gee/nitrogen_0_5cm_mean.tif")
  nitrogen_5   <- rast("./data_inputs/raster/via_gee/nitrogen_5_15cm_mean.tif")
  nitrogen_15  <- rast("./data_inputs/raster/via_gee/nitrogen_15_30cm_mean.tif")
  nitrogen_30  <- rast("./data_inputs/raster/via_gee/nitrogen_30_60cm_mean.tif")
  nitrogen_60  <- rast("./data_inputs/raster/via_gee/nitrogen_60_100cm_mean.tif")
  nitrogen_100 <- rast("./data_inputs/raster/via_gee/nitrogen_100_200cm_mean.tif")
  
  # ── 3c. Define cropland mask ──────────────────────────────────────────────────
  # Pixels with > 70% cropland cover (Copernicus) are used to obtain a globally
  # representative sample of current cropland SOC — sufficient coverage while
  # excluding transitional or mixed-use areas.
  cropland <- rast("./data_inputs/raster/via_gee/cropland_copernicus.tif")
  cropland <- cropland > 70
  plot(cropland)
  
  # ── 3d. Assemble regression data frame ───────────────────────────────────────
  # Depth midpoints match the SoilGrids layer boundaries.  The response is
  # log OCD corrected for coarse fragments, divided by 100 to convert units.
  raster_df <- tibble(
    c_dens   = log(as.vector(ocd_0[cropland]) / 100),
    mat      = (as.vector(mat[cropland] / 10) - sc_df["mat",      "center"]) / sc_df["mat",      "scale"],
    map      = (log(as.vector(map[cropland])) - sc_df["map",      "center"]) / sc_df["map",      "scale"],
    clay     = (as.vector(clay_0[cropland])   - sc_df["clay",     "center"]) / sc_df["clay",     "scale"],
    nitrogen = (log(as.vector(nitrogen_0[cropland])) - sc_df["nitrogen", "center"]) / sc_df["nitrogen", "scale"],
    depth    = (log(2.5)                      - sc_df["depth",    "center"]) / sc_df["depth",    "scale"]
  ) %>%
    bind_rows(tibble(
      c_dens   = log(as.vector(ocd_5[cropland]) / 100),
      mat      = (as.vector(mat[cropland] / 10) - sc_df["mat",      "center"]) / sc_df["mat",      "scale"],
      map      = (log(as.vector(map[cropland])) - sc_df["map",      "center"]) / sc_df["map",      "scale"],
      clay     = (as.vector(clay_5[cropland])   - sc_df["clay",     "center"]) / sc_df["clay",     "scale"],
      nitrogen = (log(as.vector(nitrogen_5[cropland])) - sc_df["nitrogen", "center"]) / sc_df["nitrogen", "scale"],
      depth    = (log(10)                       - sc_df["depth",    "center"]) / sc_df["depth",    "scale"]
    )) %>%
    bind_rows(tibble(
      c_dens   = log(as.vector(ocd_15[cropland]) / 100),
      mat      = (as.vector(mat[cropland] / 10)  - sc_df["mat",      "center"]) / sc_df["mat",      "scale"],
      map      = (log(as.vector(map[cropland]))  - sc_df["map",      "center"]) / sc_df["map",      "scale"],
      clay     = (as.vector(clay_15[cropland])   - sc_df["clay",     "center"]) / sc_df["clay",     "scale"],
      nitrogen = (log(as.vector(nitrogen_15[cropland])) - sc_df["nitrogen", "center"]) / sc_df["nitrogen", "scale"],
      depth    = (log((15 + 30) / 2)             - sc_df["depth",    "center"]) / sc_df["depth",    "scale"]
    )) %>%
    bind_rows(tibble(
      c_dens   = log(as.vector(ocd_30[cropland]) / 100),
      mat      = (as.vector(mat[cropland] / 10)  - sc_df["mat",      "center"]) / sc_df["mat",      "scale"],
      map      = (log(as.vector(map[cropland]))  - sc_df["map",      "center"]) / sc_df["map",      "scale"],
      clay     = (as.vector(clay_30[cropland])   - sc_df["clay",     "center"]) / sc_df["clay",     "scale"],
      nitrogen = (log(as.vector(nitrogen_30[cropland])) - sc_df["nitrogen", "center"]) / sc_df["nitrogen", "scale"],
      depth    = (log((30 + 60) / 2)             - sc_df["depth",    "center"]) / sc_df["depth",    "scale"]
    )) %>%
    bind_rows(tibble(
      c_dens   = log(as.vector(ocd_60[cropland]) / 100),
      mat      = (as.vector(mat[cropland] / 10)  - sc_df["mat",      "center"]) / sc_df["mat",      "scale"],
      map      = (log(as.vector(map[cropland]))  - sc_df["map",      "center"]) / sc_df["map",      "scale"],
      clay     = (as.vector(clay_60[cropland])   - sc_df["clay",     "center"]) / sc_df["clay",     "scale"],
      nitrogen = (log(as.vector(nitrogen_60[cropland])) - sc_df["nitrogen", "center"]) / sc_df["nitrogen", "scale"],
      depth    = (log((60 + 100) / 2)            - sc_df["depth",    "center"]) / sc_df["depth",    "scale"]
    )) %>%
    bind_rows(tibble(
      c_dens   = log(as.vector(ocd_100[cropland]) / 100),
      mat      = (as.vector(mat[cropland] / 10)   - sc_df["mat",      "center"]) / sc_df["mat",      "scale"],
      map      = (log(as.vector(map[cropland]))   - sc_df["map",      "center"]) / sc_df["map",      "scale"],
      clay     = (as.vector(clay_100[cropland])   - sc_df["clay",     "center"]) / sc_df["clay",     "scale"],
      nitrogen = (log(as.vector(nitrogen_100[cropland])) - sc_df["nitrogen", "center"]) / sc_df["nitrogen", "scale"],
      depth    = (log(150)                         - sc_df["depth",    "center"]) / sc_df["depth",    "scale"]
    )) %>%
    drop_na() %>%
    filter(is.finite(c_dens)) %>%
    filter(map > -6)   # exclude hyper-arid pixels
  
  # ── 3e. Fit regression and save ───────────────────────────────────────────────
  cropland_mod <- lm(c_dens ~ depth + map * mat + mat + I(mat^2) +
                       map + I(map^2) + clay + nitrogen,
                     data = raster_df)
  
  cropland_coefs <- summary(cropland_mod)
  save(cropland_coefs, file = "./data_output/cropland_coefs.RData")
  rm(cropland_coefs, cropland_mod, no_lu_coefs, no_lu_mod)
  gc()
  
  
  # =============================================================================
  # SECTION 4 — PASTURE MODEL (informs beta_init_* priors, pasture [2])
  # =============================================================================
  # Identical pipeline to Section 3, applied to pixels with > 80% pasture cover
  # (SEDAC).  A more conservative threshold is used here because pasture cover
  # estimates are noisier than cropland cover, and a higher threshold reduces
  # contamination from mixed or transitional land uses.
  
  # ── 4a. Define pasture mask ───────────────────────────────────────────────────
  pasture <- rast("./data_inputs/raster/pasture_sedac.tif")
  pasture <- pasture > 0.80
  plot(pasture)
  
  # ── 4b. Assemble regression data frame ───────────────────────────────────────
  raster_df <- tibble(
    c_dens   = log(as.vector(ocd_0[pasture]) / 100),
    mat      = (as.vector(mat[pasture] / 10) - sc_df["mat",      "center"]) / sc_df["mat",      "scale"],
    map      = (log(as.vector(map[pasture])) - sc_df["map",      "center"]) / sc_df["map",      "scale"],
    clay     = (as.vector(clay_0[pasture])   - sc_df["clay",     "center"]) / sc_df["clay",     "scale"],
    nitrogen = (log(as.vector(nitrogen_0[pasture])) - sc_df["nitrogen", "center"]) / sc_df["nitrogen", "scale"],
    depth    = (log(2.5)                     - sc_df["depth",    "center"]) / sc_df["depth",    "scale"]
  ) %>%
    bind_rows(tibble(
      c_dens   = log(as.vector(ocd_5[pasture]) / 100),
      mat      = (as.vector(mat[pasture] / 10) - sc_df["mat",      "center"]) / sc_df["mat",      "scale"],
      map      = (log(as.vector(map[pasture])) - sc_df["map",      "center"]) / sc_df["map",      "scale"],
      clay     = (as.vector(clay_5[pasture])   - sc_df["clay",     "center"]) / sc_df["clay",     "scale"],
      nitrogen = (log(as.vector(nitrogen_5[pasture])) - sc_df["nitrogen", "center"]) / sc_df["nitrogen", "scale"],
      depth    = (log(10)                      - sc_df["depth",    "center"]) / sc_df["depth",    "scale"]
    )) %>%
    bind_rows(tibble(
      c_dens   = log(as.vector(ocd_15[pasture]) / 100),
      mat      = (as.vector(mat[pasture] / 10)  - sc_df["mat",      "center"]) / sc_df["mat",      "scale"],
      map      = (log(as.vector(map[pasture]))  - sc_df["map",      "center"]) / sc_df["map",      "scale"],
      clay     = (as.vector(clay_15[pasture])   - sc_df["clay",     "center"]) / sc_df["clay",     "scale"],
      nitrogen = (log(as.vector(nitrogen_15[pasture])) - sc_df["nitrogen", "center"]) / sc_df["nitrogen", "scale"],
      depth    = (log((15 + 30) / 2)            - sc_df["depth",    "center"]) / sc_df["depth",    "scale"]
    )) %>%
    bind_rows(tibble(
      c_dens   = log(as.vector(ocd_30[pasture]) / 100),
      mat      = (as.vector(mat[pasture] / 10)  - sc_df["mat",      "center"]) / sc_df["mat",      "scale"],
      map      = (log(as.vector(map[pasture]))  - sc_df["map",      "center"]) / sc_df["map",      "scale"],
      clay     = (as.vector(clay_30[pasture])   - sc_df["clay",     "center"]) / sc_df["clay",     "scale"],
      nitrogen = (log(as.vector(nitrogen_30[pasture])) - sc_df["nitrogen", "center"]) / sc_df["nitrogen", "scale"],
      depth    = (log((30 + 60) / 2)            - sc_df["depth",    "center"]) / sc_df["depth",    "scale"]
    )) %>%
    bind_rows(tibble(
      c_dens   = log(as.vector(ocd_60[pasture]) / 100),
      mat      = (as.vector(mat[pasture] / 10)  - sc_df["mat",      "center"]) / sc_df["mat",      "scale"],
      map      = (log(as.vector(map[pasture]))  - sc_df["map",      "center"]) / sc_df["map",      "scale"],
      clay     = (as.vector(clay_60[pasture])   - sc_df["clay",     "center"]) / sc_df["clay",     "scale"],
      nitrogen = (log(as.vector(nitrogen_60[pasture])) - sc_df["nitrogen", "center"]) / sc_df["nitrogen", "scale"],
      depth    = (log((60 + 100) / 2)           - sc_df["depth",    "center"]) / sc_df["depth",    "scale"]
    )) %>%
    bind_rows(tibble(
      c_dens   = log(as.vector(ocd_100[pasture]) / 100),
      mat      = (as.vector(mat[pasture] / 10)   - sc_df["mat",      "center"]) / sc_df["mat",      "scale"],
      map      = (log(as.vector(map[pasture]))   - sc_df["map",      "center"]) / sc_df["map",      "scale"],
      clay     = (as.vector(clay_100[pasture])   - sc_df["clay",     "center"]) / sc_df["clay",     "scale"],
      nitrogen = (log(as.vector(nitrogen_100[pasture])) - sc_df["nitrogen", "center"]) / sc_df["nitrogen", "scale"],
      depth    = (log(150)                        - sc_df["depth",    "center"]) / sc_df["depth",    "scale"]
    )) %>%
    drop_na() %>%
    filter(is.finite(c_dens)) %>%
    filter(map > -6)   # exclude hyper-arid pixels
  
  # ── 4c. Fit regression and save ───────────────────────────────────────────────
  pasture_mod <- lm(c_dens ~ depth + map * mat + mat + I(mat^2) +
                      map + I(map^2) + clay + nitrogen,
                    data = raster_df)
  
  pasture_coefs <- summary(pasture_mod)
  save(pasture_coefs, file = "./data_output/pasture_coefs.RData")
  
} # end of refit / cache-miss block