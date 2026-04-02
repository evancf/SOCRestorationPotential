# =============================================================================
# 06_scenarios_estimates_maps.R
#
# PURPOSE
# -------
# Combines spatially explicit SOC stock maps (from 04_bayes_results.R) with
# aboveground carbon estimates from three literature products (CCI-BM, SBM,
# Robinson et al. 2025) and land-use area fraction rasters for four restoration
# scenarios to compute global carbon sequestration estimates (Pg C) over time.
# Credible intervals are constructed from per-sample stock maps produced by
# 05_bayes_uncertainty.R.  All manuscript scenario figures and supplementary
# tables are produced here.
#
# Sections 9–11 additionally load outputs from the SSP3-7.0 2041–2060 climate
# scenario runs of scripts 04 and 05, and produce SOC-only comparisons between
# contemporary and projected future climates.  AGC/BGB estimates are excluded
# from the future climate sections because the aboveground carbon products are
# parameterised under current climate only.
#
# SCRIPT DEPENDENCIES (must be run first)
# ----------------------------------------
#   04_bayes_results.R      — run twice: climate_scenario = "contemporary" and
#                             "ssp370_2050"; produces stock maps in
#                             ./data_output/stock_map/<scenario>/
#   05_bayes_uncertainty.R  — run twice: same two climate_scenario values;
#                             produces per-sample maps in
#                             ./data_output/uncertainty/<scenario>/stock_map_<b>/
#
# INPUTS
# ------
#   ./data_output/stock_map/contemporary/   — posterior-median SOC stock maps
#   ./data_output/stock_map/ssp370_2050/    — future posterior-median SOC maps
#   ./data_output/uncertainty/contemporary/ — per-sample SOC stock maps
#   ./data_output/uncertainty/ssp370_2050/  — future per-sample SOC maps
#   ./data_inputs/raster/                   — scenario fraction rasters, HYDE
#                                             masks, biome masks, and AGC rasters
#   ./data_inputs/raster/belowground/       — root mass fraction rasters
#
# OUTPUTS
# -------
#   ./figures/                                    — all contemporary figures
#   ./figures/delta_ci_0_100.csv
#   ./figures/delta_ci_0_30.csv
#   ./figures/supplementary_table_SOC_change.csv
#   ./figures/supplementary_table_AGC_change.csv
#   ./figures/ssp370_2050/                        — all future climate outputs
#   ./figures/ssp370_2050/delta_ci_0_100.csv
#   ./figures/ssp370_2050/delta_ci_0_30.csv
#   ./figures/ssp370_2050/supplementary_table_SOC_change_ssp370_2050.csv
#   ./figures/ssp370_2050/supplementary_table_SOC_pct_reduction_ssp370_2050.csv
#   ./figures/ssp370_2050/soil_average_rates_50_ssp370_2050.png
#   ./figures/ssp370_2050/soil_rates_diff_ssp370_2050_vs_contemporary.png
#   ./figures/ssp370_2050/maps_ssp370_2050_and_difference.png
#   ./figures/ssp370_2050/soc_contemporary_vs_ssp370_2050_bars.png
#
# KEY DESIGN NOTES
# ----------------
#   Four land-use scenarios are evaluated:
#     crop_aban       — current abandoned cropland (Potapov et al.)
#     crop_sparing    — cropland sparing (Beyer et al.)
#     crop_relocation — cropland relocation (Beyer et al.)
#     pasture_sparing — pasture sparing hotspots (Hayek et al.)
#   A combined "both_sparing" scenario sums crop and pasture sparing.
#
#   Aboveground carbon is estimated from three products (contemporary only):
#     CCI-BM  — Chapman-Richards curves from Cook-Patton et al.
#     SBM     — Chapman-Richards curves from Spawn et al.
#     Robinson 2025 — Chapman-Richards curves with delta-method SEs
#   The CCI-BM and SBM estimates are gap-filled at low tree cover using
#   treecover-scaled nearest-neighbour propagation and kNN regression.
#
#   Posterior uncertainty uses every 10th MCMC sample (agg_fact = 25 in
#   05_bayes_uncertainty_parallel.R).  Posterior-median estimates use
#   agg_fact = 10 (04_bayes_results.R).  When combining the two,
#   resample() is used rather than aggregate() to handle the non-integer
#   ratio (25/10 = 2.5).
#
# =============================================================================


# ── Packages ──────────────────────────────────────────────────────────────────
# ipak() installs and loads packages in one call
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(c("tidyverse", "terra", "ggplot2", "sf", "FNN",
       "geodata", "maps", "patchwork", "cowplot",
       "scales", "stringr"))
sf::sf_use_s2(FALSE)



# =============================================================================
# OUTPUT PATH CONFIGURATION
# =============================================================================
# The primary analysis (Sections 1–8) always uses contemporary climate inputs.
# Sections 9–11 load outputs from both the contemporary and ssp370_2050 runs
# of scripts 04 and 05 to produce the future climate sensitivity figures and
# tables.  Both sets of outputs must exist before running Sections 9–11.

stock_map_dir   <- "./data_output/stock_map/contemporary"
uncertainty_dir <- "./data_output/uncertainty/contemporary"
figures_dir     <- "./figures"

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# SECTION 1 — RASTER SETUP AND SCENARIO FRACTION LAYERS
# =============================================================================
# Loads and caches aggregated versions of all scenario fraction rasters,
# HYDE cropland/pasture masks, biome masks, and a cell-area raster.
# get_agg_0 fills NA → 0 before aggregating (used for fraction rasters);
# get_agg preserves NA (used for masks).
#
# agg_fact  — spatial aggregation factor for all rasters in this script.
# redo_rasts — set TRUE to force reprocessing of all cached rasters.

agg_fact <- 10

redo_rasts <- FALSE

cache_dir <- "./data_inputs/raster/agg_cache"
if (!dir.exists(cache_dir)) {
  dir.create(cache_dir, recursive = TRUE)
}

cache_path <- function(name) {
  file.path(cache_dir, sprintf("%s_f%d.tif", name, agg_fact))
}

get_agg_0 <- function(src_path, cache_name) {
  out <- cache_path(cache_name)
  if (!file.exists(out) || redo_rasts) {
    r <- rast(src_path)
    r[is.na(r)] <- 0
    r <- aggregate(r, fact = agg_fact, na.rm = TRUE)
    writeRaster(r, out, overwrite = TRUE)
  } else {
    r <- rast(out)
  }
  r
}

# helper: load, aggregate, cache (no NA→0 step)
get_agg <- function(src_path, cache_name) {
  out <- cache_path(cache_name)
  if (!file.exists(out) || redo_rasts) {
    r <- rast(src_path)
    r <- aggregate(r, fact = agg_fact, na.rm = TRUE)
    writeRaster(r, out, overwrite = TRUE)
  } else {
    r <- rast(out)
  }
  r
}

# Load or build each scenario fraction and mask raster

crop_aban_frac <- get_agg_0("./data_inputs/raster/abandon_frac_1km.tif", 
                            "abandon_frac")
crop_sparing_frac <- get_agg_0("./data_inputs/raster/spaArea_frac_per_grid_1km.tif", 
                               "crop_sparing_frac")
# Align to the common grid
crop_sparing_frac <- resample(crop_sparing_frac, crop_aban_frac, method = "bilinear")

crop_relocation_frac <- get_agg_0("./data_inputs/raster/spaArea_frac_grid_beyerR_1km.tif", 
                                  "crop_relocation_frac")



if(!"pasture_sparing_hotspot_frac_Hayek.tif" %in% list.files("./data_inputs/raster/")){
  writeRaster(rast("./data_inputs/raster/pasture_sparing_hotspot_area_Hayek.tif") / cellSize(rast("./data_inputs/raster/pasture_sparing_hotspot_area_Hayek.tif"), unit = "ha"),
              filename = "./data_inputs/raster/pasture_sparing_hotspot_frac_Hayek.tif")
  
}
pasture_sparing_frac <- rast("./data_inputs/raster/pasture_sparing_hotspot_frac_Hayek.tif")
# This doesn't have perfect alignment so resample to make minute corrections
pasture_sparing_frac[is.na(pasture_sparing_frac)] <- 0
pasture_sparing_frac <- resample(pasture_sparing_frac, crop_aban_frac, method = "bilinear")

# HYDE land-use masks (cells with < 5 % cover set to NA)
crop_mask <- get_agg("./data_inputs/raster/cropland_frac_2010_HYDE.tif",
                     "crop_mask")
crop_mask[crop_mask <= 0.05] <- NA
crop_mask[crop_mask > 0.05] <- 1

pasture_mask <- get_agg("./data_inputs/raster/pasture_frac_2010_HYDE.tif",
                        "pasture_mask")
pasture_mask[pasture_mask <= 0.05] <- NA
pasture_mask[pasture_mask > 0.05] <- 1

forest_mask <- get_agg("./data_inputs/raster/forest_biomes.tif",
                       "forest_mask")
nonforest_mask <- get_agg("./data_inputs/raster/nonforest_biomes.tif",
                          "nonforest_mask")

# Cell area raster (ha) — used for all Pg C calculations
area_raster <- cellSize(crop_aban_frac, unit = "ha")

# # quick plots to visualize the above
# par(mfrow = c(2,2), mar = c(1,1,2,1))
# plot(crop_aban_frac, main = "Abandoned cropland fraction")
# plot(crop_sparing_frac, main = "Cropland sparing fraction")
# plot(crop_relocation_frac,main = "Cropland relocation fraction")
# plot(pasture_sparing_frac,main = "Pasture sparing fraction")

# Report scenario areas (Mha) as a quick sanity check
mha <- function(r) sum(values(r) * values(area_raster), na.rm = TRUE)/1e6

crop_aban_frac_Mha <- mha(crop_aban_frac)
crop_sparing_frac_Mha <- mha(crop_sparing_frac)
crop_relocation_frac_Mha <- mha(crop_relocation_frac)
pasture_sparing_frac_Mha <- mha(pasture_sparing_frac)
both_sparing_frac_Mha <- (mha(crop_sparing_frac) + mha(pasture_sparing_frac))

tibble(
  crop_aban_frac_Mha = crop_aban_frac_Mha,
  crop_sparing_frac_Mha = crop_sparing_frac_Mha,
  crop_relocation_frac_Mha = crop_relocation_frac_Mha,
  pasture_sparing_frac_Mha = pasture_sparing_frac_Mha,
  both_sparing_frac_Mha = both_sparing_frac_Mha
)


# =============================================================================
# SECTION 2 — MAP HELPER FUNCTIONS
# =============================================================================
# prep_for_map applies a crop or pasture mask and combines forest/open
# transition rasters using biome masks, returning a long data frame for ggplot.
# prep_for_map_diff returns the difference relative to a reference raster.


prep_for_map <- function(forest_path,
                         mask_type = c("crop", "pasture"),
                         test_flag = FALSE) {
  # Select the appropriate land-use mask
  mask_type <- match.arg(mask_type)
  mask_obj <- switch(mask_type,
                     crop = crop_mask,
                     pasture = pasture_mask)
  
  # simple case: no transitions requested (or test_flag)
  if (!any(c(grepl("_CFN_", forest_path),
             grepl("_CON_", forest_path),
             grepl("_PFN_", forest_path),
             grepl("_PON_", forest_path),
             grepl("_CFA_", forest_path),
             grepl("_COA_", forest_path),
             grepl("_PFA_", forest_path),
             grepl("_POA_", forest_path))) || test_flag) {
    my_raster <- mask_obj * rast(forest_path)
    df <- as.data.frame(my_raster, xy = TRUE, na.rm = TRUE) %>%
      gather(key = "variable", value = "value", -x, -y)
    return(df)
  }
  
  # 1) infer the corresponding "open" file:
  open_path <- forest_path %>%
    { if (grepl("_CFN_", .)) sub("_CFN_", "_CON_", .)
      else if (grepl("_PFN_", .)) sub("_PFN_", "_PON_", .)
      else stop("forest_path must contain either '_CFN_' or '_PFN_'") }
  
  # 2) load rasters
  r_forest <- rast(forest_path)
  r_open <- rast(open_path)
  
  # 3) apply chosen mask to both transitions and biome masks
  r_forest <- mask_obj * r_forest
  r_open <- mask_obj * r_open
  mask_f <- mask_obj * forest_mask
  mask_nf <- mask_obj * nonforest_mask
  
  # 4) binarize
  bin_f <- mask_f >= 1
  bin_nf <- mask_nf >= 1
  
  # 5) piece together forest vs nonforest values
  combined <- ifel(
    bin_f, 
    r_forest,
    ifel(bin_nf, r_open, NA)
  )
  
  # 6) to data.frame for ggplot
  df <- as.data.frame(combined, xy = TRUE, na.rm = TRUE) %>%
    gather(key = "variable", value = "value", -x, -y)
  
  return(df)
}


prep_for_map_diff <- function(x, y = file.path(stock_map_dir, "stock_map_adj_init_crop_0_100.tif")) {
  raster_df <- prep_for_map(x)
  raster_df$value <- raster_df$value - (prep_for_map(y) %>% pull(value))
  return(raster_df)
}





# =============================================================================
# SECTION 3 — ABOVEGROUND CARBON ESTIMATES
# =============================================================================
# Applies Chapman-Richards growth curves at each time horizon using three
# literature products: CCI-BM (Cook-Patton et al.), SBM (Spawn et al.), and
# Robinson et al. 2025.  Estimates are gap-filled at low tree cover via
# treecover-scaled nearest-neighbour propagation followed by kNN regression.

# Chapman-Richards functions for each time horizon (output: Mg C ha⁻¹)
agc_5 <- function(max, # asymptotic accumulation
                  tau, # growth rate
                  c # a shape parameter
){
  t <- 5
  biomass_or_carbon <- "carbon"
  
  agb <- max*(1-exp(-t/tau))^c
  if(biomass_or_carbon == "biomass"){
    return(agb)
  }
  if(biomass_or_carbon == "carbon"){
    return(agb * 0.47)
  }
}

agc_15 <- function(max, # asymptotic accumulation
                   tau, # growth rate
                   c # a shape parameter
){
  t <- 15
  biomass_or_carbon <- "carbon"
  
  agb <- max*(1-exp(-t/tau))^c
  if(biomass_or_carbon == "biomass"){
    return(agb)
  }
  if(biomass_or_carbon == "carbon"){
    return(agb * 0.47)
  }
}

agc_30 <- function(max, # asymptotic accumulation
                   tau, # growth rate
                   c # a shape parameter
){
  t <- 30
  biomass_or_carbon <- "carbon"
  
  agb <- max*(1-exp(-t/tau))^c
  if(biomass_or_carbon == "biomass"){
    return(agb)
  }
  if(biomass_or_carbon == "carbon"){
    return(agb * 0.47)
  }
}

agc_50 <- function(max, # asymptotic accumulation
                   tau, # growth rate
                   c # a shape parameter
){
  t <- 50
  biomass_or_carbon <- "carbon"
  
  agb <- max*(1-exp(-t/tau))^c
  if(biomass_or_carbon == "biomass"){
    return(agb)
  }
  if(biomass_or_carbon == "carbon"){
    return(agb * 0.47)
  }
}

agc_100 <- function(max, # asymptotic accumulation
                    tau, # growth rate
                    c # a shape parameter
){
  t <- 100
  biomass_or_carbon <- "carbon"
  
  agb <- max*(1-exp(-t/tau))^c
  if(biomass_or_carbon == "biomass"){
    return(agb)
  }
  if(biomass_or_carbon == "carbon"){
    return(agb * 0.47)
  }
}

agc_200 <- function(max, # asymptotic accumulation
                    tau, # growth rate
                    c # a shape parameter
){
  t <- 200
  biomass_or_carbon <- "carbon"
  
  agb <- max*(1-exp(-t/tau))^c
  if(biomass_or_carbon == "biomass"){
    return(agb)
  }
  if(biomass_or_carbon == "carbon"){
    return(agb * 0.47)
  }
}

# Load CCI-BM Chapman-Richards parameter rasters
cci_max <- "./data_inputs/raster/CCIBM/CCIBM_max_1deg_tha.tif" %>%
  rast()
cci_tau <- "./data_inputs/raster/CCIBM/CCIBM_tau_1deg.tif" %>%
  rast()
cci_c <- "./data_inputs/raster/CCIBM/CCIBM_c_1deg.tif" %>%
  rast()



# Apply the functions across these rasters of parameter values
cci_5 <- terra::lapp(x = c(cci_max, cci_tau, cci_c), fun = agc_5)
cci_15 <- terra::lapp(x = c(cci_max, cci_tau, cci_c), fun = agc_15)
cci_30 <- terra::lapp(x = c(cci_max, cci_tau, cci_c), fun = agc_30)
cci_50 <- terra::lapp(x = c(cci_max, cci_tau, cci_c), fun = agc_50)
cci_100 <- terra::lapp(x = c(cci_max, cci_tau, cci_c), fun = agc_100)
cci_200 <- terra::lapp(x = c(cci_max, cci_tau, cci_c), fun = agc_200)


# SBM (Spawn et al.) — identical Chapman-Richards setup
sbm_max <- "./data_inputs/raster/SBM/SBM_max_1deg_tha.tif" %>%
  rast()
sbm_tau <- "./data_inputs/raster/SBM/SBM_tau_1deg.tif" %>%
  rast()
sbm_c <- "./data_inputs/raster/SBM/SBM_c_1deg.tif" %>%
  rast()

# Apply the functions across these rasters of parameter values
sbm_5 <- terra::lapp(x = c(sbm_max, sbm_tau, sbm_c), fun = agc_5)
sbm_15 <- terra::lapp(x = c(sbm_max, sbm_tau, sbm_c), fun = agc_15)
sbm_30 <- terra::lapp(x = c(sbm_max, sbm_tau, sbm_c), fun = agc_30)
sbm_50 <- terra::lapp(x = c(sbm_max, sbm_tau, sbm_c), fun = agc_50)
sbm_100 <- terra::lapp(x = c(sbm_max, sbm_tau, sbm_c), fun = agc_100)
sbm_200 <- terra::lapp(x = c(sbm_max, sbm_tau, sbm_c), fun = agc_200)



# Gap-fill at low tree cover: propagate values from the ~15 % treecover
# band scaled by local treecover fraction, then kNN-fill remaining gaps.
if(!"pot_treecover.tif" %in% list.files("./data_inputs/raster/")){
  pot_treecover <- rast("./data_inputs/raster/Total_potential.tif") %>% 
    project("EPSG:4326") %>% 
    resample(cci_max, method = "bilinear")
  
  fill_self_from_coarse <- function(
    r,
    factors = c(2, 4, 8), # coarsening pyramid
    agg_fun = mean,
    verbose = TRUE
  ) {
    stopifnot(nlyr(r) == 1)
    out <- r
    
    n_na <- function(z) sum(is.na(values(z)))
    if (verbose) message("Initial NA cells: ", n_na(out))
    
    for (f in factors) {
      if (verbose) message(" • Factor ", f, ": aggregate → resample → cover")
      # Aggregate to coarser grid (ignore NA within blocks)
      rc <- aggregate(r, fact = f, fun = agg_fun, na.rm = TRUE)
      # Bring coarse grid back to original with NEAREST (choose a single coarse value)
      rc_back <- resample(rc, r, method = "near")
      # Fill gaps in the working raster
      out_new <- cover(out, rc_back)
      filled_now <- n_na(out) - n_na(out_new)
      out <- out_new
      if (verbose) message(" filled ", filled_now, " cells | remaining NAs: ", n_na(out))
      if (all(!is.na(values(out)))) break
    }
    
    names(out) <- names(r)
    out
  }
  
  # Mask to the WorldClim land extent using a pre-downloaded bio raster.
  # wc_align is not yet defined here, so we load a WorldClim layer directly
  # and resample it to the pot_treecover grid before masking.
  wc_mask <- rast("./data_inputs/worldclim/climate/wc2.1_10m/wc2.1_10m_bio_1.tif") %>%
    resample(cci_max, method = "near")
  pot_treecover <- fill_self_from_coarse(pot_treecover, 
                                         factors = 2:40) %>% 
    mask(wc_mask)
  
  writeRaster(pot_treecover, 
              filename = "./data_inputs/raster/pot_treecover.tif",
              overwrite=T)
} else{
  pot_treecover <- rast("./data_inputs/raster/pot_treecover.tif")
}

pot_tc_less_than_15 <- pot_treecover < 15
dtc <- 3 # Define a percent threshold for the pixels (with data) that
# represent low density treecover.
pot_tc_15 <- ((pot_treecover < (15+dtc)) & (pot_treecover > (15-dtc))) & !is.na(cci_max)



# --- core helper: propagate ~15% values to <15% region, scaled by treecover ---
propagate_from_band <- function(
    src, # SpatRaster (1 layer): e.g., cci_5
    band_mask, # SpatRaster logical: TRUE in ~15% band (15% ± dtc)
    lt15_mask, # SpatRaster logical: TRUE where <15%
    pot_treecover, # SpatRaster: % treecover on same grid as src
    fill_into_source = TRUE
){
  stopifnot(nlyr(src) == 1)
  
  # Align masks and treecover to src grid
  if (!compareGeom(src, band_mask, stopOnError = FALSE)) {
    band_mask <- resample(band_mask, src, method = "near")
  }
  if (!compareGeom(src, lt15_mask, stopOnError = FALSE)) {
    lt15_mask <- resample(lt15_mask, src, method = "near")
  }
  if (!compareGeom(src, pot_treecover, stopOnError = FALSE)) {
    pot_treecover <- resample(pot_treecover, src, method = "bilinear")
  }
  
  # Keep only src values inside ~15% band
  src_band <- mask(src, band_mask, maskvalue = 0)
  
  # SOURCE points (band)
  src_df <- as.data.frame(c(src_band, pot_treecover), xy = TRUE, na.rm = TRUE)
  if (nrow(src_df) == 0) stop("propagate_from_band: no non-NA band values found.")
  names(src_df)[3:4] <- c("val", "tc_donor")
  
  # TARGET coords = all <15% pixels
  tgt_df <- as.data.frame(c(lt15_mask, pot_treecover), xy = TRUE, na.rm = FALSE)
  names(tgt_df)[3:4] <- c("lt", "tc_focal")
  tgt_df <- tgt_df[!is.na(tgt_df$lt) & tgt_df$lt, c("x","y","tc_focal")]
  if (nrow(tgt_df) == 0) {
    return(if (fill_into_source) src else mask(src, lt15_mask))
  }
  
  # Nearest donor (by XY) for each target
  nbr <- FNN::get.knnx(
    data = as.matrix(src_df[, c("x","y")]),
    query = as.matrix(tgt_df[, c("x","y")]),
    k = 1
  )
  donor_vals <- src_df$val[nbr$nn.index[,1]]
  donor_tc <- src_df$tc_donor[nbr$nn.index[,1]]
  focal_tc <- tgt_df$tc_focal
  
  # Scale by treecover ratio (guard against 0)
  donor_tc_safe <- pmax(donor_tc, 1e-6)
  scaled_vals <- donor_vals * (focal_tc / donor_tc_safe)
  
  # Build raster of propagated (scaled) values
  tgt_df$val <- scaled_vals
  prop <- rast(tgt_df[, c("x","y","val")], type = "xyz", crs = crs(src))
  prop <- resample(prop, src, method = "near")
  prop <- mask(prop, lt15_mask)
  
  if (fill_into_source) {
    cover(src, prop) # fill only where src is NA
  } else {
    prop # just the propagated layer
  }
}

# --- batch runner over a named list of rasters ---
fill_all <- function(
    rlist, # named list of SpatRasters (1 layer each)
    band_mask, # ~15% band mask
    lt15_mask, # <15% mask
    pot_treecover, # % treecover raster
    fill_into_source = TRUE,
    write_dir = NULL, # optional folder to write outputs
    file_suffix = "_filled_lt15_scaled"
){
  stopifnot(is.list(rlist), length(rlist) > 0, !is.null(names(rlist)))
  
  if (!is.null(write_dir)) dir.create(write_dir, showWarnings = FALSE, recursive = TRUE)
  
  out <- lapply(names(rlist), function(nm) {
    r <- rlist[[nm]]
    filled <- propagate_from_band(
      src = r,
      band_mask = band_mask,
      lt15_mask = lt15_mask,
      pot_treecover = pot_treecover,
      fill_into_source = fill_into_source
    )
    if (!is.null(write_dir)) {
      writeRaster(filled,
                  filename = file.path(write_dir, paste0(nm, file_suffix, ".tif")),
                  overwrite = TRUE)
    }
    filled
  })
  setNames(out, names(rlist))
}


cci_list <- list(
  cci_5 = cci_5,
  cci_15 = cci_15,
  cci_30 = cci_30,
  cci_50 = cci_50,
  cci_100 = cci_100,
  cci_200 = cci_200
)

sbm_list <- list(
  sbm_5 = sbm_5,
  sbm_15 = sbm_15,
  sbm_30 = sbm_30,
  sbm_50 = sbm_50,
  sbm_100 = sbm_100,
  sbm_200 = sbm_200
)

cci_filled <- fill_all(
  rlist = cci_list,
  band_mask = pot_tc_15,
  lt15_mask = pot_tc_less_than_15,
  pot_treecover = pot_treecover
)

sbm_filled <- fill_all(
  rlist = sbm_list,
  band_mask = pot_tc_15,
  lt15_mask = pot_tc_less_than_15,
  pot_treecover = pot_treecover
)

# Overwrite originals with treecover-filled versions
for (nm in names(cci_filled)) {
  assign(nm, cci_filled[[nm]], envir = .GlobalEnv)
}

for (nm in names(sbm_filled)) {
  assign(nm, sbm_filled[[nm]], envir = .GlobalEnv)
}




# kNN gap-fill any remaining NAs using WorldClim climate predictors

# 1) download & align WorldClim only once
# Download 19 Bioclim layers at 10' resolution (cached to disk after first run)
wc <- worldclim_global(var = "bio", res = 10, path = "./data_inputs/worldclim")

# Subset to temperature, precipitation, and their seasonality indices (bio1, bio4, bio12, bio15)
vars4 <- c("wc2.1_10m_bio_1",
           "wc2.1_10m_bio_4",
           "wc2.1_10m_bio_12",
           "wc2.1_10m_bio_15"
)
names(wc)
#wc <- wc[[vars4]]

wc_align <- resample(wc, cci_5, method = "bilinear")

# Collect CCI rasters into a named list for batch gap-filling
ages <- c(5, 15, 30, 50, 100, 200)
cci_names <- paste0("cci_", ages)
cci_list <- mget(cci_names)

# kNN regression gap-fill: predicts missing AGC values from WorldClim climate predictors
gap_fill_knn <- function(cci, age, clim_stack, k = 5) {
  # 1) pull out raster values + climate into one df
  df_cci <- as.data.frame(cci, xy = TRUE, na.rm = FALSE)
  df_clim <- as.data.frame(clim_stack, xy = TRUE, na.rm = FALSE)
  resp <- paste0("agb", age)
  names(df_cci)[3] <- resp
  
  dat <- inner_join(df_clim, df_cci, by = c("x","y"))
  
  # 2) find complete‐cases for predictors
  vars4 <- names(clim_stack)
  ok_preds <- complete.cases(dat[, vars4])
  
  # 3) define train / test indices
  train_idx <- which(!is.na(dat[[resp]]) & ok_preds)
  test_idx <- which(is.na(dat[[resp]]) & ok_preds)
  
  # 4) extract train/test matrices
  X_train <- dat[train_idx, vars4, drop = FALSE]
  y_train <- dat[train_idx, resp]
  X_test <- dat[test_idx, vars4, drop = FALSE]
  
  # 5) run kNN regression
  knn_out <- knn.reg(train = X_train,
                     test = X_test,
                     y = y_train,
                     k = k)
  
  # 6) fill in predictions
  dat[[resp]][test_idx] <- knn_out$pred
  
  # 7) rebuild raster (retain original grid)
  filled <- rast(dat[, c("x","y", resp)], type = "xyz")
  crs(filled) <- crs(cci)
  filled <- resample(filled, cci, method = "near")
  return(filled)
}

# Apply gap-fill across all CCI time horizons
cci_filled_knn <- Map(function(cci, age) {
  gap_fill_knn(cci, age, wc_align, k = 5)
}, cci_list, ages)

names(cci_filled_knn) <- paste0(cci_names, "_knn_filled")
list2env(cci_filled_knn, envir = .GlobalEnv)


# Overwrite originals with the kNN gap-filled versions, resampled to the common analysis grid
cci_5 <- cci_5_knn_filled %>%
  resample(crop_mask, method = "bilinear")
cci_15 <- cci_15_knn_filled %>%
  resample(crop_mask, method = "bilinear")
cci_30 <- cci_30_knn_filled %>%
  resample(crop_mask, method = "bilinear")
cci_50 <- cci_50_knn_filled %>%
  resample(crop_mask, method = "bilinear")
cci_100 <- cci_100_knn_filled %>%
  resample(crop_mask, method = "bilinear")
cci_200 <- cci_200_knn_filled %>%
  resample(crop_mask, method = "bilinear")

rm(list = names(cci_filled_knn))





# Apply the same kNN gap-fill to SBM

sbm_names <- paste0("sbm_", ages)
sbm_list <- mget(sbm_names)
# Apply gap-fill across all SBM time horizons
sbm_filled_knn <- Map(function(sbm, age) {
  gap_fill_knn(sbm, age, wc_align, k = 5)
}, sbm_list, ages)

names(sbm_filled_knn) <- paste0(sbm_names, "_knn_filled")
list2env(sbm_filled_knn, envir = .GlobalEnv)

# Overwrite originals with the kNN gap-filled versions, resampled to the common analysis grid
sbm_5 <- sbm_5_knn_filled %>%
  resample(crop_mask, method = "bilinear")
sbm_15 <- sbm_15_knn_filled %>%
  resample(crop_mask, method = "bilinear")
sbm_30 <- sbm_30_knn_filled %>%
  resample(crop_mask, method = "bilinear")
sbm_50 <- sbm_50_knn_filled %>%
  resample(crop_mask, method = "bilinear")
sbm_100 <- sbm_100_knn_filled %>%
  resample(crop_mask, method = "bilinear")
sbm_200 <- sbm_200_knn_filled %>%
  resample(crop_mask, method = "bilinear")

rm(list = names(sbm_filled_knn))


# =============================================================================
# GAP-FILL COVERAGE SUMMARY
# =============================================================================
# Reports the percentage of pixels with >10% agricultural cover (cropland +
# pasture combined, HYDE 2010) whose AGC value was supplied by the kNN
# gap-fill step.  Reported separately for CCI-BM and SBM since their NA
# patterns differ.
#
# The kNN gap-fill fills pixels that are NA in the original CR parameter
# rasters (cci_max / sbm_max), so is.na() on those rasters identifies exactly
# the pixels that were kNN-filled (as opposed to the earlier treecover-scaling
# step, which is not counted here).

# Reconstruct the NA mask from the original CR parameter rasters.
# These are still in memory and have not been overwritten.
cci_knn_mask <- is.na(cci_max)   # TRUE = pixel was kNN gap-filled for CCI-BM
sbm_knn_mask <- is.na(sbm_max)   # TRUE = pixel was kNN gap-filled for SBM

# Load HYDE agricultural fraction rasters and resample to the cci_max grid
ag_crop_r <- rast("./data_inputs/raster/cropland_frac_2010_HYDE.tif") %>%
  resample(cci_max, method = "bilinear")
ag_past_r <- rast("./data_inputs/raster/pasture_frac_2010_HYDE.tif") %>%
  resample(cci_max, method = "bilinear")
ag_total_r <- ag_crop_r + ag_past_r

# Cell area on the cci_max grid (ha)
area_cci <- cellSize(cci_max, unit = "ha")

# Area-weighted percentage — share of total agricultural area (crop + pasture,
# weighted by fractional cover × cell area) that falls on kNN gap-filled pixels.
pct_area_filled <- function(knn_mask, ag_frac, area) {
  total_ag_area  <- sum(values(ag_frac * area),              na.rm = TRUE)
  filled_ag_area <- sum(values(ifel(knn_mask, ag_frac * area, 0)), na.rm = TRUE)
  round(100 * filled_ag_area / total_ag_area, 2)
}

cci_pct_area <- pct_area_filled(cci_knn_mask, ag_total_r, area_cci)
sbm_pct_area <- pct_area_filled(sbm_knn_mask, ag_total_r, area_cci)

gap_fill_summary <- tibble(
  product                = c("CCI-BM", "SBM", "Average"),
  pct_ag_area_knn_filled = c(
    cci_pct_area, sbm_pct_area,
    round(mean(c(cci_pct_area, sbm_pct_area)), 2)
  )
)

write.csv(
  gap_fill_summary,
  file.path(figures_dir, "agc_knn_gap_fill_coverage.csv"),
  row.names = FALSE
)
message("Gap-fill coverage table written to ",
        file.path(figures_dir, "agc_knn_gap_fill_coverage.csv"))


# Robinson et al. 2025 Chapman-Richards curves with delta-method standard errors
# (doi:10.1038/s41558-025-02355-5).  CIs reflect parameter fitting uncertainty only.

# Parameters 
fact <- 10 # aggregation factor - use values > 1 for quicker testing
times <- c(0, 5, 15, 30, 50, 100, 200)
toPg <- 1e-9 # Mg → Pg

# Filepaths
base <- "./data_inputs/raster/Robinson et al. 2025"
A_hi_fp <- file.path(base, "cr_a.tif")
b_hi_fp <- file.path(base, "cr_b.tif")
k_hi_fp <- file.path(base, "cr_k.tif")
seA_hi_fp <- file.path(base, "cr_a_error.tif")
seB_hi_fp <- file.path(base, "cr_b_error.tif")
seK_hi_fp <- file.path(base, "cr_k_error.tif")

raw_fps <- list(
  crop_aban = "./data_inputs/raster/abandon_frac_1km.tif",
  crop_sparing = "./data_inputs/raster/spaArea_frac_per_grid_1km.tif",
  crop_relocation = "./data_inputs/raster/spaArea_frac_grid_beyerR_1km.tif",
  pasture_sparing = "./data_inputs/raster/pasture_sparing_hotspot_frac_Hayek.tif"
)
raw_fps$both_sparing <- c(raw_fps$crop_sparing, raw_fps$pasture_sparing)

# Load & aggregate CR rasters 
agg_one <- function(r) if(fact != 1) aggregate(r, fact, fun="mean", na.rm=TRUE) else(r)
A_ag <- agg_one(rast(A_hi_fp))
b_ag <- agg_one(rast(b_hi_fp))
k_ag <- agg_one(rast(k_hi_fp))
seA_ag <- agg_one(rast(seA_hi_fp))
seB_ag <- agg_one(rast(seB_hi_fp))
seK_ag <- agg_one(rast(seK_hi_fp))

area_ag <- cellSize(A_ag, unit="ha")

# Save rasters of the estimates for given years and their standard errors
m <- 0.67 # Defined in the study
exponent <- 1 / (1 - m)

cr_with_se <- function(A, b, k, seA, seB, seK, t) {
  u <- 1 - b * exp(-k * t)
  
  # mean in Mg/ha
  mean <- A * u^exponent
  
  # partial derivatives
  dA <- u^exponent
  db <- A * exponent * u^(exponent - 1) * (-exp(-k * t))
  dk <- A * exponent * u^(exponent - 1) * ( b * t * exp(-k * t))
  
  # delta-method variance → se
  var_y <- (dA * seA)^2 + (db * seB)^2 + (dk * seK)^2
  se_y <- sqrt(var_y)
  
  c(mean, se_y)
}

for (ti in times) {
  out <- terra::lapp(
    x = c(A_ag, b_ag, k_ag, seA_ag, seB_ag, seK_ag),
    fun = function(A, b, k, seA, seB, seK) cr_with_se(A, b, k, seA, seB, seK, ti)
  )
  assign(paste0("rob_", ti), out[[1]]) # mean (Mg/ha)
  assign(paste0("rob_", ti, "_se"), out[[2]]) # se (Mg/ha)
}



# Root mass fraction (woody components only) — gap-filled; note this is
# time-invariant and thus does not capture dynamic below-ground accumulation.
rmf_forest_r <- rast("./data_inputs/raster/belowground/forest_rmf_map_mean_20200517-0000000000-0000000000.tif") %>%
  merge(rast("./data_inputs/raster/belowground/forest_rmf_map_mean_20200517-0000000000-0000032768.tif")) %>%
  resample(A_ag, method="bilinear") %>%
  `/`(100) # convert %→fraction

rmf_shrub_r <- rast("./data_inputs/raster/belowground/shrub_rmf_map_mean_20200517-0000000000-0000000000.tif") %>%
  merge(rast("./data_inputs/raster/belowground/shrub_rmf_map_mean_20200517-0000000000-0000032768.tif")) %>%
  resample(A_ag, method="bilinear") %>%
  `/`(100) # convert %→fraction

# Sum shrub and forest rasters, treating NAs as zero
rmf_r <- mean(rmf_forest_r, rmf_shrub_r, na.rm = TRUE)

# Use aggregated uncertainty estimates
rmf_boot_r <- rast("./data_inputs/raster/belowground/rmf_boot_aggregate_vc_202311170000000000-0000000000.tif") %>%
  merge(rast("./data_inputs/raster/belowground/rmf_boot_aggregate_vc_202311170000000000-0000023296.tif")) %>%
  resample(A_ag, method="bilinear") %>%
  `/`(100)

# set all zeros (true absence of data) to NA
rmf_r[rmf_r == 0] <- NA
rmf_boot_r[rmf_boot_r == 0] <- NA

# similar gap‐fill function
gap_fill_knn <- function(x, clim_stack, k = 5) {
  # x: a single‐layer SpatRaster with NAs to fill
  df_x <- as.data.frame(x, xy=TRUE, na.rm=FALSE)
  df_clim <- as.data.frame(clim_stack, xy=TRUE, na.rm=FALSE)
  names(df_x)[3] <- "resp"
  
  # training = non‐NA resp & complete climate
  ok_preds <- complete.cases(df_clim)
  train_idx <- which(!is.na(df_x$resp) & ok_preds)
  test_idx <- which( is.na(df_x$resp) & ok_preds)
  
  X_train <- df_clim[train_idx, , drop=FALSE]
  y_train <- df_x$resp[train_idx]
  X_test <- df_clim[test_idx, , drop=FALSE]
  
  knn_out <- knn.reg(train = X_train, test = X_test, y = y_train, k = k)
  df_x$resp[test_idx] <- knn_out$pred
  
  # rebuild raster
  filled <- rast(df_x, type="xyz")
  crs(filled) <- crs(x)
  resample(filled, x, method="near")
}

# Realign WorldClim predictors to the RMF raster grid for gap-filling
wc_align <- resample(wc, rmf_r, method = "bilinear")

# Gap-fill RMF mean and bootstrap SD using kNN in climate space
rmf_filled <- gap_fill_knn(rmf_r, wc_align, k=5)
rmf_boot_filled <- gap_fill_knn(rmf_boot_r, wc_align, k=5)


# Note: gap-filling extends values into areas without agriculture (e.g., Antarctica,
# Sahara); these cells are masked out in all downstream calculations and do not
# affect results.
# plot(rmf_filled, main="RMF (filled)")
# plot(rmf_boot_filled, main="RMF bootstraps (filled)")

# Replace originals with gap-filled versions
rmf_r <- rmf_filled
rmf_boot_r <- rmf_boot_filled
rm(rmf_filled, rmf_boot_filled)




# =============================================================================
# SECTION 4 — SCENARIO CARBON ESTIMATES (POSTERIOR MEDIAN)
# =============================================================================
# Computes Pg C sequestered under each scenario at each time horizon using
# posterior-median SOC stock maps (0–100 cm and 0–30 cm) and the three AGC
# products.  pg_fun handles the forest/open biome-mask blending;
# pg_fun_single_rast is used for time-independent (init, max) stock maps.

pg_fun_single_rast <- function(stock_tif,
                               frac){
  stock_rast <- rast(stock_tif)
  mat <- as.matrix(stock_rast * frac * area_raster)
  return(sum(mat, na.rm = T) / (10^9))
}


# prep_for_map without external crop/pasture masks
prep_for_map_no_mask <- function(forest_path) {
  # 1) infer "open" path from the "forest" transition file
  open_path <- if (grepl("_CFN_", forest_path)) {
    sub("_CFN_", "_CON_", forest_path)
  } else if (grepl("_PFN_", forest_path)) {
    sub("_PFN_", "_PON_", forest_path)
  } else if (grepl("_CFA_", forest_path)) {
    sub("_CFA_", "_COA_", forest_path)
  } else if (grepl("_PFA_", forest_path)) {
    sub("_PFA_", "_POA_", forest_path)
  } else {
    stop("Filename must contain one of '_CFN_', '_PFN_', '_CFA_' or '_PFA_'")
  }
  
  # 2) load rasters
  r_forest <- rast(forest_path)
  r_open <- rast(open_path)
  mask_f <- forest_mask
  mask_nf <- nonforest_mask
  
  # 3) align masks to the transition raster grid
  mask_f <- resample(mask_f, r_forest, method = "near")
  mask_nf <- resample(mask_nf, r_forest, method = "near")
  
  # 4) binarize biome masks
  bin_f <- mask_f >= 1
  bin_nf <- mask_nf >= 1
  
  # 5) build combined raster: use "forest" values on forest pixels, "open" on nonforest
  combined <- ifel(
    bin_f,
    r_forest,
    ifel(bin_nf, r_open, NA)
  )
  
  # 6) convert to long data.frame for ggplot
  df <- as.data.frame(combined, xy = TRUE, na.rm = TRUE) %>%
    gather(key = "variable", value = "value", -x, -y)
  
  return(df)
}

pg_fun <- function(stock_tif,
                   frac_rast,
                   area_raster_rast = area_raster) {
  # 1) infer mask type from filename
  mask_type <- if (grepl("_CFN_", stock_tif) || grepl("_CFA_", stock_tif)) {
    "crop"
  } else if (grepl("_PFN_", stock_tif) || grepl("_PFA_", stock_tif)) {
    "pasture"
  } else {
    stop("stock_tif must contain either '_CFN_'/'_CFA_' or '_PFN_'/'_PFA_'")
  }
  
  # 2) build the long df of (x, y, value) for that transition
  df <- prep_for_map_no_mask(stock_tif) %>%
    select(x, y, value)
  
  # 3) pull in frac and cell area as data.frames
  frac_df <- as.data.frame(frac_rast, xy = TRUE, na.rm = FALSE) %>%
    rename(frac = 3)
  area_df <- as.data.frame(area_raster_rast, xy = TRUE, na.rm = FALSE) %>%
    rename(area_raster = 3)
  
  # 4) join them all and compute Pg C
  result <- df %>%
    left_join(frac_df, by = c("x", "y")) %>%
    left_join(area_df, by = c("x", "y")) %>%
    summarize(
      PgC = sum(value * frac * area_raster, na.rm = TRUE) / 1e9
    ) %>%
    pull(PgC)
  
  return(result)
}

# 0–100 cm depth range

out_df_0_100 <- tibble(age = c(5,15,30,50,100,200),
                       crop_aban_C_N_soil_PgC = c(pg_fun(file.path(stock_map_dir, "stock_map_adj_5_years_CFN_0_100.tif"), crop_aban_frac),
                                                  pg_fun(file.path(stock_map_dir, "stock_map_adj_15_years_CFN_0_100.tif"), crop_aban_frac),
                                                  pg_fun(file.path(stock_map_dir, "stock_map_adj_30_years_CFN_0_100.tif"), crop_aban_frac),
                                                  pg_fun(file.path(stock_map_dir, "stock_map_adj_50_years_CFN_0_100.tif"), crop_aban_frac),
                                                  pg_fun(file.path(stock_map_dir, "stock_map_adj_100_years_CFN_0_100.tif"), crop_aban_frac),
                                                  pg_fun(file.path(stock_map_dir, "stock_map_adj_200_years_CFN_0_100.tif"), crop_aban_frac)),
                       crop_sparing_C_N_soil_PgC = c(pg_fun(file.path(stock_map_dir, "stock_map_adj_5_years_CFN_0_100.tif"), crop_sparing_frac),
                                                     pg_fun(file.path(stock_map_dir, "stock_map_adj_15_years_CFN_0_100.tif"), crop_sparing_frac),
                                                     pg_fun(file.path(stock_map_dir, "stock_map_adj_30_years_CFN_0_100.tif"), crop_sparing_frac),
                                                     pg_fun(file.path(stock_map_dir, "stock_map_adj_50_years_CFN_0_100.tif"), crop_sparing_frac),
                                                     pg_fun(file.path(stock_map_dir, "stock_map_adj_100_years_CFN_0_100.tif"), crop_sparing_frac),
                                                     pg_fun(file.path(stock_map_dir, "stock_map_adj_200_years_CFN_0_100.tif"), crop_sparing_frac)),
                       crop_relocation_C_N_soil_PgC = c(pg_fun(file.path(stock_map_dir, "stock_map_adj_5_years_CFN_0_100.tif"), crop_relocation_frac),
                                                        pg_fun(file.path(stock_map_dir, "stock_map_adj_15_years_CFN_0_100.tif"), crop_relocation_frac),
                                                        pg_fun(file.path(stock_map_dir, "stock_map_adj_30_years_CFN_0_100.tif"), crop_relocation_frac),
                                                        pg_fun(file.path(stock_map_dir, "stock_map_adj_50_years_CFN_0_100.tif"), crop_relocation_frac),
                                                        pg_fun(file.path(stock_map_dir, "stock_map_adj_100_years_CFN_0_100.tif"), crop_relocation_frac),
                                                        pg_fun(file.path(stock_map_dir, "stock_map_adj_200_years_CFN_0_100.tif"), crop_relocation_frac)),
                       pasture_sparing_P_N_soil_PgC = c(pg_fun(file.path(stock_map_dir, "stock_map_adj_5_years_PFN_0_100.tif"), pasture_sparing_frac),
                                                        pg_fun(file.path(stock_map_dir, "stock_map_adj_15_years_PFN_0_100.tif"), pasture_sparing_frac),
                                                        pg_fun(file.path(stock_map_dir, "stock_map_adj_30_years_PFN_0_100.tif"), pasture_sparing_frac),
                                                        pg_fun(file.path(stock_map_dir, "stock_map_adj_50_years_PFN_0_100.tif"), pasture_sparing_frac),
                                                        pg_fun(file.path(stock_map_dir, "stock_map_adj_100_years_PFN_0_100.tif"), pasture_sparing_frac),
                                                        pg_fun(file.path(stock_map_dir, "stock_map_adj_200_years_PFN_0_100.tif"), pasture_sparing_frac)),
                       
                       
                       
                       crop_aban_Crop_soil_PgC = pg_fun_single_rast(file.path(stock_map_dir, "stock_map_adj_init_crop_0_100.tif"), crop_aban_frac),
                       crop_sparing_Crop_soil_PgC = pg_fun_single_rast(file.path(stock_map_dir, "stock_map_adj_init_crop_0_100.tif"), crop_sparing_frac),
                       crop_relocation_Crop_soil_PgC = pg_fun_single_rast(file.path(stock_map_dir, "stock_map_adj_init_crop_0_100.tif"), crop_relocation_frac),
                       pasture_sparing_Pasture_soil_PgC = pg_fun_single_rast(file.path(stock_map_dir, "stock_map_adj_init_pasture_0_100.tif"), pasture_sparing_frac)
                       
                       ,
                       crop_aban_agc_PgC = c(sum(as.matrix(cci_5 * crop_aban_frac * area_raster), na.rm = T)/ (10^9),
                                             sum(as.matrix(cci_15 * crop_aban_frac * area_raster), na.rm = T)/ (10^9),
                                             sum(as.matrix(cci_30 * crop_aban_frac * area_raster), na.rm = T)/ (10^9),
                                             sum(as.matrix(cci_50 * crop_aban_frac * area_raster), na.rm = T)/ (10^9),
                                             sum(as.matrix(cci_100 * crop_aban_frac * area_raster), na.rm = T)/ (10^9),
                                             sum(as.matrix(cci_200 * crop_aban_frac * area_raster), na.rm = T)/ (10^9)),
                       crop_sparing_agc_PgC = c(sum(as.matrix(cci_5 * crop_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                sum(as.matrix(cci_15 * crop_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                sum(as.matrix(cci_30 * crop_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                sum(as.matrix(cci_50 * crop_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                sum(as.matrix(cci_100 * crop_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                sum(as.matrix(cci_200 * crop_sparing_frac * area_raster), na.rm = T)/ (10^9)),
                       crop_relocation_agc_PgC = c(sum(as.matrix(cci_5 * crop_relocation_frac * area_raster), na.rm = T)/ (10^9),
                                                   sum(as.matrix(cci_15 * crop_relocation_frac * area_raster), na.rm = T)/ (10^9),
                                                   sum(as.matrix(cci_30 * crop_relocation_frac * area_raster), na.rm = T)/ (10^9),
                                                   sum(as.matrix(cci_50 * crop_relocation_frac * area_raster), na.rm = T)/ (10^9),
                                                   sum(as.matrix(cci_100 * crop_relocation_frac * area_raster), na.rm = T)/ (10^9),
                                                   sum(as.matrix(cci_200 * crop_relocation_frac * area_raster), na.rm = T)/ (10^9)),
                       pasture_sparing_agc_PgC = c(sum(as.matrix(cci_5 * pasture_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                   sum(as.matrix(cci_15 * pasture_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                   sum(as.matrix(cci_30 * pasture_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                   sum(as.matrix(cci_50 * pasture_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                   sum(as.matrix(cci_100 * pasture_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                   sum(as.matrix(cci_200 * pasture_sparing_frac * area_raster), na.rm = T)/ (10^9))
                       
                       ,
                       crop_aban_agc_sbm_PgC = c(sum(as.matrix(sbm_5 * crop_aban_frac * area_raster), na.rm = T)/ (10^9),
                                                 sum(as.matrix(sbm_15 * crop_aban_frac * area_raster), na.rm = T)/ (10^9),
                                                 sum(as.matrix(sbm_30 * crop_aban_frac * area_raster), na.rm = T)/ (10^9),
                                                 sum(as.matrix(sbm_50 * crop_aban_frac * area_raster), na.rm = T)/ (10^9),
                                                 sum(as.matrix(sbm_100 * crop_aban_frac * area_raster), na.rm = T)/ (10^9),
                                                 sum(as.matrix(sbm_200 * crop_aban_frac * area_raster), na.rm = T)/ (10^9)),
                       crop_sparing_agc_sbm_PgC = c(sum(as.matrix(sbm_5 * crop_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                    sum(as.matrix(sbm_15 * crop_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                    sum(as.matrix(sbm_30 * crop_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                    sum(as.matrix(sbm_50 * crop_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                    sum(as.matrix(sbm_100 * crop_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                    sum(as.matrix(sbm_200 * crop_sparing_frac * area_raster), na.rm = T)/ (10^9)),
                       crop_relocation_agc_sbm_PgC = c(sum(as.matrix(sbm_5 * crop_relocation_frac * area_raster), na.rm = T)/ (10^9),
                                                       sum(as.matrix(sbm_15 * crop_relocation_frac * area_raster), na.rm = T)/ (10^9),
                                                       sum(as.matrix(sbm_30 * crop_relocation_frac * area_raster), na.rm = T)/ (10^9),
                                                       sum(as.matrix(sbm_50 * crop_relocation_frac * area_raster), na.rm = T)/ (10^9),
                                                       sum(as.matrix(sbm_100 * crop_relocation_frac * area_raster), na.rm = T)/ (10^9),
                                                       sum(as.matrix(sbm_200 * crop_relocation_frac * area_raster), na.rm = T)/ (10^9)),
                       pasture_sparing_agc_sbm_PgC = c(sum(as.matrix(sbm_5 * pasture_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                       sum(as.matrix(sbm_15 * pasture_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                       sum(as.matrix(sbm_30 * pasture_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                       sum(as.matrix(sbm_50 * pasture_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                       sum(as.matrix(sbm_100 * pasture_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                       sum(as.matrix(sbm_200 * pasture_sparing_frac * area_raster), na.rm = T)/ (10^9))
)
gc()














# 0 to 30 cm version

out_df_0_30 <- tibble(age = c(5,15,30,50,100,200),
                      crop_aban_C_N_soil_PgC = c(pg_fun(file.path(stock_map_dir, "stock_map_adj_5_years_CFN_0_30.tif"), crop_aban_frac),
                                                 pg_fun(file.path(stock_map_dir, "stock_map_adj_15_years_CFN_0_30.tif"), crop_aban_frac),
                                                 pg_fun(file.path(stock_map_dir, "stock_map_adj_30_years_CFN_0_30.tif"), crop_aban_frac),
                                                 pg_fun(file.path(stock_map_dir, "stock_map_adj_50_years_CFN_0_30.tif"), crop_aban_frac),
                                                 pg_fun(file.path(stock_map_dir, "stock_map_adj_100_years_CFN_0_30.tif"), crop_aban_frac),
                                                 pg_fun(file.path(stock_map_dir, "stock_map_adj_200_years_CFN_0_30.tif"), crop_aban_frac)),
                      crop_sparing_C_N_soil_PgC = c(pg_fun(file.path(stock_map_dir, "stock_map_adj_5_years_CFN_0_30.tif"), crop_sparing_frac),
                                                    pg_fun(file.path(stock_map_dir, "stock_map_adj_15_years_CFN_0_30.tif"), crop_sparing_frac),
                                                    pg_fun(file.path(stock_map_dir, "stock_map_adj_30_years_CFN_0_30.tif"), crop_sparing_frac),
                                                    pg_fun(file.path(stock_map_dir, "stock_map_adj_50_years_CFN_0_30.tif"), crop_sparing_frac),
                                                    pg_fun(file.path(stock_map_dir, "stock_map_adj_100_years_CFN_0_30.tif"), crop_sparing_frac),
                                                    pg_fun(file.path(stock_map_dir, "stock_map_adj_200_years_CFN_0_30.tif"), crop_sparing_frac)),
                      crop_relocation_C_N_soil_PgC = c(pg_fun(file.path(stock_map_dir, "stock_map_adj_5_years_CFN_0_30.tif"), crop_relocation_frac),
                                                       pg_fun(file.path(stock_map_dir, "stock_map_adj_15_years_CFN_0_30.tif"), crop_relocation_frac),
                                                       pg_fun(file.path(stock_map_dir, "stock_map_adj_30_years_CFN_0_30.tif"), crop_relocation_frac),
                                                       pg_fun(file.path(stock_map_dir, "stock_map_adj_50_years_CFN_0_30.tif"), crop_relocation_frac),
                                                       pg_fun(file.path(stock_map_dir, "stock_map_adj_100_years_CFN_0_30.tif"), crop_relocation_frac),
                                                       pg_fun(file.path(stock_map_dir, "stock_map_adj_200_years_CFN_0_30.tif"), crop_relocation_frac)),
                      pasture_sparing_P_N_soil_PgC = c(pg_fun(file.path(stock_map_dir, "stock_map_adj_5_years_PFN_0_30.tif"), pasture_sparing_frac),
                                                       pg_fun(file.path(stock_map_dir, "stock_map_adj_15_years_PFN_0_30.tif"), pasture_sparing_frac),
                                                       pg_fun(file.path(stock_map_dir, "stock_map_adj_30_years_PFN_0_30.tif"), pasture_sparing_frac),
                                                       pg_fun(file.path(stock_map_dir, "stock_map_adj_50_years_PFN_0_30.tif"), pasture_sparing_frac),
                                                       pg_fun(file.path(stock_map_dir, "stock_map_adj_100_years_PFN_0_30.tif"), pasture_sparing_frac),
                                                       pg_fun(file.path(stock_map_dir, "stock_map_adj_200_years_PFN_0_30.tif"), pasture_sparing_frac))
                      
                      # ,
                      # crop_aban_Crop_soil_PgC = pg_fun_single_rast(file.path(stock_map_dir, "stock_map_adj_init_crop_0_30.tif"), crop_aban_frac),
                      # crop_sparing_Crop_soil_PgC = pg_fun_single_rast(file.path(stock_map_dir, "stock_map_adj_init_crop_0_30.tif"), crop_sparing_frac),
                      # crop_relocation_Crop_soil_PgC = pg_fun_single_rast(file.path(stock_map_dir, "stock_map_adj_init_crop_0_30.tif"), crop_relocation_frac),
                      # pasture_sparing_Pasture_soil_PgC = pg_fun_single_rast(file.path(stock_map_dir, "stock_map_adj_init_pasture_0_30.tif"), pasture_sparing_frac)
                      # Alternative remote sensing AGC approach (commented out)
                      ,
                      crop_aban_agc_PgC = c(sum(as.matrix(cci_5 * crop_aban_frac * area_raster), na.rm = T)/ (10^9),
                                            sum(as.matrix(cci_15 * crop_aban_frac * area_raster), na.rm = T)/ (10^9),
                                            sum(as.matrix(cci_30 * crop_aban_frac * area_raster), na.rm = T)/ (10^9),
                                            sum(as.matrix(cci_50 * crop_aban_frac * area_raster), na.rm = T)/ (10^9),
                                            sum(as.matrix(cci_100 * crop_aban_frac * area_raster), na.rm = T)/ (10^9),
                                            sum(as.matrix(cci_200 * crop_aban_frac * area_raster), na.rm = T)/ (10^9)),
                      crop_sparing_agc_PgC = c(sum(as.matrix(cci_5 * crop_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                               sum(as.matrix(cci_15 * crop_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                               sum(as.matrix(cci_30 * crop_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                               sum(as.matrix(cci_50 * crop_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                               sum(as.matrix(cci_100 * crop_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                               sum(as.matrix(cci_200 * crop_sparing_frac * area_raster), na.rm = T)/ (10^9)),
                      crop_relocation_agc_PgC = c(sum(as.matrix(cci_5 * crop_relocation_frac * area_raster), na.rm = T)/ (10^9),
                                                  sum(as.matrix(cci_15 * crop_relocation_frac * area_raster), na.rm = T)/ (10^9),
                                                  sum(as.matrix(cci_30 * crop_relocation_frac * area_raster), na.rm = T)/ (10^9),
                                                  sum(as.matrix(cci_50 * crop_relocation_frac * area_raster), na.rm = T)/ (10^9),
                                                  sum(as.matrix(cci_100 * crop_relocation_frac * area_raster), na.rm = T)/ (10^9),
                                                  sum(as.matrix(cci_200 * crop_relocation_frac * area_raster), na.rm = T)/ (10^9)),
                      pasture_sparing_agc_PgC = c(sum(as.matrix(cci_5 * pasture_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                  sum(as.matrix(cci_15 * pasture_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                  sum(as.matrix(cci_30 * pasture_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                  sum(as.matrix(cci_50 * pasture_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                  sum(as.matrix(cci_100 * pasture_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                  sum(as.matrix(cci_200 * pasture_sparing_frac * area_raster), na.rm = T)/ (10^9))
                      ,
                      crop_aban_agc_sbm_PgC = c(sum(as.matrix(sbm_5 * crop_aban_frac * area_raster), na.rm = T)/ (10^9),
                                                sum(as.matrix(sbm_15 * crop_aban_frac * area_raster), na.rm = T)/ (10^9),
                                                sum(as.matrix(sbm_30 * crop_aban_frac * area_raster), na.rm = T)/ (10^9),
                                                sum(as.matrix(sbm_50 * crop_aban_frac * area_raster), na.rm = T)/ (10^9),
                                                sum(as.matrix(sbm_100 * crop_aban_frac * area_raster), na.rm = T)/ (10^9),
                                                sum(as.matrix(sbm_200 * crop_aban_frac * area_raster), na.rm = T)/ (10^9)),
                      crop_sparing_agc_sbm_PgC = c(sum(as.matrix(sbm_5 * crop_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                   sum(as.matrix(sbm_15 * crop_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                   sum(as.matrix(sbm_30 * crop_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                   sum(as.matrix(sbm_50 * crop_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                   sum(as.matrix(sbm_100 * crop_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                   sum(as.matrix(sbm_200 * crop_sparing_frac * area_raster), na.rm = T)/ (10^9)),
                      crop_relocation_agc_sbm_PgC = c(sum(as.matrix(sbm_5 * crop_relocation_frac * area_raster), na.rm = T)/ (10^9),
                                                      sum(as.matrix(sbm_15 * crop_relocation_frac * area_raster), na.rm = T)/ (10^9),
                                                      sum(as.matrix(sbm_30 * crop_relocation_frac * area_raster), na.rm = T)/ (10^9),
                                                      sum(as.matrix(sbm_50 * crop_relocation_frac * area_raster), na.rm = T)/ (10^9),
                                                      sum(as.matrix(sbm_100 * crop_relocation_frac * area_raster), na.rm = T)/ (10^9),
                                                      sum(as.matrix(sbm_200 * crop_relocation_frac * area_raster), na.rm = T)/ (10^9)),
                      pasture_sparing_agc_sbm_PgC = c(sum(as.matrix(sbm_5 * pasture_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                      sum(as.matrix(sbm_15 * pasture_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                      sum(as.matrix(sbm_30 * pasture_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                      sum(as.matrix(sbm_50 * pasture_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                      sum(as.matrix(sbm_100 * pasture_sparing_frac * area_raster), na.rm = T)/ (10^9),
                                                      sum(as.matrix(sbm_200 * pasture_sparing_frac * area_raster), na.rm = T)/ (10^9))
)



# =============================================================================
# SECTION 5 — CREDIBLE INTERVALS FROM POSTERIOR SAMPLES
# =============================================================================
# Iterates over per-sample stock maps from 05_bayes_uncertainty.R to build
# a distribution of scenario Pg C deltas, then summarises into 50 % and 95 %
# credible intervals.  The uncertainty rasters may be at a non-integer multiple
# of the posterior-median grid (e.g. agg_fact = 25 vs 10), so resample() is
# used throughout rather than aggregate() to handle non-integer ratios.
# compute_delta_ci() handles both 0–100 cm and 0–30 cm depth ranges.


# ——————————————————————————————————————————————
# 0) parameters & sample dirs
# ——————————————————————————————————————————————
# Load one uncertainty raster to use as the target grid template for all
# resample() calls.  This ensures every coarsened raster aligns exactly
# to the uncertainty raster grid regardless of the agg_fact ratio.
.unc_template <- rast(
  list.files(
    list.dirs(uncertainty_dir, recursive = FALSE)[
      grepl("stock_map_", list.dirs(uncertainty_dir, recursive = FALSE))
    ][1],
    pattern = "\\.tif$", full.names = TRUE
  )[1]
)

# Specify folders to access
sample_dirs <- list.dirs(uncertainty_dir, recursive = FALSE) %>%
  keep(~ grepl("stock_map", .)) %>%
  keep(~ as.numeric(gsub(paste0(uncertainty_dir, "/stock_map_"), "", .)) %in% seq(1, 3000, by = 10)) # This is the particular iterations used (thinned by 10 across all chains). Would need to update under alternate number of samples!!!!!

# scenarios: fraction rasters resampled to the uncertainty grid
scenarios <- list(
  crop_aban = list(code = "CFN", frac10 = resample(crop_aban_frac,      .unc_template, method = "bilinear")),
  crop_sparing = list(code = "CFN", frac10 = resample(crop_sparing_frac,  .unc_template, method = "bilinear")),
  crop_relocation = list(code = "CFN", frac10 = resample(crop_relocation_frac, .unc_template, method = "bilinear")),
  pasture_sparing = list(code = "PFN", frac10 = resample(pasture_sparing_frac, .unc_template, method = "bilinear"))
)

# =============================================================================
# PRE-COMPUTE SCENARIO WEIGHTS (runs once, not per sample)
# =============================================================================
# The uncertainty stock rasters are already on the uncertainty grid
# (.unc_template) so they need no resampling.  We only resample the fraction
# rasters and biome masks once, then pre-compute frac × area weight vectors
# used in every Pg C calculation — avoiding repeated resample() calls inside
# per-sample loops.

# Biome masks resampled to uncertainty grid once
.mask_f  <- resample(forest_mask,    .unc_template, method = "near") >= 1
.mask_nf <- resample(nonforest_mask, .unc_template, method = "near") >= 1

# Cell area at uncertainty resolution (ha)
.area_unc <- cellSize(.unc_template, unit = "ha")
.area_vec  <- values(.area_unc)[, 1]

# Per-scenario: resample frac to uncertainty grid once and pre-compute
# frac × area weight vectors (in Pg C units per Mg/ha)
for (.scn in names(scenarios)) {
  .frac_unc <- resample(scenarios[[.scn]]$frac10, .unc_template, method = "bilinear")
  scenarios[[.scn]]$weight_vec <- values(.frac_unc)[, 1] * .area_vec / 1e9
}
rm(.scn, .frac_unc)

# Fast Pg C from an uncertainty stock raster (already on the uncertainty grid).
# For CFN/PFN rasters, blends forest and open biome values using pre-computed
# biome masks.  For init/max rasters (no transition type in name), no blending.
pg_unc <- function(stock_tif, weight_vec) {
  r <- rast(stock_tif)
  if (grepl("_CFN_", stock_tif) || grepl("_PFN_", stock_tif)) {
    open_path <- stock_tif %>%
      sub("_CFN_", "_CON_", .) %>%
      sub("_PFN_", "_PON_", .)
    r_o <- rast(open_path)
    v   <- values(ifel(.mask_f, r, ifel(.mask_nf, r_o, NA)))[, 1]
  } else {
    v <- values(r)[, 1]
  }
  sum(v * weight_vec, na.rm = TRUE)
}

# Compute CI tables for BOTH depth suffixes in a single pass over sample_dirs.
# This halves the number of disk reads compared to calling compute_delta_ci()
# twice with separate depth suffixes.
compute_delta_ci_both <- function(sample_dirs) {
  
  raw <- map_dfr(sample_dirs, function(dir) {
    scn_names <- names(scenarios)
    
    # Init stocks (time-independent) — read once per scenario per depth
    init_pgc <- setNames(
      lapply(scn_names, function(scn) {
        lu <- ifelse(grepl("^pasture", scn), "pasture", "crop")
        w  <- scenarios[[scn]]$weight_vec
        list(
          d100 = pg_unc(file.path(dir, sprintf("stock_map_adj_init_%s_0_100.tif", lu)), w),
          d30  = pg_unc(file.path(dir, sprintf("stock_map_adj_init_%s_0_30.tif",  lu)), w)
        )
      }),
      scn_names
    )
    
    # Age stocks — both depths in one inner loop
    map_dfr(ages, function(age) {
      map_dfr(scn_names, function(scn) {
        info <- scenarios[[scn]]
        w    <- info$weight_vec
        tibble(
          sample    = basename(dir),
          age       = age,
          scenario  = scn,
          delta_100 = pg_unc(file.path(dir, sprintf("stock_map_adj_%d_years_%s_0_100.tif", age, info$code)), w) - init_pgc[[scn]]$d100,
          delta_30  = pg_unc(file.path(dir, sprintf("stock_map_adj_%d_years_%s_0_30.tif",  age, info$code)), w) - init_pgc[[scn]]$d30
        )
      })
    })
  })
  
  # Add both_sparing = crop_sparing + pasture_sparing
  both <- raw %>%
    filter(scenario %in% c("crop_sparing", "pasture_sparing")) %>%
    group_by(sample, age) %>%
    summarize(delta_100 = sum(delta_100),
              delta_30  = sum(delta_30),
              scenario  = "both_sparing",
              .groups   = "drop")
  
  raw <- bind_rows(raw, both)
  
  # Summarise to quantiles then pivot wide
  summarise_ci <- function(delta_col) {
    raw %>%
      group_by(age, scenario) %>%
      summarize(
        lower  = quantile(.data[[delta_col]], 0.025),
        lo50   = quantile(.data[[delta_col]], 0.25),
        median = quantile(.data[[delta_col]], 0.50),
        hi50   = quantile(.data[[delta_col]], 0.75),
        upper  = quantile(.data[[delta_col]], 0.975),
        .groups = "drop"
      ) %>%
      pivot_wider(
        names_from  = scenario,
        values_from = c(lower, lo50, median, hi50, upper),
        names_glue  = "{scenario}_{.value}"
      ) %>%
      arrange(age)
  }
  
  list(ci_0_100 = summarise_ci("delta_100"),
       ci_0_30  = summarise_ci("delta_30"))
}

# Compute credible intervals for both depth ranges in a single pass
message("Computing credible intervals ...")
ci_results    <- compute_delta_ci_both(sample_dirs)
ci_wide_0_100 <- ci_results$ci_0_100
ci_wide_0_30  <- ci_results$ci_0_30

# Write credible interval tables
write.csv(ci_wide_0_100, file.path(figures_dir, "delta_ci_0_100.csv"), row.names = FALSE)
write.csv(ci_wide_0_30, file.path(figures_dir, "delta_ci_0_30.csv"), row.names = FALSE)







# =============================================================================
# SECTION 6 — SUPPLEMENTARY TABLES
# =============================================================================
# Formats credible interval outputs into manuscript-ready supplementary tables
# for both SOC (build_supp_table) and AGC (agc_supp_table), then writes to
# ./figures/ as CSV files.


build_supp_table <- function(ci_wide, depth_label) {
  ci_wide %>%
    # create both_sparing columns by summing the two scenarios
    mutate(
      both_sparing_lower = crop_sparing_lower + pasture_sparing_lower,
      both_sparing_median = crop_sparing_median + pasture_sparing_median,
      both_sparing_upper = crop_sparing_upper + pasture_sparing_upper
    ) %>%
    # Want a per hectare per year column for maximum too.
    mutate(both_sparing_median_per_ha = both_sparing_median/0.617/age,
           both_sparing_lower_per_ha = both_sparing_lower/0.617/age,
           both_sparing_upper_per_ha = both_sparing_upper/0.617/age) %>% 
    transmute(
      depth_interval = depth_label,
      age,
      crop_aban = sprintf("%0.1f [%0.1f, %0.1f]",
                          crop_aban_median,
                          crop_aban_lower,
                          crop_aban_upper),
      crop_sparing = sprintf("%0.1f [%0.1f, %0.1f]",
                             crop_sparing_median,
                             crop_sparing_lower,
                             crop_sparing_upper),
      crop_relocation = sprintf("%0.1f [%0.1f, %0.1f]",
                                crop_relocation_median,
                                crop_relocation_lower,
                                crop_relocation_upper),
      pasture_sparing = sprintf("%0.1f [%0.1f, %0.1f]",
                                pasture_sparing_median,
                                pasture_sparing_lower,
                                pasture_sparing_upper),
      both_sparing = sprintf("%0.1f [%0.1f, %0.1f]",
                             both_sparing_median,
                             both_sparing_lower,
                             both_sparing_upper),
      both_sparing_per_ha = sprintf("%0.2f [%0.2f, %0.2f]",
                                    both_sparing_median_per_ha,
                                    both_sparing_lower_per_ha,
                                    both_sparing_upper_per_ha)
    )
}

supp100 <- build_supp_table(ci_wide_0_100, "0–100 cm")
supp30 <- build_supp_table(ci_wide_0_30, "0–30 cm")
supp_table <- bind_rows(supp100, supp30)
supp_table


write.csv(
  supp_table,
  file.path(figures_dir, "supplementary_table_SOC_change.csv"),
  row.names = FALSE
)









# =============================================================================
# SECTION 7 — SCENARIO FIGURES
# =============================================================================
# Combines SOC credible intervals with mean AGC estimates (averaged across
# CCI-BM and SBM) to produce the main manuscript scenario panels and
# composite multi-panel figures.

# Average AGC across the two remote sensing products
# 1. Add both_sparing for each method
out2 <- out_df_0_100 %>%
  mutate(
    both_sparing_agc_PgC = crop_sparing_agc_PgC + pasture_sparing_agc_PgC,
    both_sparing_agc_sbm_PgC = crop_sparing_agc_sbm_PgC + pasture_sparing_agc_sbm_PgC
  )

# 2. Pivot longer across both methods, then summarise mean & sd
agc_stats <- out2 %>%
  pivot_longer(
    cols = -age,
    names_to = c("scenario", "method"),
    names_pattern = "(.*)_agc(?:_)?(sbm)?_PgC",
    values_to = "agc"
  ) %>%
  mutate(
    # method == NA → remote‐sensing, method == "sbm" → sbm
    method = ifelse(is.na(method), "rs", "sbm")
  ) %>%
  group_by(age, scenario) %>%
  summarise(
    agc_mean = mean(agc),
    agc_sd = sd(agc),
    agc_lower = agc_mean - agc_sd,
    agc_upper = agc_mean + agc_sd,
    .groups = "drop"
  )

# Similarly reshape the soil‐C output (ci_wide_0_100)
soil_long <- ci_wide_0_100 %>%
  #filter(age < 101) %>% # Can decide to stop at 100 years
  pivot_longer(
    cols = -age,
    names_to = c("scenario","stat"),
    names_pattern = "(.*)_(.*)",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  )
# soil_long now has: age, scenario, lower, lo50, median, hi50, upper

# 3. Join to soil_long and re‐factor the scenarios
df_plot <- soil_long %>%
  left_join(agc_stats, by = c("age", "scenario")) %>%
  mutate(
    scenario = factor(
      scenario,
      levels = c(
        "crop_aban", "crop_sparing", "crop_relocation",
        "pasture_sparing", "both_sparing"
      ),
      labels = c(
        "Current abandoned\ncropland",
        "Cropland\nsparing",
        "Cropland\nrelocation",
        "Pasture\nsparing",
        "Maximum\ncropland + pasture sparing"
      )
    )
  )

# Write the output table
agc_supp_table <- agc_stats %>%
  mutate(
    # main AGC columns → 1 decimal place
    agc_fmt = sprintf(
      "%.1f [%.1f, %.1f]",
      agc_mean, agc_lower, agc_upper
    ),
    # calculate per-hectare per-year rates just for both_sparing
    rate_per_ha_yr = if_else(
      scenario == "both_sparing",
      agc_mean / (0.617 * age),
      NA_real_
    ),
    rate_lower = if_else(
      scenario == "both_sparing",
      agc_lower / (0.617 * age),
      NA_real_
    ),
    rate_upper = if_else(
      scenario == "both_sparing",
      agc_upper / (0.617 * age),
      NA_real_
    ),
    # rate columns → 2 decimal places
    rate_fmt = if_else(
      scenario == "both_sparing",
      sprintf("%.2f [%.2f, %.2f]",
              rate_per_ha_yr, rate_lower, rate_upper),
      NA_character_
    )
  ) %>%
  select(age, agc_fmt, rate_fmt, scenario) %>% 
  pivot_wider(
    names_from = scenario,
    values_from = c(agc_fmt, rate_fmt)
  ) %>%
  select(age, 
         agc_fmt_crop_aban,
         agc_fmt_crop_sparing,
         agc_fmt_crop_relocation,
         agc_fmt_pasture_sparing,
         agc_fmt_both_sparing,
         rate_fmt_both_sparing)


# Write out
write.csv(
  agc_supp_table,
  file.path(figures_dir, "supplementary_table_AGC_change.csv"),
  row.names = FALSE
)

# Colour palette for soil C and AGC
agc_col <- "#228B22"
soil_col <- "#8B5A2B"

alt_make_panel <- function(dat) {
  scen <- unique(dat$scenario)
  ggplot(dat, aes(x = age)) +
    # Soil carbon ribbons & line
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = soil_col, alpha = 0.2) +
    geom_ribbon(aes(ymin = lo50, ymax = hi50 ), fill = soil_col, alpha = 0.2) +
    geom_line( aes(y = median), color = soil_col, size = 1) +
    geom_point( aes(y = median), color = soil_col) +
    
    # --- NEW: Aboveground ribbon ---
    geom_ribbon(aes(ymin = agc_lower, ymax = agc_upper),
                fill = agc_col, alpha = 0.2) +
    
    # AGC mean line & points
    geom_line( aes(y = agc_mean), color = agc_col,
               linetype = "dashed", size = 1) +
    geom_point( aes(y = agc_mean), color = agc_col) +
    
    labs(
      x = "Years since restoration",
      y = if (scen %in% 
              c("Current abandoned\ncropland", 
                "Maximum\ncropland + pasture sparing"))
        "Carbon change (PgC)" else NULL,
      title = scen
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "italic", hjust = 0.5, size = 10),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 7)
    )
}




# Build one panel per scenario
alt_panels <- df_plot %>% 
  split(.$scenario) %>% 
  purrr::map(alt_make_panel)

# Panel tag theme
tag_theme <- theme(
  plot.tag.position = c(0.02, 0.98), # x,y in npc
  plot.tag = element_text(
    face = "plain",
    size = 12,
  )
)


# Inset maps showing spatial distribution of each scenario fraction
frac_list <- list(
  "Current abandoned\ncropland" = crop_aban_frac,
  "Cropland\nsparing" = crop_sparing_frac,
  "Cropland\nrelocation" = crop_relocation_frac,
  "Pasture\nsparing" = pasture_sparing_frac,
  "Maximum\ncropland + pasture sparing" = (crop_sparing_frac + pasture_sparing_frac)
)
# Continental outlines for inset maps

#Defining a general CRS
mycrs <- "+proj=longlat +datum=WGS84 +no_defs"

#Downloading world map (maps package) and converting 'map' object to 'sf' object
world <- maps::map("world", fill=TRUE, plot = FALSE) 
world <- sf::st_as_sf(world)
world <- sf::st_transform(world, mycrs) #optional
world <- world %>%
  filter(ID != "Antarctica")

#The resulting map has invalid geometries. Solving issue
if(isFALSE(all(sf::st_is_valid(world)))){
  world <- sf::st_make_valid(world)
}

world_spat <- vect(world)
inset_map <- function(scen_label) {
  # 1) coarsen by 5×
  r5 <- aggregate(frac_list[[scen_label]], fact = 5, fun = mean, na.rm = TRUE)
  
  # 2) mask *ocean* (i.e. outside the world polygon) → these become NA,
  # but zeros *inside* the outline stay 0
  r5_land <- mask(r5, world_spat)
  
  # 3) to data.frame, keeping zeros on land but dropping ocean NAs
  df <- as.data.frame(r5_land, xy = TRUE, na.rm = FALSE)
  names(df)[3] <- "frac"
  df <- filter(df, !is.na(frac)) # now only land
  
  # 4) plot with viridis and transparent background
  ggplot(df, aes(x = x, y = y, fill = frac)) +
    geom_tile() +
    geom_sf(
      data = world, 
      inherit.aes = FALSE,
      fill = NA,
      color = NA,
      size = 0.2
    ) +
    scale_fill_viridis_c(
      option = "rocket",
      direction = -1,
      limits = c(0,1),
      trans = "sqrt",
      na.value = "transparent"
    ) +
    coord_sf(expand = FALSE) +
    theme_void() +
    theme(legend.position = "none") 
}

l1 <- 0.02
b1 <- 0.66
r1 <- 0.6
t1 <- 0.98

alt_p1 <- alt_panels[[1]] + coord_cartesian(ylim = c(0, 50)) + labs(tag = "A") + tag_theme # abandoned
alt_p2 <- alt_panels[[2]] + coord_cartesian(ylim = c(0, 50)) + labs(tag = "B") + tag_theme # cropland sparing
alt_p3 <- alt_panels[[3]] + coord_cartesian(ylim = c(0, 50)) + labs(tag = "C") + tag_theme # relocation
alt_p4 <- alt_panels[[4]] + coord_cartesian(ylim = c(0, 50)) + labs(tag = "D") + tag_theme # pasture sparing
alt_p5 <- alt_panels[[5]] + coord_cartesian(ylim = c(0, 80)) + labs(tag = "E") + tag_theme # both sparing


alt_p1a <- alt_p1 +
  # Soil carbon: solid line
  geom_segment(
    x = 2, xend = 48, y = 28, yend = 28,
    color = soil_col,
    size = 1,
    inherit.aes = FALSE
  ) +
  annotate(
    "text",
    x = 52, y = 28,
    label = "SOC",
    hjust = 0,
    size = 2.75
  ) +
  
  # Aboveground carbon: dashed line
  geom_segment(
    x = 2, xend = 48, y = 24, yend = 24,
    color = agc_col,
    linetype = 2,
    size = 1,
    inherit.aes = FALSE
  ) +
  annotate(
    "text",
    x = 52, y = 24,
    label = "AGC",
    hjust = 0,
    size = 2.75
  )

# Add little inset maps
alt_p1a_with_inset <- alt_p1a +
  inset_element(
    inset_map("Current abandoned\ncropland"),
    left = l1, 
    bottom = b1, 
    right = r1, 
    top = t1
  )

alt_p2_with_inset <- alt_p2 +
  inset_element(
    inset_map("Cropland\nsparing"),
    left = l1, 
    bottom = b1, 
    right = r1, 
    top = t1
  )


alt_p3_with_inset <- alt_p3 +
  inset_element(
    inset_map("Cropland\nrelocation"),
    left = l1, 
    bottom = b1, 
    right = r1, 
    top = t1
  )


alt_p4_with_inset <- alt_p4 +
  inset_element(
    inset_map("Pasture\nsparing"),
    left = l1, 
    bottom = b1, 
    right = r1, 
    top = t1
  )


alt_p5_with_inset <- alt_p5 +
  inset_element(
    inset_map("Maximum\ncropland + pasture sparing"),
    left = l1, 
    bottom = b1, 
    right = r1, 
    top = t1
  )

alt_layout <- "
ABX
CDE
"

alt_combined <- alt_p1a_with_inset + alt_p2_with_inset + alt_p3_with_inset + alt_p4_with_inset + alt_p5_with_inset +
  plot_layout(design = alt_layout)

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png")
    png(file.path(figures_dir, "five scenario no map_mean_ci_agc.png"), width = 6.5, height = 5, unit = "in", res = 440)
  else
    pdf(file.path(figures_dir, "five scenario no map_mean_ci_agc.pdf"), width = 6.5, height = 5)
  print(alt_combined)
  dev.off()
}

























# =============================================================================
# SECTION 7 — SOC LOSS MAP AND SANDERMAN COMPARISON
# =============================================================================
# Estimates pixel-level SOC losses attributable to cropland and pasture by
# comparing modelled stocks under current land use to a 100-year regrowth
# baseline.  The 100-year horizon is used as the natural-stock reference because
# <2% of training data come from ages >100 years, avoiding extrapolation.
# Calculations are restricted to 0–100 cm depth for the same reason.
# HYDE 2010 land-use fractions are used to match the Sanderman et al. comparison.

cfn100 <- rast(file.path(stock_map_dir, "stock_map_adj_100_years_CFN_0_100.tif"))
con100 <- rast(file.path(stock_map_dir, "stock_map_adj_100_years_CON_0_100.tif"))
pfn100 <- rast(file.path(stock_map_dir, "stock_map_adj_100_years_PFN_0_100.tif"))
pon100 <- rast(file.path(stock_map_dir, "stock_map_adj_100_years_PON_0_100.tif"))

# Average forest and pasture transitions within each biome type to produce
# a single "natural baseline" stock map at year 100
forest_100 <- (forest_mask * cfn100 + forest_mask * pfn100)/2
nonforest_100 <- (nonforest_mask * con100 + nonforest_mask * pon100)/2
both_100 <- mean(forest_100, nonforest_100, na.rm = T)

# Global natural SOC stock under 100-year regrowth baseline (Pg C)
sum(as.matrix(both_100 * area_raster), na.rm = T) / (10^9)
# About 1738 at 100 cm

# Current modelled SOC stocks under cropland and pasture land use
init_crop <- rast(file.path(stock_map_dir, "stock_map_adj_init_crop_0_100.tif"))
init_pasture<- rast(file.path(stock_map_dir, "stock_map_adj_init_pasture_0_100.tif"))

# HYDE 2010 land-use fractions — used to match the Sanderman et al. comparison year
crop_frac <- rast("./data_inputs/raster/cropland_frac_2010_HYDE.tif")
past_frac <- rast("./data_inputs/raster/pasture_frac_2010_HYDE.tif")
crop_frac <- resample(crop_frac, cfn100, method = "bilinear")
past_frac <- resample(past_frac, cfn100, method = "bilinear")

lu_frac <- crop_frac + past_frac
hist(as.matrix(lu_frac))

# Pixel-level SOC losses due to cropland and pasture relative to the natural baseline
loss_crop100 <- (both_100 - init_crop) * crop_frac
loss_past100 <- (both_100 - init_pasture) * past_frac

percent_losses <- (loss_crop100+loss_past100)/both_100 * 100
# Clamp extreme negative values (spurious apparent gains) to -50%
percent_losses <- ifel(percent_losses < -50, -50, percent_losses)

# Total global SOC loss (Pg C) — compare with Sanderman et al.
sum(as.matrix(loss_crop100 * area_raster + loss_past100 * area_raster), na.rm = T) / (10^9)



# Convert rasters to a data frame
df <- as.data.frame(c(percent_losses, lu_frac), xy = FALSE, na.rm = TRUE)
colnames(df) <- c("percent_loss", "lu_frac")

# Classify into low/high land-use fraction
df <- df %>%
  mutate(lu_class = ifelse(lu_frac <= 0.1, "≤10%", ">10%"))

ggplot(df, aes(x = percent_loss, fill = lu_class)) +
  geom_histogram(position = "stack", bins = 50, color = "black", linewidth = 0.1) +
  scale_fill_manual(values = c("≤10%" = "gray70", ">10%" = "firebrick")) +
  labs(
    x = "SOC loss (% of natural stock)",
    y = "Pixel count",
    fill = "Crop + pasture\nlandcover fraction"
  ) +
  theme_minimal()


loss_both_100 <- mean(c(loss_crop100, loss_past100), na.rm = T)


# Convert to dataframe (include coordinates)
r_df <- as.data.frame(c(loss_both_100, lu_frac), xy = TRUE, na.rm = TRUE)
colnames(r_df)[3] <- "value" # Rename value column
colnames(r_df)[4] <- "lu_frac" # Rename value column

r_df <- r_df %>% 
  mutate(lu_class = ifelse(lu_frac <= 0.1, "≤10%", ">10%"))





hist_plot <- ggplot(r_df, aes(x = value, fill = lu_class)) +
  geom_histogram(position = "stack", bins = 20, size = 0.1) +
  scale_fill_manual(values = c("≤10%" = "gray60", ">10%" = "firebrick")) +
  labs(x = "SOC loss (% of baseline)", y = NULL, fill = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.y=element_blank(),
    legend.position = "none",
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    panel.grid = element_blank(),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7)
  )

# Convert histogram plot to grob
hist_grob <- ggplotGrob(hist_plot)


# Dummy plot to extract the legend
legend_plot <- ggplot(r_df, aes(x = value, fill = lu_class)) +
  geom_histogram() +
  scale_fill_manual(values = c("≤10%" = "gray60", ">10%" = "firebrick"),
                    name = "Crop + pasture\nlandcover fraction") +
  #theme_minimal() +
  theme_minimal(base_size = 7) + # Base font size
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "lines") # Size of fill boxes
  )

# Extract just the legend as a grob
legend_grob <- get_legend(legend_plot)

# A version of theme_void that avoids some other formatting issues
theme_blank <- function(){
  theme(line = element_blank(),
        rect = element_blank(),
        
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        
        axis.text.x = element_blank(),
        axis.text.y = element_blank())
}

# Main map plot

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png")
    png(filename = file.path(figures_dir, "SOC loss map and histogram.png"), res = 440, units = "in", width = 8, height = 3)
  else
    pdf(file = file.path(figures_dir, "SOC loss map and histogram.pdf"), width = 8, height = 3)
  print(ggplot(r_df, aes(x = x, y = y, fill = value)) +
          geom_tile() +
          scale_fill_gradient2(
            low = "blue",
            mid = "grey90",
            high = "red",
            midpoint = 0,
            name = "SOC loss\n(Mg C/ha)"
          ) +
          coord_equal() +
          theme_minimal() +
          theme_blank() +
          annotation_custom(
            grob = hist_grob,
            xmin = -195, xmax = -90, # Adjust for position (lower left)
            ymin = -60, ymax = 20
          ) +
          annotation_custom(
            grob = legend_grob,
            xmin = -145, xmax = -105,
            ymin = -30, ymax = 10
          ))
  dev.off()
}








b_samp <- seq(1, 3000, by = 10)

# Set up vectors to recieve the estimates
global_stock <- c()
lost_stock <- c()

# Resample masks and fraction rasters to the uncertainty grid
forest_mask10    <- resample(forest_mask,    .unc_template, method = "near")
nonforest_mask10 <- resample(nonforest_mask, .unc_template, method = "near")
area_raster10    <- cellSize(forest_mask10, unit = "ha")
crop_frac10      <- resample(crop_frac, .unc_template, method = "bilinear")
past_frac10      <- resample(past_frac, .unc_template, method = "bilinear")

for(i in 1:length(b_samp)){
  cfn100 <- rast(file.path(uncertainty_dir, paste0("stock_map_", b_samp[i]), "stock_map_adj_100_years_CFN_0_100.tif"))
  con100 <- rast(file.path(uncertainty_dir, paste0("stock_map_", b_samp[i]), "stock_map_adj_100_years_CON_0_100.tif"))
  pfn100 <- rast(file.path(uncertainty_dir, paste0("stock_map_", b_samp[i]), "stock_map_adj_100_years_PFN_0_100.tif"))
  pon100 <- rast(file.path(uncertainty_dir, paste0("stock_map_", b_samp[i]), "stock_map_adj_100_years_PON_0_100.tif"))
  
  # Apply these estimates to the appropriate forest/non-forest biomes
  forest_100 <- (forest_mask10 * cfn100 + forest_mask10 * pfn100)/2
  nonforest_100 <- (nonforest_mask10 * con100 + nonforest_mask10 * pon100)/2
  
  # Add together to represent this 100 year "maximum" or "natural baseline" 
  # estimate
  both_100 <- mean(forest_100, nonforest_100, na.rm = T)
  
  # Global natural SOC stock for this posterior sample (Pg C)
  global_stock[i] <- sum(as.matrix(both_100 * area_raster10), na.rm = T) / (10^9)
  
  # Current modelled SOC stocks under cropland and pasture for this sample
  init_crop <- rast(file.path(uncertainty_dir, paste0("stock_map_", b_samp[i]), "stock_map_adj_init_crop_0_100.tif"))
  init_pasture<- rast(file.path(uncertainty_dir, paste0("stock_map_", b_samp[i]), "stock_map_adj_init_pasture_0_100.tif"))
  
  # Pixel-level SOC losses due to cropland and pasture relative to the natural baseline
  loss_crop100 <- (both_100 - init_crop) * crop_frac10
  loss_past100 <- (both_100 - init_pasture) * past_frac10
  
  # Total global SOC loss for this sample (Pg C)
  lost_stock[i] <- sum(as.matrix(loss_crop100 * area_raster10 + loss_past100 * area_raster10), na.rm = T) / (10^9)
}

par(mfrow=c(2,1))
hist(global_stock)
hist(lost_stock)

# Summarise posterior distributions and write to disk.
# All three quantiles use the sample loop for consistency.
sanderman_comparison <- tibble(
  metric = c("global_natural_stock_PgC", "lost_stock_PgC"),
  lower  = c(quantile(global_stock, 0.025), quantile(lost_stock, 0.025)),
  median = c(quantile(global_stock, 0.500), quantile(lost_stock, 0.500)),
  upper  = c(quantile(global_stock, 0.975), quantile(lost_stock, 0.975))
)

write.csv(
  sanderman_comparison,
  file.path(figures_dir, "sanderman_comparison.csv"),
  row.names = FALSE
)











# =============================================================================
# SECTION 7 (CONTINUED) — AGC/SOC RATIO MAP AND BIOME-LEVEL RATE COMPARISONS
# =============================================================================
# Constructs per-pixel maps comparing AGC and SOC accrual rates at year 50,
# and builds biome-level bar charts of accumulation rates for both pools.

# AGC at year 50: mean of CCI-BM and SBM products
delta_agc_px <- (cci_50 + sbm_50)/2

# SOC change at year 50 relative to initial stocks, per scenario
soil_crop_aban <- prep_for_map_diff(file.path(stock_map_dir, "stock_map_adj_50_years_CFN_0_100.tif"))
soil_crop_sparing <- prep_for_map_diff(file.path(stock_map_dir, "stock_map_adj_50_years_CFN_0_100.tif"))
soil_crop_relocation <- prep_for_map_diff(file.path(stock_map_dir, "stock_map_adj_50_years_CFN_0_100.tif"))
soil_pasture_sparing <- prep_for_map_diff(file.path(stock_map_dir, "stock_map_adj_50_years_PFN_0_100.tif"))

# rename the soil column
crop_aban_df <- soil_crop_aban %>% rename(soil_diff = value)
crop_sparing_df <- soil_crop_sparing %>% rename(soil_diff = value)
crop_relocation_df <- soil_crop_relocation %>% rename(soil_diff = value)
pasture_sparing_df <- soil_pasture_sparing %>% rename(soil_diff = value)

# Join the delta_agc_px into each scenario df
crop_aban_df <- crop_aban_df %>%
  left_join(
    as.data.frame(delta_agc_px, xy = TRUE, na.rm = TRUE) %>%
      #rename(agc_diff = cr_a),
      rename(agc_diff = agb50),
    by = c("x","y")
  )

crop_sparing_df <- crop_sparing_df %>%
  left_join(
    as.data.frame(delta_agc_px, xy = TRUE, na.rm = TRUE) %>%
      #rename(agc_diff = cr_a),
      rename(agc_diff = agb50),
    by = c("x","y")
  )

crop_relocation_df <- crop_relocation_df %>%
  left_join(
    as.data.frame(delta_agc_px, xy = TRUE, na.rm = TRUE) %>%
      #rename(agc_diff = cr_a),
      rename(agc_diff = agb50),
    by = c("x","y")
  )

pasture_sparing_df <- pasture_sparing_df %>%
  left_join(
    as.data.frame(delta_agc_px, xy = TRUE, na.rm = TRUE) %>%
      #rename(agc_diff = cr_a),
      rename(agc_diff = agb50),
    by = c("x","y")
  )

# 3) Build the “both_sparing” map by summing the two scenarios
both_sparing_df <- full_join(
  crop_sparing_df %>% select(x,y,soil_diff, agc_diff) %>% rename(soil1 = soil_diff, agc1 = agc_diff),
  pasture_sparing_df %>% select(x,y,soil_diff, agc_diff) %>% rename(soil2 = soil_diff, agc2 = agc_diff),
  by = c("x","y")
) %>%
  transmute(
    x, y,
    soil_diff = rowSums(cbind(soil1, soil2), na.rm = TRUE),
    agc_diff = rowSums(cbind(agc1, agc2 ), na.rm = TRUE)
  )


# 1) compute the raw ratio, and drop any zeros to avoid Inf/NaN
ratio_df <- both_sparing_df %>%
  filter(soil_diff > 0, agc_diff > 0) %>%
  mutate(ratio = agc_diff / soil_diff)

# 2) plot with a diverging brown–white–green scale, log‐transformed
# Prep continent outline
limitless_world <- sf::st_union(world)

p_ratio_map <- ggplot(ratio_df, aes(x = x, y = y, fill = ratio)) +
  geom_raster() +
  scale_fill_gradient2(
    low = soil_col, # brown
    mid = "grey",
    high = agc_col, # green
    midpoint = , # white at ratio = 1
    trans = "log10", # log10 scale
    limits = c(0.2, 5), # clamp to [0.5, 2]
    breaks = c(5, 2, 1, 0.5, 0.2),
    labels = c(
      ">5x AGC/SOC", # at ratio = 4
      "2x", # at ratio = 2
      "1:1 SOC/AGC", # at ratio = 1
      "2x", # at ratio = 0.5
      ">5x SOC/AGC" # at ratio = 0.25
    ),
    oob = squish, # squish anything outside
    name = 
      expression(atop("Accumulation",
                      "rate ratio (50 yr)"))
  ) +
  coord_equal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA)
  ) +
  geom_sf(data = limitless_world, inherit.aes = FALSE, fill = 'NA', color = 'grey', linewidth = 0.1)



# Crop/pasture-fraction-weighted mean SOC accrual rate map at year 50

# SOC change at year 50 relative to initial stocks, for cropland and pasture transitions
soil_crop <- prep_for_map_diff(file.path(stock_map_dir, "stock_map_adj_50_years_CFN_0_100.tif"))
soil_past <- prep_for_map_diff(file.path(stock_map_dir, "stock_map_adj_50_years_PFN_0_100.tif"))

# rename the soil column
crop_df <- soil_crop %>% rename(crop_soil_diff = value)
past_df <- soil_past %>% rename(past_soil_diff = value)

# Get equivalent dataframes for the crop and pasture fraction to get
# a weighted average

crop_frac_df <- as.data.frame(crop_frac, xy = TRUE, na.rm = F) %>%
  rename(crop_frac = cropland_113)
past_frac_df <- as.data.frame(past_frac, xy = TRUE, na.rm = F) %>%
  rename(past_frac = pasture_113)

both_df <- crop_df %>%
  mutate(past_soil_diff = past_df$past_soil_diff) %>%
  left_join(crop_frac_df) %>%
  left_join(past_frac_df)

# Now weight by crop or pasture fraction to get a mean value for the pixel
both_df <- both_df %>%
  mutate(
    sum_frac = crop_frac + past_frac,
    weighted_soil_diff = if_else(
      sum_frac > 0,
      (crop_frac * crop_soil_diff + past_frac * past_soil_diff) / sum_frac,
      NA_real_
    )
  ) %>%
  select(-sum_frac) # drop the helper column if you like

# Prep continent outline
limitless_world <- sf::st_union(world)

soil_average_rates_50 <- ggplot(both_df, aes(x = x, y = y, fill = weighted_soil_diff/50)) +
  geom_raster() +
  scale_fill_gradient2(
    low = "#74A5D4", # brown
    mid = "grey",
    high = "#8B5A2B", # green
    midpoint = 0,
    limits = c(-0.25, 2.25),
    oob = squish,
    name = expression(
      atop("SOC rate",
           "(Mg C " * ha^-1 * yr^-1 * ")")
    )
  ) +
  coord_equal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.title = element_text(hjust = 0) # left justify
  ) +
  geom_sf(data = limitless_world, inherit.aes = FALSE, fill = 'NA', color = 'grey', linewidth = 0.1)




# Panel B showing mean rates by biome ------------------------------------------

# ————————————————————————————————————————————————————————————————
# PARAMETERS
# ————————————————————————————————————————————————————————————————
years <- 50
# all of the Monte Carlo sample dirs
sample_dirs <- list.dirs(uncertainty_dir, recursive = FALSE) %>%
  keep(~ grepl("stock_map_", .)) %>%
  keep(~ as.numeric(sub(".*stock_map_","",.)) %in% seq(1,3000,by=10))

# ————————————————————————————————————————————————————————————————
# LOAD ECOREGION GROUPS & MAKE MASKS
# ————————————————————————————————————————————————————————————————
eco_v <- vect("./data_inputs/vector/Ecoregions2017/Ecoregions2017.shp") %>% makeValid()

# pull out names
biomes <- eco_v$BIOME_NAME

# 1) climate_zone
tropicals <- c(
  "Tropical & Subtropical Moist Broadleaf Forests",
  "Tropical & Subtropical Grasslands, Savannas & Shrublands",
  "Tropical & Subtropical Dry Broadleaf Forests",
  "Tropical & Subtropical Coniferous Forests",
  "Mangroves"
)
temperates <- c(
  "Temperate Conifer Forests",
  "Temperate Broadleaf & Mixed Forests",
  "Temperate Grasslands, Savannas & Shrublands",
  "Mediterranean Forests, Woodlands & Scrub"
)

eco_v$climate_zone <- NA_character_
eco_v$climate_zone[biomes %in% tropicals] <- "Tropical"
eco_v$climate_zone[biomes %in% temperates] <- "Temperate"

# 2) ecosystem_type
eco_v$ecosystem_type <- NA_character_
# everything with “Forest” in the name plus Mangroves
forest_idx <- grepl("Forest", biomes) | biomes == "Mangroves"
eco_v$ecosystem_type[forest_idx] <- "Forest"

# the explicit “Open” list
opens <- c(
  "Tropical & Subtropical Grasslands, Savannas & Shrublands",
  "Temperate Grasslands, Savannas & Shrublands",
  "Deserts & Xeric Shrublands",
  "Montane Grasslands & Shrublands",
  "Flooded Grasslands & Savannas",
  "Tundra"
)
eco_v$ecosystem_type[biomes %in% opens] <- "Open"

# 3) drop anything missing either one
eco_v <- eco_v[!is.na(eco_v$climate_zone) & !is.na(eco_v$ecosystem_type), ]

# make a coarse template from one of the frac10 rasters
template <- rast((sample_dirs[1] %>% 
                    list.files(pattern="stock_map_adj_init_crop_0_100.tif", full.names=TRUE))[1]) 

# rasterize each group to a binary mask
groups <- unique(paste(eco_v$climate_zone, eco_v$ecosystem_type, sep="_"))
group_masks <- set_names(
  lapply(groups, function(gr) {
    sv <- eco_v[paste(eco_v$climate_zone, eco_v$ecosystem_type, sep = "_") == gr, ]
    # use NA background so off-group cells don't drag down means
    m0 <- rasterize(sv, template, field = 1, background = NA_real_)
    # convert any covered cell to exactly 1
    m0[!is.na(m0)] <- 1
    m0
  }),
  groups
)

# A finer biome grouping (forest, savanna, shrubland/grassland) is used when
# Temperate_Savanna is present in the ecoregion data; the if/else blocks
# below handle both the coarser (Forest/Open) and finer grouping cases.

# ---- 0) load the potveg raster and pick a layer (change [[1]] if needed) ----
potveg <- rast("./data_inputs/raster/potveg_nc/vegtype_5min.nc")
potveg_layer <- potveg[[1]] # <-- change this if a different layer is desired

# ---- 1) define category IDs from potveg ----
forest_ids <- c(1,2,3,4,5,6,7,8)
savanna_ids <- c(9)
shrublandgrassland_ids <- c(10,11,12)

# if missing, set CRS (most likely lon/lat WGS84)
if (crs(potveg_layer) == "") {
  crs(potveg_layer) <- "EPSG:4326"
}

# ---- 3) project/align potveg to the template grid (nearest for categorical) ----
potveg_t <- project(potveg_layer, template, method = "near")

# helper to turn a set of integer classes into a 1/NA mask
make_type_mask <- function(r, ids) {
  classify(r, cbind(ids, rep(1, length(ids))), others = NA)
}

type_masks <- list(
  Forest = make_type_mask(potveg_t, forest_ids),
  Savanna = make_type_mask(potveg_t, savanna_ids),
  ShrublandGrassland = make_type_mask(potveg_t, shrublandgrassland_ids)
)


# ---- 4) climate masks from eco_v (Tropical / Temperate) ----
climate_groups <- c("Tropical", "Temperate")

climate_masks <- set_names(
  lapply(climate_groups, function(z) {
    sv <- eco_v[eco_v$climate_zone == z, ]
    m <- rasterize(sv, template, field = 1, background = NA_real_)
    m[!is.na(m)] <- 1
    m
  }),
  climate_groups
)

# ---- 5) build the 8 group names and masks by intersecting climate x type ----
groups <- as.vector(outer(climate_groups, names(type_masks), paste, sep = "_"))

group_masks <- set_names(
  lapply(groups, function(gr) {
    parts <- str_split(gr, "_")[[1]]
    clim <- parts[1]
    etype <- parts[2]
    # elementwise multiply 1/NA masks → 1 inside intersection, NA elsewhere
    m <- climate_masks[[clim]] * type_masks[[etype]]
    m
  }),
  groups
)

# --- Build a union of the 6 masks (Trop/Temp × Forest/Savanna/Shrub&Grassland) ---
gm <- rast(group_masks) # the 6 layers
union_mask <- app(gm, function(v) ifelse(any(v == 1, na.rm = TRUE), 1, NA_real_))
union_mask <- resample(union_mask, template, method = "near")

# --- Forest vs Open (Savanna + Shrub&Grassland) within the union ---
forest_global <- resample(type_masks$Forest, template, method = "near") * union_mask
open_global <- resample(
  ifel(!is.na(type_masks$Savanna) | !is.na(type_masks$ShrublandGrassland), 1, NA),
  template, method = "near"
) * union_mask


crop_frac_100 <- resample(crop_frac, 
                          rast(file.path(uncertainty_dir, "stock_map_2991", "stock_map_adj_50_years_CON_0_100.tif")),
                          method = "bilinear")
past_frac_100 <- resample(past_frac, 
                          rast(file.path(uncertainty_dir, "stock_map_2991", "stock_map_adj_50_years_CON_0_100.tif")),
                          method = "bilinear")
soil_rates <- map_dfr(groups, function(gr) {
  eco <- str_split(gr, "_")[[1]][2] # "Forest" or "Open"
  mask <- group_masks[[gr]] # a binary (0/1) SpatRaster
  
  # --- compute group‐level weights from the frac rasters ----------
  cf_vals <- values(crop_frac_100 * mask, mat = FALSE)
  pf_vals <- values(past_frac_100 * mask, mat = FALSE)
  
  sum_cf <- sum(cf_vals, na.rm = TRUE)
  sum_pf <- sum(pf_vals, na.rm = TRUE)
  total <- sum_cf + sum_pf
  
  # avoid divide‐by‐zero and pixels that have no ag lands
  if (total < 0.0001) { #(total < 0.001) {
    p_crop <- p_past <- NA_real_
  } else {
    p_crop <- sum_cf / total
    p_past <- sum_pf / total
  }
  
  # --- Monte Carlo draw of raw rates under each scenario ----------
  rates_crop <- map_dbl(sample_dirs, function(dir) {
    r0 <- rast(file.path(dir, "stock_map_adj_init_crop_0_100.tif"))
    r50 <- rast(file.path(dir, sprintf(
      "stock_map_adj_50_years_%s_0_100.tif",
      ifelse(eco=="Forest","CFN","CON")
    )))
    mean( values((r50 - r0) / years * mask), na.rm = TRUE )
  })
  
  rates_pasture <- map_dbl(sample_dirs, function(dir) {
    r0 <- rast(file.path(dir, "stock_map_adj_init_pasture_0_100.tif"))
    r50 <- rast(file.path(dir, sprintf(
      "stock_map_adj_50_years_%s_0_100.tif",
      ifelse(eco=="Forest","PFN","PON")
    )))
    mean( values((r50 - r0) / years * mask), na.rm = TRUE )
  })
  
  # --- build weighted draw vector ---
  weighted_rates <- p_crop * rates_crop + p_past * rates_pasture
  
  # --- assemble the tibble with three “pools” -------------
  tibble(
    group = gr,
    pool = c("SOC after cropland",
             "SOC after pasture",
             "SOC"),
    mean = c(mean(rates_crop),
             mean(rates_pasture),
             mean(weighted_rates)),
    lower = c(quantile(rates_crop, 0.025),
              quantile(rates_pasture,0.025),
              quantile(weighted_rates,0.025)),
    upper = c(quantile(rates_crop, 0.975),
              quantile(rates_pasture,0.975),
              quantile(weighted_rates,0.975))
  )
})

# AGC rates restricted to the same any-ag presence mask used for SOC rates above

# Binary any-ag presence mask (cells with any cropland or pasture cover)
ag_frac_mask <- (crop_frac_100 + past_frac_100) > 0.0001#0.001

# Annual AGC rate and inter-product variance (CCI-BM vs SBM) at year 50
mean_agc <- (cci_50 + sbm_50) / 2 / years
var_pix <- (cci_50 - sbm_50)^2 / 2 / years^2

# Area-weighted mean AGC rate per group, restricted to presence-masked pixels
agc_rates_any_ag <- imap_dfr(group_masks, function(mask, gr) {
  # Combine the group mask with the ag-presence mask and align to mean_agc grid
  combined_mask <- mask * ag_frac_mask
  cm_aligned <- resample(combined_mask, mean_agc, method = "near")
  
  # Group mean and SE weighted by cell area
  grp_mean <- global(mean_agc * cm_aligned, fun = "mean", na.rm = TRUE)[1,1]
  sum_var <- global(var_pix * cm_aligned, fun = "sum", na.rm = TRUE)[1,1]
  n_pix <- global(cm_aligned, fun = "sum", na.rm = TRUE)[1,1]
  grp_se <- sqrt(sum_var) / n_pix
  
  tibble(
    group = gr,
    pool = "AGC",
    mean = grp_mean,
    lower = grp_mean - 1.96 * grp_se,
    upper = grp_mean + 1.96 * grp_se
  )
})

# 4. Inspect
agc_rates_any_ag




# 1) bind agc & soil, split group names, keep trop/temps
if(!"Temperate_Savanna" %in% groups){
  
  plot_df3 <- bind_rows(soil_rates, agc_rates_any_ag) %>%
    separate(group, into = c("climate_zone","ecosystem_type"), sep = "_") %>%
    filter(climate_zone %in% c("Tropical","Temperate")) %>%
    
    # 2) drop the two pure‐SOC pools, keep only AGC & weighted‐SOC
    filter(pool %in% c("AGC", "SOC")) %>%
    
    # 3) factor‐order for facet & legend
    mutate(
      climate_zone = factor(climate_zone, levels = c("Tropical","Temperate")),
      ecosystem_type = factor(ecosystem_type, levels = c("Forest","Open")),
      pool = factor(pool, levels = c("AGC","SOC"))
    )
  
  # 4) plot
  p_mean_bars <- ggplot(plot_df3, aes(x = ecosystem_type, y = mean, fill = pool)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.6) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.75), width = 0.2) +
    scale_fill_manual(
      values = c(
        "AGC" = "#228B22",
        "SOC" = "#8B5A2B"
      ),
      name = NULL
    ) +
    ylim(c(0, 2)) +
    facet_wrap(~ climate_zone, nrow = 1) +
    labs(
      x = NULL,
      y = expression(
        atop("Biome-wide accumulation",
             "rate (Mg C " * ha^-1 * yr^-1 * ")"
        )
      )
    ) +
    guides(
      fill = guide_legend(
        direction = "horizontal", title = NULL,
        nrow = 1, byrow = TRUE
      )
    ) +
    theme_classic(base_size = 14) +
    theme(
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_text(margin = margin(r = 10)),
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key.width = unit(1.5, "lines")
    )
  
} else{
  # 1) build plotting df (now with 4 ecosystem types)
  plot_df3 <- bind_rows(soil_rates, agc_rates_any_ag) %>%
    separate(group, into = c("climate_zone","ecosystem_type"), sep = "_", remove = FALSE) %>%
    filter(climate_zone %in% c("Tropical","Temperate")) %>%
    filter(pool %in% c("AGC", "SOC")) %>%
    mutate(ecosystem_type = gsub("ShrublandGrassland", "Shrub & Grassland", ecosystem_type)) %>% 
    mutate(
      climate_zone = factor(climate_zone, levels = c("Temperate","Tropical")),
      ecosystem_type = factor(
        ecosystem_type,
        levels = c("Forest","Savanna","Shrub & Grassland")
      ),
      pool = factor(pool, levels = c("AGC","SOC"))
    ) %>%
    # drop rows that didn't match the four ecosystem types (just in case)
    filter(!is.na(ecosystem_type))
  
  # 2) choose a sensible y limit from the data (keeps things from clipping)
  ymax <- max(plot_df3$upper, na.rm = TRUE)
  ymax <- ceiling(ymax * 1.05 * 10) / 10 # small headroom, round up
  
  # 3) plot
  p_mean_bars <- ggplot(plot_df3, aes(x = ecosystem_type, y = mean, fill = pool)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.6) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.75), width = 0.2) +
    scale_fill_manual(
      values = c("AGC" = "#228B22", "SOC" = "#8B5A2B"),
      name = NULL
    ) +
    coord_cartesian(ylim = c(0, ymax)) +
    facet_wrap(~ climate_zone, nrow = 1) +
    labs(
      x = NULL,
      y = expression(
        atop("Biome-wide accumulation",
             "rate (Mg C " * ha^-1 * yr^-1 * ")")
      )
    ) +
    guides(fill = guide_legend(direction = "horizontal", nrow = 1, byrow = TRUE)) +
    theme_classic(base_size = 14) +
    theme(
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_text(margin = margin(r = 10)),
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key.width = unit(1.5, "lines")
    )
  
}

# display it
p_mean_bars


# This version of the panel compares a biome-wide rate with average
# rates weighted across the maximum scenario.


# Area-weighted mean over an intersection of a group mask and a presence mask (both binary 0/1).
# No fraction weighting here—just presence/absence × area.
aw_mean_presence <- function(rate_vals, weight_vec) {
  # rate_vals: numeric vector of rate values (already on the uncertainty grid)
  # weight_vec: pre-computed group_mask * presence_mask * area vector
  W <- sum(weight_vec, na.rm = TRUE)
  if (is.na(W) || W == 0) return(NA_real_)
  sum(rate_vals * weight_vec, na.rm = TRUE) / W
}

aw_mean_fraction <- function(rate_vals, weight_vec) {
  # weight_vec: pre-computed group_mask * frac * area vector
  W <- sum(weight_vec, na.rm = TRUE)
  if (is.na(W) || W == 0) return(NA_real_)
  sum(rate_vals * weight_vec, na.rm = TRUE) / W
}

# Per-draw SOC rate rasters (Mg C ha^-1 yr^-1) for CROPLAND / PASTURE transitions.
# Biome masks are pre-resampled outside this function and passed in as mF/mO
# to avoid repeating the resample() call for every sample.
soc_rate_rasters <- function(dir, T_years, depth_suffix, mF, mO) {
  r0c  <- rast(file.path(dir, sprintf("stock_map_adj_init_crop_%s.tif",         depth_suffix)))
  rCFN <- rast(file.path(dir, sprintf("stock_map_adj_%d_years_CFN_%s.tif", T_years, depth_suffix)))
  rCON <- rast(file.path(dir, sprintf("stock_map_adj_%d_years_CON_%s.tif", T_years, depth_suffix)))
  r0p  <- rast(file.path(dir, sprintf("stock_map_adj_init_pasture_%s.tif",      depth_suffix)))
  rPFN <- rast(file.path(dir, sprintf("stock_map_adj_%d_years_PFN_%s.tif", T_years, depth_suffix)))
  rPON <- rast(file.path(dir, sprintf("stock_map_adj_%d_years_PON_%s.tif", T_years, depth_suffix)))
  
  rate_crop <- ifel(mF, (rCFN - r0c)/T_years, ifel(mO, (rCON - r0c)/T_years, NA))
  rate_past <- ifel(mF, (rPFN - r0p)/T_years, ifel(mO, (rPON - r0p)/T_years, NA))
  
  list(rate_crop = rate_crop, rate_past = rate_past)
}

# ------------------------ ANY-AG PRESENCE MASKS ------------------------
crop_any  <- (crop_frac_100 > 0) * union_mask
past_any  <- (past_frac_100 > 0) * union_mask
any_ag_any <- ((crop_frac_100 + past_frac_100) > 0) * union_mask

# ── Pre-compute sample-invariant quantities for soc_rate_rasters ─────────────
# Biome masks resampled to the uncertainty grid once
.rt_unc <- rast(file.path(sample_dirs[[1]], "stock_map_adj_init_crop_0_100.tif"))
.mF_unc <- resample(forest_mask,    .rt_unc, method = "near") >= 1
.mO_unc <- resample(nonforest_mask, .rt_unc, method = "near") >= 1
.A_unc  <- values(cellSize(.rt_unc, unit = "ha"))[, 1]

# Fraction rasters resampled to uncertainty grid once
.cf_unc       <- values(resample(crop_frac_100,  .rt_unc, method = "bilinear"))[, 1]
.pf_unc       <- values(resample(past_frac_100,  .rt_unc, method = "bilinear"))[, 1]
.den_frac_unc <- .cf_unc + .pf_unc

# Restoration weights resampled to uncertainty grid once
.cs_unc  <- values(resample(scenarios$crop_sparing$frac10  * union_mask, .rt_unc, method = "bilinear"))[, 1]
.ps_unc  <- values(resample(scenarios$pasture_sparing$frac10 * union_mask, .rt_unc, method = "bilinear"))[, 1]
.den_rest_unc <- .cs_unc + .ps_unc

# Presence masks resampled to uncertainty grid once
.crop_any_unc  <- values(resample(crop_any,  .rt_unc, method = "near"))[, 1]
.past_any_unc  <- values(resample(past_any,  .rt_unc, method = "near"))[, 1]
.any_ag_unc    <- values(resample(any_ag_any, .rt_unc, method = "near"))[, 1]
.union_unc     <- values(resample(union_mask, .rt_unc, method = "near"))[, 1]

# Per-group weight vectors: group_mask * presence_mask * area
# Pre-compute for each group × pool combination used in soil_rates_any_ag
.gw <- lapply(group_masks, function(gm) {
  gv <- values(resample(gm, .rt_unc, method = "near"))[, 1]
  gv[is.na(gv)] <- 0
  list(
    crop   = gv * replace(.crop_any_unc,  is.na(.crop_any_unc),  0) * .A_unc,
    past   = gv * replace(.past_any_unc,  is.na(.past_any_unc),  0) * .A_unc,
    any_ag = gv * replace(.any_ag_unc,    is.na(.any_ag_unc),    0) * .A_unc,
    cs     = gv * replace(.cs_unc,        is.na(.cs_unc),        0) * .A_unc,
    ps     = gv * replace(.ps_unc,        is.na(.ps_unc),        0) * .A_unc,
    all    = gv * replace(.den_rest_unc,  is.na(.den_rest_unc),  0) * .A_unc
  )
})

# ---------------- SOC rates: any-ag presence ----------------
# Single pass over sample_dirs per group: read each sample once,
# compute all three rates (crop, pasture, mosaic) simultaneously.
soil_rates_any_ag <- map_dfr(groups, function(gr) {
  gw <- .gw[[gr]]
  
  per_sample <- map_dfr(sample_dirs, function(dir) {
    rr <- soc_rate_rasters(dir, T_years = years, depth_suffix = "0_100",
                           mF = .mF_unc, mO = .mO_unc)
    vc <- values(rr$rate_crop)[, 1]
    vp <- values(rr$rate_past)[, 1]
    
    # Fraction-blended mosaic rate
    vm <- ifelse(is.na(.den_frac_unc) | .den_frac_unc <= 0, NA_real_,
                 (vc * .cf_unc + vp * .pf_unc) / .den_frac_unc)
    
    tibble(
      rc = aw_mean_presence(vc, gw$crop),
      rp = aw_mean_presence(vp, gw$past),
      rs = aw_mean_presence(vm, gw$any_ag)
    )
  })
  
  tibble(
    group = gr,
    pool  = c("SOC after cropland", "SOC after pasture", "SOC"),
    mean  = c(mean(per_sample$rc, na.rm = TRUE),
              mean(per_sample$rp, na.rm = TRUE),
              mean(per_sample$rs, na.rm = TRUE)),
    lower = c(quantile(per_sample$rc, 0.025, na.rm = TRUE),
              quantile(per_sample$rp, 0.025, na.rm = TRUE),
              quantile(per_sample$rs, 0.025, na.rm = TRUE)),
    upper = c(quantile(per_sample$rc, 0.975, na.rm = TRUE),
              quantile(per_sample$rp, 0.975, na.rm = TRUE),
              quantile(per_sample$rs, 0.975, na.rm = TRUE))
  )
})


# AGC rates (presence-masked) — mirrors the SOC any-ag calculation above
agc_rates_any_ag <- imap_dfr(group_masks, function(mask_g, gr) {
  m_group <- resample(mask_g, mean_agc, method = "near")
  m_any <- resample(any_ag_any, mean_agc, method = "near")
  A <- cellSize(mean_agc, unit = "ha")
  
  W <- global(m_group * m_any * A, fun = "sum", na.rm = TRUE)[1,1]
  if (is.na(W) || W == 0) {
    return(tibble(group = gr, pool = "AGC",
                  mean = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  
  M <- global(mean_agc * m_group * m_any * A, fun = "sum", na.rm = TRUE)[1,1]
  Vw <- global(var_pix * (m_group * m_any * A)^2, fun = "sum", na.rm = TRUE)[1,1]
  
  mu <- M / W
  se <- sqrt(Vw) / W
  
  tibble(group = gr, pool = "AGC",
         mean = mu, lower = mu - 1.96*se, upper = mu + 1.96*se)
})

# --------------------------------------------------------------------
# 2) SOC (maximum_weighted) — restoration-weighted, with per-pixel
# blending by (cs, ps): (rc*cs + rp*ps)/(cs+ps)
# --------------------------------------------------------------------
soil_rates_maximum_weighted <- map_dfr(groups, function(gr) {
  gw <- .gw[[gr]]
  
  per_sample <- map_dfr(sample_dirs, function(dir) {
    rr <- soc_rate_rasters(dir, T_years = years, depth_suffix = "0_100",
                           mF = .mF_unc, mO = .mO_unc)
    vc <- values(rr$rate_crop)[, 1]
    vp <- values(rr$rate_past)[, 1]
    
    # Restoration-weighted blend
    vm <- ifelse(is.na(.den_rest_unc) | .den_rest_unc <= 0, NA_real_,
                 (vc * .cs_unc + vp * .ps_unc) / .den_rest_unc)
    
    tibble(
      rc = aw_mean_fraction(vc, gw$cs),
      rp = aw_mean_fraction(vp, gw$ps),
      rs = aw_mean_fraction(vm, gw$all)
    )
  })
  
  tibble(
    group = gr,
    pool  = c("SOC after cropland", "SOC after pasture", "SOC (restoration areas)"),
    mean  = c(mean(per_sample$rc, na.rm = TRUE),
              mean(per_sample$rp, na.rm = TRUE),
              mean(per_sample$rs, na.rm = TRUE)),
    lower = c(quantile(per_sample$rc, 0.025, na.rm = TRUE),
              quantile(per_sample$rp, 0.025, na.rm = TRUE),
              quantile(per_sample$rs, 0.025, na.rm = TRUE)),
    upper = c(quantile(per_sample$rc, 0.975, na.rm = TRUE),
              quantile(per_sample$rp, 0.975, na.rm = TRUE),
              quantile(per_sample$rs, 0.975, na.rm = TRUE))
  )
})

# Per-pixel AGC mean & variance (per-ha per-yr) — restoration-weighted, paralleling SOC above
mean_agc <- (cci_50 + sbm_50) / 2 / years
var_pix <- (cci_50 - sbm_50)^2 / 2 / years^2

# --- AGC over restoration areas (max sparing): area×fraction–weighted means by group ---
agc_rates_maximum_weighted <- purrr::imap_dfr(group_masks, function(mask_g, gr) {
  # restoration weights (crop + pasture sparing fractions)
  w_all <- scenarios$crop_sparing$frac10 + scenarios$pasture_sparing$frac10
  
  # align everything to the AGC grid
  m <- terra::resample(mask_g, mean_agc, method = "near")
  w <- terra::resample(w_all, mean_agc, method = "bilinear")
  
  # Restrict to the analysis union mask if present (intersection of all biome groups)
  if (exists("union_mask")) {
    u <- terra::resample(union_mask, mean_agc, method = "near")
    m <- m * u
  }
  
  # cell areas (ha) on this grid
  A <- terra::cellSize(mean_agc, unit = "ha")
  
  # denominator = total area-weighted restoration weight
  W <- terra::global(w * m * A, fun = "sum", na.rm = TRUE)[1,1]
  if (is.na(W) || W == 0) {
    return(tibble::tibble(
      group = gr, pool = "AGC (restoration areas)",
      mean = NA_real_, lower = NA_real_, upper = NA_real_
    ))
  }
  
  # numerator = area×fraction–weighted sum of per-ha rates
  M <- terra::global(mean_agc * w * m * A, fun = "sum", na.rm = TRUE)[1,1]
  
  # variance of the weighted mean (add variances with squared weights)
  Vw <- terra::global(var_pix * (w * m * A)^2, fun = "sum", na.rm = TRUE)[1,1]
  
  mu <- M / W
  se <- sqrt(Vw) / W
  
  tibble::tibble(
    group = gr,
    pool = "AGC (restoration areas)",
    mean = mu,
    lower = mu - 1.96 * se,
    upper = mu + 1.96 * se
  )
})

# --------------------------------------------------------------------
# 3) GLOBAL means (SOC + AGC)
# - SOC any_ag: presence‑masked; per‑pixel blend by crop/past *fractions*
# - SOC max: restoration‑weighted; per‑pixel blend by cs/ps
# - AGC any_ag: presence‑masked (no fraction weighting)
# - AGC max: restoration‑weighted (cs+ps)
# --------------------------------------------------------------------

# Global weight vectors (union_mask only, no group subdivision)
# These reuse the pre-computed vectors already on the uncertainty grid.
.w_union_crop   <- replace(.union_unc, is.na(.union_unc), 0) *
  replace(.crop_any_unc,  is.na(.crop_any_unc),  0) * .A_unc
.w_union_past   <- replace(.union_unc, is.na(.union_unc), 0) *
  replace(.past_any_unc,  is.na(.past_any_unc),  0) * .A_unc
.w_union_any_ag <- replace(.union_unc, is.na(.union_unc), 0) *
  replace(.any_ag_unc,    is.na(.any_ag_unc),    0) * .A_unc
.w_union_cs     <- replace(.union_unc, is.na(.union_unc), 0) *
  replace(.cs_unc, is.na(.cs_unc), 0) * .A_unc
.w_union_ps     <- replace(.union_unc, is.na(.union_unc), 0) *
  replace(.ps_unc, is.na(.ps_unc), 0) * .A_unc
.w_union_all    <- replace(.union_unc, is.na(.union_unc), 0) *
  replace(.den_rest_unc, is.na(.den_rest_unc), 0) * .A_unc

z2na <- function(x) if (is.na(x) || x == 0) NA_real_ else x

soc_draws_global <- function(dir) {
  rr <- soc_rate_rasters(dir, T_years = years, depth_suffix = "0_100",
                         mF = .mF_unc, mO = .mO_unc)
  vc <- values(rr$rate_crop)[, 1]
  vp <- values(rr$rate_past)[, 1]
  
  # any-ag fraction-blended rate
  vm_any <- ifelse(is.na(.den_frac_unc) | .den_frac_unc <= 0, NA_real_,
                   (vc * .cf_unc + vp * .pf_unc) / .den_frac_unc)
  
  # restoration-weighted blend
  vm_max <- ifelse(is.na(.den_rest_unc) | .den_rest_unc <= 0, NA_real_,
                   (vc * .cs_unc + vp * .ps_unc) / .den_rest_unc)
  
  tibble(
    crop_any = aw_mean_presence(vc,     .w_union_crop),
    past_any = aw_mean_presence(vp,     .w_union_past),
    soc_any  = aw_mean_presence(vm_any, .w_union_any_ag),
    crop_max = aw_mean_fraction(vc,     .w_union_cs),
    past_max = aw_mean_fraction(vp,     .w_union_ps),
    soc_max  = aw_mean_fraction(vm_max, .w_union_all)
  )
}

soc_draws <- purrr::map_dfr(sample_dirs, soc_draws_global)

soil_rates_any_ag_global <- tibble::tibble(
  group = "Global_All",
  pool = c("SOC after cropland","SOC after pasture","SOC"),
  mean = c(mean(soc_draws$crop_any, na.rm=TRUE),
           mean(soc_draws$past_any, na.rm=TRUE),
           mean(soc_draws$soc_any, na.rm=TRUE)),
  lower = c(quantile(soc_draws$crop_any, 0.025, na.rm=TRUE),
            quantile(soc_draws$past_any, 0.025, na.rm=TRUE),
            quantile(soc_draws$soc_any, 0.025, na.rm=TRUE)),
  upper = c(quantile(soc_draws$crop_any, 0.975, na.rm=TRUE),
            quantile(soc_draws$past_any, 0.975, na.rm=TRUE),
            quantile(soc_draws$soc_any, 0.975, na.rm=TRUE))
)

soil_rates_maximum_weighted_global <- tibble::tibble(
  group = "Global_All",
  pool = c("SOC after cropland","SOC after pasture","SOC (restoration areas)"),
  mean = c(mean(soc_draws$crop_max, na.rm=TRUE),
           mean(soc_draws$past_max, na.rm=TRUE),
           mean(soc_draws$soc_max, na.rm=TRUE)),
  lower = c(quantile(soc_draws$crop_max, 0.025, na.rm=TRUE),
            quantile(soc_draws$past_max, 0.025, na.rm=TRUE),
            quantile(soc_draws$soc_max, 0.025, na.rm=TRUE)),
  upper = c(quantile(soc_draws$crop_max, 0.975, na.rm=TRUE),
            quantile(soc_draws$past_max, 0.975, na.rm=TRUE),
            quantile(soc_draws$soc_max, 0.975, na.rm=TRUE))
)

# ---------------- AGC global ----------------
agc_rate <- (cci_50 + sbm_50) / 2 / years
# any‑ag presence (binary), not fractions:
m_any_global <- resample(any_ag_any, agc_rate, method = "near")
Aha_agc <- cellSize(agc_rate, unit = "ha")

W_any <- global(m_any_global * Aha_agc, fun="sum", na.rm=TRUE)[1,1]
M_any <- global(agc_rate * m_any_global * Aha_agc, fun="sum", na.rm=TRUE)[1,1]
V_any <- global(var_pix * (m_any_global * Aha_agc)^2, fun="sum", na.rm=TRUE)[1,1]
mu_any <- M_any / z2na(W_any)
se_any <- sqrt(V_any) / z2na(W_any)

# Restoration weight rasters (SpatRaster form, needed for AGC grid resampling)
cs_rt <- scenarios$crop_sparing$frac10  * union_mask
ps_rt <- scenarios$pasture_sparing$frac10 * union_mask

# restoration areas (cs+ps)
w_all_max_agc <- resample(cs_rt + ps_rt, agc_rate, method = "bilinear")
W_max <- global(w_all_max_agc * Aha_agc, fun="sum", na.rm=TRUE)[1,1]
M_max <- global(agc_rate * w_all_max_agc * Aha_agc, fun="sum", na.rm=TRUE)[1,1]
V_max <- global(var_pix * (w_all_max_agc * Aha_agc)^2, fun="sum", na.rm=TRUE)[1,1]
mu_max <- M_max / z2na(W_max)
se_max <- sqrt(V_max) / z2na(W_max)

agc_any_ag_global <- tibble::tibble(
  group = "Global_All", pool = "AGC",
  mean = mu_any, lower = mu_any - 1.96*se_any, upper = mu_any + 1.96*se_any
)

agc_max_weighted_global <- tibble::tibble(
  group = "Global_All", pool = "AGC (restoration areas)",
  mean = mu_max, lower = mu_max - 1.96*se_max, upper = mu_max + 1.96*se_max
)

# clean old Global_All rows and append new ones
soil_rates_any_ag <- soil_rates_any_ag %>% filter(group != "Global_All")
soil_rates_maximum_weighted <- soil_rates_maximum_weighted %>% filter(group != "Global_All")
agc_rates_any_ag <- agc_rates_any_ag %>% filter(group != "Global_All")
agc_rates_maximum_weighted <- agc_rates_maximum_weighted %>% filter(group != "Global_All")

soil_rates_any_ag <- bind_rows(soil_rates_any_ag, soil_rates_any_ag_global)
soil_rates_maximum_weighted <- bind_rows(soil_rates_maximum_weighted, soil_rates_maximum_weighted_global)
agc_rates_any_ag <- bind_rows(agc_rates_any_ag, agc_any_ag_global)
agc_rates_maximum_weighted <- bind_rows(agc_rates_maximum_weighted, agc_max_weighted_global)


# Build plotting data frame: branch on whether the finer biome grouping
# (Forest/Savanna/Shrubland+Grassland) is present, or only Forest/Open.
if(!"Temperate_Savanna" %in% groups){
  plot_df4 <- bind_rows(
    soil_rates_any_ag %>% filter(pool == "SOC") %>% mutate(pool = "SOC (biome-wide)"),
    soil_rates_maximum_weighted %>% filter(pool == "SOC (restoration areas)"),
    agc_rates_any_ag %>% mutate(pool = "AGC (biome-wide)"),
    agc_rates_maximum_weighted
  ) %>%
    separate(group, into = c("climate_zone","ecosystem_type"), sep = "_") %>%
    filter(climate_zone %in% c("Temperate","Tropical")) %>%
    mutate(
      climate_zone = factor(climate_zone, levels = c("Tropical","Temperate")),
      ecosystem_type = factor(ecosystem_type, levels = c("Forest","Open")),
      pool = factor(pool, levels = c(
        "AGC (biome-wide)",
        "AGC (restoration areas)",
        "SOC (biome-wide)",
        "SOC (restoration areas)"
      ))
    )
  
  # 1) turn ecosystem_type into a numeric x-axis
  plot_df5 <- plot_df4 %>%
    mutate(
      x_base = as.numeric(ecosystem_type), # 1 = Forest, 2 = Open
      
      # 2) assign an offset so AGC lives on the left cluster,
      # SOC on the right, with a visible gap in between
      x_off = case_when(
        pool == "AGC (biome-wide)" ~ -0.35,
        pool == "AGC (restoration areas)" ~ -0.15,
        pool == "SOC (biome-wide)" ~ +0.15,
        pool == "SOC (restoration areas)" ~ +0.35,
      ),
      
      x = x_base + x_off
    )
  
  # 3) plot “identity” style on those manual x's:
  p_mean_bars <- ggplot(plot_df5, aes(x = x, y = mean, fill = pool)) +
    geom_col(width = 0.19) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  width = 0.05) +
    
    # 4) put the real labels at 1 & 2
    scale_x_continuous(
      breaks = unique(plot_df5$x_base),
      labels = levels(plot_df5$ecosystem_type)
    ) +
    
    # 5) the facets + colours + theme as before
    facet_wrap(~ climate_zone, nrow = 1) +
    scale_fill_manual(
      values = c(
        "AGC (biome-wide)" = "#74a7a3",
        "AGC (restoration areas)" = "#74a7a380",
        "SOC (biome-wide)" = "#976020",
        "SOC (restoration areas)" = "#97602080"
      ),
      name = NULL
    ) +
    guides(
      fill = guide_legend(
        ncol = 2,
        byrow = F
      )
    ) +
    ylim(0,2) +
    labs(
      x = NULL,
      y = expression(
        atop("Accumulation rate",
             "(Mg C " * ha^-1 * yr^-1 * ")")
      )
    ) +
    theme_classic(base_size = 14) +
    theme(
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_text(margin = margin(r = 10)),
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key.width = unit(1.5, "lines")
    )
} else{
  plot_df4 <- bind_rows(
    soil_rates_any_ag %>% filter(pool == "SOC") %>% mutate(pool = "SOC (biome-wide)"),
    soil_rates_maximum_weighted %>% filter(pool == "SOC (restoration areas)"),
    agc_rates_any_ag %>% mutate(pool = "AGC (biome-wide)"),
    agc_rates_maximum_weighted
  ) %>%
    separate(group, into = c("climate_zone","ecosystem_type"), sep = "_") %>%
    # now allow Global too
    filter(climate_zone %in% c("Temperate","Tropical","Global")) %>% 
    # Change to "both" from "global"
    mutate(climate_zone = gsub("Global", "Both", climate_zone)) %>% 
    mutate(ecosystem_type = gsub("ShrublandGrassland", "Shrubland &\nGrassland", ecosystem_type)) %>% 
    mutate(ecosystem_type = gsub("All", "All\nBiomes", ecosystem_type)) %>% 
    mutate(
      climate_zone = factor(climate_zone, levels = c("Temperate","Tropical","Both")),
      ecosystem_type = factor(
        ecosystem_type,
        levels = c("Forest","Savanna","Shrubland &\nGrassland","All\nBiomes") # add "All" for Global facet
      ),
      pool = factor(pool, levels = c(
        "AGC (biome-wide)",
        "AGC (restoration areas)",
        "SOC (biome-wide)",
        "SOC (restoration areas)"
      ))
    ) %>%
    filter(!is.na(ecosystem_type)) %>% 
    droplevels()
  
  # --- compute a data-driven y-limit with a little headroom ---
  ymax <- max(plot_df4$upper, na.rm = TRUE)
  ymax <- ceiling((ymax * 1.25) * 10) / 10
  ymin <- min(plot_df4$lower, na.rm = TRUE)
  
  # --- numeric x + offsets (same spacing) ---
  plot_df5 <- plot_df4 %>%
    mutate(
      x_base = as.numeric(ecosystem_type),
      x_off = dplyr::case_when(
        pool == "AGC (biome-wide)" ~ -0.26,
        pool == "AGC (restoration areas)" ~ -0.14,
        pool == "SOC (biome-wide)" ~ +0.14,
        pool == "SOC (restoration areas)" ~ +0.26,
        TRUE ~ 0
      ),
      x = x_base + x_off
    )
  
  p_mean_bars <- ggplot(plot_df5, aes(x = x, y = mean, fill = pool)) +
    geom_col(width = 0.16) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.05) +
    scale_x_continuous(
      breaks = sort(unique(plot_df5$x_base)),
      labels = levels(plot_df5$ecosystem_type)
    ) +
    facet_grid(~ climate_zone, scales = "free", space = 'free') + # Tropical | Temperate | Global
    scale_fill_manual(
      values = c(
        "AGC (biome-wide)" = "#228B22",
        "AGC (restoration areas)" = "#228B2280",
        "SOC (biome-wide)" = "#8B5A2B",
        "SOC (restoration areas)" = "#8B5A2B80"
      ),
      name = NULL
    ) +
    guides(fill = guide_legend(
      ncol = 2,
      byrow = FALSE,
      keywidth = unit(1.2, "lines"), # width of legend keys
      keyheight = unit(0.9, "lines"), # height of legend keys
      title.theme = element_text(size = 10),
      label.theme = element_text(size = 10) # text size in legend
    )) +
    coord_cartesian(ylim = c(-0.05, ymax)) +
    labs(
      x = NULL,
      y = expression(atop("Accumulation rate", "(Mg C " * ha^-1 * yr^-1 * ")"))
    ) +
    theme_classic(base_size = 14) +
    theme(
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_text(margin = margin(r = 10)),
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = NA, color = NA),
      legend.key.width = unit(1.5, "lines"),
      legend.key.size = unit(0.9, "lines"), # overall key size
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10)
    )
  
  
  p_mean_bars
}





# Now we'll also include a very similar thing as a table with some bells and 
# whistles - including belowground biomass.

# ───────────────────────────────────────────────────────────────────────────────
# PREP: Align RMF to the two working grids once (rate grid for SOC, AGC grid)
# Assumes we have already created & gap-filled: rmf_r (fraction, 0–1)
# ───────────────────────────────────────────────────────────────────────────────
rmf_sd <- rmf_boot_r # <- replace if this is SE, see note above
var_rmf <- rmf_sd^2 # variance of RMF

# Align to AGC grid
var_rmf_agc <- resample(var_rmf, agc_rate, method = "bilinear")

# Guardrails for r
rmf_agc <- clamp(resample(rmf_r, agc_rate, method = "bilinear"), lower=0, upper=0.99)

# Mean: AGC+BGC rate
agcbgc_rate <- agc_rate / (1 - rmf_agc)

# Variance components
varA_term <- var_pix / (1 - rmf_agc)^2
varR_term <- (agc_rate^2) * var_rmf_agc / (1 - rmf_agc)^4
var_pix_agcbgc <- varA_term + varR_term

# ───────────────────────────────────────────────────────────────────────────────
# GROUP-LEVEL: AGC+BGC (presence-masked and restoration-weighted)
# ───────────────────────────────────────────────────────────────────────────────

# Any-ag presence (binary) on AGC grid (reuse the any_ag_any logic)
m_any_agc <- resample(any_ag_any, agcbgc_rate, method = "near")
Aha_agc <- cellSize(agcbgc_rate, unit = "ha")

agcbgc_rates_any_ag <- imap_dfr(group_masks, function(mask_g, gr) {
  m_group <- resample(mask_g, agcbgc_rate, method = "near")
  W <- global(m_group * m_any_agc * Aha_agc, fun = "sum", na.rm = TRUE)[1,1]
  if (is.na(W) || W == 0) {
    return(tibble(group = gr, pool = "AGC+BGC (biome-wide)",
                  mean = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  M <- global(agcbgc_rate * m_group * m_any_agc * Aha_agc, fun = "sum", na.rm = TRUE)[1,1]
  Vw <- global(var_pix_agcbgc * (m_group * m_any_agc * Aha_agc)^2, fun = "sum", na.rm = TRUE)[1,1]
  mu <- M / W
  se <- sqrt(Vw) / W
  tibble(group = gr, pool = "AGC+BGC (biome-wide)",
         mean = mu, lower = mu - 1.96*se, upper = mu + 1.96*se)
})

# Restoration-weighted version (cs+ps on AGC grid)
w_all_agc <- resample(cs_rt + ps_rt, agcbgc_rate, method = "bilinear")
agcbgc_rates_maximum_weighted <- imap_dfr(group_masks, function(mask_g, gr) {
  m_group <- resample(mask_g, agcbgc_rate, method = "near")
  W <- global(w_all_agc * m_group * Aha_agc, fun = "sum", na.rm = TRUE)[1,1]
  if (is.na(W) || W == 0) {
    return(tibble(group = gr, pool = "AGC+BGC (restoration areas)",
                  mean = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  M <- global(agcbgc_rate * w_all_agc * m_group * Aha_agc, fun = "sum", na.rm = TRUE)[1,1]
  Vw <- global(var_pix_agcbgc * (w_all_agc * m_group * Aha_agc)^2, fun = "sum", na.rm = TRUE)[1,1]
  mu <- M / W
  se <- sqrt(Vw) / W
  tibble(group = gr, pool = "AGC+BGC (restoration areas)",
         mean = mu, lower = mu - 1.96*se, upper = mu + 1.96*se)
})

# ───────────────────────────────────────────────────────────────────────────────
# GLOBAL: AGC+BGC (presence-masked + restoration-weighted)
# ───────────────────────────────────────────────────────────────────────────────
m_any_global_agcbgc <- resample(any_ag_any, agcbgc_rate, method = "near")
W_any <- global(m_any_global_agcbgc * Aha_agc, fun = "sum", na.rm = TRUE)[1,1]
M_any <- global(agcbgc_rate * m_any_global_agcbgc * Aha_agc, fun = "sum", na.rm = TRUE)[1,1]
V_any <- global(var_pix_agcbgc * (m_any_global_agcbgc * Aha_agc)^2, fun = "sum", na.rm = TRUE)[1,1]
mu_any <- M_any / W_any; se_any <- sqrt(V_any) / W_any

W_max <- global(w_all_agc * Aha_agc, fun = "sum", na.rm = TRUE)[1,1]
M_max <- global(agcbgc_rate * w_all_agc * Aha_agc, fun = "sum", na.rm = TRUE)[1,1]
V_max <- global(var_pix_agcbgc * (w_all_agc * Aha_agc)^2, fun = "sum", na.rm = TRUE)[1,1]
mu_max <- M_max / W_max; se_max <- sqrt(V_max) / W_max

agcbgc_any_ag_global <- tibble::tibble(
  group = "Global_All", pool = "AGC+BGC (biome-wide)",
  mean = mu_any, lower = mu_any - 1.96*se_any, upper = mu_any + 1.96*se_any
)
agcbgc_max_weighted_global <- tibble::tibble(
  group = "Global_All", pool = "AGC+BGC (restoration areas)",
  mean = mu_max, lower = mu_max - 1.96*se_max, upper = mu_max + 1.96*se_max
)

# Bind into the existing data frames (after removing any old Global_All rows)
agcbgc_rates_any_ag <- agcbgc_rates_any_ag %>% filter(group != "Global_All")
agcbgc_rates_maximum_weighted <- agcbgc_rates_maximum_weighted %>% filter(group != "Global_All")

agcbgc_rates_any_ag <- bind_rows(agcbgc_rates_any_ag, agcbgc_any_ag_global)
agcbgc_rates_maximum_weighted <- bind_rows(agcbgc_rates_maximum_weighted, agcbgc_max_weighted_global)




# ───────────────────────────────────────────────────────────────────────────────
# ROBINSON (forest-only): propagate uncertainty and aggregate with weights
# ───────────────────────────────────────────────────────────────────────────────

# 1) Annualized AGC rate and its variance from Robinson 50-yr rasters
rob_rate <- rob_50 / 50
rob_var <- (rob_50_se / 50)^2 # variance for AGC (rate)

# 2) RMF (mean and variance) aligned to Robinson grid
rmf_rob <- clamp(resample(rmf_r, rob_rate, method="bilinear"), 0, 0.99)
var_rmf_rob <- resample(var_rmf, rob_rate, method="bilinear")

# 3) AGC+BGC via delta method (assume Cov(AGC, RMF) = 0)
# f = A / (1-r), Var[f] ≈ VarA/(1-r)^2 + A^2*VarR/(1-r)^4
rob_agcbgc <- rob_rate / (1 - rmf_rob)
rob_agcbgc_var <- rob_var / (1 - rmf_rob)^2 + (rob_rate^2) * var_rmf_rob / (1 - rmf_rob)^4

# 4) Weights on the Robinson grid
any_ag_rob <- resample(any_ag_any, rob_rate, method = "near")
cs_rob <- resample(cs_rt, rob_rate, method = "bilinear")
ps_rob <- resample(ps_rt, rob_rate, method = "bilinear")
w_all_rob <- cs_rob + ps_rob
Aha_rob <- cellSize(rob_rate, unit = "ha")

# 5) Aggregators (biome-wide uses presence mask; restoration uses cs+ps)
agg_presence <- function(mean_r, var_r, mask_g) {
  mg <- resample(mask_g, mean_r, method = "near")
  W <- global(mg * any_ag_rob * Aha_rob, fun="sum", na.rm=TRUE)[1,1]
  if (is.na(W) || W == 0) return(list(mu=NA_real_, se=NA_real_))
  M <- global(mean_r * mg * any_ag_rob * Aha_rob, fun="sum", na.rm=TRUE)[1,1]
  Vw <- global(var_r * (mg * any_ag_rob * Aha_rob)^2, fun="sum", na.rm=TRUE)[1,1]
  mu <- M / W
  se <- sqrt(Vw) / W
  list(mu=mu, se=se)
}
agg_restoration <- function(mean_r, var_r, mask_g) {
  mg <- resample(mask_g, mean_r, method = "near")
  W <- global(w_all_rob * mg * Aha_rob, fun="sum", na.rm=TRUE)[1,1]
  if (is.na(W) || W == 0) return(list(mu=NA_real_, se=NA_real_))
  M <- global(mean_r * w_all_rob * mg * Aha_rob, fun="sum", na.rm=TRUE)[1,1]
  Vw <- global(var_r * (w_all_rob * mg * Aha_rob)^2, fun="sum", na.rm=TRUE)[1,1]
  mu <- M / W
  se <- sqrt(Vw) / W
  list(mu=mu, se=se)
}

forest_groups <- c("Temperate_Forest","Tropical_Forest")

robinson_forest_stats <- purrr::map_dfr(forest_groups, function(gr){
  gmask <- group_masks[[gr]]
  # AGC
  pres_agc <- agg_presence(rob_rate, rob_var, gmask)
  rest_agc <- agg_restoration(rob_rate, rob_var, gmask)
  # AGC+BGC
  pres_ab <- agg_presence(rob_agcbgc, rob_agcbgc_var, gmask)
  rest_ab <- agg_restoration(rob_agcbgc,rob_agcbgc_var, gmask)
  tibble::tibble(
    group = gr,
    AGC_Robinson_biome_wide_mean = pres_agc$mu,
    AGC_Robinson_biome_wide_se = pres_agc$se,
    AGC_Robinson_rest_mean = rest_agc$mu,
    AGC_Robinson_rest_se = rest_agc$se,
    AGCBGC_Robinson_biome_wide_mean = pres_ab$mu,
    AGCBGC_Robinson_biome_wide_se = pres_ab$se,
    AGCBGC_Robinson_rest_mean = rest_ab$mu,
    AGCBGC_Robinson_rest_se = rest_ab$se
  )
})




# ───────────────────────────────────────────────────────────────────────────────
# WIDE SUMMARY TABLE with CI strings, one row per group
# ───────────────────────────────────────────────────────────────────────────────

fmt_ci <- function(mu, lo, up) {
  ifelse(is.na(mu), NA_character_, sprintf("%.2f [%.2f, %.2f]", mu, lo, up))
}

# 1) Reduce the existing long frames to wide-by-pool CI strings
pick_ci <- function(df, pool_label, out_col) {
  df %>%
    dplyr::filter(pool == pool_label) %>%
    dplyr::transmute(group, !!out_col := fmt_ci(mean, lower, upper))
}

soc_bw <- pick_ci(soil_rates_any_ag, "SOC", "SOC_biome_wide")
soc_rest <- pick_ci(soil_rates_maximum_weighted, "SOC (restoration areas)", "SOC_restoration_areas")

agc_bw <- pick_ci(agc_rates_any_ag, "AGC", "AGC_biome_wide") %>%
  { if (!"AGC_biome_wide" %in% names(.)) dplyr::rename(., AGC_biome_wide = pool) else . }
agc_rest <- pick_ci(agc_rates_maximum_weighted, "AGC (restoration areas)", "AGC_restoration_areas")

ab_bw <- pick_ci(agcbgc_rates_any_ag, "AGC+BGC (biome-wide)", "AGCBGC_biome_wide")
ab_rest <- pick_ci(agcbgc_rates_maximum_weighted,"AGC+BGC (restoration areas)","AGCBGC_restoration_areas")

# 2) Robinson forest-only → CI strings; NA for non-forest after join
rob_wide <- robinson_forest_stats %>%
  dplyr::mutate(
    AGC_Robinson_biome_wide = fmt_ci(AGC_Robinson_biome_wide_mean,
                                     AGC_Robinson_biome_wide_mean - 1.96*AGC_Robinson_biome_wide_se,
                                     AGC_Robinson_biome_wide_mean + 1.96*AGC_Robinson_biome_wide_se),
    AGC_Robinson_restoration_weighted = fmt_ci(AGC_Robinson_rest_mean,
                                               AGC_Robinson_rest_mean - 1.96*AGC_Robinson_rest_se,
                                               AGC_Robinson_rest_mean + 1.96*AGC_Robinson_rest_se),
    AGCBGC_Robinson_biome_wide = fmt_ci(AGCBGC_Robinson_biome_wide_mean,
                                        AGCBGC_Robinson_biome_wide_mean - 1.96*AGCBGC_Robinson_biome_wide_se,
                                        AGCBGC_Robinson_biome_wide_mean + 1.96*AGCBGC_Robinson_biome_wide_se),
    AGCBGC_Robinson_restoration_weighted = fmt_ci(AGCBGC_Robinson_rest_mean,
                                                  AGCBGC_Robinson_rest_mean - 1.96*AGCBGC_Robinson_rest_se,
                                                  AGCBGC_Robinson_rest_mean + 1.96*AGCBGC_Robinson_rest_se)
  ) %>%
  dplyr::select(
    group,
    AGC_Robinson_biome_wide,
    AGC_Robinson_restoration_weighted,
    AGCBGC_Robinson_biome_wide,
    AGCBGC_Robinson_restoration_weighted
  )

# 3) Assemble: start from all groups present in the SOC frame (one row per group)
all_groups <- soil_rates_any_ag %>% dplyr::distinct(group)

summary_tbl <- all_groups %>%
  dplyr::left_join(soc_bw, by="group") %>%
  dplyr::left_join(soc_rest, by="group") %>%
  dplyr::left_join(agc_bw, by="group") %>%
  dplyr::left_join(agc_rest, by="group") %>%
  dplyr::left_join(rob_wide, by="group") %>% # fills only forest rows; others NA
  dplyr::left_join(ab_bw, by="group") %>%
  dplyr::left_join(ab_rest, by="group") %>%
  # column order matches the summary table layout
  dplyr::select(
    group,
    `SOC (biome_wide)` = SOC_biome_wide,
    `SOC (restoration areas)` = SOC_restoration_areas,
    `AGB (biome wide)` = AGC_biome_wide,
    `AGB (restoration areas)` = AGC_restoration_areas,
    `AGB_Robinson (biome wide)` = AGC_Robinson_biome_wide,
    `AGB_Robinson (restoration_weighted)` = AGC_Robinson_restoration_weighted,
    `AGB+BGB (biome wide)` = AGCBGC_biome_wide,
    `AGB+BGB (restoration areas)` = AGCBGC_restoration_areas,
    `AGB+BGB_Robinson (biome wide)` = AGCBGC_Robinson_biome_wide,
    `AGB+BGB_Robinson (restoration_weighted)` = AGCBGC_Robinson_restoration_weighted
  ) %>%
  dplyr::arrange(group)

# View
summary_tbl

write.csv(summary_tbl, file = file.path(figures_dir, "soc_agb_bgb_table_50_year_mean.csv"))



# Assemble plotting data frame including AGC+BGC pools and the Global_All group
plot_df4_with_bgc <- bind_rows(
  soil_rates_any_ag %>% filter(pool == "SOC") %>% mutate(pool = "SOC (biome-wide)"),
  soil_rates_maximum_weighted %>% filter(pool == "SOC (restoration areas)"),
  agcbgc_rates_any_ag %>% mutate(pool = "AGC+BGC (biome-wide)"),
  agcbgc_rates_maximum_weighted %>% mutate(pool = "AGC+BGC (restoration areas)")
) %>%
  separate(group, into = c("climate_zone","ecosystem_type"), sep = "_", fill = "right") %>%
  filter(climate_zone %in% c("Temperate","Tropical","Global")) %>% 
  mutate(
    climate_zone = gsub("Global", "Both", climate_zone),
    ecosystem_type = gsub("ShrublandGrassland", "Shrubland &\nGrassland", ecosystem_type),
    ecosystem_type = gsub("All", "All\nBiomes", ecosystem_type),
    climate_zone = factor(climate_zone, levels = c("Temperate","Tropical","Both")),
    ecosystem_type = factor(
      ecosystem_type,
      levels = c("Forest","Savanna","Shrubland &\nGrassland","All\nBiomes")
    ),
    pool = factor(pool, levels = c(
      "AGC+BGC (biome-wide)",
      "AGC+BGC (restoration areas)",
      "SOC (biome-wide)",
      "SOC (restoration areas)"
    ))
  ) %>%
  filter(!is.na(ecosystem_type)) %>%
  droplevels()

# --- compute a data-driven y-limit with a little headroom ---
ymax <- max(plot_df4_with_bgc$upper, na.rm = TRUE)
ymax <- ceiling((ymax * 1.25) * 10) / 10
ymin <- min(plot_df4_with_bgc$lower, na.rm = TRUE)

# --- numeric x + offsets ---
plot_df5_with_bgc <- plot_df4_with_bgc %>%
  mutate(
    x_base = as.numeric(ecosystem_type),
    x_off = dplyr::case_when(
      pool == "AGC+BGC (biome-wide)" ~ -0.26,
      pool == "AGC+BGC (restoration areas)" ~ -0.14,
      pool == "SOC (biome-wide)" ~ +0.14,
      pool == "SOC (restoration areas)" ~ +0.26,
      TRUE ~ 0
    ),
    x = x_base + x_off
  )

p_mean_bars_with_bgc <- ggplot(plot_df5_with_bgc, aes(x = x, y = mean, fill = pool)) +
  geom_col(width = 0.16) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.05) +
  scale_x_continuous(
    breaks = sort(unique(plot_df5$x_base)),
    labels = levels(plot_df5$ecosystem_type)
  ) +
  facet_grid(~ climate_zone, scales = "free", space = 'free') +
  scale_fill_manual(
    values = c(
      "AGC+BGC (biome-wide)" = "#1f7a1f", # deeper green than AGC
      "AGC+BGC (restoration areas)" = "#1f7a1f80",
      "SOC (biome-wide)" = "#8B5A2B",
      "SOC (restoration areas)" = "#8B5A2B80"
    ),
    name = NULL
  ) +
  guides(fill = guide_legend(
    ncol = 2, byrow = FALSE,
    keywidth = unit(1.2, "lines"),
    keyheight = unit(0.9, "lines"),
    title.theme = element_text(size = 10),
    label.theme = element_text(size = 10)
  )) +
  coord_cartesian(ylim = c(-0.05, ymax)) +
  labs(
    x = NULL,
    y = expression(atop("Accumulation rate", "(Mg C " * ha^-1 * yr^-1 * ")"))
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = NA, color = NA),
    legend.key.width = unit(1.5, "lines"),
    legend.key.size = unit(0.9, "lines"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  )

p_mean_bars_with_bgc

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png")
    png(file.path(figures_dir, "p_mean_bars_with_bgc.png"), width = 8.5, height = 3, unit = "in", res = 440)
  else
    pdf(file.path(figures_dir, "p_mean_bars_with_bgc.pdf"), width = 8.5, height = 3)
  print(p_mean_bars_with_bgc)
  dev.off()
}



# And lastly one where it's broken up between agc and bgc

# -- AGC base
agc_rows <- plot_df5 %>%
  filter(grepl("^AGC", pool)) %>%
  mutate(
    mode = ifelse(grepl("restoration", pool), "restoration areas", "biome-wide"),
    fill_group = paste("AGC", mode)
  ) %>%
  select(climate_zone, ecosystem_type, x_base, x,
         agc_mean = mean, agc_lower = lower, agc_upper = upper, fill_group, mode)

# -- AGC+BGC totals
tot_rows <- plot_df5_with_bgc %>%
  filter(grepl("^AGC\\+BGC", pool)) %>%
  mutate(
    mode = ifelse(grepl("restoration", pool), "restoration areas", "biome-wide"),
    fill_group = paste("BGC", mode)
  ) %>%
  mutate(pool = gsub("^AGC\\+BGC", "AGC", pool)) %>%
  select(climate_zone, ecosystem_type, x_base, x,
         tot_mean = mean, tot_lower = lower, tot_upper = upper, fill_group, mode)

# -- Join for stacking
agc_join <- left_join(agc_rows, tot_rows,
                      by = c("climate_zone", "ecosystem_type", "x_base", "x", "mode"))

stack_df <- bind_rows(
  agc_join %>%
    mutate(component = "AGC", ymin = 0, ymax = agc_mean, fill_group = paste("AGC", mode)),
  agc_join %>%
    mutate(component = "BGC",
           ymin = agc_mean,
           ymax = pmax(tot_mean, agc_mean), # ensure positive bar height
           fill_group = paste("BGC", mode))
)

# -- Total error bars (for AGC+BGC)
tot_err_df <- agc_join %>%
  mutate(fill_group = paste("BGC", mode)) %>%
  distinct(climate_zone, ecosystem_type, x, tot_lower, tot_upper, fill_group, mode)

# -- AGC error bars
agc_err_df <- agc_rows

# -- SOC data
soc_rows <- plot_df5 %>%
  filter(grepl("^SOC", pool)) %>%
  mutate(
    mode = ifelse(grepl("restoration", pool), "restoration areas", "biome-wide"),
    fill_group = paste("SOC", mode)
  )

# -- Legend colors and transparency
fill_colors <- c(
  "AGC biome-wide" = "#228B2280",
  "AGC restoration areas" = "#228B22",
  "BGC biome-wide" = "#B6CBCD80",
  "BGC restoration areas" = "#B6CBCD",
  "SOC biome-wide" = "#8B5A2B80",
  "SOC restoration areas" = "#8B5A2B"
)

alpha_vals <- c(
  "biome-wide" = 0.5,
  "restoration areas" = 1
)

# -- Y axis limit
ymax2 <- max(c(plot_df5$upper, plot_df5_with_bgc$upper), na.rm = TRUE)
ymax2 <- ceiling((ymax2 * 1.3) * 10) / 10

# -- Final plot
p_stacked_clean <- ggplot() +
  
  # AGC + BGC stacked bars
  geom_rect(
    data = stack_df,
    aes(xmin = x - 0.08, xmax = x + 0.08, ymin = ymin, ymax = ymax,
        fill = fill_group, alpha = mode)
  ) +
  
  # AGC error bars
  geom_errorbar(
    data = agc_err_df,
    aes(x = x, ymin = agc_lower, ymax = agc_upper),
    width = 0.05, color = "grey20", linewidth = 0.4
  ) +
  
  # Total AGC+BGC error bars
  geom_errorbar(
    data = tot_err_df,
    aes(x = x, ymin = tot_lower, ymax = tot_upper),
    width = 0.05, color = "grey20", linewidth = 0.4
  ) +
  
  # SOC bars
  geom_col(
    data = soc_rows,
    aes(x = x, y = mean, fill = fill_group, alpha = mode),
    width = 0.16
  ) +
  
  # SOC error bars
  geom_errorbar(
    data = soc_rows,
    aes(x = x, ymin = lower, ymax = upper),
    width = 0.05, color = "grey20", linewidth = 0.4
  ) +
  
  scale_x_continuous(
    breaks = sort(unique(plot_df5$x_base)),
    labels = levels(plot_df5$ecosystem_type)
  ) +
  facet_grid(~ climate_zone, scales = "free", space = 'free') +
  scale_fill_manual(values = fill_colors, name = NULL) +
  scale_alpha_manual(values = alpha_vals, guide = "none") +
  coord_cartesian(ylim = c(-0.05, ymax2)) +
  labs(
    x = NULL,
    y = expression(atop("Accumulation rate", "(Mg C " * ha^-1 * yr^-1 * ")"))
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = NA, color = NA),
    legend.key.width = unit(1.5, "lines"),
    legend.key.size = unit(0.9, "lines")
  ) +
  guides(fill = guide_legend(nrow = 2, ncol = 3, byrow = F))

p_stacked_clean

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png")
    png(file.path(figures_dir, "p_mean_bars_agc_vs_bgc.png"), width = 8.5, height = 3, unit = "in", res = 440)
  else
    pdf(file.path(figures_dir, "p_mean_bars_agc_vs_bgc.pdf"), width = 8.5, height = 3)
  print(p_stacked_clean)
  dev.off()
}







# ---------------------------
# Stacked bar variant: AGC (bottom) + BGC (top) + SOC, with additional
# per-land-use SOC bars (after cropland / after pasture) at the right of each cluster
# ---------------------------

# 1) Build AGC rows (from original plot_df5) and AGC+BGC totals (from plot_df5_with_bgc)
agc_rows <- plot_df5 %>%
  dplyr::filter(grepl("^AGC", pool)) %>%
  dplyr::mutate(
    mode = ifelse(grepl("restoration", pool), "restoration areas", "biome-wide"),
    fill_group = paste("AGC", mode)
  ) %>%
  dplyr::select(climate_zone, ecosystem_type, x_base, x,
                agc_mean = mean, agc_lower = lower, agc_upper = upper, fill_group, mode)

tot_rows <- plot_df5_with_bgc %>%
  dplyr::filter(grepl("^AGC\\+BGC", pool)) %>%
  dplyr::mutate(
    mode = ifelse(grepl("restoration", pool), "restoration areas", "biome-wide"),
    fill_group = paste("BGC", mode),
    pool = gsub("^AGC\\+BGC", "AGC", pool) # align pool naming to join with AGC rows
  ) %>%
  dplyr::select(climate_zone, ecosystem_type, x_base, x,
                tot_mean = mean, tot_lower = lower, tot_upper = upper, fill_group, mode)

# 2) Join at identical x positions (so BGC stacks on top of the AGC bar)
agc_join <- dplyr::left_join(
  agc_rows, tot_rows,
  by = c("climate_zone","ecosystem_type","x_base","x","mode")
)

# 3) Rects for stacked bars: AGC from 0→agc_mean; BGC from agc_mean→tot_mean
stack_df <- dplyr::bind_rows(
  agc_join %>% dplyr::transmute(
    climate_zone, ecosystem_type, x_base, x, mode,
    fill_group = paste("AGC", mode),
    ymin = 0, ymax = agc_mean
  ),
  agc_join %>% dplyr::transmute(
    climate_zone, ecosystem_type, x_base, x, mode,
    fill_group = paste("BGC", mode),
    ymin = agc_mean, ymax = pmax(tot_mean, agc_mean) # guard against negative slice
  )
)

# 4) Errorbars data frames (black)
agc_err_df <- agc_rows %>% dplyr::select(climate_zone, ecosystem_type, x, agc_lower, agc_upper)
tot_err_df <- agc_join %>% dplyr::select(climate_zone, ecosystem_type, x, tot_lower, tot_upper) %>%
  dplyr::distinct()

# 5) SOC bars at existing offsets + the two extra SOC bars at the far right of "All\nBiomes"
soc_rows <- plot_df5_with_bgc %>%
  dplyr::filter(grepl("^SOC", pool)) %>%
  dplyr::mutate(
    mode = ifelse(grepl("restoration", pool), "restoration areas", "biome-wide"),
    fill_group = paste("SOC", mode)
  )

# find x_base for "All\nBiomes" (last cluster) and add two new SOC bars there
all_biomes_x <- which(levels(plot_df5_with_bgc$ecosystem_type) == "All\nBiomes")
stopifnot(length(all_biomes_x) == 1)

extra_soc <- soil_rates_maximum_weighted_global %>%
  dplyr::filter(pool %in% c("SOC after cropland","SOC after pasture")) %>%
  dplyr::transmute(
    climate_zone = factor("Both", levels = c("Temperate","Tropical","Both")),
    ecosystem_type = factor("All\nBiomes", levels = levels(plot_df5_with_bgc$ecosystem_type)),
    pool, mean, lower, upper,
    x_base = all_biomes_x,
    x = all_biomes_x + ifelse(pool == "SOC after cropland", 0.54, 0.82),
    fill_group = pool,
    mode = "restoration areas"
  )

soc_rows_ext <- dplyr::bind_rows(soc_rows, extra_soc)

# 6) Legend: 8 entries (AGC/BGC by mode + SOC by mode + two extra SOC bars)
fill_colors <- c(
  "AGC biome-wide" = "#228B2280",
  "AGC restoration areas" = "#228B22",
  "BGC biome-wide" = "#B6CBCD80", # lighter blue-green for roots
  "BGC restoration areas" = "#B6CBCD",
  "SOC biome-wide" = "#8B5A2B80",
  "SOC restoration areas" = "#8B5A2B",
  "SOC after cropland" = "#d95f02",
  "SOC after pasture" = "#7570b3"
)

soc_rows_ext$fill_group <- factor(
  soc_rows_ext$fill_group,
  levels = c(
    "AGC biome-wide","AGC restoration areas",
    "BGC biome-wide","BGC restoration areas",
    "SOC biome-wide","SOC restoration areas",
    "SOC after cropland","SOC after pasture"
  )
)

stack_df$fill_group <- factor(
  stack_df$fill_group,
  levels = levels(soc_rows_ext$fill_group)
)

alpha_vals <- c("biome-wide" = 0.5, "restoration areas" = 1)

# 7) Y-axis limit with headroom
ymax2 <- max(
  c(plot_df5_with_bgc$upper, agc_rows$agc_upper, tot_rows$tot_upper, extra_soc$upper),
  na.rm = TRUE
)
ymax2 <- ceiling((ymax2 * 1.3) * 10) / 10

#
# midpoint of all x's drawn for "All\nBiomes" (stacked AGC/BGC + SOC, Both facet)
all_x_vals <- c(
  soc_rows_ext %>%
    dplyr::filter(climate_zone == "Both", ecosystem_type == "All\nBiomes") %>%
    dplyr::pull(x),
  stack_df %>%
    dplyr::filter(climate_zone == "Both", ecosystem_type == "All\nBiomes") %>%
    dplyr::pull(x)
)
all_biomes_mid <- mean(range(all_x_vals, na.rm = TRUE))

# original breaks, but replace the "All\nBiomes" entry with the midpoint
breaks_vec <- sort(unique(plot_df5_with_bgc$x_base))
breaks_vec[all_biomes_x] <- all_biomes_mid

# 8) Plot — keeps original spacing & ticks
p_stacked_clean_with_past_crop_bars <- ggplot() +
  # Stacked AGC+BGC at the AGC x (AGC bottom, BGC top)
  geom_rect(
    data = stack_df,
    aes(xmin = x - 0.08, xmax = x + 0.08, ymin = ymin, ymax = ymax,
        fill = fill_group, alpha = mode)
  ) +
  # AGC (bottom slice) error bars
  geom_errorbar(
    data = agc_err_df,
    aes(x = x, ymin = agc_lower, ymax = agc_upper),
    width = 0.05, color = "black", linewidth = 0.4
  ) +
  # Total (AGC+BGC) error bars
  geom_errorbar(
    data = tot_err_df,
    aes(x = x, ymin = tot_lower, ymax = tot_upper),
    width = 0.05, color = "black", linewidth = 0.4
  ) +
  # SOC bars (existing four + the two extra at right)
  geom_col(
    data = soc_rows_ext,
    aes(x = x, y = mean, fill = fill_group, alpha = mode),
    width = 0.16
  ) +
  # SOC error bars (black)
  geom_errorbar(
    data = soc_rows_ext,
    aes(x = x, ymin = lower, ymax = upper),
    width = 0.05, color = "black", linewidth = 0.4
  ) +
  # x scale: reuse the ORIGINAL tick positions/labels so “All\nBiomes” stays singular
  scale_x_continuous(
    breaks = breaks_vec,
    labels = levels(plot_df5_with_bgc$ecosystem_type)
  ) +
  facet_grid(~ climate_zone, scales = "free", space = 'free') +
  scale_fill_manual(values = fill_colors, name = NULL) +
  scale_alpha_manual(values = alpha_vals, guide = "none") +
  coord_cartesian(ylim = c(-0.05, ymax2)) +
  labs(
    x = NULL,
    y = expression(atop("Accumulation rate", "(Mg C " * ha^-1 * yr^-1 * ")"))
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = NA, color = NA),
    legend.text = element_text(size = 10),
    legend.key.width = unit(1, "lines"),
    legend.key.size = unit(0.9, "lines")
  ) +
  guides(fill = guide_legend(nrow = 2, ncol = 4, byrow = F))

p_stacked_clean_with_past_crop_bars

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png")
    png(file.path(figures_dir, "p_mean_bars_agc_vs_bgc_with_past_crop_bars.png"), width = 8.5, height = 3, unit = "in", res = 440)
  else
    pdf(file.path(figures_dir, "p_mean_bars_agc_vs_bgc_with_past_crop_bars.pdf"), width = 8.5, height = 3)
  print(p_stacked_clean_with_past_crop_bars)
  dev.off()
}
















# Restoration-areas-only variant: biome-wide bars removed; AGC/BGC stacked, SOC alongside

fill_colors <- c(
  "AGC" = "#228B22",
  "BGC" = "#B6CBCD",
  "SOC" = "#8B5A2B",
  "SOC after cropland" = "#d95f02",
  "SOC after pasture" = "#7570b3"
)

# relabel fill_group to drop "restoration areas"
soc_rows_ext$fill_group <- forcats::fct_recode(
  soc_rows_ext$fill_group,
  "AGC" = "AGC restoration areas",
  "BGC" = "BGC restoration areas",
  "SOC" = "SOC restoration areas"
)

stack_df$fill_group <- forcats::fct_recode(
  stack_df$fill_group,
  "AGC" = "AGC restoration areas",
  "BGC" = "BGC restoration areas"
)

# relevel so the legend order is nice
soc_rows_ext$fill_group <- factor(
  soc_rows_ext$fill_group,
  levels = c("AGC","BGC","SOC","SOC after cropland","SOC after pasture")
)
stack_df$fill_group <- factor(
  stack_df$fill_group,
  levels = levels(soc_rows_ext$fill_group)
)

# since we dropped biome-wide, alpha isn't needed anymore
alpha_vals <- c("restoration areas" = 1)


ymax2 <- 1.9


# ---------------------------
# Re-center bar positions
# ---------------------------

# Lane centers relative to the tick (cluster center)
left_lane <- -0.175 # AGC+BGC
right_lane <- 0.175 # SOC family (SOC, SOC after cropland/pasture)
bar_w <- 0.18

# 1) Center AGC/BGC stack at the left lane of each cluster
stack_df_centered <- stack_df %>%
  dplyr::filter(mode != "biome-wide") %>%
  dplyr::mutate(
    x_center = x_base + left_lane
  ) %>% 
  mutate( xmin = x_center - bar_w/2,
          xmax = x_center + bar_w/2)

# Match error bars for AGC and total (AGC+BGC)
agc_err_centered <- agc_rows %>%
  dplyr::filter(mode != "biome-wide") %>%
  dplyr::mutate(x_center = x_base + left_lane)

tot_err_centered <- agc_join %>%
  dplyr::filter(mode != "biome-wide") %>%
  dplyr::mutate(x_center = x_base + left_lane)

# 2) Evenly space SOC bars at the right lane of each cluster
# - For clusters with 1 SOC bar: centered at right_lane
# - With 2 bars: +/- 0.06 around right_lane
# - With 3 bars: -0.08, 0, +0.08 around right_lane (order: SOC, SOC after cropland, SOC after pasture)

soc_order <- c("SOC", "SOC after cropland", "SOC after pasture")

soc_rows_centered <- soc_rows_ext %>%
  dplyr::filter(mode != "biome-wide") %>%
  dplyr::mutate(fill_group = factor(fill_group, levels = soc_order)) %>%
  dplyr::group_by(climate_zone, ecosystem_type) %>%
  dplyr::arrange(fill_group, .by_group = TRUE) %>%
  dplyr::mutate(
    k = dplyr::n(),
    r = dplyr::row_number(),
    base_center = x_base + right_lane,
    offset = dplyr::case_when(
      k == 1 ~ 0,
      k == 2 ~ ifelse(r == 1, -0.06, +0.06),
      k >= 3 ~ c(0, 0.25, 0.5)[pmin(r, 3)]
    ),
    x_soc = base_center + offset
  ) %>%
  dplyr::ungroup()

# ---------------------------
# Plot with centered positions
# ---------------------------
p_stacked_clean_without_biome_wide_bars <- ggplot() +
  # Stacked AGC+BGC (left lane)
  geom_rect(
    data = stack_df_centered,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = fill_group),
    alpha = 1
  ) +
  # AGC (bottom slice) error bars — left lane
  geom_errorbar(
    data = agc_err_centered,
    aes(x = x_center, ymin = agc_lower, ymax = agc_upper),
    width = 0.05, color = "black", linewidth = 0.4
  ) +
  # Total (AGC+BGC) error bars — left lane
  geom_errorbar(
    data = tot_err_centered,
    aes(x = x_center, ymin = tot_lower, ymax = tot_upper),
    width = 0.05, color = "black", linewidth = 0.4
  ) +
  # SOC bars (right lane, evenly spaced)
  geom_col(
    data = soc_rows_centered,
    aes(x = x_soc, y = mean, fill = fill_group),
    width = bar_w
  ) +
  # SOC error bars — right lane
  geom_errorbar(
    data = soc_rows_centered,
    aes(x = x_soc, ymin = lower, ymax = upper),
    width = 0.05, color = "black", linewidth = 0.4
  ) +
  # Ticks at original cluster centers
  scale_x_continuous(
    breaks = sort(unique(plot_df5_with_bgc$x_base)),
    labels = levels(plot_df5_with_bgc$ecosystem_type)
  ) +
  facet_grid(~ climate_zone, scales = "free", space = 'free') +
  scale_fill_manual(values = fill_colors, name = NULL) +
  coord_cartesian(ylim = c(-0.05, ymax2)) +
  labs(
    x = NULL,
    y = expression(atop("Accumulation rate", "(Mg C " * ha^-1 * yr^-1 * ")"))
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = NA, color = NA),
    legend.text = element_text(size = 10),
    legend.key.width = unit(1, "lines"),
    legend.key.size = unit(0.9, "lines")
  ) +
  guides(fill = guide_legend(nrow = 1, ncol = 5, byrow = FALSE))

p_stacked_clean_without_biome_wide_bars



for (dev_type in c("png", "pdf")) {
  if (dev_type == "png")
    png(file.path(figures_dir, "p_stacked_clean_without_biome_wide_bars.png"), width = 8.5, height = 3, unit = "in", res = 440)
  else
    pdf(file.path(figures_dir, "p_stacked_clean_without_biome_wide_bars.pdf"), width = 8.5, height = 3)
  print(p_stacked_clean_without_biome_wide_bars)
  dev.off()
}







# =============================================================================
# SECTION 8 — COMPARISON FIGURE (SANDERMAN MAP)
# =============================================================================
# Reproduces the Sanderman et al. SOC loss map as a reference panel and
# assembles multi-panel composite figures.

p_sanderman_map <- ggplot(r_df, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = expression(
      atop("SOC loss",
           "(Mg C " * ha^-1 * ")")
    )
  ) +
  coord_equal() +
  theme_minimal() +
  theme_blank() +
  annotation_custom(
    grob = hist_grob,
    xmin = -195, xmax = -90, # Adjust for position (lower left)
    ymin = -60, ymax = 20
  ) +
  annotation_custom(
    grob = legend_grob,
    xmin = -145, xmax = -105,
    ymin = -30, ymax = 10
  ) +
  geom_sf(data = limitless_world, inherit.aes = FALSE, fill = 'NA', color = 'grey', linewidth = 0.1)




# Assemble multi-panel composite figure

leg_soc_rate <- get_legend(
  soil_average_rates_50 +
    guides(fill = guide_colorbar(
      barwidth = unit(1.5, "lines"),
      title.position = "top",
      title.hjust = 0,
      label.hjust = 0
    )) +
    theme(
      legend.position = "right",
      legend.justification = c(0, 0.5), # anchor at left+middle
      legend.box.just = "left", # box hugs left
      legend.margin = margin(0,0,0,0),
      legend.box.margin = margin(0,0,0,0),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12,
                                 hjust = 0),
      legend.text.align = 0
    )
)

legratio_map <- get_legend(
  p_ratio_map +
    guides(fill = guide_colorbar(
      barwidth = unit(1.5, "lines"),
      title.position = "top",
      title.hjust = 0,
      label.hjust = 0
    )) +
    theme(
      legend.position = "right",
      legend.justification = c(0, 0.5), # anchor at left+middle
      legend.box.just = "left", # box hugs left
      legend.margin = margin(0,0,0,0),
      legend.box.margin = margin(0,0,0,0),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12,
                                 hjust = 0),
      legend.text.align = 0
    )
)

leg_sanderman <- get_legend(
  p_sanderman_map +
    guides(fill = guide_colorbar(
      barwidth = unit(1.5, "lines"),
      title.position = "top",
      title.hjust = 0,
      label.hjust = 0
    )) +
    theme(
      legend.position = "right",
      legend.justification = c(0, 0.5),
      legend.box.just = "left",
      legend.margin = margin(0,0,0,0),
      legend.box.margin = margin(0,0,0,0),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12,
                                 hjust = 0),
      legend.text.align = 0
    )
)

# 2) Strip legends from the map plots
psoc_rate <- soil_average_rates_50 + theme(plot.margin = unit(c(0,0,0,0), "pt"), legend.position="none")
psanderman <- p_sanderman_map + theme(plot.margin = unit(c(0,0,0,0), "pt"), legend.position="none")
pratio_map <- p_ratio_map + theme(plot.margin = unit(c(0,0,0,0), "pt"), legend.position="none")


# 3) Rebuild each panel as Plot + Legend, forcing the map to be the same relative width
panelsoc_rate <- plot_grid(
  psoc_rate, leg_soc_rate,
  ncol = 2,
  rel_widths = c(1, 0.21), # map gets 1 unit, legend 0.2 units
  align = "h"
)
panelsanderman <- plot_grid(
  psanderman, leg_sanderman,
  ncol = 2,
  rel_widths = c(1, 0.21),
  align = "h"
)

panelratio_map <- plot_grid(
  pratio_map, legratio_map,
  ncol = 2,
  rel_widths = c(1, 0.21), # map gets 1 unit, legend 0.2 units
  align = "h"
)

# 4) Stitch everything together: A | B | C
# final_fig <- plot_grid(
# panelsoc_rate,
# p_mean_bars,
# panelsanderman,
# nrow = 3,
# rel_heights = c(1, 1, 1), 
# align = "hv"
# )

final_fig <- plot_grid(
  panelsoc_rate,
  p_mean_bars,
  panelratio_map,
  nrow = 3, # three across
  rel_heights = c(1, 1, 1), # equal map‐widths
  align = "h",
  axis = "l",
  labels = c("A", "B", "C"), # add labels
  label_size = 16, # pick the size
  #label_fontface = "bold",
  label_x = c(0, 0, 0), # left‐align each label
  label_y = c(1, 1, 1) # top‐align each label
)

# And draw it
for (dev_type in c("png", "pdf")) {
  if (dev_type == "png")
    png(file.path(figures_dir, "figure 5_new_with_average_rates.png"), width = 8.5, height = 8.75, unit = "in", res = 440)
  else
    pdf(file.path(figures_dir, "figure 5_new_with_average_rates.pdf"), width = 8.5, height = 8.75)
  print(final_fig)
  dev.off()
}




# 4-panel composite: SOC rate map / biome bars (AGC+SOC) / AGC:SOC ratio map / Sanderman loss map
final_fig <- plot_grid(
  panelsoc_rate,
  p_mean_bars,
  panelratio_map,
  panelsanderman,
  nrow = 4, # three across
  rel_heights = c(1, 1, 1), # equal map‐widths
  align = "h",
  axis = "l",
  labels = c("A", "B", "C"), # add labels
  label_size = 16, # pick the size
  #label_fontface = "bold",
  label_x = c(0, 0, 0), # left‐align each label
  label_y = c(1, 1, 1) # top‐align each label
)

# And draw it
for (dev_type in c("png", "pdf")) {
  if (dev_type == "png")
    png(file.path(figures_dir, "figure 5_new_four_panel.png"), width = 8.5, height = 11.75, unit = "in", res = 440)
  else
    pdf(file.path(figures_dir, "figure 5_new_four_panel.pdf"), width = 8.5, height = 11.75)
  print(final_fig)
  dev.off()
}


# 4-panel composite: SOC rate map / stacked AGC+BGC+SOC bars / AGC:SOC ratio map / Sanderman loss map
final_fig <- plot_grid(
  panelsoc_rate,
  p_stacked_clean,
  panelratio_map,
  panelsanderman,
  nrow = 4, # three across
  rel_heights = c(1, 1, 1), # equal map‐widths
  align = "h",
  axis = "l",
  labels = c("A", "B", "C"), # add labels
  label_size = 16, # pick the size
  #label_fontface = "bold",
  label_x = c(0, 0, 0), # left‐align each label
  label_y = c(1, 1, 1) # top‐align each label
)

# And draw it
for (dev_type in c("png", "pdf")) {
  if (dev_type == "png")
    png(file.path(figures_dir, "figure 5_new_four_panel_with_bgc.png"), width = 8.5, height = 11.75, unit = "in", res = 440)
  else
    pdf(file.path(figures_dir, "figure 5_new_four_panel_with_bgc.pdf"), width = 8.5, height = 11.75)
  print(final_fig)
  dev.off()
}


# 4-panel composite: SOC rate map / stacked bars with per-land-use SOC / AGC:SOC ratio map / Sanderman loss map
final_fig <- plot_grid(
  panelsoc_rate,
  p_stacked_clean_with_past_crop_bars,
  panelratio_map,
  panelsanderman,
  nrow = 4, # three across
  rel_heights = c(1, 1, 1), # equal map‐widths
  align = "h",
  axis = "l",
  labels = c("A", "B", "C", "D"), # add labels
  label_size = 16, # pick the size
  #label_fontface = "bold",
  label_x = c(0, 0, 0), # left‐align each label
  label_y = c(1, 1, 1) # top‐align each label
)

# And draw it
for (dev_type in c("png", "pdf")) {
  if (dev_type == "png")
    png(file.path(figures_dir, "figure 5_new_four_panel_with_bgc_with_past_crop_bars.png"), width = 8.5, height = 11.75, unit = "in", res = 440)
  else
    pdf(file.path(figures_dir, "figure 5_new_four_panel_with_bgc_with_past_crop_bars.pdf"), width = 8.5, height = 11.75)
  print(final_fig)
  dev.off()
}









# 4-panel composite: restoration-areas-only stacked bars (no biome-wide)
final_fig <- plot_grid(
  panelsoc_rate,
  p_stacked_clean_without_biome_wide_bars,
  panelratio_map,
  panelsanderman,
  nrow = 4, # three across
  rel_heights = c(1, 1, 1), # equal map‐widths
  align = "h",
  axis = "l",
  labels = c("A", "B", "C", "D"), # add labels
  label_size = 16, # pick the size
  #label_fontface = "bold",
  label_x = c(0, 0, 0), # left‐align each label
  label_y = c(1, 1, 1) # top‐align each label
)

# And draw it
for (dev_type in c("png", "pdf")) {
  if (dev_type == "png")
    png(file.path(figures_dir, "figure 5_new_four_panel_without_bione_wide_bars.png"), width = 8.5, height = 11.75, unit = "in", res = 440)
  else
    pdf(file.path(figures_dir, "figure 5_new_four_panel_without_bione_wide_bars.pdf"), width = 8.5, height = 11.75)
  print(final_fig)
  dev.off()
}


# =============================================================================
# SECTION 9 — SSP3-7.0 2041–2060 SOC ESTIMATES
# =============================================================================
# Loads posterior-median and per-sample SOC stock maps produced by running
# 04_bayes_results.R and 05_bayes_uncertainty.R with climate_scenario set to
# "ssp370_2050".  Computes credible-interval Pg C deltas and a supplementary
# table in the same format as the contemporary results.
#
# NOTE: AGC/BGB estimates are not included here because the aboveground carbon
# products (CCI-BM, SBM, Robinson et al. 2025) are parameterised under current
# climate and cannot be projected under future MAT/MAP without additional
# modelling.  All comparisons in Sections 10–11 are therefore SOC-only.

stock_map_dir_fut   <- "./data_output/stock_map/ssp370_2050"
uncertainty_dir_fut <- "./data_output/uncertainty/ssp370_2050"
figures_dir_fut     <- "./figures/ssp370_2050"

dir.create(figures_dir_fut, recursive = TRUE, showWarnings = FALSE)

# Future sample directories (same thinning as contemporary)
sample_dirs_fut <- list.dirs(uncertainty_dir_fut, recursive = FALSE) %>%
  keep(~ grepl("stock_map_", .)) %>%
  keep(~ as.numeric(sub(".*stock_map_", "", .)) %in% seq(1, 3000, by = 10))

# Raster cache must be cleared between contemporary and future runs so that
# get_rast() and get_masks() do not serve stale contemporary rasters
rm(.rast_cache, .mask_cache)
.rast_cache <- new.env(parent = emptyenv())
.mask_cache <- new.env(parent = emptyenv())

# Compute credible intervals for the future scenario
ci_results_fut    <- compute_delta_ci_both(sample_dirs_fut)
ci_wide_0_100_fut <- ci_results_fut$ci_0_100
ci_wide_0_30_fut  <- ci_results_fut$ci_0_30

write.csv(ci_wide_0_100_fut,
          file.path(figures_dir_fut, "delta_ci_0_100.csv"),
          row.names = FALSE)
write.csv(ci_wide_0_30_fut,
          file.path(figures_dir_fut, "delta_ci_0_30.csv"),
          row.names = FALSE)

# Supplementary table: SOC change under ssp370_2050 (same format as contemporary)
supp100_fut <- build_supp_table(ci_wide_0_100_fut, "0–100 cm")
supp30_fut  <- build_supp_table(ci_wide_0_30_fut,  "0–30 cm")
supp_table_fut <- bind_rows(supp100_fut, supp30_fut)

write.csv(
  supp_table_fut,
  file.path(figures_dir_fut, "supplementary_table_SOC_change_ssp370_2050.csv"),
  row.names = FALSE
)


# =============================================================================
# SECTION 10 — SSP3-7.0 MAPS: MEDIAN SOC ACCRUAL AND DIFFERENCE FROM
#              CONTEMPORARY
# =============================================================================
# Produces two global maps:
#   (a) Median SOC accrual rate under ssp370_2050 at year 50 (same style as
#       soil_average_rates_50 in Section 7).
#   (b) Difference map: ssp370_2050 minus contemporary accrual rate at year
#       50, showing where future climate reduces or increases accrual potential.

# ── Helper: reset caches before reading future rasters ─────────────────────
rm(.rast_cache, .mask_cache)
.rast_cache <- new.env(parent = emptyenv())
.mask_cache <- new.env(parent = emptyenv())

# ── (a) Future median accrual rate map ─────────────────────────────────────

prep_for_map_diff_dir <- function(stock_tif, ref_dir = stock_map_dir) {
  # Like prep_for_map_diff but takes explicit dirs so we can mix scenarios
  ref_tif <- file.path(ref_dir, "stock_map_adj_init_crop_0_100.tif")
  raster_df <- prep_for_map(stock_tif)
  raster_df$value <- raster_df$value -
    (prep_for_map(ref_tif) %>% pull(value))
  raster_df
}

soil_crop_fut  <- prep_for_map_diff_dir(
  file.path(stock_map_dir_fut, "stock_map_adj_50_years_CFN_0_100.tif"),
  ref_dir = stock_map_dir_fut
)
soil_past_fut  <- prep_for_map_diff_dir(
  file.path(stock_map_dir_fut, "stock_map_adj_50_years_PFN_0_100.tif"),
  ref_dir = stock_map_dir_fut
)

crop_df_fut  <- soil_crop_fut  %>% rename(crop_soil_diff = value)
past_df_fut  <- soil_past_fut  %>% rename(past_soil_diff = value)

both_df_fut <- crop_df_fut %>%
  mutate(past_soil_diff = past_df_fut$past_soil_diff) %>%
  left_join(crop_frac_df, by = c("x", "y")) %>%
  left_join(past_frac_df, by = c("x", "y")) %>%
  mutate(
    sum_frac = crop_frac + past_frac,
    weighted_soil_diff = if_else(
      sum_frac > 0,
      (crop_frac * crop_soil_diff + past_frac * past_soil_diff) / sum_frac,
      NA_real_
    )
  ) %>%
  select(-sum_frac)

soil_average_rates_50_fut <- ggplot(
  both_df_fut, aes(x = x, y = y, fill = weighted_soil_diff / 50)
) +
  geom_raster() +
  scale_fill_gradient2(
    low     = "#74A5D4",
    mid     = "grey",
    high    = "#8B5A2B",
    midpoint = 0,
    limits  = c(-0.25, 2.25),
    oob     = squish,
    name    = expression(atop("SOC rate (SSP3-7.0)",
                              "(Mg C " * ha^-1 * yr^-1 * ")"))
  ) +
  coord_equal() +
  theme(
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    axis.title      = element_blank(),
    panel.grid      = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA),
    legend.title     = element_text(hjust = 0)
  ) +
  geom_sf(data = limitless_world, inherit.aes = FALSE,
          fill = "NA", color = "grey", linewidth = 0.1)

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png")
    png(file.path(figures_dir_fut, "soil_average_rates_50_ssp370_2050.png"), width = 8.5, height = 4, unit = "in", res = 440)
  else
    pdf(file.path(figures_dir_fut, "soil_average_rates_50_ssp370_2050.pdf"), width = 8.5, height = 4)
  print(soil_average_rates_50_fut)
  dev.off()
}

# ── (b) Difference map: ssp370_2050 minus contemporary ─────────────────────
# Reset caches again before reading contemporary rasters
rm(.rast_cache, .mask_cache)
.rast_cache <- new.env(parent = emptyenv())
.mask_cache <- new.env(parent = emptyenv())

soil_crop_con <- prep_for_map_diff(
  file.path(stock_map_dir, "stock_map_adj_50_years_CFN_0_100.tif")
)
soil_past_con <- prep_for_map_diff(
  file.path(stock_map_dir, "stock_map_adj_50_years_PFN_0_100.tif")
)

# Join on coordinates rather than assuming row alignment — the future and
# contemporary rasters may have slightly different spatial coverage.
both_df_diff <- both_df_fut %>%
  rename(weighted_fut = weighted_soil_diff) %>%
  left_join(
    both_df %>% select(x, y, weighted_con = weighted_soil_diff),
    by = c("x", "y")
  ) %>%
  mutate(diff_per_yr = (weighted_fut - weighted_con) / 50)

soil_rates_diff_map <- ggplot(
  both_df_diff, aes(x = x, y = y, fill = diff_per_yr)
) +
  geom_raster() +
  scale_fill_gradient2(
    low      = "#74A5D4",
    mid      = "grey",
    high     = "#8B5A2B",
    midpoint = 0,
    limits   = c(-0.5, 0.5),
    oob      = squish,
    name     = expression(atop("SOC rate change",
                               "(Mg C " * ha^-1 * yr^-1 * ")"))
  ) +
  coord_equal() +
  theme(
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    axis.title       = element_blank(),
    panel.grid       = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA),
    legend.title     = element_text(hjust = 0)
  ) +
  geom_sf(data = limitless_world, inherit.aes = FALSE,
          fill = "NA", color = "grey", linewidth = 0.1)

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png")
    png(file.path(figures_dir_fut, "soil_rates_diff_ssp370_2050_vs_contemporary.png"), width = 8.5, height = 4, unit = "in", res = 440)
  else
    pdf(file.path(figures_dir_fut, "soil_rates_diff_ssp370_2050_vs_contemporary.pdf"), width = 8.5, height = 4)
  print(soil_rates_diff_map)
  dev.off()
}

# ── Two-panel composite: future map + difference map ───────────────────────
p_map_composite <- plot_grid(
  soil_average_rates_50_fut + theme(legend.position = "right"),
  soil_rates_diff_map        + theme(legend.position = "right"),
  nrow   = 2,
  labels = c("A", "B"),
  label_size = 14,
  align  = "v",
  axis   = "l"
)

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png")
    png(file.path(figures_dir_fut, "maps_ssp370_2050_and_difference.png"), width = 8.5, height = 8, unit = "in", res = 440)
  else
    pdf(file.path(figures_dir_fut, "maps_ssp370_2050_and_difference.pdf"), width = 8.5, height = 8)
  print(p_map_composite)
  dev.off()
}


# =============================================================================
# SECTION 11 — SOC SCENARIO BARS: CONTEMPORARY VS SSP3-7.0 2041–2060
# =============================================================================
# Side-by-side per-hectare bar charts comparing SOC accrual under contemporary
# and ssp370_2050 climates for each scenario and time horizon.  Uses the same
# credible interval structure as the main manuscript bar figures.
# A supplementary table reports the percentage reduction in SOC accrual
# potential (with 95 % CIs) attributable to projected climate change.

# ── Reshape both ci_wide objects into a plottable long format ──────────────

reshape_ci <- function(ci_wide, label) {
  ci_wide %>%
    pivot_longer(
      cols      = -age,
      names_to  = c("scenario", "stat"),
      names_pattern = "(.*)_(lower|lo50|median|hi50|upper)",
      values_to = "value"
    ) %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    mutate(climate = label)
}

soil_long_con <- reshape_ci(ci_wide_0_100,     "Contemporary")
soil_long_fut <- reshape_ci(ci_wide_0_100_fut, "SSP3-7.0 2041–2060")

soil_both <- bind_rows(soil_long_con, soil_long_fut) %>%
  filter(scenario %in% c("crop_aban", "crop_sparing",
                         "crop_relocation", "pasture_sparing",
                         "both_sparing")) %>%
  mutate(
    scenario = factor(scenario,
                      levels = c("crop_aban", "crop_sparing", "crop_relocation",
                                 "pasture_sparing", "both_sparing"),
                      labels = c("Current abandoned\ncropland", "Cropland\nsparing",
                                 "Cropland\nrelocation",       "Pasture\nsparing",
                                 "Maximum\ncropland +\npasture sparing")
    ),
    climate = factor(climate, levels = c("Contemporary", "SSP3-7.0 2041–2060"))
  )

# Per-scenario area denominators (Mha) — each scenario is normalised by its
# own total area so the y-axis reflects Mg C ha⁻¹ yr⁻¹ for that scenario.
# both_sparing sums crop and pasture sparing areas.
scenario_area_Mha <- tibble(
  scenario = c("Current abandoned
cropland",
               "Cropland
sparing",
               "Cropland
relocation",
               "Pasture
sparing",
               "Maximum
cropland +
pasture sparing"),
  area_Mha = c(crop_aban_frac_Mha,
               crop_sparing_frac_Mha,
               crop_relocation_frac_Mha,
               pasture_sparing_frac_Mha,
               crop_sparing_frac_Mha + pasture_sparing_frac_Mha)
)

soil_both_rate <- soil_both %>%
  left_join(scenario_area_Mha, by = "scenario") %>%
  mutate(across(c(lower, lo50, median, hi50, upper), ~ . / area_Mha / age))

# Build panel function matching the style of alt_make_panel
make_climate_panel <- function(dat) {
  scen <- as.character(unique(dat$scenario))
  ggplot(dat, aes(x = age, y = median, colour = climate, fill = climate)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
    geom_ribbon(aes(ymin = lo50,  ymax = hi50),  alpha = 0.20, colour = NA) +
    geom_line(size = 1) +
    geom_point(size = 1.5) +
    scale_colour_manual(
      values = c("Contemporary" = soil_col, "SSP3-7.0 2041–2060" = "#4E84C4"),
      name   = NULL
    ) +
    scale_fill_manual(
      values = c("Contemporary" = soil_col, "SSP3-7.0 2041–2060" = "#4E84C4"),
      name   = NULL
    ) +
    labs(
      x     = "Years since restoration",
      y     = if (scen %in% c("Current abandoned\ncropland",
                              "Maximum\ncropland +\npasture sparing"))
        "SOC change (Mg C ha\u207b\u00b9 yr\u207b\u00b9)" else NULL,
      title = scen
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title   = element_text(face = "italic", hjust = 0.5, size = 10),
      axis.title   = element_text(size = 10),
      axis.text    = element_text(size = 7),
      legend.position = if (scen == "Current abandoned\ncropland") "top" else "none",
      legend.text  = element_text(size = 8)
    )
}

climate_panels <- soil_both_rate %>%
  split(.$scenario) %>%
  purrr::map(make_climate_panel)

p_climate_comparison <- plot_grid(
  plotlist  = climate_panels,
  nrow      = 1,
  labels    = LETTERS[seq_along(climate_panels)],
  label_size = 12,
  align     = "hv",
  axis      = "tblr"
)

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png")
    png(file.path(figures_dir_fut, "soc_contemporary_vs_ssp370_2050_bars.png"), width = 14, height = 3.5, unit = "in", res = 440)
  else
    pdf(file.path(figures_dir_fut, "soc_contemporary_vs_ssp370_2050_bars.pdf"), width = 14, height = 3.5)
  print(p_climate_comparison)
  dev.off()
}

# ── Supplementary table: percent reduction in SOC potential ────────────────
# For each scenario × age × depth, reports:
#   contemporary estimate, ssp370_2050 estimate, and % reduction
#   (contemporary − future) / contemporary × 100, with 95 % CIs propagated
#   by applying the ratio to the contemporary credible bounds.

build_pct_reduction_table <- function(ci_con, ci_fut, depth_label) {
  # Pull median and CI columns for each scenario from both ci_wide objects
  scenarios_vec <- c("crop_aban", "crop_sparing", "crop_relocation",
                     "pasture_sparing", "both_sparing")
  
  map_dfr(scenarios_vec, function(scn) {
    med_con   <- ci_con[[paste0(scn, "_median")]]
    lower_con <- ci_con[[paste0(scn, "_lower")]]
    upper_con <- ci_con[[paste0(scn, "_upper")]]
    med_fut   <- ci_fut[[paste0(scn, "_median")]]
    lower_fut <- ci_fut[[paste0(scn, "_lower")]]
    upper_fut <- ci_fut[[paste0(scn, "_upper")]]
    
    pct_med   <- 100 * (med_con   - med_fut)   / abs(med_con)
    pct_lower <- 100 * (lower_con - upper_fut) / abs(lower_con)  # conservative CI
    pct_upper <- 100 * (upper_con - lower_fut) / abs(upper_con)  # conservative CI
    
    tibble(
      depth_interval   = depth_label,
      age              = ci_con$age,
      scenario         = scn,
      contemporary     = sprintf("%0.1f [%0.1f, %0.1f]",
                                 med_con, lower_con, upper_con),
      ssp370_2050      = sprintf("%0.1f [%0.1f, %0.1f]",
                                 med_fut, lower_fut, upper_fut),
      pct_reduction    = sprintf("%0.1f [%0.1f, %0.1f]",
                                 pct_med, pct_lower, pct_upper)
    )
  })
}

# Add both_sparing to each ci_wide before passing
add_both_sparing <- function(ci_wide) {
  ci_wide %>%
    mutate(
      both_sparing_lower  = crop_sparing_lower  + pasture_sparing_lower,
      both_sparing_median = crop_sparing_median + pasture_sparing_median,
      both_sparing_upper  = crop_sparing_upper  + pasture_sparing_upper
    )
}

ci_con_both <- add_both_sparing(ci_wide_0_100)
ci_fut_both <- add_both_sparing(ci_wide_0_100_fut)

pct_table_100 <- build_pct_reduction_table(ci_con_both, ci_fut_both, "0–100 cm")

ci_con_both30 <- add_both_sparing(ci_wide_0_30)
ci_fut_both30 <- add_both_sparing(ci_wide_0_30_fut)

pct_table_30  <- build_pct_reduction_table(ci_con_both30, ci_fut_both30, "0–30 cm")

pct_table <- bind_rows(pct_table_100, pct_table_30)

write.csv(
  pct_table,
  file.path(figures_dir_fut,
            "supplementary_table_SOC_pct_reduction_ssp370_2050.csv"),
  row.names = FALSE
)