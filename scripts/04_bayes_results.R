# =============================================================================
# 04_bayes_results.R
#
# PURPOSE
# -------
# Applies posterior median coefficients from the hierarchical Bayesian model
# to global raster layers (WorldClim, SoilGrids) to produce spatially explicit
# SOC stock maps for each land-use transition type, depth interval, and time
# horizon.  Depth-integrated stocks are then corrected for coarse fragment
# volume to yield final 0–30 cm and 0–100 cm stock maps.
#
# SCRIPT DEPENDENCIES (must be run first)
# ----------------------------------------
#   03_bayes_analysis.R  — produces bayes_output.RData and sc_df.csv
#
# INPUTS
# ------
#   ./data_output/bayes_output.RData
#   ./data_inputs/sc_df.csv
#   ./data_inputs/raster/via_gee/  — WorldClim (MAP, MAT) and SoilGrids
#                                    (clay, nitrogen, cfvo) rasters by depth
#   ./data_inputs/raster/gmia_plus_0.1.tif  — irrigation raster
#
# OUTPUTS
# -------
#   ./data_inputs/raster/via_gee/*_agg<N>.tif  — aggregated/scaled covariate
#                                                 rasters (cached on first run)
#   ./data_output/rast_component/  — per-interval raster components of the
#                                    linear predictor (init, r, max; one file
#                                    per transition type × depth interval)
#   ./data_output/stock_map/       — per-year × transition × depth stock maps,
#                                    plus coarse-fraction-adjusted composites
#                                    for 0–30 cm and 0–100 cm depth ranges
#
# KEY DESIGN NOTES
# ----------------
#   Transition type codes follow the convention [starting LU][final LU][MGMT]:
#     C = Cropland, P = Pasture, F = Forest, O = Other, A = Active, N = None
#   Index mapping: CFA=1, CFN=2, COA=3, CON=4, PFA=5, PFN=6, POA=7, PON=8
#
#   redo_rasts controls whether raster components are recomputed or read from
#   cached files.  Set to FALSE to skip the computationally expensive
#   integration loop when stock maps already exist.
#
#   agg_fact controls spatial aggregation of input rasters (default 10×).
#   Reduce for higher-resolution outputs at the cost of compute time.
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

ipak(c("tidyverse", "sf", "stars", "MCMCvis"))


# =============================================================================
# CLIMATE SCENARIO CONFIGURATION
# =============================================================================
# Set climate_scenario to "contemporary" for the primary analysis, or to a
# future scenario identifier (e.g. "ssp370_2050") to run with projected
# MAT/MAP rasters.  All raster component, stock map, and figure outputs are
# written to scenario-specific subdirectories so runs never overwrite each other.
#
# For "contemporary":  MAT/MAP sourced from ./data_inputs/raster/via_gee/
# For "ssp370_2050":   MAT/MAP sourced from ./data_inputs/raster/cmip6/ssp370_2050/

climate_scenario <- "contemporary"  # "contemporary" | "ssp370_2050"

# Resolve all scenario-dependent paths from the flag above
climate_raster_dir <- if (climate_scenario == "contemporary") {
  "./data_inputs/raster/via_gee"
} else {
  file.path("./data_inputs/raster/cmip6", climate_scenario)
}

rast_component_dir <- file.path("./data_output/rast_component", climate_scenario)
stock_map_dir      <- file.path("./data_output/stock_map",      climate_scenario)
figures_dir        <- if (climate_scenario == "contemporary") {
  "./figures"
} else {
  file.path("./figures", climate_scenario)
}

# SoilGrids variables (clay, nitrogen, cfvo) are climate-independent and always
# read from via_gee regardless of climate_scenario.  Only map (tap.tif) and
# mat (amt.tif) vary with climate scenario and are sourced from climate_raster_dir.
soilgrids_dir <- "./data_inputs/raster/via_gee"

for (d in c(rast_component_dir, stock_map_dir, figures_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# =============================================================================
# SECTION 1 — HELPER FUNCTIONS
# =============================================================================
# sc_fun converts a covariate from its natural scale to the centred/scaled
# form used in the model.  Pass natural_val for point predictions or sd_val
# to supply a pre-computed z-score directly.

sc_fun <- function(var,
                   natural_val = NA,
                   sd_val = NA){
  if(!is.na(natural_val)[1]){
    natural_val <- ifelse(sc_df[var, "log"], 
                          log(natural_val), 
                          natural_val)
    out_val <- (natural_val - sc_df[var, "center"]) / sc_df[var, "scale"]
    
  } else{
    out_val <- sd_val
  }
  return(out_val)
}


# =============================================================================
# SECTION 2 — LOAD POSTERIOR SAMPLES
# =============================================================================
# Loads the primary MCMC output from 03_bayes_analysis.R.  Posterior medians
# (bsm) and 95 % credible interval bounds (bs025, bs975) are computed across
# all chains for use in point-estimate raster predictions.

load('./data_output/bayes_output.RData')

# Summarise posterior across chains
mod_samp_df <- do.call(rbind, mod_samp)
colMedians <- function(x, na.rm = FALSE) {
  # coerce to matrix if it isn’t one already
  m <- as.matrix(x)
  # apply median over the *columns* (MARGIN = 2)
  apply(m, 2, median, na.rm = na.rm)
}
bsm <- colMedians(mod_samp_df)
bs975 <-apply(mod_samp_df, 2, function(x) quantile(x, 0.975))
bs025 <-apply(mod_samp_df, 2, function(x) quantile(x, 0.025))




# =============================================================================
# SECTION 3 — RASTER SETUP AND CONFIGURATION
# =============================================================================
# Defines transition type codes, SoilGrids depth intervals, and the
# load_scaled_agg() helper that applies optional transforms, log-scaling,
# centering, and spatial aggregation to each covariate raster.  Aggregated
# rasters are cached to disk so subsequent runs skip reprocessing.
#
# redo_rasts — set TRUE to (re)compute all raster components and stock maps;
#              set FALSE to use previously cached outputs.
# agg_fact   — spatial aggregation factor applied to all input rasters.

ipak("terra")

# Land cover types
lu_start <- c("crop", "pasture")
final_type <- c("CFA","CFN","COA","CON",
                "PFA","PFN","POA","PON")

# # A tibble: 8 × 2
# num char 
# <int> <chr>
# 1     1 C F A
# 2     2 C F N
# 3     3 C O A
# 4     4 C O N
# 5     5 P F A
# 6     6 P F N
# 7     7 P O A
# 8     8 P O N


# Set up the intervals used in soilgrids
sg_intervals <- c("0_5", "5_15", "15_30", "30_60", "60_100", "100_200")
sg_interval_lo <- c(0.1, 5, 15, 30, 60, 100) # Note 0.1 cm rather than 0 cm because of log scale
sg_interval_hi <- c(5, 15, 30, 60, 100, 200)

# Decide whether we want to remake the rast_components
redo_rasts <- T # This determines whether the plots should be rewritten

# For testing, give an aggregation factor
agg_fact <- 10

# Helper function to load (or write on first usage) the aggregated raster
load_scaled_agg <- function(varname,
                            interval,
                            sc_df,
                            agg_fact,
                            indir  = climate_raster_dir,
                            outdir = climate_raster_dir,
                            transforms = list(
                              mat = function(x) x / 10   # e.g. divide amt by 10
                            ),
                            snap_to = NULL) {
  # map the variable to the actual input filename
  special_files <- list(
    map         = "tap.tif",
    mat         = "amt.tif",
    irrigation  = "gmia_plus_0.1.tif"
  )
  
  if (varname %in% names(special_files)) {
    infile_name  <- special_files[[varname]]
    outfile_name <- sprintf("%s_agg%d.tif", varname, agg_fact)
  } else {
    infile_name  <- sprintf("%s_%scm_mean.tif", varname, interval)
    outfile_name <- sprintf("%s_%scm_mean_agg%d.tif",
                            varname, interval, agg_fact)
  }
  
  infile  <- file.path(indir,  infile_name)
  outfile <- file.path(outdir, outfile_name)
  
  if (file.exists(outfile)) {
    message("Reading cached: ", outfile)
    return(rast(outfile))
  }
  
  message("Processing & caching: ", outfile)
  r0 <- rast(infile)
  
  # 1) optional var-specific pre-transform
  if (!is.null(transforms[[varname]])) {
    r0 <- transforms[[varname]](r0)
  }
  
  # 2) log if requested
  if (isTRUE(sc_df[varname, "log"])) {
    r0 <- log(r0)
  }
  
  # 3) center & scale
  cen <- sc_df[varname, "center"]
  sc  <- sc_df[varname, "scale"]
  r1  <- (r0 - cen) / sc
  
  # 4) aggregate
  r2  <- aggregate(r1, fact = agg_fact, fun = mean, na.rm = TRUE)
  
  # 5) if a reference raster is supplied, resample to match its grid exactly.
  #    This is needed when climate rasters (map, mat) have a slightly different
  #    extent or origin than the SoilGrids rasters they will be combined with.
  if (!is.null(snap_to)) {
    r2 <- resample(r2, snap_to, method = "bilinear")
  }
  
  # 6) write & return
  writeRaster(r2, filename = outfile, overwrite = TRUE)
  return(r2)
}


# =============================================================================
# SECTION 4 — COMPUTE PER-INTERVAL RASTER COMPONENTS
# =============================================================================
# For each SoilGrids depth interval, evaluates the raster-valued portion of
# the linear predictor for init (crop and pasture), r (all eight transition
# types), and max.  Each component raster is written to
# rast_component_dir and freed from memory immediately.
# The depth coefficient is excluded here and applied analytically during
# integration in Section 5.

if(redo_rasts == T){
  
  # Load a SoilGrids reference raster to use as the snap_to grid for climate
  # rasters.  CMIP6 and via_gee rasters may have slightly different extents
  # after aggregation; resampling map/mat to match clay ensures pixel alignment.
  grid_ref <- load_scaled_agg("clay", "0_5", sc_df, agg_fact,
                              indir = soilgrids_dir, outdir = soilgrids_dir)
  
  # The worldclim rasters, centered and scaled as appropriate
  map <- load_scaled_agg("map", sg_intervals[i], sc_df, agg_fact,
                         snap_to = grid_ref)
  mat <- load_scaled_agg("mat", sg_intervals[i], sc_df, agg_fact,
                         snap_to = grid_ref)
  
  # Also get the irrigation raster, standardized
  irrigation <- load_scaled_agg("irrigation", sg_intervals[i], sc_df, agg_fact,
                                indir = soilgrids_dir, outdir = soilgrids_dir)
  
  # Loop over these intervals to save the rast components at each interval
  for(i in 1:length(sg_intervals)){
    
    print(sg_intervals[i])
    
    # Make filenames
    
    # Init
    rast_component_init_crop_filename <- paste0(rast_component_dir, "/rast_component_init_crop_",
                                                sg_intervals[i],".tif")
    
    rast_component_init_pasture_filename <- paste0(rast_component_dir, "/rast_component_init_pasture_",
                                                   sg_intervals[i],".tif")
    
    # R val "CFA","CFN","COA","CON","PFA","PFN","POA","PON"
    rast_component_r_CFA_filename <- paste0(rast_component_dir, "/rast_component_r_CFA_",
                                            sg_intervals[i],".tif")
    
    rast_component_r_CFN_filename <- paste0(rast_component_dir, "/rast_component_r_CFN_",
                                            sg_intervals[i],".tif")
    
    rast_component_r_COA_filename <- paste0(rast_component_dir, "/rast_component_r_COA_",
                                            sg_intervals[i],".tif")
    
    rast_component_r_CON_filename <- paste0(rast_component_dir, "/rast_component_r_CON_",
                                            sg_intervals[i],".tif")
    
    rast_component_r_PFA_filename <- paste0(rast_component_dir, "/rast_component_r_PFA_",
                                            sg_intervals[i],".tif")
    
    rast_component_r_PFN_filename <- paste0(rast_component_dir, "/rast_component_r_PFN_",
                                            sg_intervals[i],".tif")
    
    rast_component_r_POA_filename <- paste0(rast_component_dir, "/rast_component_r_POA_",
                                            sg_intervals[i],".tif")
    
    rast_component_r_PON_filename <- paste0(rast_component_dir, "/rast_component_r_PON_",
                                            sg_intervals[i],".tif")
    
    # Max
    rast_component_max_filename <- paste0(rast_component_dir, "/rast_component_max_",
                                          sg_intervals[i],".tif")
    
    
    
    
    # Get the appropriate clay and nitrogen rasters for this interval
    clay <- load_scaled_agg("clay", sg_intervals[i], sc_df, agg_fact,
                            indir = soilgrids_dir, outdir = soilgrids_dir) 
    
    nitrogen <- load_scaled_agg("nitrogen", sg_intervals[i], sc_df, agg_fact,
                                indir = soilgrids_dir, outdir = soilgrids_dir)
    
    # Init rast component
    rast_component_init_crop <- bsm["beta_init_intercept[1]"] +
      
      bsm["beta_init_clay[1]"] * clay + 
      
      bsm["beta_init_nitrogen[1]"] * nitrogen + 
      
      bsm["beta_init_map_mat[1]"] * map * mat +
      
      bsm["beta_init_irrigation[1]"] * irrigation + 
      
      bsm["beta_init_map[1]"] * map + 
      
      bsm["beta_init_mat[1]"] * mat
    
    writeRaster(rast_component_init_crop,
                overwrite = T,
                filename = rast_component_init_crop_filename)
    rm(rast_component_init_crop)
    gc()
    
    
    rast_component_init_pasture <- bsm["beta_init_intercept[2]"] +
      
      bsm["beta_init_clay[2]"] * clay + 
      
      bsm["beta_init_nitrogen[2]"] * nitrogen + 
      
      bsm["beta_init_map_mat[2]"] * map * mat +
      
      bsm["beta_init_irrigation[2]"] * irrigation + 
      
      bsm["beta_init_map[2]"] * map + 
      
      bsm["beta_init_mat[2]"] * mat
    
    writeRaster(rast_component_init_pasture,
                overwrite = T,
                filename = rast_component_init_pasture_filename)
    rm(rast_component_init_pasture)
    gc()
    
    # R rast component "CFA","CFN","COA","CON","PFA","PFN","POA","PON"
    # The first four are for initially cropland sites
    
    # This chunk the same for final types 1-4, so will calculate only once
    rast_subcomponent_r_C <- bsm["beta_r_clay[1]"] * clay + # This and the next rows' indices represents the starting land use, crop (1) or pasture (2)
      
      bsm["beta_r_nitrogen[1]"] * nitrogen +
      
      bsm["beta_r_map_mat[1]"] * map * mat +
      
      bsm["beta_r_map[1]"] * map +
      
      bsm["beta_r_mat[1]"] * mat
    
    rast_component_r_CFA <- bsm["beta_r_intercept[1]"] + # This index represents the final type (1-8)
      rast_subcomponent_r_C
    writeRaster(rast_component_r_CFA,
                overwrite = T,
                filename = rast_component_r_CFA_filename)
    rm(rast_component_r_CFA)
    gc()
    
    rast_component_r_CFN <- bsm["beta_r_intercept[2]"] + # This index represents the final type (1-8)
      rast_subcomponent_r_C
    writeRaster(rast_component_r_CFN,
                overwrite = T,
                filename = rast_component_r_CFN_filename)
    rm(rast_component_r_CFN)
    gc()
    
    rast_component_r_COA <- bsm["beta_r_intercept[3]"] + # This index represents the final type (1-8)
      rast_subcomponent_r_C
    writeRaster(rast_component_r_COA,
                overwrite = T,
                filename = rast_component_r_COA_filename)
    rm(rast_component_r_COA)
    gc()
    
    rast_component_r_CON <- bsm["beta_r_intercept[4]"] + # This index represents the final type (1-8)
      rast_subcomponent_r_C
    writeRaster(rast_component_r_CON,
                overwrite = T,
                filename = rast_component_r_CON_filename)
    rm(rast_component_r_CON)
    gc()
    
    rm(rast_subcomponent_r_C)
    gc()
    
    # The next four are for initially pasture sites   
    # This chunk the same for final types 1-4, so will calculate only once
    rast_subcomponent_r_P <- bsm["beta_r_clay[2]"] * clay + # This and the next rows' indices represents the starting land use, crop (1) or pasture (2)
      
      bsm["beta_r_nitrogen[2]"] * nitrogen +
      
      bsm["beta_r_map_mat[2]"] * map * mat +
      
      bsm["beta_r_map[2]"] * map +
      
      bsm["beta_r_mat[2]"] * mat
    
    rast_component_r_PFA <- bsm["beta_r_intercept[5]"] + # This index represents the final type (1-8)
      rast_subcomponent_r_P
    writeRaster(rast_component_r_PFA,
                overwrite = T,
                filename = rast_component_r_PFA_filename)
    rm(rast_component_r_PFA)
    gc()
    
    rast_component_r_PFN <- bsm["beta_r_intercept[6]"] + # This index represents the final type (1-8)
      rast_subcomponent_r_P
    writeRaster(rast_component_r_PFN,
                overwrite = T,
                filename = rast_component_r_PFN_filename)
    rm(rast_component_r_PFN)
    gc()
    
    
    rast_component_r_POA <- bsm["beta_r_intercept[7]"] + # This index represents the final type (1-8)
      rast_subcomponent_r_P
    writeRaster(rast_component_r_POA,
                overwrite = T,
                filename = rast_component_r_POA_filename)
    rm(rast_component_r_POA)
    gc()
    
    rast_component_r_PON <- bsm["beta_r_intercept[8]"] + # This index represents the final type (1-8)
      rast_subcomponent_r_P
    writeRaster(rast_component_r_PON,
                overwrite = T,
                filename = rast_component_r_PON_filename)
    rm(rast_component_r_PON)
    gc()
    
    rm(rast_subcomponent_r_P)
    gc()
    
    # Max
    rast_component_max <- bsm["beta_max_intercept"] +
      bsm["beta_max_clay"] * clay +
      
      bsm["beta_max_nitrogen"] * nitrogen +
      
      bsm["beta_max_map_mat"] * map * mat +
      
      bsm["beta_max_map"] * map +
      
      bsm["beta_max_mat"] * mat
    
    writeRaster(rast_component_max,
                overwrite = T,
                filename = rast_component_max_filename)
    rm(rast_component_max)
    gc()
    
  }
}


# =============================================================================
# SECTION 5 — DEPTH INTEGRATION AND STOCK MAP GENERATION
# =============================================================================
# Completes the linear predictor by adding the depth coefficient analytically
# and numerically integrates the resulting SOC density function over each
# depth interval using pracma::quadv.  Produces one stock map raster per
# transition type × depth interval × time point and writes to stock_map_dir.

# Time points (years since abandonment) for which stock maps are produced
time_points <- c(5, 15, 30, 50, 100, 200)


time1 <- Sys.time()
if(redo_rasts){
  for(yr in time_points){
    print(paste("starting on year", yr, "at", Sys.time()))
    for(i in 1:length(sg_intervals)){
      print(paste("--",sg_intervals[i], "at", Sys.time()))
      
      # Get rast component filenames
      # Init
      rast_component_init_crop_filename <- paste0(rast_component_dir, "/rast_component_init_crop_",
                                                  sg_intervals[i],".tif")
      
      rast_component_init_pasture_filename <- paste0(rast_component_dir, "/rast_component_init_pasture_",
                                                     sg_intervals[i],".tif")
      
      # R val "CFA","CFN","COA","CON","PFA","PFN","POA","PON"
      rast_component_r_CFA_filename <- paste0(rast_component_dir, "/rast_component_r_CFA_",
                                              sg_intervals[i],".tif")
      
      rast_component_r_CFN_filename <- paste0(rast_component_dir, "/rast_component_r_CFN_",
                                              sg_intervals[i],".tif")
      
      rast_component_r_COA_filename <- paste0(rast_component_dir, "/rast_component_r_COA_",
                                              sg_intervals[i],".tif")
      
      rast_component_r_CON_filename <- paste0(rast_component_dir, "/rast_component_r_CON_",
                                              sg_intervals[i],".tif")
      
      rast_component_r_PFA_filename <- paste0(rast_component_dir, "/rast_component_r_PFA_",
                                              sg_intervals[i],".tif")
      
      rast_component_r_PFN_filename <- paste0(rast_component_dir, "/rast_component_r_PFN_",
                                              sg_intervals[i],".tif")
      
      rast_component_r_POA_filename <- paste0(rast_component_dir, "/rast_component_r_POA_",
                                              sg_intervals[i],".tif")
      
      rast_component_r_PON_filename <- paste0(rast_component_dir, "/rast_component_r_PON_",
                                              sg_intervals[i],".tif")
      
      # Max
      rast_component_max_filename <- paste0(rast_component_dir, "/rast_component_max_",
                                            sg_intervals[i],".tif")
      
      
      
      final_type_filenames <- c(
        rast_component_r_CFA_filename,rast_component_r_CFN_filename,
        rast_component_r_COA_filename,rast_component_r_CON_filename,
        rast_component_r_PFA_filename,rast_component_r_PFN_filename,
        rast_component_r_POA_filename,rast_component_r_PON_filename)
      
      
      
      
      
      # Now perform the integration to get to stock
      # Cropland
      # Read in rast_component_init_crop_filename
      rast_component_init_crop <- rast(rast_component_init_crop_filename)
      # Define a function to integrate over with respect to depth (x) for cropland
      depth_func_init_crop <- function(x, init_crop_rast) {
        return(exp(bsm["beta_init_depth[1]"] * ((log(x) - sc_df["depth", "center"]) / sc_df["depth", "scale"]) + init_crop_rast))
      }
      # This wont work with any NAs, so just replace them with zero for now, and
      # later, we will put the NAs back in.
      na_inds <- is.na(rast_component_init_crop) | is.infinite(rast_component_init_crop)
      rast_component_init_crop[na_inds] <- 0
      crop_stock <- app(rast_component_init_crop, 
                        fun = function(pixel_value) {
                          pracma::quadv(depth_func_init_crop, 
                                        a = sg_interval_lo[i], 
                                        b = sg_interval_hi[i], 
                                        init_crop_rast = pixel_value)$Q})
      crop_stock[na_inds] <- NA
      writeRaster(crop_stock, 
                  filename = gsub("rast_component", "stock_map", 
                                  rast_component_init_crop_filename),
                  overwrite = T)
      
      # Pasture
      # Read in rast_component_init_pasture_filename
      rast_component_init_pasture <- rast(rast_component_init_pasture_filename)
      # Define a function to integrate over with respect to depth (x) for pastureland
      depth_func_init_pasture <- function(x, init_pasture_rast) {
        return(exp(bsm["beta_init_depth[2]"] * ((log(x) - sc_df["depth", "center"]) / sc_df["depth", "scale"]) + init_pasture_rast))
      }
      na_inds_pasture <- is.na(rast_component_init_pasture) | is.infinite(rast_component_init_pasture)
      rast_component_init_pasture[na_inds_pasture] <- 0
      pasture_stock <- app(rast_component_init_pasture, 
                           fun = function(pixel_value) {
                             pracma::quadv(depth_func_init_pasture, 
                                           a = sg_interval_lo[i], 
                                           b = sg_interval_hi[i], 
                                           init_pasture_rast = pixel_value)$Q})
      pasture_stock[na_inds_pasture] <- NA
      writeRaster(pasture_stock, 
                  filename = gsub("rast_component", "stock_map", 
                                  rast_component_init_pasture_filename),
                  overwrite = T)
      
      # Max
      # Read in rast_component_max
      rast_component_max <- rast(rast_component_max_filename)
      # Define a function to integrate over with respect to depth (x) for max
      depth_func_max <- function(x, max_rast) {
        return(exp(bsm["beta_max_depth"] * ((log(x) - sc_df["depth", "center"]) / sc_df["depth", "scale"]) + max_rast))
      }
      na_inds_max <- is.na(rast_component_max) | is.infinite(rast_component_max)
      rast_component_max[na_inds_max] <- 0
      max_stock <- app(rast_component_max, 
                       fun = function(pixel_value) {
                         pracma::quadv(depth_func_max, 
                                       a = sg_interval_lo[i], 
                                       b = sg_interval_hi[i], 
                                       max_rast = pixel_value)$Q})
      max_stock[na_inds_max] <- NA
      writeRaster(max_stock, 
                  filename = gsub("rast_component", "stock_map", 
                                  rast_component_max_filename),
                  overwrite = T)
      
      
      
      # Now for each of the eight final types.
      # When there are multiple rasters used together, need to format the integration
      # in a somewhat different way.
      
      compute_stock_func <- function(init_rast, r_rast, max_rast, final_type_num) {
        
        pracma::quadv(function(x, 
                               init_rast = init_rast,
                               r_rast = r_rast,
                               max_rast = max_rast,
                               final_type_num = final_type_num){
          
          # First deal with indexing used for the init, r_val, and m_val
          init_string <- ifelse(final_type_num < 5, 
                                "beta_init_depth[1]", 
                                "beta_init_depth[2]")
          
          r_string <- ifelse(final_type_num < 5, 
                             "beta_r_depth[1]", 
                             "beta_r_depth[2]")
          
          m_string <- paste0("m_val[", final_type_num, "]")
          
          depth_val <- ((log(x) - sc_df["depth", "center"]) / sc_df["depth", "scale"])
          
          init_val <- bsm[init_string] * depth_val + init_rast
          
          r_val <- bsm[r_string] * depth_val + r_rast
          
          max_val <- bsm["beta_max_depth"] * depth_val + max_rast
          
          m_val <- bsm[m_string]
          
          # Now calculate the output
          exp(init_val) + 
            (exp(max_val) - exp(init_val)) * 
            (1 - exp(-yr/exp(r_val))) ^ exp(m_val)
          
        },
        
        a = sg_interval_lo[i],
        b = sg_interval_hi[i],
        init_rast = init_rast,
        r_rast = r_rast,
        max_rast = max_rast,
        final_type_num = final_type_num)$Q
        
      }
      
      # Now loop over the final types to compute their stocks
      for(j in 1:length(final_type_filenames)){
        init_filename <- ifelse(j < 5, 
                                rast_component_init_crop_filename,
                                rast_component_init_pasture_filename)
        
        rast_final_type <- rast(init_filename) %>% 
          c(rast(final_type_filenames[j])) %>% 
          c(rast(rast_component_max_filename))
        
        # Mask out any NA or Inf values across all three layers before integration
        na_inds_j <- is.na(rast_final_type[[1]]) | is.infinite(rast_final_type[[1]]) |
          is.na(rast_final_type[[2]]) | is.infinite(rast_final_type[[2]]) |
          is.na(rast_final_type[[3]]) | is.infinite(rast_final_type[[3]])
        rast_final_type[[1]][na_inds_j] <- 0
        rast_final_type[[2]][na_inds_j] <- 0
        rast_final_type[[3]][na_inds_j] <- 0
        
        stock_map_final_type <- lapp(rast_final_type, compute_stock_func, final_type_num = j)
        stock_map_final_type[na_inds_j] <- NA
        
        writeRaster(stock_map_final_type, 
                    filename = gsub("rast_component", "stock_map", 
                                    final_type_filenames[j]) %>% 
                      gsub("map_", paste0("map_", yr, "_"), .) %>% 
                      gsub("_r_", "_years_", .),
                    overwrite = T)
        
      }
    }
  }
}
time2 <- Sys.time()

# Clean up raster objects no longer needed after integration
rm(map, mat, clay, nitrogen, irrigation, grid_ref)
gc()




# =============================================================================
# SECTION 6 — DEPTH AGGREGATION AND COARSE FRACTION CORRECTION
# =============================================================================
# Combines per-interval stock maps into 0–30 cm and 0–100 cm column totals,
# weighting each layer by (1 – coarse fragment volume fraction) derived from
# SoilGrids cfvo rasters.  Produces one adjusted composite raster per
# transition type × time point × depth range and writes to
# stock_map_dir with the "adj" prefix.

# Load and aggregate the coarse-fraction correction factors by depth interval
cfc_0_5 <- (1 - rast("./data_inputs/raster/via_gee/cfvo_0_5cm_mean.tif") %>% 
              aggregate(fact = agg_fact, fun = mean, na.rm = T)/1000)
cfc_5_15 <- (1 - rast("./data_inputs/raster/via_gee/cfvo_5_15cm_mean.tif") %>% 
               aggregate(fact = agg_fact, fun = mean, na.rm = T)/1000)
cfc_15_30 <- (1 - rast("./data_inputs/raster/via_gee/cfvo_15_30cm_mean.tif") %>% 
                aggregate(fact = agg_fact, fun = mean, na.rm = T)/1000)
cfc_30_60 <- (1 - rast("./data_inputs/raster/via_gee/cfvo_30_60cm_mean.tif") %>% 
                aggregate(fact = agg_fact, fun = mean, na.rm = T)/1000)
cfc_60_100 <- (1 - rast("./data_inputs/raster/via_gee/cfvo_60_100cm_mean.tif") %>% 
                 aggregate(fact = agg_fact, fun = mean, na.rm = T)/1000)
cfc_100_200 <- (1 - rast("./data_inputs/raster/via_gee/cfvo_100_200cm_mean.tif") %>% 
                  aggregate(fact = agg_fact, fun = mean, na.rm = T)/1000)


# 0–100 cm composites: init (time-independent) and max ────────────────────────

# Cropland initial stock
stock_map_adj_init_crop_0_100 <- rast(file.path(stock_map_dir, "stock_map_init_crop_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_init_crop_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_init_crop_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_init_crop_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_init_crop_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_init_crop_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_init_crop_0_100.tif"),
            overwrite = T)

# Pasture initial stock
stock_map_adj_init_pasture_0_100 <- rast(file.path(stock_map_dir, "stock_map_init_pasture_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_init_pasture_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_init_pasture_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_init_pasture_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_init_pasture_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_init_pasture_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_init_pasture_0_100.tif"),
            overwrite = T)

# Maximum attainable stock (time-independent)
stock_map_adj_max_0_100 <- rast(file.path(stock_map_dir, "stock_map_max_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_max_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_max_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_max_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_max_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_max_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_max_0_100.tif"),
            overwrite = T)



# Now the other types for each year "CFA","CFN","COA","CON","PFA","PFN","POA","PON"

# 5_years
stock_map_adj_5_years_CFA_0_100 <- rast(file.path(stock_map_dir, "stock_map_5_years_CFA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_5_years_CFA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_5_years_CFA_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_5_years_CFA_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_5_years_CFA_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_5_years_CFA_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_5_years_CFA_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_5_years_CFA_0_100)

stock_map_adj_5_years_CFN_0_100 <- rast(file.path(stock_map_dir, "stock_map_5_years_CFN_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_5_years_CFN_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_5_years_CFN_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_5_years_CFN_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_5_years_CFN_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_5_years_CFN_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_5_years_CFN_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_5_years_CFN_0_100)

stock_map_adj_5_years_COA_0_100 <- rast(file.path(stock_map_dir, "stock_map_5_years_COA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_5_years_COA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_5_years_COA_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_5_years_COA_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_5_years_COA_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_5_years_COA_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_5_years_COA_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_5_years_COA_0_100)

stock_map_adj_5_years_CON_0_100 <- rast(file.path(stock_map_dir, "stock_map_5_years_CON_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_5_years_CON_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_5_years_CON_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_5_years_CON_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_5_years_CON_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_5_years_CON_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_5_years_CON_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_5_years_CON_0_100)

stock_map_adj_5_years_PFA_0_100 <- rast(file.path(stock_map_dir, "stock_map_5_years_PFA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_5_years_PFA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_5_years_PFA_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_5_years_PFA_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_5_years_PFA_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_5_years_PFA_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_5_years_PFA_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_5_years_PFA_0_100)

stock_map_adj_5_years_PFN_0_100 <- rast(file.path(stock_map_dir, "stock_map_5_years_PFN_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_5_years_PFN_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_5_years_PFN_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_5_years_PFN_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_5_years_PFN_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_5_years_PFN_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_5_years_PFN_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_5_years_PFN_0_100)

stock_map_adj_5_years_POA_0_100 <- rast(file.path(stock_map_dir, "stock_map_5_years_POA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_5_years_POA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_5_years_POA_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_5_years_POA_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_5_years_POA_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_5_years_POA_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_5_years_POA_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_5_years_POA_0_100)

stock_map_adj_5_years_PON_0_100 <- rast(file.path(stock_map_dir, "stock_map_5_years_PON_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_5_years_PON_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_5_years_PON_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_5_years_PON_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_5_years_PON_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_5_years_PON_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_5_years_PON_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_5_years_PON_0_100)


# 15_years
stock_map_adj_15_years_CFA_0_100 <- rast(file.path(stock_map_dir, "stock_map_15_years_CFA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_15_years_CFA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_15_years_CFA_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_15_years_CFA_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_15_years_CFA_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_15_years_CFA_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_15_years_CFA_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_15_years_CFA_0_100)

stock_map_adj_15_years_CFN_0_100 <- rast(file.path(stock_map_dir, "stock_map_15_years_CFN_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_15_years_CFN_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_15_years_CFN_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_15_years_CFN_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_15_years_CFN_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_15_years_CFN_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_15_years_CFN_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_15_years_CFN_0_100)

stock_map_adj_15_years_COA_0_100 <- rast(file.path(stock_map_dir, "stock_map_15_years_COA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_15_years_COA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_15_years_COA_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_15_years_COA_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_15_years_COA_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_15_years_COA_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_15_years_COA_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_15_years_COA_0_100)

stock_map_adj_15_years_CON_0_100 <- rast(file.path(stock_map_dir, "stock_map_15_years_CON_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_15_years_CON_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_15_years_CON_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_15_years_CON_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_15_years_CON_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_15_years_CON_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_15_years_CON_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_15_years_CON_0_100)

stock_map_adj_15_years_PFA_0_100 <- rast(file.path(stock_map_dir, "stock_map_15_years_PFA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_15_years_PFA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_15_years_PFA_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_15_years_PFA_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_15_years_PFA_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_15_years_PFA_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_15_years_PFA_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_15_years_PFA_0_100)

stock_map_adj_15_years_PFN_0_100 <- rast(file.path(stock_map_dir, "stock_map_15_years_PFN_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_15_years_PFN_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_15_years_PFN_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_15_years_PFN_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_15_years_PFN_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_15_years_PFN_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_15_years_PFN_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_15_years_PFN_0_100)

stock_map_adj_15_years_POA_0_100 <- rast(file.path(stock_map_dir, "stock_map_15_years_POA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_15_years_POA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_15_years_POA_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_15_years_POA_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_15_years_POA_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_15_years_POA_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_15_years_POA_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_15_years_POA_0_100)

stock_map_adj_15_years_PON_0_100 <- rast(file.path(stock_map_dir, "stock_map_15_years_PON_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_15_years_PON_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_15_years_PON_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_15_years_PON_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_15_years_PON_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_15_years_PON_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_15_years_PON_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_15_years_PON_0_100)


# 30_years
stock_map_adj_30_years_CFA_0_100 <- rast(file.path(stock_map_dir, "stock_map_30_years_CFA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_30_years_CFA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_30_years_CFA_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_30_years_CFA_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_30_years_CFA_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_30_years_CFA_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_30_years_CFA_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_30_years_CFA_0_100)

stock_map_adj_30_years_CFN_0_100 <- rast(file.path(stock_map_dir, "stock_map_30_years_CFN_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_30_years_CFN_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_30_years_CFN_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_30_years_CFN_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_30_years_CFN_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_30_years_CFN_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_30_years_CFN_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_30_years_CFN_0_100)

stock_map_adj_30_years_COA_0_100 <- rast(file.path(stock_map_dir, "stock_map_30_years_COA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_30_years_COA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_30_years_COA_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_30_years_COA_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_30_years_COA_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_30_years_COA_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_30_years_COA_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_30_years_COA_0_100)

stock_map_adj_30_years_CON_0_100 <- rast(file.path(stock_map_dir, "stock_map_30_years_CON_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_30_years_CON_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_30_years_CON_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_30_years_CON_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_30_years_CON_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_30_years_CON_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_30_years_CON_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_30_years_CON_0_100)

stock_map_adj_30_years_PFA_0_100 <- rast(file.path(stock_map_dir, "stock_map_30_years_PFA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_30_years_PFA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_30_years_PFA_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_30_years_PFA_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_30_years_PFA_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_30_years_PFA_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_30_years_PFA_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_30_years_PFA_0_100)

stock_map_adj_30_years_PFN_0_100 <- rast(file.path(stock_map_dir, "stock_map_30_years_PFN_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_30_years_PFN_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_30_years_PFN_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_30_years_PFN_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_30_years_PFN_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_30_years_PFN_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_30_years_PFN_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_30_years_PFN_0_100)

stock_map_adj_30_years_POA_0_100 <- rast(file.path(stock_map_dir, "stock_map_30_years_POA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_30_years_POA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_30_years_POA_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_30_years_POA_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_30_years_POA_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_30_years_POA_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_30_years_POA_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_30_years_POA_0_100)

stock_map_adj_30_years_PON_0_100 <- rast(file.path(stock_map_dir, "stock_map_30_years_PON_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_30_years_PON_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_30_years_PON_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_30_years_PON_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_30_years_PON_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_30_years_PON_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_30_years_PON_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_30_years_PON_0_100)

# 50_years
stock_map_adj_50_years_CFA_0_100 <- rast(file.path(stock_map_dir, "stock_map_50_years_CFA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_50_years_CFA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_50_years_CFA_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_50_years_CFA_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_50_years_CFA_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_50_years_CFA_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_50_years_CFA_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_50_years_CFA_0_100)

stock_map_adj_50_years_CFN_0_100 <- rast(file.path(stock_map_dir, "stock_map_50_years_CFN_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_50_years_CFN_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_50_years_CFN_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_50_years_CFN_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_50_years_CFN_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_50_years_CFN_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_50_years_CFN_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_50_years_CFN_0_100)

stock_map_adj_50_years_COA_0_100 <- rast(file.path(stock_map_dir, "stock_map_50_years_COA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_50_years_COA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_50_years_COA_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_50_years_COA_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_50_years_COA_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_50_years_COA_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_50_years_COA_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_50_years_COA_0_100)

stock_map_adj_50_years_CON_0_100 <- rast(file.path(stock_map_dir, "stock_map_50_years_CON_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_50_years_CON_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_50_years_CON_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_50_years_CON_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_50_years_CON_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_50_years_CON_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_50_years_CON_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_50_years_CON_0_100)

stock_map_adj_50_years_PFA_0_100 <- rast(file.path(stock_map_dir, "stock_map_50_years_PFA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_50_years_PFA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_50_years_PFA_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_50_years_PFA_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_50_years_PFA_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_50_years_PFA_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_50_years_PFA_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_50_years_PFA_0_100)

stock_map_adj_50_years_PFN_0_100 <- rast(file.path(stock_map_dir, "stock_map_50_years_PFN_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_50_years_PFN_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_50_years_PFN_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_50_years_PFN_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_50_years_PFN_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_50_years_PFN_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_50_years_PFN_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_50_years_PFN_0_100)

stock_map_adj_50_years_POA_0_100 <- rast(file.path(stock_map_dir, "stock_map_50_years_POA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_50_years_POA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_50_years_POA_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_50_years_POA_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_50_years_POA_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_50_years_POA_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_50_years_POA_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_50_years_POA_0_100)

stock_map_adj_50_years_PON_0_100 <- rast(file.path(stock_map_dir, "stock_map_50_years_PON_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_50_years_PON_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_50_years_PON_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_50_years_PON_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_50_years_PON_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_50_years_PON_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_50_years_PON_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_50_years_PON_0_100)


# 100_years
stock_map_adj_100_years_CFA_0_100 <- rast(file.path(stock_map_dir, "stock_map_100_years_CFA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_100_years_CFA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_100_years_CFA_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_100_years_CFA_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_100_years_CFA_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_100_years_CFA_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_100_years_CFA_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_100_years_CFA_0_100)

stock_map_adj_100_years_CFN_0_100 <- rast(file.path(stock_map_dir, "stock_map_100_years_CFN_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_100_years_CFN_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_100_years_CFN_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_100_years_CFN_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_100_years_CFN_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_100_years_CFN_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_100_years_CFN_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_100_years_CFN_0_100)

stock_map_adj_100_years_COA_0_100 <- rast(file.path(stock_map_dir, "stock_map_100_years_COA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_100_years_COA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_100_years_COA_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_100_years_COA_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_100_years_COA_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_100_years_COA_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_100_years_COA_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_100_years_COA_0_100)

stock_map_adj_100_years_CON_0_100 <- rast(file.path(stock_map_dir, "stock_map_100_years_CON_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_100_years_CON_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_100_years_CON_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_100_years_CON_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_100_years_CON_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_100_years_CON_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_100_years_CON_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_100_years_CON_0_100)

stock_map_adj_100_years_PFA_0_100 <- rast(file.path(stock_map_dir, "stock_map_100_years_PFA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_100_years_PFA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_100_years_PFA_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_100_years_PFA_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_100_years_PFA_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_100_years_PFA_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_100_years_PFA_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_100_years_PFA_0_100)

stock_map_adj_100_years_PFN_0_100 <- rast(file.path(stock_map_dir, "stock_map_100_years_PFN_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_100_years_PFN_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_100_years_PFN_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_100_years_PFN_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_100_years_PFN_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_100_years_PFN_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_100_years_PFN_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_100_years_PFN_0_100)

stock_map_adj_100_years_POA_0_100 <- rast(file.path(stock_map_dir, "stock_map_100_years_POA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_100_years_POA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_100_years_POA_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_100_years_POA_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_100_years_POA_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_100_years_POA_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_100_years_POA_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_100_years_POA_0_100)

stock_map_adj_100_years_PON_0_100 <- rast(file.path(stock_map_dir, "stock_map_100_years_PON_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_100_years_PON_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_100_years_PON_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_100_years_PON_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_100_years_PON_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_100_years_PON_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_100_years_PON_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_100_years_PON_0_100)


# 200_years
stock_map_adj_200_years_CFA_0_100 <- rast(file.path(stock_map_dir, "stock_map_200_years_CFA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_200_years_CFA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_200_years_CFA_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_200_years_CFA_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_200_years_CFA_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_200_years_CFA_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_200_years_CFA_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_200_years_CFA_0_100)

stock_map_adj_200_years_CFN_0_100 <- rast(file.path(stock_map_dir, "stock_map_200_years_CFN_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_200_years_CFN_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_200_years_CFN_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_200_years_CFN_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_200_years_CFN_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_200_years_CFN_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_200_years_CFN_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_200_years_CFN_0_100)

stock_map_adj_200_years_COA_0_100 <- rast(file.path(stock_map_dir, "stock_map_200_years_COA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_200_years_COA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_200_years_COA_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_200_years_COA_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_200_years_COA_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_200_years_COA_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_200_years_COA_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_200_years_COA_0_100)

stock_map_adj_200_years_CON_0_100 <- rast(file.path(stock_map_dir, "stock_map_200_years_CON_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_200_years_CON_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_200_years_CON_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_200_years_CON_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_200_years_CON_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_200_years_CON_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_200_years_CON_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_200_years_CON_0_100)

stock_map_adj_200_years_PFA_0_100 <- rast(file.path(stock_map_dir, "stock_map_200_years_PFA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_200_years_PFA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_200_years_PFA_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_200_years_PFA_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_200_years_PFA_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_200_years_PFA_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_200_years_PFA_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_200_years_PFA_0_100)

stock_map_adj_200_years_PFN_0_100 <- rast(file.path(stock_map_dir, "stock_map_200_years_PFN_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_200_years_PFN_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_200_years_PFN_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_200_years_PFN_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_200_years_PFN_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_200_years_PFN_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_200_years_PFN_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_200_years_PFN_0_100)

stock_map_adj_200_years_POA_0_100 <- rast(file.path(stock_map_dir, "stock_map_200_years_POA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_200_years_POA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_200_years_POA_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_200_years_POA_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_200_years_POA_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_200_years_POA_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_200_years_POA_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_200_years_POA_0_100)

stock_map_adj_200_years_PON_0_100 <- rast(file.path(stock_map_dir, "stock_map_200_years_PON_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_200_years_PON_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_200_years_PON_15_30.tif")) * cfc_15_30 +
  rast(file.path(stock_map_dir, "stock_map_200_years_PON_30_60.tif")) * cfc_30_60 +
  rast(file.path(stock_map_dir, "stock_map_200_years_PON_60_100.tif")) * cfc_60_100
writeRaster(stock_map_adj_200_years_PON_0_100,
            file = file.path(stock_map_dir, "stock_map_adj_200_years_PON_0_100.tif"),
            overwrite = T)
rm(stock_map_adj_200_years_PON_0_100)



# 0 to 30 cm version

# crop (which is independent of year)
stock_map_adj_init_crop_0_30 <- rast(file.path(stock_map_dir, "stock_map_init_crop_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_init_crop_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_init_crop_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_init_crop_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_init_crop_0_30.tif"),
            overwrite = T)

# pasture (which is independent of year)
stock_map_adj_init_pasture_0_30 <- rast(file.path(stock_map_dir, "stock_map_init_pasture_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_init_pasture_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_init_pasture_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_init_pasture_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_init_pasture_0_30.tif"),
            overwrite = T)

# max (which is independent of year)
stock_map_adj_max_0_30 <- rast(file.path(stock_map_dir, "stock_map_max_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_max_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_max_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_max_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_max_0_30.tif"),
            overwrite = T)



# Now the other types for each year "CFA","CFN","COA","CON","PFA","PFN","POA","PON"

# 5_years
stock_map_adj_5_years_CFA_0_30 <- rast(file.path(stock_map_dir, "stock_map_5_years_CFA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_5_years_CFA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_5_years_CFA_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_5_years_CFA_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_5_years_CFA_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_5_years_CFA_0_30)

stock_map_adj_5_years_CFN_0_30 <- rast(file.path(stock_map_dir, "stock_map_5_years_CFN_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_5_years_CFN_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_5_years_CFN_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_5_years_CFN_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_5_years_CFN_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_5_years_CFN_0_30)

stock_map_adj_5_years_COA_0_30 <- rast(file.path(stock_map_dir, "stock_map_5_years_COA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_5_years_COA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_5_years_COA_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_5_years_COA_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_5_years_COA_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_5_years_COA_0_30)

stock_map_adj_5_years_CON_0_30 <- rast(file.path(stock_map_dir, "stock_map_5_years_CON_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_5_years_CON_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_5_years_CON_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_5_years_CON_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_5_years_CON_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_5_years_CON_0_30)

stock_map_adj_5_years_PFA_0_30 <- rast(file.path(stock_map_dir, "stock_map_5_years_PFA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_5_years_PFA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_5_years_PFA_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_5_years_PFA_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_5_years_PFA_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_5_years_PFA_0_30)

stock_map_adj_5_years_PFN_0_30 <- rast(file.path(stock_map_dir, "stock_map_5_years_PFN_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_5_years_PFN_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_5_years_PFN_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_5_years_PFN_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_5_years_PFN_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_5_years_PFN_0_30)

stock_map_adj_5_years_POA_0_30 <- rast(file.path(stock_map_dir, "stock_map_5_years_POA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_5_years_POA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_5_years_POA_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_5_years_POA_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_5_years_POA_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_5_years_POA_0_30)

stock_map_adj_5_years_PON_0_30 <- rast(file.path(stock_map_dir, "stock_map_5_years_PON_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_5_years_PON_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_5_years_PON_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_5_years_PON_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_5_years_PON_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_5_years_PON_0_30)


# 15_years
stock_map_adj_15_years_CFA_0_30 <- rast(file.path(stock_map_dir, "stock_map_15_years_CFA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_15_years_CFA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_15_years_CFA_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_15_years_CFA_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_15_years_CFA_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_15_years_CFA_0_30)

stock_map_adj_15_years_CFN_0_30 <- rast(file.path(stock_map_dir, "stock_map_15_years_CFN_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_15_years_CFN_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_15_years_CFN_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_15_years_CFN_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_15_years_CFN_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_15_years_CFN_0_30)

stock_map_adj_15_years_COA_0_30 <- rast(file.path(stock_map_dir, "stock_map_15_years_COA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_15_years_COA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_15_years_COA_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_15_years_COA_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_15_years_COA_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_15_years_COA_0_30)

stock_map_adj_15_years_CON_0_30 <- rast(file.path(stock_map_dir, "stock_map_15_years_CON_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_15_years_CON_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_15_years_CON_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_15_years_CON_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_15_years_CON_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_15_years_CON_0_30)

stock_map_adj_15_years_PFA_0_30 <- rast(file.path(stock_map_dir, "stock_map_15_years_PFA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_15_years_PFA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_15_years_PFA_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_15_years_PFA_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_15_years_PFA_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_15_years_PFA_0_30)

stock_map_adj_15_years_PFN_0_30 <- rast(file.path(stock_map_dir, "stock_map_15_years_PFN_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_15_years_PFN_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_15_years_PFN_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_15_years_PFN_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_15_years_PFN_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_15_years_PFN_0_30)

stock_map_adj_15_years_POA_0_30 <- rast(file.path(stock_map_dir, "stock_map_15_years_POA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_15_years_POA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_15_years_POA_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_15_years_POA_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_15_years_POA_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_15_years_POA_0_30)

stock_map_adj_15_years_PON_0_30 <- rast(file.path(stock_map_dir, "stock_map_15_years_PON_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_15_years_PON_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_15_years_PON_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_15_years_PON_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_15_years_PON_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_15_years_PON_0_30)


# 30_years
stock_map_adj_30_years_CFA_0_30 <- rast(file.path(stock_map_dir, "stock_map_30_years_CFA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_30_years_CFA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_30_years_CFA_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_30_years_CFA_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_30_years_CFA_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_30_years_CFA_0_30)

stock_map_adj_30_years_CFN_0_30 <- rast(file.path(stock_map_dir, "stock_map_30_years_CFN_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_30_years_CFN_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_30_years_CFN_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_30_years_CFN_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_30_years_CFN_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_30_years_CFN_0_30)

stock_map_adj_30_years_COA_0_30 <- rast(file.path(stock_map_dir, "stock_map_30_years_COA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_30_years_COA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_30_years_COA_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_30_years_COA_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_30_years_COA_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_30_years_COA_0_30)

stock_map_adj_30_years_CON_0_30 <- rast(file.path(stock_map_dir, "stock_map_30_years_CON_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_30_years_CON_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_30_years_CON_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_30_years_CON_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_30_years_CON_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_30_years_CON_0_30)

stock_map_adj_30_years_PFA_0_30 <- rast(file.path(stock_map_dir, "stock_map_30_years_PFA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_30_years_PFA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_30_years_PFA_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_30_years_PFA_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_30_years_PFA_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_30_years_PFA_0_30)

stock_map_adj_30_years_PFN_0_30 <- rast(file.path(stock_map_dir, "stock_map_30_years_PFN_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_30_years_PFN_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_30_years_PFN_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_30_years_PFN_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_30_years_PFN_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_30_years_PFN_0_30)

stock_map_adj_30_years_POA_0_30 <- rast(file.path(stock_map_dir, "stock_map_30_years_POA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_30_years_POA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_30_years_POA_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_30_years_POA_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_30_years_POA_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_30_years_POA_0_30)

stock_map_adj_30_years_PON_0_30 <- rast(file.path(stock_map_dir, "stock_map_30_years_PON_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_30_years_PON_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_30_years_PON_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_30_years_PON_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_30_years_PON_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_30_years_PON_0_30)

# 50_years
stock_map_adj_50_years_CFA_0_30 <- rast(file.path(stock_map_dir, "stock_map_50_years_CFA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_50_years_CFA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_50_years_CFA_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_50_years_CFA_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_50_years_CFA_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_50_years_CFA_0_30)

stock_map_adj_50_years_CFN_0_30 <- rast(file.path(stock_map_dir, "stock_map_50_years_CFN_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_50_years_CFN_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_50_years_CFN_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_50_years_CFN_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_50_years_CFN_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_50_years_CFN_0_30)

stock_map_adj_50_years_COA_0_30 <- rast(file.path(stock_map_dir, "stock_map_50_years_COA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_50_years_COA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_50_years_COA_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_50_years_COA_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_50_years_COA_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_50_years_COA_0_30)

stock_map_adj_50_years_CON_0_30 <- rast(file.path(stock_map_dir, "stock_map_50_years_CON_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_50_years_CON_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_50_years_CON_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_50_years_CON_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_50_years_CON_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_50_years_CON_0_30)

stock_map_adj_50_years_PFA_0_30 <- rast(file.path(stock_map_dir, "stock_map_50_years_PFA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_50_years_PFA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_50_years_PFA_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_50_years_PFA_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_50_years_PFA_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_50_years_PFA_0_30)

stock_map_adj_50_years_PFN_0_30 <- rast(file.path(stock_map_dir, "stock_map_50_years_PFN_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_50_years_PFN_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_50_years_PFN_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_50_years_PFN_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_50_years_PFN_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_50_years_PFN_0_30)

stock_map_adj_50_years_POA_0_30 <- rast(file.path(stock_map_dir, "stock_map_50_years_POA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_50_years_POA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_50_years_POA_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_50_years_POA_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_50_years_POA_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_50_years_POA_0_30)

stock_map_adj_50_years_PON_0_30 <- rast(file.path(stock_map_dir, "stock_map_50_years_PON_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_50_years_PON_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_50_years_PON_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_50_years_PON_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_50_years_PON_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_50_years_PON_0_30)


# 100_years
stock_map_adj_100_years_CFA_0_30 <- rast(file.path(stock_map_dir, "stock_map_100_years_CFA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_100_years_CFA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_100_years_CFA_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_100_years_CFA_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_100_years_CFA_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_100_years_CFA_0_30)

stock_map_adj_100_years_CFN_0_30 <- rast(file.path(stock_map_dir, "stock_map_100_years_CFN_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_100_years_CFN_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_100_years_CFN_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_100_years_CFN_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_100_years_CFN_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_100_years_CFN_0_30)

stock_map_adj_100_years_COA_0_30 <- rast(file.path(stock_map_dir, "stock_map_100_years_COA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_100_years_COA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_100_years_COA_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_100_years_COA_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_100_years_COA_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_100_years_COA_0_30)

stock_map_adj_100_years_CON_0_30 <- rast(file.path(stock_map_dir, "stock_map_100_years_CON_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_100_years_CON_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_100_years_CON_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_100_years_CON_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_100_years_CON_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_100_years_CON_0_30)

stock_map_adj_100_years_PFA_0_30 <- rast(file.path(stock_map_dir, "stock_map_100_years_PFA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_100_years_PFA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_100_years_PFA_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_100_years_PFA_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_100_years_PFA_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_100_years_PFA_0_30)

stock_map_adj_100_years_PFN_0_30 <- rast(file.path(stock_map_dir, "stock_map_100_years_PFN_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_100_years_PFN_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_100_years_PFN_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_100_years_PFN_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_100_years_PFN_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_100_years_PFN_0_30)

stock_map_adj_100_years_POA_0_30 <- rast(file.path(stock_map_dir, "stock_map_100_years_POA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_100_years_POA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_100_years_POA_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_100_years_POA_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_100_years_POA_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_100_years_POA_0_30)

stock_map_adj_100_years_PON_0_30 <- rast(file.path(stock_map_dir, "stock_map_100_years_PON_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_100_years_PON_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_100_years_PON_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_100_years_PON_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_100_years_PON_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_100_years_PON_0_30)


# 200_years
stock_map_adj_200_years_CFA_0_30 <- rast(file.path(stock_map_dir, "stock_map_200_years_CFA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_200_years_CFA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_200_years_CFA_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_200_years_CFA_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_200_years_CFA_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_200_years_CFA_0_30)

stock_map_adj_200_years_CFN_0_30 <- rast(file.path(stock_map_dir, "stock_map_200_years_CFN_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_200_years_CFN_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_200_years_CFN_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_200_years_CFN_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_200_years_CFN_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_200_years_CFN_0_30)

stock_map_adj_200_years_COA_0_30 <- rast(file.path(stock_map_dir, "stock_map_200_years_COA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_200_years_COA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_200_years_COA_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_200_years_COA_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_200_years_COA_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_200_years_COA_0_30)

stock_map_adj_200_years_CON_0_30 <- rast(file.path(stock_map_dir, "stock_map_200_years_CON_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_200_years_CON_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_200_years_CON_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_200_years_CON_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_200_years_CON_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_200_years_CON_0_30)

stock_map_adj_200_years_PFA_0_30 <- rast(file.path(stock_map_dir, "stock_map_200_years_PFA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_200_years_PFA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_200_years_PFA_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_200_years_PFA_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_200_years_PFA_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_200_years_PFA_0_30)

stock_map_adj_200_years_PFN_0_30 <- rast(file.path(stock_map_dir, "stock_map_200_years_PFN_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_200_years_PFN_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_200_years_PFN_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_200_years_PFN_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_200_years_PFN_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_200_years_PFN_0_30)

stock_map_adj_200_years_POA_0_30 <- rast(file.path(stock_map_dir, "stock_map_200_years_POA_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_200_years_POA_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_200_years_POA_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_200_years_POA_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_200_years_POA_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_200_years_POA_0_30)

stock_map_adj_200_years_PON_0_30 <- rast(file.path(stock_map_dir, "stock_map_200_years_PON_0_5.tif")) * cfc_0_5 +
  rast(file.path(stock_map_dir, "stock_map_200_years_PON_5_15.tif")) * cfc_5_15 +
  rast(file.path(stock_map_dir, "stock_map_200_years_PON_15_30.tif")) * cfc_15_30
writeRaster(stock_map_adj_200_years_PON_0_30,
            file = file.path(stock_map_dir, "stock_map_adj_200_years_PON_0_30.tif"),
            overwrite = T)
rm(stock_map_adj_200_years_PON_0_30)