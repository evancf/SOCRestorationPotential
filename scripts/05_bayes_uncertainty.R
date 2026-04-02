# =============================================================================
# 05_bayes_uncertainty_parallel.R
#
# PURPOSE
# -------
# Parallelised version of 05_bayes_uncertainty.R.  Runs the SOC uncertainty
# pipeline for all posterior samples across both climate scenarios
# (contemporary and ssp370_2050) using mclapply (fork-based, Mac/Linux).
#
# PRE-REQUISITES
# --------------
# bayes_output.RData must exist (produced by 03_bayes_analysis.R).
# Source rasters (SoilGrids, WorldClim, CMIP6) must be present.
# Covariate rasters are aggregated and cached at the script's own agg_fact
# in Section 2 before workers launch — no manual pre-caching required.
#
# CONFIGURATION
# -------------
# Set n_workers and agg_fact at the top.  agg_fact is independent of the
# value used in 04_bayes_results.R — each script caches at its own resolution.
# =============================================================================


# ── Packages ──────────────────────────────────────────────────────────────────
# ipak() installs and loads packages in one call
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(c("tidyverse", "terra", "parallel"))


# =============================================================================
# CONFIGURATION
# =============================================================================

n_workers <- 10    # parallel workers; keep n_workers × ~2.5 GB < available RAM
agg_fact  <- 25   # independent of agg_fact in 04_bayes_results.R

climate_scenarios <- c("contemporary", "ssp370_2050")


# =============================================================================
# SECTION 1 — SHARED SETUP (runs once in main process)
# =============================================================================
# sc_fun converts a covariate from its natural scale to the centred/scaled
# form used in the model.  See 04_bayes_results.R for full documentation.

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

load('./data_output/bayes_output.RData')
mod_samp_df <- do.call(rbind, mod_samp)


# =============================================================================
# SECTION 2 — PRE-CACHE COVARIATE RASTERS AT THIS SCRIPT'S agg_fact
# =============================================================================
# Builds and caches all aggregated covariate rasters at the agg_fact set
# above.  Safe to re-run: existing cache files are read without reprocessing.
# This runs once in the main process before workers launch, so there is no
# risk of multiple workers writing the same file simultaneously.
#
# load_scaled_agg mirrors the version in 04_bayes_results.R exactly, so
# both scripts produce identically-named cache files at their own agg_fact.

load_scaled_agg <- function(varname,
                            interval,
                            sc_df,
                            agg_fact,
                            indir  = "./data_inputs/raster/via_gee",
                            outdir = "./data_inputs/raster/via_gee",
                            transforms = list(
                              mat = function(x) x / 10
                            ),
                            snap_to = NULL) {
  special_files <- list(
    map        = "tap.tif",
    mat        = "amt.tif",
    irrigation = "gmia_plus_0.1.tif"
  )
  if (varname %in% names(special_files)) {
    infile_name  <- special_files[[varname]]
    outfile_name <- sprintf("%s_agg%d.tif", varname, agg_fact)
  } else {
    infile_name  <- sprintf("%s_%scm_mean.tif", varname, interval)
    outfile_name <- sprintf("%s_%scm_mean_agg%d.tif", varname, interval, agg_fact)
  }
  infile  <- file.path(indir,  infile_name)
  outfile <- file.path(outdir, outfile_name)
  
  if (file.exists(outfile)) {
    message("Reading cached: ", outfile)
    return(rast(outfile))
  }
  
  message("Caching: ", outfile)
  r0 <- rast(infile)
  if (!is.null(transforms[[varname]])) r0 <- transforms[[varname]](r0)
  if (isTRUE(sc_df[varname, "log"]))   r0 <- log(r0)
  cen <- sc_df[varname, "center"]
  sc  <- sc_df[varname, "scale"]
  r1  <- (r0 - cen) / sc
  r2  <- aggregate(r1, fact = agg_fact, fun = mean, na.rm = TRUE)
  if (!is.null(snap_to)) r2 <- resample(r2, snap_to, method = "bilinear")
  writeRaster(r2, filename = outfile, overwrite = TRUE)
  return(r2)
}

# Pre-cache coarse fraction correction rasters (cfvo, SoilGrids)
load_cfc_agg <- function(interval, agg_fact,
                         indir  = "./data_inputs/raster/via_gee",
                         outdir = "./data_inputs/raster/via_gee") {
  infile  <- file.path(indir,  sprintf("cfvo_%scm_mean.tif",   interval))
  outfile <- file.path(outdir, sprintf("cfc_%s_agg%d.tif", interval, agg_fact))
  if (file.exists(outfile)) {
    message("Reading cached: ", outfile)
    return(rast(outfile))
  }
  message("Caching: ", outfile)
  r0    <- rast(infile)
  r1    <- aggregate(r0, fact = agg_fact, fun = mean, na.rm = TRUE)
  r_cfc <- 1 - (r1 / 1000)
  writeRaster(r_cfc, filename = outfile, overwrite = TRUE)
  return(r_cfc)
}

sg_intervals <- c("0_5", "5_15", "15_30", "30_60", "60_100", "100_200")
soilgrids_dir <- "./data_inputs/raster/via_gee"

message("Pre-caching covariate rasters at agg_fact = ", agg_fact, " ...")

# Coarse fraction correction factors (climate-independent)
for (iv in sg_intervals) {
  load_cfc_agg(iv, agg_fact)
}

# SoilGrids variables (climate-independent, all depth intervals)
for (iv in sg_intervals) {
  load_scaled_agg("clay",     iv, sc_df, agg_fact,
                  indir = soilgrids_dir, outdir = soilgrids_dir)
  load_scaled_agg("nitrogen", iv, sc_df, agg_fact,
                  indir = soilgrids_dir, outdir = soilgrids_dir)
}

# Climate-dependent variables — one entry per scenario
for (scen in climate_scenarios) {
  clim_dir <- if (scen == "contemporary") {
    soilgrids_dir
  } else {
    file.path("./data_inputs/raster/cmip6", scen)
  }
  
  # Use clay as the snap_to reference so climate rasters align to SoilGrids grid
  grid_ref <- rast(file.path(soilgrids_dir,
                             sprintf("clay_0_5cm_mean_agg%d.tif", agg_fact)))
  
  load_scaled_agg("map",        "0_5", sc_df, agg_fact,
                  indir = clim_dir, outdir = clim_dir, snap_to = grid_ref)
  load_scaled_agg("mat",        "0_5", sc_df, agg_fact,
                  indir = clim_dir, outdir = clim_dir, snap_to = grid_ref)
  load_scaled_agg("irrigation", "0_5", sc_df, agg_fact,
                  indir = soilgrids_dir, outdir = soilgrids_dir)
}

message("Pre-caching complete.")


# =============================================================================
# SECTION 3 — BUILD JOB LIST
# =============================================================================
# All (b_samp, climate_scenario) pairs that don't yet have output.
# Restartable: already-completed stock_map_<b_samp> directories are skipped.

all_samples <- seq(1, nrow(mod_samp_df), by = 10)

jobs <- expand.grid(
  b_samp           = all_samples,
  climate_scenario = climate_scenarios,
  stringsAsFactors = FALSE
)

jobs <- jobs[!mapply(function(b, scen) {
  dir.exists(file.path("./data_output/uncertainty", scen,
                       paste0("stock_map_", b)))
}, jobs$b_samp, jobs$climate_scenario), ]

message(sprintf("%d jobs remaining (%d samples x %d scenarios)",
                nrow(jobs), length(all_samples), length(climate_scenarios)))


# =============================================================================
# SECTION 4 — WORKER FUNCTION
# =============================================================================

run_one_job <- function(b_samp, climate_scenario, agg_fact,
                        mod_samp_df, sc_df) {
  
  ipak(c("terra", "tidyverse", "pracma"))
  
  # Prevent terra spawning sub-threads inside a forked worker
  terraOptions(parallel = FALSE)
  
  # ── Paths ──────────────────────────────────────────────────────────────────
  soilgrids_dir <- "./data_inputs/raster/via_gee"
  climate_raster_dir <- if (climate_scenario == "contemporary") {
    soilgrids_dir
  } else {
    file.path("./data_inputs/raster/cmip6", climate_scenario)
  }
  uncertainty_dir <- file.path("./data_output/uncertainty", climate_scenario)
  dir.create(uncertainty_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ── Skip if already done ───────────────────────────────────────────────────
  if (dir.exists(file.path(uncertainty_dir, paste0("stock_map_", b_samp)))) {
    return(invisible(NULL))
  }
  
  message(sprintf("[%s / b_samp=%d] Starting at %s",
                  climate_scenario, b_samp, format(Sys.time(), "%H:%M:%S")))
  
  # ── Read-only helper: load pre-cached aggregated raster ───────────────────
  # Never writes — all caching is done in the main process / prior run.
  load_scaled_agg_cached <- function(varname, interval,
                                     indir = climate_raster_dir) {
    special_files <- list(map = "tap.tif", mat = "amt.tif",
                          irrigation = "gmia_plus_0.1.tif")
    outfile_name <- if (varname %in% names(special_files)) {
      sprintf("%s_agg%d.tif", varname, agg_fact)
    } else {
      sprintf("%s_%scm_mean_agg%d.tif", varname, interval, agg_fact)
    }
    rast(file.path(indir, outfile_name))
  }
  
  # ── Coarse fraction correction rasters ────────────────────────────────────
  sg_intervals <- c("0_5", "5_15", "15_30", "30_60", "60_100", "100_200")
  cfc_list <- setNames(
    lapply(sg_intervals, function(iv) {
      rast(file.path(soilgrids_dir,
                     sprintf("cfc_%s_agg%d.tif", iv, agg_fact)))
    }),
    paste0("cfc_", sg_intervals)
  )
  
  sg_interval_lo <- c(0.1, 5, 15, 30, 60, 100)
  sg_interval_hi <- c(5, 15, 30, 60, 100, 200)
  time_points    <- c(5, 15, 30, 50, 100, 200)
  redo_rasts     <- TRUE
  
  # ── Posterior sample ───────────────────────────────────────────────────────
  bsm <- mod_samp_df[b_samp, ]
  
  # ── Output directories ─────────────────────────────────────────────────────
  dir.create(file.path(uncertainty_dir, paste0("rast_component_", b_samp)),
             showWarnings = FALSE)
  dir.create(file.path(uncertainty_dir, paste0("stock_map_", b_samp)),
             showWarnings = FALSE)
  
  # ── Climate + irrigation rasters (depth-invariant) ────────────────────────
  grid_ref   <- load_scaled_agg_cached("clay",       "0_5",
                                       indir = soilgrids_dir)
  map        <- load_scaled_agg_cached("map",        "0_5")
  mat        <- load_scaled_agg_cached("mat",        "0_5")
  irrigation <- load_scaled_agg_cached("irrigation", "0_5",
                                       indir = soilgrids_dir)
  
  # ── SECTION A: Build raster components ──────────────────────────────────
  if(redo_rasts == T){
    # Loop over these intervals to save the rast components at each interval
    for(i in 1:length(sg_intervals)){
      
      print(sg_intervals[i])
      
      # Make filenames
      
      # Init
      rast_component_init_crop_filename <- file.path(uncertainty_dir, paste0("rast_component_", b_samp), paste0("rast_component_init_crop_", sg_intervals[i], ".tif"))
      
      rast_component_init_pasture_filename <- file.path(uncertainty_dir, paste0("rast_component_", b_samp), paste0("rast_component_init_pasture_", sg_intervals[i], ".tif"))
      
      # R val "CFA","CFN","COA","CON","PFA","PFN","POA","PON"
      rast_component_r_CFA_filename <- file.path(uncertainty_dir, paste0("rast_component_", b_samp), paste0("rast_component_r_CFA_", sg_intervals[i], ".tif"))
      
      rast_component_r_CFN_filename <- file.path(uncertainty_dir, paste0("rast_component_", b_samp), paste0("rast_component_r_CFN_", sg_intervals[i], ".tif"))
      
      rast_component_r_COA_filename <- file.path(uncertainty_dir, paste0("rast_component_", b_samp), paste0("rast_component_r_COA_", sg_intervals[i], ".tif"))
      
      rast_component_r_CON_filename <- file.path(uncertainty_dir, paste0("rast_component_", b_samp), paste0("rast_component_r_CON_", sg_intervals[i], ".tif"))
      
      rast_component_r_PFA_filename <- file.path(uncertainty_dir, paste0("rast_component_", b_samp), paste0("rast_component_r_PFA_", sg_intervals[i], ".tif"))
      
      rast_component_r_PFN_filename <- file.path(uncertainty_dir, paste0("rast_component_", b_samp), paste0("rast_component_r_PFN_", sg_intervals[i], ".tif"))
      
      rast_component_r_POA_filename <- file.path(uncertainty_dir, paste0("rast_component_", b_samp), paste0("rast_component_r_POA_", sg_intervals[i], ".tif"))
      
      rast_component_r_PON_filename <- file.path(uncertainty_dir, paste0("rast_component_", b_samp), paste0("rast_component_r_PON_", sg_intervals[i], ".tif"))
      
      # Max
      rast_component_max_filename <- file.path(uncertainty_dir, paste0("rast_component_", b_samp), paste0("rast_component_max_", sg_intervals[i], ".tif"))
      
      
      
      
      # Get the appropriate clay and nitrogen rasters for this interval
      clay <- load_scaled_agg_cached("clay", sg_intervals[i], indir = soilgrids_dir)
      
      nitrogen <- load_scaled_agg_cached("nitrogen", sg_intervals[i], indir = soilgrids_dir)
      
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
      
      # rast_component_r_CFA <- bsm["beta_r_intercept[1]"] + # This index represents the final type (1-8)
      #   rast_subcomponent_r_C
      # writeRaster(rast_component_r_CFA,
      #             overwrite = T,
      #             filename = rast_component_r_CFA_filename)
      # rm(rast_component_r_CFA)
      # gc()
      
      rast_component_r_CFN <- bsm["beta_r_intercept[2]"] + # This index represents the final type (1-8)
        rast_subcomponent_r_C
      writeRaster(rast_component_r_CFN,
                  overwrite = T,
                  filename = rast_component_r_CFN_filename)
      rm(rast_component_r_CFN)
      gc()
      
      # rast_component_r_COA <- bsm["beta_r_intercept[3]"] + # This index represents the final type (1-8)
      #   rast_subcomponent_r_C
      # writeRaster(rast_component_r_COA,
      #             overwrite = T,
      #             filename = rast_component_r_COA_filename)
      # rm(rast_component_r_COA)
      # gc()
      
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
      
      # rast_component_r_PFA <- bsm["beta_r_intercept[5]"] + # This index represents the final type (1-8)
      #   rast_subcomponent_r_P
      # writeRaster(rast_component_r_PFA,
      #             overwrite = T,
      #             filename = rast_component_r_PFA_filename)
      # rm(rast_component_r_PFA)
      # gc()
      
      rast_component_r_PFN <- bsm["beta_r_intercept[6]"] + # This index represents the final type (1-8)
        rast_subcomponent_r_P
      writeRaster(rast_component_r_PFN,
                  overwrite = T,
                  filename = rast_component_r_PFN_filename)
      rm(rast_component_r_PFN)
      gc()
      
      
      # rast_component_r_POA <- bsm["beta_r_intercept[7]"] + # This index represents the final type (1-8)
      #   rast_subcomponent_r_P
      # writeRaster(rast_component_r_POA,
      #             overwrite = T,
      #             filename = rast_component_r_POA_filename)
      # rm(rast_component_r_POA)
      # gc()
      
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
  } # end redo_rasts
  
  # ── SECTION B: Depth integration ─────────────────────────────────────────
  
  
  time1 <- Sys.time()
  if(redo_rasts){
    for(yr in time_points){
      print(paste("starting on year", yr, "at", Sys.time()))
      for(i in 1:length(sg_intervals)){
        print(paste("--",sg_intervals[i], "at", Sys.time()))
        
        # Get rast component filenames
        # Init
        rast_component_init_crop_filename <- file.path(uncertainty_dir, paste0("rast_component_", b_samp), paste0("rast_component_init_crop_", sg_intervals[i], ".tif"))
        
        rast_component_init_pasture_filename <- file.path(uncertainty_dir, paste0("rast_component_", b_samp), paste0("rast_component_init_pasture_", sg_intervals[i], ".tif"))
        
        # R val "CFA","CFN","COA","CON","PFA","PFN","POA","PON"
        rast_component_r_CFA_filename <- file.path(uncertainty_dir, paste0("rast_component_", b_samp), paste0("rast_component_r_CFA_", sg_intervals[i], ".tif"))
        
        rast_component_r_CFN_filename <- file.path(uncertainty_dir, paste0("rast_component_", b_samp), paste0("rast_component_r_CFN_", sg_intervals[i], ".tif"))
        
        rast_component_r_COA_filename <- file.path(uncertainty_dir, paste0("rast_component_", b_samp), paste0("rast_component_r_COA_", sg_intervals[i], ".tif"))
        
        rast_component_r_CON_filename <- file.path(uncertainty_dir, paste0("rast_component_", b_samp), paste0("rast_component_r_CON_", sg_intervals[i], ".tif"))
        
        rast_component_r_PFA_filename <- file.path(uncertainty_dir, paste0("rast_component_", b_samp), paste0("rast_component_r_PFA_", sg_intervals[i], ".tif"))
        
        rast_component_r_PFN_filename <- file.path(uncertainty_dir, paste0("rast_component_", b_samp), paste0("rast_component_r_PFN_", sg_intervals[i], ".tif"))
        
        rast_component_r_POA_filename <- file.path(uncertainty_dir, paste0("rast_component_", b_samp), paste0("rast_component_r_POA_", sg_intervals[i], ".tif"))
        
        rast_component_r_PON_filename <- file.path(uncertainty_dir, paste0("rast_component_", b_samp), paste0("rast_component_r_PON_", sg_intervals[i], ".tif"))
        
        # Max
        rast_component_max_filename <- file.path(uncertainty_dir, paste0("rast_component_", b_samp), paste0("rast_component_max_", sg_intervals[i], ".tif"))
        
        
        # Don't comment out the active ones so that indexing works later on
        final_type_filenames <- c(
          rast_component_r_CFA_filename,
          rast_component_r_CFN_filename,
          rast_component_r_COA_filename,
          rast_component_r_CON_filename,
          rast_component_r_PFA_filename,
          rast_component_r_PFN_filename,
          rast_component_r_POA_filename,
          rast_component_r_PON_filename)
        
        
        
        
        
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
        for(j in c(2,4,6,8)){
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
                        gsub("map_r", paste0("map_", yr, "_years"), .),
                      overwrite = T)
          
        }
      }
    }
  }
  time2 <- Sys.time()
  
  
  # ── SECTION C: Coarse fraction adjustment ────────────────────────────────
  # Stock summed across depth and adjusted by coarse fraction --------------------
  
  
  
  
  # 0 to 100 cm version
  
  # crop (which is independent of year)
  stock_map_adj_init_crop_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_init_crop_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_init_crop_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_init_crop_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_init_crop_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_init_crop_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_init_crop_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_init_crop_0_100.tif"),
              overwrite = T)
  
  # pasture (which is independent of year)
  stock_map_adj_init_pasture_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_init_pasture_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_init_pasture_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_init_pasture_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_init_pasture_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_init_pasture_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_init_pasture_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_init_pasture_0_100.tif"),
              overwrite = T)
  
  # max (which is independent of year)
  stock_map_adj_max_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_max_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_max_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_max_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_max_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_max_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_max_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_max_0_100.tif"),
              overwrite = T)
  
  
  
  # Now the other types for each year "CFA","CFN","COA","CON","PFA","PFN","POA","PON"
  
  # 5_years
  # stock_map_adj_5_years_CFA_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_CFA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_CFA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_CFA_15_30.tif")) * cfc_list[["cfc_15_30"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_CFA_30_60.tif")) * cfc_list[["cfc_30_60"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_CFA_60_100.tif")) * cfc_list[["cfc_60_100"]]
  # writeRaster(stock_map_adj_5_years_CFA_0_100,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_5_years_CFA_0_100.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_5_years_CFA_0_100)
  
  stock_map_adj_5_years_CFN_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_CFN_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_CFN_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_CFN_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_CFN_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_CFN_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_5_years_CFN_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_5_years_CFN_0_100.tif"),
              overwrite = T)
  rm(stock_map_adj_5_years_CFN_0_100)
  
  # stock_map_adj_5_years_COA_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_COA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_COA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_COA_15_30.tif")) * cfc_list[["cfc_15_30"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_COA_30_60.tif")) * cfc_list[["cfc_30_60"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_COA_60_100.tif")) * cfc_list[["cfc_60_100"]]
  # writeRaster(stock_map_adj_5_years_COA_0_100,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_5_years_COA_0_100.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_5_years_COA_0_100)
  
  stock_map_adj_5_years_CON_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_CON_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_CON_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_CON_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_CON_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_CON_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_5_years_CON_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_5_years_CON_0_100.tif"),
              overwrite = T)
  rm(stock_map_adj_5_years_CON_0_100)
  
  # stock_map_adj_5_years_PFA_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_PFA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_PFA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_PFA_15_30.tif")) * cfc_list[["cfc_15_30"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_PFA_30_60.tif")) * cfc_list[["cfc_30_60"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_PFA_60_100.tif")) * cfc_list[["cfc_60_100"]]
  # writeRaster(stock_map_adj_5_years_PFA_0_100,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_5_years_PFA_0_100.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_5_years_PFA_0_100)
  
  stock_map_adj_5_years_PFN_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_PFN_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_PFN_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_PFN_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_PFN_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_PFN_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_5_years_PFN_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_5_years_PFN_0_100.tif"),
              overwrite = T)
  rm(stock_map_adj_5_years_PFN_0_100)
  
  # stock_map_adj_5_years_POA_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_POA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_POA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_POA_15_30.tif")) * cfc_list[["cfc_15_30"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_POA_30_60.tif")) * cfc_list[["cfc_30_60"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_POA_60_100.tif")) * cfc_list[["cfc_60_100"]]
  # writeRaster(stock_map_adj_5_years_POA_0_100,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_5_years_POA_0_100.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_5_years_POA_0_100)
  
  stock_map_adj_5_years_PON_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_PON_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_PON_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_PON_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_PON_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_PON_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_5_years_PON_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_5_years_PON_0_100.tif"),
              overwrite = T)
  rm(stock_map_adj_5_years_PON_0_100)
  
  
  # 15_years
  # stock_map_adj_15_years_CFA_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_CFA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_CFA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_CFA_15_30.tif")) * cfc_list[["cfc_15_30"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_CFA_30_60.tif")) * cfc_list[["cfc_30_60"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_CFA_60_100.tif")) * cfc_list[["cfc_60_100"]]
  # writeRaster(stock_map_adj_15_years_CFA_0_100,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_15_years_CFA_0_100.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_15_years_CFA_0_100)
  
  stock_map_adj_15_years_CFN_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_CFN_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_CFN_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_CFN_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_CFN_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_CFN_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_15_years_CFN_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_15_years_CFN_0_100.tif"),
              overwrite = T)
  rm(stock_map_adj_15_years_CFN_0_100)
  
  # stock_map_adj_15_years_COA_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_COA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_COA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_COA_15_30.tif")) * cfc_list[["cfc_15_30"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_COA_30_60.tif")) * cfc_list[["cfc_30_60"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_COA_60_100.tif")) * cfc_list[["cfc_60_100"]]
  # writeRaster(stock_map_adj_15_years_COA_0_100,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_15_years_COA_0_100.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_15_years_COA_0_100)
  
  stock_map_adj_15_years_CON_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_CON_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_CON_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_CON_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_CON_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_CON_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_15_years_CON_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_15_years_CON_0_100.tif"),
              overwrite = T)
  rm(stock_map_adj_15_years_CON_0_100)
  
  # stock_map_adj_15_years_PFA_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_PFA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_PFA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_PFA_15_30.tif")) * cfc_list[["cfc_15_30"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_PFA_30_60.tif")) * cfc_list[["cfc_30_60"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_PFA_60_100.tif")) * cfc_list[["cfc_60_100"]]
  # writeRaster(stock_map_adj_15_years_PFA_0_100,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_15_years_PFA_0_100.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_15_years_PFA_0_100)
  
  stock_map_adj_15_years_PFN_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_PFN_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_PFN_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_PFN_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_PFN_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_PFN_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_15_years_PFN_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_15_years_PFN_0_100.tif"),
              overwrite = T)
  rm(stock_map_adj_15_years_PFN_0_100)
  
  # stock_map_adj_15_years_POA_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_POA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_POA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_POA_15_30.tif")) * cfc_list[["cfc_15_30"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_POA_30_60.tif")) * cfc_list[["cfc_30_60"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_POA_60_100.tif")) * cfc_list[["cfc_60_100"]]
  # writeRaster(stock_map_adj_15_years_POA_0_100,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_15_years_POA_0_100.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_15_years_POA_0_100)
  
  stock_map_adj_15_years_PON_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_PON_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_PON_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_PON_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_PON_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_PON_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_15_years_PON_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_15_years_PON_0_100.tif"),
              overwrite = T)
  rm(stock_map_adj_15_years_PON_0_100)
  
  
  # 30_years
  # stock_map_adj_30_years_CFA_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_CFA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_CFA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_CFA_15_30.tif")) * cfc_list[["cfc_15_30"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_CFA_30_60.tif")) * cfc_list[["cfc_30_60"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_CFA_60_100.tif")) * cfc_list[["cfc_60_100"]]
  # writeRaster(stock_map_adj_30_years_CFA_0_100,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_30_years_CFA_0_100.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_30_years_CFA_0_100)
  
  stock_map_adj_30_years_CFN_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_CFN_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_CFN_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_CFN_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_CFN_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_CFN_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_30_years_CFN_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_30_years_CFN_0_100.tif"),
              overwrite = T)
  rm(stock_map_adj_30_years_CFN_0_100)
  
  # stock_map_adj_30_years_COA_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_COA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_COA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_COA_15_30.tif")) * cfc_list[["cfc_15_30"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_COA_30_60.tif")) * cfc_list[["cfc_30_60"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_COA_60_100.tif")) * cfc_list[["cfc_60_100"]]
  # writeRaster(stock_map_adj_30_years_COA_0_100,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_30_years_COA_0_100.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_30_years_COA_0_100)
  
  stock_map_adj_30_years_CON_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_CON_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_CON_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_CON_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_CON_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_CON_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_30_years_CON_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_30_years_CON_0_100.tif"),
              overwrite = T)
  rm(stock_map_adj_30_years_CON_0_100)
  
  # stock_map_adj_30_years_PFA_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_PFA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_PFA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_PFA_15_30.tif")) * cfc_list[["cfc_15_30"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_PFA_30_60.tif")) * cfc_list[["cfc_30_60"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_PFA_60_100.tif")) * cfc_list[["cfc_60_100"]]
  # writeRaster(stock_map_adj_30_years_PFA_0_100,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_30_years_PFA_0_100.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_30_years_PFA_0_100)
  
  stock_map_adj_30_years_PFN_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_PFN_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_PFN_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_PFN_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_PFN_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_PFN_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_30_years_PFN_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_30_years_PFN_0_100.tif"),
              overwrite = T)
  rm(stock_map_adj_30_years_PFN_0_100)
  
  # stock_map_adj_30_years_POA_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_POA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_POA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_POA_15_30.tif")) * cfc_list[["cfc_15_30"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_POA_30_60.tif")) * cfc_list[["cfc_30_60"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_POA_60_100.tif")) * cfc_list[["cfc_60_100"]]
  # writeRaster(stock_map_adj_30_years_POA_0_100,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_30_years_POA_0_100.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_30_years_POA_0_100)
  
  stock_map_adj_30_years_PON_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_PON_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_PON_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_PON_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_PON_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_PON_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_30_years_PON_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_30_years_PON_0_100.tif"),
              overwrite = T)
  rm(stock_map_adj_30_years_PON_0_100)
  
  # 50_years
  # stock_map_adj_50_years_CFA_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_CFA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_CFA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_CFA_15_30.tif")) * cfc_list[["cfc_15_30"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_CFA_30_60.tif")) * cfc_list[["cfc_30_60"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_CFA_60_100.tif")) * cfc_list[["cfc_60_100"]]
  # writeRaster(stock_map_adj_50_years_CFA_0_100,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_50_years_CFA_0_100.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_50_years_CFA_0_100)
  
  stock_map_adj_50_years_CFN_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_CFN_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_CFN_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_CFN_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_CFN_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_CFN_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_50_years_CFN_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_50_years_CFN_0_100.tif"),
              overwrite = T)
  rm(stock_map_adj_50_years_CFN_0_100)
  
  # stock_map_adj_50_years_COA_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_COA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_COA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_COA_15_30.tif")) * cfc_list[["cfc_15_30"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_COA_30_60.tif")) * cfc_list[["cfc_30_60"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_COA_60_100.tif")) * cfc_list[["cfc_60_100"]]
  # writeRaster(stock_map_adj_50_years_COA_0_100,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_50_years_COA_0_100.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_50_years_COA_0_100)
  
  stock_map_adj_50_years_CON_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_CON_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_CON_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_CON_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_CON_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_CON_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_50_years_CON_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_50_years_CON_0_100.tif"),
              overwrite = T)
  rm(stock_map_adj_50_years_CON_0_100)
  
  # stock_map_adj_50_years_PFA_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_PFA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_PFA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_PFA_15_30.tif")) * cfc_list[["cfc_15_30"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_PFA_30_60.tif")) * cfc_list[["cfc_30_60"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_PFA_60_100.tif")) * cfc_list[["cfc_60_100"]]
  # writeRaster(stock_map_adj_50_years_PFA_0_100,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_50_years_PFA_0_100.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_50_years_PFA_0_100)
  
  stock_map_adj_50_years_PFN_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_PFN_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_PFN_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_PFN_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_PFN_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_PFN_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_50_years_PFN_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_50_years_PFN_0_100.tif"),
              overwrite = T)
  rm(stock_map_adj_50_years_PFN_0_100)
  
  # stock_map_adj_50_years_POA_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_POA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_POA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_POA_15_30.tif")) * cfc_list[["cfc_15_30"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_POA_30_60.tif")) * cfc_list[["cfc_30_60"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_POA_60_100.tif")) * cfc_list[["cfc_60_100"]]
  # writeRaster(stock_map_adj_50_years_POA_0_100,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_50_years_POA_0_100.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_50_years_POA_0_100)
  
  stock_map_adj_50_years_PON_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_PON_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_PON_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_PON_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_PON_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_PON_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_50_years_PON_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_50_years_PON_0_100.tif"),
              overwrite = T)
  rm(stock_map_adj_50_years_PON_0_100)
  
  
  # 100_years
  # stock_map_adj_100_years_CFA_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_CFA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_CFA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_CFA_15_30.tif")) * cfc_list[["cfc_15_30"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_CFA_30_60.tif")) * cfc_list[["cfc_30_60"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_CFA_60_100.tif")) * cfc_list[["cfc_60_100"]]
  # writeRaster(stock_map_adj_100_years_CFA_0_100,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_100_years_CFA_0_100.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_100_years_CFA_0_100)
  
  stock_map_adj_100_years_CFN_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_CFN_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_CFN_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_CFN_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_CFN_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_CFN_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_100_years_CFN_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_100_years_CFN_0_100.tif"),
              overwrite = T)
  rm(stock_map_adj_100_years_CFN_0_100)
  
  # stock_map_adj_100_years_COA_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_COA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_COA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_COA_15_30.tif")) * cfc_list[["cfc_15_30"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_COA_30_60.tif")) * cfc_list[["cfc_30_60"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_COA_60_100.tif")) * cfc_list[["cfc_60_100"]]
  # writeRaster(stock_map_adj_100_years_COA_0_100,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_100_years_COA_0_100.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_100_years_COA_0_100)
  
  stock_map_adj_100_years_CON_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_CON_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_CON_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_CON_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_CON_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_CON_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_100_years_CON_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_100_years_CON_0_100.tif"),
              overwrite = T)
  rm(stock_map_adj_100_years_CON_0_100)
  
  # stock_map_adj_100_years_PFA_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_PFA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_PFA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_PFA_15_30.tif")) * cfc_list[["cfc_15_30"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_PFA_30_60.tif")) * cfc_list[["cfc_30_60"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_PFA_60_100.tif")) * cfc_list[["cfc_60_100"]]
  # writeRaster(stock_map_adj_100_years_PFA_0_100,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_100_years_PFA_0_100.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_100_years_PFA_0_100)
  
  stock_map_adj_100_years_PFN_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_PFN_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_PFN_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_PFN_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_PFN_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_PFN_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_100_years_PFN_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_100_years_PFN_0_100.tif"),
              overwrite = T)
  rm(stock_map_adj_100_years_PFN_0_100)
  
  # stock_map_adj_100_years_POA_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_POA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_POA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_POA_15_30.tif")) * cfc_list[["cfc_15_30"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_POA_30_60.tif")) * cfc_list[["cfc_30_60"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_POA_60_100.tif")) * cfc_list[["cfc_60_100"]]
  # writeRaster(stock_map_adj_100_years_POA_0_100,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_100_years_POA_0_100.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_100_years_POA_0_100)
  
  stock_map_adj_100_years_PON_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_PON_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_PON_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_PON_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_PON_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_PON_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_100_years_PON_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_100_years_PON_0_100.tif"),
              overwrite = T)
  rm(stock_map_adj_100_years_PON_0_100)
  
  
  # 200_years
  # stock_map_adj_200_years_CFA_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_CFA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_CFA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_CFA_15_30.tif")) * cfc_list[["cfc_15_30"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_CFA_30_60.tif")) * cfc_list[["cfc_30_60"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_CFA_60_100.tif")) * cfc_list[["cfc_60_100"]]
  # writeRaster(stock_map_adj_200_years_CFA_0_100,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_200_years_CFA_0_100.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_200_years_CFA_0_100)
  
  stock_map_adj_200_years_CFN_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_CFN_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_CFN_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_CFN_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_CFN_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_CFN_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_200_years_CFN_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_200_years_CFN_0_100.tif"),
              overwrite = T)
  rm(stock_map_adj_200_years_CFN_0_100)
  
  # stock_map_adj_200_years_COA_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_COA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_COA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_COA_15_30.tif")) * cfc_list[["cfc_15_30"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_COA_30_60.tif")) * cfc_list[["cfc_30_60"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_COA_60_100.tif")) * cfc_list[["cfc_60_100"]]
  # writeRaster(stock_map_adj_200_years_COA_0_100,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_200_years_COA_0_100.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_200_years_COA_0_100)
  
  stock_map_adj_200_years_CON_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_CON_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_CON_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_CON_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_CON_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_CON_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_200_years_CON_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_200_years_CON_0_100.tif"),
              overwrite = T)
  rm(stock_map_adj_200_years_CON_0_100)
  
  # stock_map_adj_200_years_PFA_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_PFA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_PFA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_PFA_15_30.tif")) * cfc_list[["cfc_15_30"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_PFA_30_60.tif")) * cfc_list[["cfc_30_60"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_PFA_60_100.tif")) * cfc_list[["cfc_60_100"]]
  # writeRaster(stock_map_adj_200_years_PFA_0_100,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_200_years_PFA_0_100.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_200_years_PFA_0_100)
  
  stock_map_adj_200_years_PFN_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_PFN_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_PFN_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_PFN_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_PFN_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_PFN_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_200_years_PFN_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_200_years_PFN_0_100.tif"),
              overwrite = T)
  rm(stock_map_adj_200_years_PFN_0_100)
  
  # stock_map_adj_200_years_POA_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_POA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_POA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_POA_15_30.tif")) * cfc_list[["cfc_15_30"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_POA_30_60.tif")) * cfc_list[["cfc_30_60"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_POA_60_100.tif")) * cfc_list[["cfc_60_100"]]
  # writeRaster(stock_map_adj_200_years_POA_0_100,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_200_years_POA_0_100.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_200_years_POA_0_100)
  
  stock_map_adj_200_years_PON_0_100 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_PON_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_PON_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_PON_15_30.tif")) * cfc_list[["cfc_15_30"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_PON_30_60.tif")) * cfc_list[["cfc_30_60"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_PON_60_100.tif")) * cfc_list[["cfc_60_100"]]
  writeRaster(stock_map_adj_200_years_PON_0_100,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_200_years_PON_0_100.tif"),
              overwrite = T)
  rm(stock_map_adj_200_years_PON_0_100)
  
  
  
  # 0 to 30 cm version
  
  # crop (which is independent of year)
  stock_map_adj_init_crop_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_init_crop_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_init_crop_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_init_crop_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_init_crop_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_init_crop_0_30.tif"),
              overwrite = T)
  
  # pasture (which is independent of year)
  stock_map_adj_init_pasture_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_init_pasture_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_init_pasture_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_init_pasture_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_init_pasture_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_init_pasture_0_30.tif"),
              overwrite = T)
  
  # max (which is independent of year)
  stock_map_adj_max_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_max_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_max_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_max_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_max_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_max_0_30.tif"),
              overwrite = T)
  
  
  
  # Now the other types for each year "CFA","CFN","COA","CON","PFA","PFN","POA","PON"
  
  # 5_years
  # stock_map_adj_5_years_CFA_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_CFA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_CFA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_CFA_15_30.tif")) * cfc_list[["cfc_15_30"]]
  # writeRaster(stock_map_adj_5_years_CFA_0_30,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_5_years_CFA_0_30.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_5_years_CFA_0_30)
  
  stock_map_adj_5_years_CFN_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_CFN_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_CFN_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_CFN_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_5_years_CFN_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_5_years_CFN_0_30.tif"),
              overwrite = T)
  rm(stock_map_adj_5_years_CFN_0_30)
  
  # stock_map_adj_5_years_COA_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_COA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_COA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_COA_15_30.tif")) * cfc_list[["cfc_15_30"]]
  # writeRaster(stock_map_adj_5_years_COA_0_30,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_5_years_COA_0_30.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_5_years_COA_0_30)
  
  stock_map_adj_5_years_CON_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_CON_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_CON_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_CON_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_5_years_CON_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_5_years_CON_0_30.tif"),
              overwrite = T)
  rm(stock_map_adj_5_years_CON_0_30)
  
  # stock_map_adj_5_years_PFA_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_PFA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_PFA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_PFA_15_30.tif")) * cfc_list[["cfc_15_30"]]
  # writeRaster(stock_map_adj_5_years_PFA_0_30,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_5_years_PFA_0_30.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_5_years_PFA_0_30)
  
  stock_map_adj_5_years_PFN_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_PFN_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_PFN_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_PFN_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_5_years_PFN_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_5_years_PFN_0_30.tif"),
              overwrite = T)
  rm(stock_map_adj_5_years_PFN_0_30)
  
  # stock_map_adj_5_years_POA_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_POA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_POA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_POA_15_30.tif")) * cfc_list[["cfc_15_30"]]
  # writeRaster(stock_map_adj_5_years_POA_0_30,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_5_years_POA_0_30.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_5_years_POA_0_30)
  
  stock_map_adj_5_years_PON_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_PON_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_PON_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_5_years_PON_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_5_years_PON_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_5_years_PON_0_30.tif"),
              overwrite = T)
  rm(stock_map_adj_5_years_PON_0_30)
  
  
  # 15_years
  # stock_map_adj_15_years_CFA_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_CFA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_CFA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_CFA_15_30.tif")) * cfc_list[["cfc_15_30"]]
  # writeRaster(stock_map_adj_15_years_CFA_0_30,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_15_years_CFA_0_30.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_15_years_CFA_0_30)
  
  stock_map_adj_15_years_CFN_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_CFN_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_CFN_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_CFN_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_15_years_CFN_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_15_years_CFN_0_30.tif"),
              overwrite = T)
  rm(stock_map_adj_15_years_CFN_0_30)
  
  # stock_map_adj_15_years_COA_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_COA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_COA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_COA_15_30.tif")) * cfc_list[["cfc_15_30"]]
  # writeRaster(stock_map_adj_15_years_COA_0_30,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_15_years_COA_0_30.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_15_years_COA_0_30)
  
  stock_map_adj_15_years_CON_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_CON_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_CON_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_CON_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_15_years_CON_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_15_years_CON_0_30.tif"),
              overwrite = T)
  rm(stock_map_adj_15_years_CON_0_30)
  
  # stock_map_adj_15_years_PFA_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_PFA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_PFA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_PFA_15_30.tif")) * cfc_list[["cfc_15_30"]]
  # writeRaster(stock_map_adj_15_years_PFA_0_30,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_15_years_PFA_0_30.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_15_years_PFA_0_30)
  
  stock_map_adj_15_years_PFN_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_PFN_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_PFN_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_PFN_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_15_years_PFN_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_15_years_PFN_0_30.tif"),
              overwrite = T)
  rm(stock_map_adj_15_years_PFN_0_30)
  
  # stock_map_adj_15_years_POA_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_POA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_POA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_POA_15_30.tif")) * cfc_list[["cfc_15_30"]]
  # writeRaster(stock_map_adj_15_years_POA_0_30,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_15_years_POA_0_30.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_15_years_POA_0_30)
  
  stock_map_adj_15_years_PON_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_PON_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_PON_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_15_years_PON_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_15_years_PON_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_15_years_PON_0_30.tif"),
              overwrite = T)
  rm(stock_map_adj_15_years_PON_0_30)
  
  
  # 30_years
  # stock_map_adj_30_years_CFA_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_CFA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_CFA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_CFA_15_30.tif")) * cfc_list[["cfc_15_30"]]
  # writeRaster(stock_map_adj_30_years_CFA_0_30,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_30_years_CFA_0_30.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_30_years_CFA_0_30)
  
  stock_map_adj_30_years_CFN_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_CFN_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_CFN_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_CFN_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_30_years_CFN_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_30_years_CFN_0_30.tif"),
              overwrite = T)
  rm(stock_map_adj_30_years_CFN_0_30)
  
  # stock_map_adj_30_years_COA_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_COA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_COA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_COA_15_30.tif")) * cfc_list[["cfc_15_30"]]
  # writeRaster(stock_map_adj_30_years_COA_0_30,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_30_years_COA_0_30.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_30_years_COA_0_30)
  
  stock_map_adj_30_years_CON_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_CON_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_CON_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_CON_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_30_years_CON_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_30_years_CON_0_30.tif"),
              overwrite = T)
  rm(stock_map_adj_30_years_CON_0_30)
  
  # stock_map_adj_30_years_PFA_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_PFA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_PFA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_PFA_15_30.tif")) * cfc_list[["cfc_15_30"]]
  # writeRaster(stock_map_adj_30_years_PFA_0_30,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_30_years_PFA_0_30.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_30_years_PFA_0_30)
  
  stock_map_adj_30_years_PFN_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_PFN_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_PFN_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_PFN_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_30_years_PFN_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_30_years_PFN_0_30.tif"),
              overwrite = T)
  rm(stock_map_adj_30_years_PFN_0_30)
  
  # stock_map_adj_30_years_POA_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_POA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_POA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_POA_15_30.tif")) * cfc_list[["cfc_15_30"]]
  # writeRaster(stock_map_adj_30_years_POA_0_30,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_30_years_POA_0_30.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_30_years_POA_0_30)
  
  stock_map_adj_30_years_PON_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_PON_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_PON_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_30_years_PON_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_30_years_PON_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_30_years_PON_0_30.tif"),
              overwrite = T)
  rm(stock_map_adj_30_years_PON_0_30)
  
  # 50_years
  # stock_map_adj_50_years_CFA_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_CFA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_CFA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_CFA_15_30.tif")) * cfc_list[["cfc_15_30"]]
  # writeRaster(stock_map_adj_50_years_CFA_0_30,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_50_years_CFA_0_30.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_50_years_CFA_0_30)
  
  stock_map_adj_50_years_CFN_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_CFN_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_CFN_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_CFN_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_50_years_CFN_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_50_years_CFN_0_30.tif"),
              overwrite = T)
  rm(stock_map_adj_50_years_CFN_0_30)
  
  # stock_map_adj_50_years_COA_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_COA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_COA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_COA_15_30.tif")) * cfc_list[["cfc_15_30"]]
  # writeRaster(stock_map_adj_50_years_COA_0_30,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_50_years_COA_0_30.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_50_years_COA_0_30)
  
  stock_map_adj_50_years_CON_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_CON_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_CON_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_CON_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_50_years_CON_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_50_years_CON_0_30.tif"),
              overwrite = T)
  rm(stock_map_adj_50_years_CON_0_30)
  
  # stock_map_adj_50_years_PFA_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_PFA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_PFA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_PFA_15_30.tif")) * cfc_list[["cfc_15_30"]]
  # writeRaster(stock_map_adj_50_years_PFA_0_30,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_50_years_PFA_0_30.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_50_years_PFA_0_30)
  
  stock_map_adj_50_years_PFN_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_PFN_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_PFN_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_PFN_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_50_years_PFN_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_50_years_PFN_0_30.tif"),
              overwrite = T)
  rm(stock_map_adj_50_years_PFN_0_30)
  
  # stock_map_adj_50_years_POA_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_POA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_POA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_POA_15_30.tif")) * cfc_list[["cfc_15_30"]]
  # writeRaster(stock_map_adj_50_years_POA_0_30,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_50_years_POA_0_30.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_50_years_POA_0_30)
  
  stock_map_adj_50_years_PON_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_PON_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_PON_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_50_years_PON_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_50_years_PON_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_50_years_PON_0_30.tif"),
              overwrite = T)
  rm(stock_map_adj_50_years_PON_0_30)
  
  
  # 100_years
  # stock_map_adj_100_years_CFA_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_CFA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_CFA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_CFA_15_30.tif")) * cfc_list[["cfc_15_30"]]
  # writeRaster(stock_map_adj_100_years_CFA_0_30,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_100_years_CFA_0_30.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_100_years_CFA_0_30)
  
  stock_map_adj_100_years_CFN_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_CFN_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_CFN_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_CFN_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_100_years_CFN_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_100_years_CFN_0_30.tif"),
              overwrite = T)
  rm(stock_map_adj_100_years_CFN_0_30)
  
  # stock_map_adj_100_years_COA_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_COA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_COA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_COA_15_30.tif")) * cfc_list[["cfc_15_30"]]
  # writeRaster(stock_map_adj_100_years_COA_0_30,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_100_years_COA_0_30.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_100_years_COA_0_30)
  
  stock_map_adj_100_years_CON_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_CON_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_CON_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_CON_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_100_years_CON_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_100_years_CON_0_30.tif"),
              overwrite = T)
  rm(stock_map_adj_100_years_CON_0_30)
  
  # stock_map_adj_100_years_PFA_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_PFA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_PFA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_PFA_15_30.tif")) * cfc_list[["cfc_15_30"]]
  # writeRaster(stock_map_adj_100_years_PFA_0_30,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_100_years_PFA_0_30.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_100_years_PFA_0_30)
  
  stock_map_adj_100_years_PFN_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_PFN_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_PFN_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_PFN_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_100_years_PFN_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_100_years_PFN_0_30.tif"),
              overwrite = T)
  rm(stock_map_adj_100_years_PFN_0_30)
  
  # stock_map_adj_100_years_POA_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_POA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_POA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_POA_15_30.tif")) * cfc_list[["cfc_15_30"]]
  # writeRaster(stock_map_adj_100_years_POA_0_30,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_100_years_POA_0_30.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_100_years_POA_0_30)
  
  stock_map_adj_100_years_PON_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_PON_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_PON_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_100_years_PON_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_100_years_PON_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_100_years_PON_0_30.tif"),
              overwrite = T)
  rm(stock_map_adj_100_years_PON_0_30)
  
  
  # 200_years
  # stock_map_adj_200_years_CFA_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_CFA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_CFA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_CFA_15_30.tif")) * cfc_list[["cfc_15_30"]]
  # writeRaster(stock_map_adj_200_years_CFA_0_30,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_200_years_CFA_0_30.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_200_years_CFA_0_30)
  
  stock_map_adj_200_years_CFN_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_CFN_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_CFN_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_CFN_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_200_years_CFN_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_200_years_CFN_0_30.tif"),
              overwrite = T)
  rm(stock_map_adj_200_years_CFN_0_30)
  
  # stock_map_adj_200_years_COA_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_COA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_COA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_COA_15_30.tif")) * cfc_list[["cfc_15_30"]]
  # writeRaster(stock_map_adj_200_years_COA_0_30,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_200_years_COA_0_30.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_200_years_COA_0_30)
  
  stock_map_adj_200_years_CON_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_CON_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_CON_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_CON_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_200_years_CON_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_200_years_CON_0_30.tif"),
              overwrite = T)
  rm(stock_map_adj_200_years_CON_0_30)
  
  # stock_map_adj_200_years_PFA_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_PFA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_PFA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_PFA_15_30.tif")) * cfc_list[["cfc_15_30"]]
  # writeRaster(stock_map_adj_200_years_PFA_0_30,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_200_years_PFA_0_30.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_200_years_PFA_0_30)
  
  stock_map_adj_200_years_PFN_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_PFN_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_PFN_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_PFN_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_200_years_PFN_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_200_years_PFN_0_30.tif"),
              overwrite = T)
  rm(stock_map_adj_200_years_PFN_0_30)
  
  # stock_map_adj_200_years_POA_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_POA_0_5.tif")) * cfc_list[["cfc_0_5"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_POA_5_15.tif")) * cfc_list[["cfc_5_15"]] +
  #   rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_POA_15_30.tif")) * cfc_list[["cfc_15_30"]]
  # writeRaster(stock_map_adj_200_years_POA_0_30,
  #             file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_200_years_POA_0_30.tif"),
  #             overwrite = T)
  # rm(stock_map_adj_200_years_POA_0_30)
  
  stock_map_adj_200_years_PON_0_30 <- rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_PON_0_5.tif")) * cfc_list[["cfc_0_5"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_PON_5_15.tif")) * cfc_list[["cfc_5_15"]] +
    rast(paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_200_years_PON_15_30.tif")) * cfc_list[["cfc_15_30"]]
  writeRaster(stock_map_adj_200_years_PON_0_30,
              file = paste0(uncertainty_dir, "/stock_map_", b_samp,"//stock_map_adj_200_years_PON_0_30.tif"),
              overwrite = T)
  rm(stock_map_adj_200_years_PON_0_30)
  
  
  # Rast_component intermediates are no longer needed; delete to free disk space
  unlink(file.path(uncertainty_dir, paste0("rast_component_", b_samp)), recursive = TRUE)
  
  
  # ── Cleanup ────────────────────────────────────────────────────────────────
  rm(map, mat, clay, nitrogen, irrigation, grid_ref)
  gc()
  terra::tmpFiles(remove = TRUE)
  
  message(sprintf("[%s / b_samp=%d] Done at %s",
                  climate_scenario, b_samp, format(Sys.time(), "%H:%M:%S")))
  return(invisible(NULL))
}


# =============================================================================
# SECTION 5 — LAUNCH
# =============================================================================

t1 <- Sys.time()

mclapply(
  seq_len(nrow(jobs)),
  function(idx) {
    run_one_job(
      b_samp           = jobs$b_samp[idx],
      climate_scenario = jobs$climate_scenario[idx],
      agg_fact         = agg_fact,
      mod_samp_df      = mod_samp_df,
      sc_df            = sc_df
    )
  },
  mc.cores       = n_workers,
  mc.preschedule = FALSE  # complete each job before starting next
)

t2 <- Sys.time()
message(sprintf("All jobs completed in %.1f hours",
                as.numeric(t2 - t1, units = "hours")))