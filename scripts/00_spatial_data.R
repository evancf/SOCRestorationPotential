# Our goal here is to get several global rasters in the same
# resolution / crs saved locally from the GEE catalogue.

# Note that this may take hours to download. Likely about 10 minutes per raster.

# Users must get set up with Google Earth Engine for this code to run.
# Follow the instructions at the following link to get set up with EE and the
# rgee package: https://cran.r-project.org/web/packages/rgee/vignettes/rgee01.html

# Toward the bottom are some other rasters that can be
# downloaded directly from the linked URLs.
# Note also that users are supplied with ./data_inputs/raster/abandon_frac_1km.tif
# that is not available elsewhere but generated for this publication.

topwd <- getwd()

# ipak() installs and loads packages in one call
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(c("tidyverse", "rgee"))
ee_Initialize(user = "USERNAME", gcs = T, drive = T)


# If the above does not work, may need to use some of the following.
# #https://github.com/r-spatial/rgee/issues/355
# # Initialize - this will connect to a project. You should always call this
# # before working with rgee. It is IMPORTANT THAT YOU SPECIFY A PROJECT using
# # the project parameter. If you forget what project IDs you have access to, find them
# # here: console.cloud.google.com/project
# ee$Authenticate(auth_mode='notebook')
# ee$Initialize(project='ee-USERNAME')  # <-- EDIT THIS FOR YOUR PROJECT
# ee_Initialize(email = NULL, drive = TRUE, gcs = FALSE)
# ee_Initialize(user = "USERNAME", gcs = T, drive = T)

# Functions
merge_project_gee_rasters <- function(x){
  t0 <- Sys.time()
  filename <- deparse(substitute(x)) %>%
    gsub("_raster", ".tif", ., fixed = T)
  
  y <- terra::merge(terra::rast(x[1]),
                    terra::rast(x[2])) %>%
    terra::project(aban_frac_terra)
  terra::writeRaster(y,
                     filename,
                     overwrite = T)
  t1 <- Sys.time()
  print("Done!")
  print(t1-t0)
}



# Will use this as the template
setwd(topwd)
aban_frac_terra <- terra::rast("./data_inputs/raster/abandon_frac_1km.tif")

# Get a bounding box
ee_bbox <- ee$Geometry$BBox(-180, -56, 180, 84)

# Want to work in this via_gee folder
setwd(topwd)
setwd("./data_inputs/raster/via_gee/")

# Temperature raster -----------------------------------------------------------
if(!"amt.tif" %in% list.files()){
  amt_raster <- ee$
    Image("WORLDCLIM/V1/BIO")$select("bio01")$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "amt.tif")
  # Now save a projected file
  merge_project_gee_rasters(amt_raster)
}


# Precipitation raster ---------------------------------------------------------
if(!"tap.tif" %in% list.files()){
  tap_raster <- ee$
    Image("WORLDCLIM/V1/BIO")$select("bio12")$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "tap.tif")
  # Now save a projected file
  merge_project_gee_rasters(tap_raster)
}

tap_raster <- ee$
  Image("WORLDCLIM/V1/BIO")$select("bio12")


# Future climate CMIP6 ---------------------------------------------------------

ipak("terra")

out_dir <- "./data_inputs/raster/cmip6/ssp370_2050"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# All nine GCMs available for SSP3-7.0 at 10-minute resolution
gcms <- c(
  "ACCESS-CM2", "BCC-CSM2-MR", "CMCC-ESM2", "EC-Earth3-Veg",
  "GISS-E2-1-G", "INM-CM5-0", "IPSL-CM6A-LR",
  "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL"
)
# Note: FIO-ESM-2-0 and HadGEM3-GC31-LL have no SSP3-7.0 entry for 2041-2060

base_url <- "https://geodata.ucdavis.edu/cmip6/10m"
period   <- "2041-2060"
ssp      <- "ssp370"

# Download all bioc files and extract BIO1 (MAT) and BIO12 (MAP)
mat_rasts <- list()
map_rasts <- list()

for (gcm in gcms) {
  url      <- sprintf("%s/%s/%s/wc2.1_10m_bioc_%s_%s_%s.tif",
                      base_url, gcm, ssp, gcm, ssp, period)
  destfile <- file.path(out_dir, sprintf("bioc_%s.tif", gcm))
  
  if (!file.exists(destfile)) {
    message("Downloading: ", gcm)
    download.file(url, destfile, mode = "wb", quiet = TRUE)
  } else {
    message("Cached: ", gcm)
  }
  
  r <- rast(destfile)
  mat_rasts[[gcm]] <- r[[1]]   # BIO1 = mean annual temperature (°C)
  map_rasts[[gcm]] <- r[[12]]  # BIO12 = annual precipitation (mm)
}

# Average across GCMs
amt_future <- mean(rast(mat_rasts))
tap_future <- mean(rast(map_rasts))

# WorldClim CMIP6 BIO1 is in plain °C; multiply by 10 to match the
# °C × 10 convention of the contemporary amt.tif (WorldClim V1 via GEE),
# which load_scaled_agg() expects to divide by 10 before centering/scaling.
amt_future <- amt_future * 10

writeRaster(amt_future, file.path(out_dir, "amt.tif"), overwrite = TRUE)
writeRaster(tap_future, file.path(out_dir, "tap.tif"), overwrite = TRUE)

message("Done. Future MAT/MAP rasters written to ", out_dir)
rm(amt_future, tap_future, map_rasts, mat_rasts)

# Soil Grids -------------------------------------------------------------------
# Info here https://gee-community-catalog.org/projects/isric/
# Bulk density  ----------------------------------------------------------------
if(!"bdod_0_5cm_mean.tif" %in% list.files()){
  bdod_0_5cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/bdod_mean")$select("bdod_0-5cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "bdod_0_5cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(bdod_0_5cm_mean_raster)
}

if(!"bdod_5_15cm_mean.tif" %in% list.files()){
  bdod_5_15cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/bdod_mean")$select("bdod_5-15cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "bdod_5_15cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(bdod_5_15cm_mean_raster)
}

if(!"bdod_15_30cm_mean.tif" %in% list.files()){
  bdod_15_30cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/bdod_mean")$select("bdod_15-30cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "bdod_15_30cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(bdod_15_30cm_mean_raster)
}

if(!"bdod_30_60cm_mean.tif" %in% list.files()){
  bdod_30_60cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/bdod_mean")$select("bdod_30-60cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "bdod_30_60cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(bdod_30_60cm_mean_raster)
}

if(!"bdod_60_100cm_mean.tif" %in% list.files()){
  bdod_60_100cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/bdod_mean")$select("bdod_60-100cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "bdod_60_100cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(bdod_60_100cm_mean_raster)
}

if(!"bdod_100_200cm_mean.tif" %in% list.files()){
  bdod_100_200cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/bdod_mean")$select("bdod_100-200cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 2000,
      via = "drive",
      dsn = "bdod_100_200cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(bdod_100_200cm_mean_raster)
}


# Coarse fragments -------------------------------------------------------------

if(!"cfvo_0_5cm_mean.tif" %in% list.files()){
  cfvo_0_5cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/cfvo_mean")$select("cfvo_0-5cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "cfvo_0_5cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(cfvo_0_5cm_mean_raster)
}

if(!"cfvo_5_15cm_mean.tif" %in% list.files()){
  cfvo_5_15cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/cfvo_mean")$select("cfvo_5-15cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "cfvo_5_15cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(cfvo_5_15cm_mean_raster)
}

if(!"cfvo_15_30cm_mean.tif" %in% list.files()){
  cfvo_15_30cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/cfvo_mean")$select("cfvo_15-30cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "cfvo_15_30cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(cfvo_15_30cm_mean_raster)
}

if(!"cfvo_30_60cm_mean.tif" %in% list.files()){
  cfvo_30_60cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/cfvo_mean")$select("cfvo_30-60cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "cfvo_30_60cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(cfvo_30_60cm_mean_raster)
}

if(!"cfvo_60_100cm_mean.tif" %in% list.files()){
  cfvo_60_100cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/cfvo_mean")$select("cfvo_60-100cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "cfvo_60_100cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(cfvo_60_100cm_mean_raster)
}

if(!"cfvo_100_200cm_mean.tif" %in% list.files()){
  cfvo_100_200cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/cfvo_mean")$select("cfvo_100-200cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 2000,
      via = "drive",
      dsn = "cfvo_100_200cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(cfvo_100_200cm_mean_raster)
}



# OCD --------------------------------------------------------------------------

if(!"ocd_0_5cm_mean.tif" %in% list.files()){
  ocd_0_5cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/ocd_mean")$select("ocd_0-5cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "ocd_0_5cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(ocd_0_5cm_mean_raster)
}

if(!"ocd_5_15cm_mean.tif" %in% list.files()){
  ocd_5_15cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/ocd_mean")$select("ocd_5-15cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "ocd_5_15cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(ocd_5_15cm_mean_raster)
}

if(!"ocd_15_30cm_mean.tif" %in% list.files()){
  ocd_15_30cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/ocd_mean")$select("ocd_15-30cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "ocd_15_30cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(ocd_15_30cm_mean_raster)
}

if(!"ocd_30_60cm_mean.tif" %in% list.files()){
  ocd_30_60cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/ocd_mean")$select("ocd_30-60cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "ocd_30_60cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(ocd_30_60cm_mean_raster)
}

if(!"ocd_60_100cm_mean.tif" %in% list.files()){
  ocd_60_100cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/ocd_mean")$select("ocd_60-100cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "ocd_60_100cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(ocd_60_100cm_mean_raster)
}


if(!"ocd_100_200cm_mean.tif" %in% list.files()){
  ocd_100_200cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/ocd_mean")$select("ocd_100-200cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 2000,
      via = "drive",
      dsn = "ocd_100_200cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(ocd_100_200cm_mean_raster)
}


# OCS --------------------------------------------------------------------------

if(!"ocs_mean.tif" %in% list.files()){
  ocs_mean_raster <- ee$
    Image("projects/soilgrids-isric/ocs_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "ocs_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(ocs_mean_raster)
}


# clay -------------------------------------------------------------------------
if(!"clay_0_5cm_mean.tif" %in% list.files()){
  clay_0_5cm_mean_raster <- ee$
      Image("projects/soilgrids-isric/clay_mean")$select("clay_0-5cm_mean")$
      clip(ee_bbox)$
      reduceResolution(
        reducer = ee$Reducer$mean(),
        maxPixels = 2000) %>%
      ee_as_raster(
        scale = 1000,
        via = "drive",
        dsn = "clay_0_5cm_mean.tif")
    # Now save a projected file
    merge_project_gee_rasters(clay_mean_raster)
}

if(!"clay_5_15cm_mean.tif" %in% list.files()){
  clay_5_15cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/clay_mean")$select("clay_5-15cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "clay_5_15cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(clay_5_15cm_mean_raster)
}

if(!"clay_15_30cm_mean.tif" %in% list.files()){
  clay_15_30cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/clay_mean")$select("clay_15-30cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "clay_15_30cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(clay_15_30cm_mean_raster)
}

if(!"clay_30_60cm_mean.tif" %in% list.files()){
  clay_30_60cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/clay_mean")$select("clay_30-60cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "clay_30_60cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(clay_30_60cm_mean_raster)
}

if(!"clay_60_100cm_mean.tif" %in% list.files()){
  clay_60_100cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/clay_mean")$select("clay_60-100cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "clay_60_100cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(clay_60_100cm_mean_raster)
}

if(!"clay_100_200cm_mean.tif" %in% list.files()){
  clay_100_200cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/clay_mean")$select("clay_100-200cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "clay_100_200cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(clay_100_200cm_mean_raster)
}

# Want to get an average weighted by depth too
# This is a little clunky to manage memory
if(!"clay_top_meter_mean.tif" %in% list.files()){
  ipak("stars")
  clay_top_meter_mean <- read_stars("clay_60_100cm_mean.tif",
                                    proxy = F) * 0.4 + 
    read_stars("clay_30_60cm_mean.tif",
               proxy = F) * 0.3
  gc()
  clay_top_meter_mean <- clay_top_meter_mean + 
    read_stars("clay_15_30cm_mean.tif",
               proxy = F) * 0.15
  gc()
  clay_top_meter_mean <- clay_top_meter_mean + 
    read_stars("clay_5_15cm_mean.tif",
               proxy = F) * 0.1
  gc()
  clay_top_meter_mean <- clay_top_meter_mean + 
    read_stars("clay_0_5cm_mean.tif",
               proxy = F) * 0.05
  gc()
  # Write this out
  clay_top_meter_mean %>% 
    write_stars(dsn = "clay_top_meter_mean.tif")
}


# Nitrogen ---------------------------------------------------------------------

if(!"nitrogen_0_5cm_mean.tif" %in% list.files()){
  nitrogen_0_5cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/nitrogen_mean")$select("nitrogen_0-5cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "nitrogen_0_5cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(nitrogen_0_5cm_mean_raster)
}

if(!"nitrogen_5_15cm_mean.tif" %in% list.files()){
  nitrogen_5_15cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/nitrogen_mean")$select("nitrogen_5-15cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "nitrogen_5_15cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(nitrogen_5_15cm_mean_raster)
}

if(!"nitrogen_15_30cm_mean.tif" %in% list.files()){
  nitrogen_15_30cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/nitrogen_mean")$select("nitrogen_15-30cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "nitrogen_15_30cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(nitrogen_15_30cm_mean_raster)
}

if(!"nitrogen_30_60cm_mean.tif" %in% list.files()){
  nitrogen_30_60cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/nitrogen_mean")$select("nitrogen_30-60cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "nitrogen_30_60cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(nitrogen_30_60cm_mean_raster)
}

if(!"nitrogen_60_100cm_mean.tif" %in% list.files()){
  nitrogen_60_100cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/nitrogen_mean")$select("nitrogen_60-100cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "nitrogen_60_100cm_mean.tif")
  # Now save a projected file
  merge_project_gee_rasters(nitrogen_60_100cm_mean_raster)
}


if(!"nitrogen_100_200cm_mean.tif" %in% list.files()){
  nitrogen_100_200cm_mean_raster <- ee$
    Image("projects/soilgrids-isric/nitrogen_mean")$select("nitrogen_100-200cm_mean")$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 2000,
      via = "drive",
      dsn = "nitrogen_100_200cm_mean.tif")
  
  
  # For some reason, this returns as single raster...
  terra::writeRaster(terra::rast(nitrogen_100_200cm_mean_raster) %>%
                       terra::project(aban_frac_terra),
                     "nitrogen_100_200cm_mean.tif",
                     overwrite = T)
  # # Now save a projected file 
  # merge_project_gee_rasters(nitrogen_100_200cm_mean_raster)
}








# Copernicus land cover --------------------------------------------------------

# Need to modify the following function
merge_project_gee_rasters_4 <- function(x){
  t0 <- Sys.time()
  filename <- deparse(substitute(x)) %>%
    gsub("_raster", ".tif", ., fixed = T)
  
  y <- terra::merge(terra::rast(x[1]),
                    terra::rast(x[2]),
                    terra::rast(x[3]),
                    terra::rast(x[4])) %>%
    terra::project(aban_frac_terra)
  terra::writeRaster(y,
                     filename,
                     overwrite = T)
  t1 <- Sys.time()
  print("Done!")
  print(t1-t0)
}

if(!"cropland_copernicus.tif" %in% list.files()){
  cropland_copernicus_raster <- ee$
    ImageCollection("COPERNICUS/Landcover/100m/Proba-V-C3/Global")$
    filter(ee$Filter$eq('system:index', '2019'))$ # Select the year 2019
    first()$
    select('crops-coverfraction')$
    clip(ee_bbox)$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "cropland_copernicus.tif",
      maxPixels = 3e9)
  # Now save a projected file
  merge_project_gee_rasters_4(cropland_copernicus_raster)
}


# SEDAC Pasture area -----------------------------------------------------------
# The product is from this link:
# https://www.earthdata.nasa.gov/data/catalog/sedac-ciesin-sedac-aglands-pas2000-1.00

setwd(topwd)
setwd("./data_inputs/raster/")
if(!"pasture_sedac.tif" %in% list.files()){
  pasture <- rast("./gl-pastures-geotif/pasture.tif")
  pasture <- resample(pasture, aban_frac_terra)
  terra::writeRaster(pasture,
                     "pasture_sedac.tif",
                     overwrite = T)
}


# Global map of irrigated areas  -----------------------------------------------
# This is FAO's 5 minute resolution map of areas equipped for irrigation
setwd(topwd)
setwd("./data_inputs/raster/")
if(!"gmia.tif" %in% list.files()){
  irrigation_url <- "https://firebasestorage.googleapis.com/v0/b/fao-aquastat.appspot.com/o/GIS%2Fgmia_v5_aei_pct_asc.zip?alt=media&token=e448ce53-296f-4756-90c1-75c87f74e569"
  download_file <- "gmia_v5_aei_pct_asc.zip"
  download.file(irrigation_url, download_file)
  unzip(download_file, exdir = "gmia_data")
  irrigation_path <- list.files("gmia_data", pattern = "\\.asc$", full.names = TRUE)
  irrigation_raster <- rast(irrigation_path)
  irrigation_resampled_raster <- resample(irrigation_raster, aban_frac_terra)
  terra::writeRaster(irrigation_resampled_raster,
                     "gmia.tif",
                     overwrite = T)
  
  # Also want a irrigation+0.1 raster to facilitate log transformation later
  terra::writeRaster(irrigation_resampled_raster + 0.1,
                     "gmia_plus_0.1.tif",
                     overwrite = T)
}



# HYDE Pasture area ------------------------------------------------------------
# The product is from this link:
# https://geo.public.data.uu.nl/vault-hyde/HYDE%203.3%5B1710493486%5D/original/hyde33_c7_base_mrt2023/NetCDF/
# Area file here: https://geo.public.data.uu.nl/vault-hyde/HYDE%203.3%5B1710493486%5D/original/
setwd(topwd)
setwd("./data_inputs/raster/")
if(!"pasture_frac_2010_HYDE.tif" %in% list.files()){
  hyde_area <- rast("./HYDE3.3/general_files/general_files/garea_cr.asc")
  pasture <- rast("./HYDE3.3/pasture.nc")
  time(pasture)
  pasture_2010 <- pasture[[which(time(pasture) == "2010-05-01")]]
  pasture_2023 <- pasture[[which(time(pasture) == "2023-05-01")]]
  
  # Make these fractional area
  pasture_2010 <- pasture_2010 / hyde_area
  pasture_2023 <- pasture_2023 / hyde_area
  
  pasture_2010 <- resample(pasture_2010, aban_frac_terra)
  pasture_2023 <- resample(pasture_2023, aban_frac_terra)
  
  terra::writeRaster(pasture_2010,
                     "pasture_frac_2010_HYDE.tif",
                     overwrite = T)
  terra::writeRaster(pasture_2023,
                     "pasture_frac_2023_HYDE.tif",
                     overwrite = T)
  rm(pasture, pasture_2010, pasture_2023, hyde_area)
}


# HYDE cropland area ------------------------------------------------------------
# The product is from this link:
# https://geo.public.data.uu.nl/vault-hyde/HYDE%203.3%5B1710493486%5D/original/hyde33_c7_base_mrt2023/NetCDF/
# Area file here: https://geo.public.data.uu.nl/vault-hyde/HYDE%203.3%5B1710493486%5D/original/
setwd(topwd)
setwd("./data_inputs/raster/")
if(!"cropland_frac_2010_HYDE.tif" %in% list.files()){
  hyde_area <- rast("./HYDE3.3/general_files/general_files/garea_cr.asc")
  cropland <- rast("./HYDE3.3/cropland.nc")
  time(cropland)
  cropland_2010 <- cropland[[which(time(cropland) == "2010-05-01")]]
  cropland_2023 <- cropland[[which(time(cropland) == "2023-05-01")]]
  
  # Make these fractional area
  cropland_2010 <- cropland_2010 / hyde_area
  cropland_2023 <- cropland_2023 / hyde_area
  
  cropland_2010 <- resample(cropland_2010, aban_frac_terra)
  cropland_2023 <- resample(cropland_2023, aban_frac_terra)
  
  terra::writeRaster(cropland_2010,
                     "cropland_frac_2010_HYDE.tif",
                     overwrite = T)
  terra::writeRaster(cropland_2023,
                     "cropland_frac_2023_HYDE.tif",
                     overwrite = T)
  rm(cropland, cropland_2010, cropland_2023, hyde_area)
}


# HYDE grazing area ------------------------------------------------------------
# The product is from this link:
# https://geo.public.data.uu.nl/vault-hyde/HYDE%203.3%5B1710493486%5D/original/hyde33_c7_base_mrt2023/NetCDF/
# Area file here: https://geo.public.data.uu.nl/vault-hyde/HYDE%203.3%5B1710493486%5D/original/
setwd(topwd)
setwd("./data_inputs/raster/")
if(!"grazing_frac_2010_HYDE.tif" %in% list.files()){
  hyde_area <- rast("./HYDE3.3/general_files/general_files/garea_cr.asc")
  grazing <- rast("./HYDE3.3/grazing_land.nc")
  time(grazing)
  grazing_2010 <- grazing[[which(time(grazing) == "2010-05-01")]]
  grazing_2023 <- grazing[[which(time(grazing) == "2023-05-01")]]
  
  # Make these fractional area
  grazing_2010 <- grazing_2010 / hyde_area
  grazing_2023 <- grazing_2023 / hyde_area
  
  grazing_2010 <- resample(grazing_2010, aban_frac_terra)
  grazing_2023 <- resample(grazing_2023, aban_frac_terra)
  
  terra::writeRaster(grazing_2010,
                     "grazing_frac_2010_HYDE.tif",
                     overwrite = T)
  terra::writeRaster(grazing_2023,
                     "grazing_frac_2023_HYDE.tif",
                     overwrite = T)
  rm(grazing, grazing_2010, grazing_2023, hyde_area)
}






# Biomes  ----------------------------------------------------------------------
# This product is available at this link (click About, then download shapefile):
# https://ecoregions.appspot.com/

ipak("terra")

# Load shapefile
ecoregions <- vect("./data_inputs/vector/Ecoregions2017/Ecoregions2017.shp")

# Define which BIOME_NAME values are forest vs non‐forest
forest_types <- c(
  "Tropical & Subtropical Moist Broadleaf Forests",
  "Tropical & Subtropical Dry Broadleaf Forests",
  "Tropical & Subtropical Coniferous Forests",
  "Temperate Conifer Forests",
  "Temperate Broadleaf & Mixed Forests",
  "Boreal Forests/Taiga",
  "Mediterranean Forests, Woodlands & Scrub"
)

nonforest_types <- c(
  "Tundra",
  "Deserts & Xeric Shrublands",
  "Temperate Grasslands, Savannas & Shrublands",
  "Montane Grasslands & Shrublands",
  "Flooded Grasslands & Savannas",
  "Tropical & Subtropical Grasslands, Savannas & Shrublands"
)



# Ensure same CRS
if (!identical(crs(ecoregions), crs(aban_frac_terra))) {
  ecoregions <- project(ecoregions, aban_frac_terra)
}


# 1. define codes: 1 = forest, 2 = non‑forest, everything else stays NA
ecoregions$cat <- NA_integer_
ecoregions$cat[ecoregions$BIOME_NAME %in% forest_types]    <- 1
ecoregions$cat[ecoregions$BIOME_NAME %in% nonforest_types] <- 2

cat_rast <- rasterize(ecoregions, aban_frac_terra, field="cat", touches=TRUE)

# logical rasters
forest_rast    <- cat_rast == 1
nonforest_rast <- cat_rast == 2

# Make 1/NA instead of TRUE/FALSE:
forest_rast    <- as.numeric(forest_rast)    # TRUE→1, FALSE→0
forest_rast[forest_rast == 0]    <- NA
nonforest_rast <- as.numeric(nonforest_rast)
nonforest_rast[nonforest_rast == 0] <- NA

# Write to disk
setwd(topwd)
setwd("./data_inputs/raster/")
writeRaster(forest_rast,    "forest_biomes.tif",    overwrite=TRUE)
writeRaster(nonforest_rast, "nonforest_biomes.tif", overwrite=TRUE)




# May want a treecover raster
if(!"treecover2000.tif" %in% list.files()){
  treecover2000_raster <- ee$
    Image("UMD/hansen/global_forest_change_2024_v1_12")$select("treecover2000")$
    reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = 2000) %>%
    ee_as_raster(
      scale = 1000,
      via = "drive",
      dsn = "treecover2000.tif")
  # Now save a projected file
  merge_project_gee_rasters(treecover2000_raster)
}








# Notes about other products used

# Pasture sparing scenario from Hayek et al. 2024
# https://www.pnas.org/doi/epdf/10.1073/pnas.2405758121
# Use data and scripts presented here: https://zenodo.org/records/12688280
# Run through line 372 in the Pastures_COI_analysis.R script and write
# out the HS.area raster. Note that this will need to be converted to a fraction
# raster to comport with the other scenario rasters.
# writeRaster(HS.area, filename = "pasture_sparing_hotspot_frac_Hayek.tif")

# Belowground ratio raster to be put in a folder called ./data_inputs/raster/belowground
# https://www.nature.com/articles/s41559-021-01485-1
# https://figshare.com/projects/The_global_distribution_and_environmental_drivers_of_aboveground_versus_belowground_plant_biomass/120897

# Sanderman et al. 2017 NoLU_10km tifs are saved in ./data_inputs/raster/Sanderman2017
# https://www.pnas.org/doi/10.1073/pnas.1706103114
# https://github.com/whrc/Soil-Carbon-Debt/tree/master/OCD

# Xu et al. 2026 in ./data_inputs/raster/CCIBM and SBM
# https://www.nature.com/articles/s41586-025-09870-7
# https://zenodo.org/records/15869647

# Robinson et al. 2025 tifs are saved in ./data_inputs/raster/Robinson et al. 2025
# https://www.nature.com/articles/s41558-025-02355-5
# https://zenodo.org/records/15090826

# And a potential treecover raster
# Available here: https://www.research-collection.ethz.ch/entities/researchdata/7234f0a1-e37d-4362-99e9-d0e6e75793dc
#"./data_inputs/raster/Total_potential.tif"