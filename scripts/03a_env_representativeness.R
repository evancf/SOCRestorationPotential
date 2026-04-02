# =============================================================================
# Environmental representativeness: comparing predictor variable distributions
# between field data and global cropland / pasture areas, separately by
# land-use type (C vs P).
#
# Cropland and pasture raster densities are weighted by fractional cover.
# Field data points are unweighted (equal weight).
#
# Outputs:
#   ./figures/env_joint_coverage_table.csv   — multivariate coverage table
#   ./figures/env_density_combined.pdf/.png  — density figure (Extended Data)
#
# This script is self-contained: it reads c_dat_analysis_ready.csv and
# reconstructs c_dat and init_covars, avoiding dplyr conflicts that arise
# when sourcing 03_bayes_analysis.R in the same session.
# =============================================================================

ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(c("terra", "tidyverse"))

# =============================================================================
# Section 0: Reconstruct c_dat and init_covars from the analysis-ready CSV
# =============================================================================

c_dat <- read.csv("./data_inputs/c_dat_analysis_ready.csv",
                  stringsAsFactors = FALSE)

c_dat <- c_dat %>%
  dplyr::mutate(coord_id = paste(LAT, LONG) %>%
                  as.factor() %>% as.numeric())

c_dat$init_id <- c_dat %>%
  dplyr::select(coord_id, depth_avg, Authors, Year) %>%
  apply(1, function(x) paste(x, collapse = "")) %>%
  as.factor() %>% as.numeric()

c_dat$final_id <- c_dat %>%
  dplyr::select(coord_id, depth_avg, Authors, Year, AGE, lu_end, MGMT) %>%
  apply(1, function(x) paste(x, collapse = "")) %>%
  as.factor() %>% as.numeric()

c_dat <- c_dat %>%
  dplyr::mutate(title_id = as.numeric(as.factor(Title)))

c_dat <- c_dat[order(c_dat$final_id), ]

init_dat <- c_dat %>%
  dplyr::select(-dplyr::starts_with("Aban"),
                -dplyr::contains("accrual"),
                -dplyr::contains("Perc_change"),
                -dplyr::contains("Cum_Aban"),
                -dplyr::starts_with("SE"),
                -dplyr::starts_with("SD"),
                -obs,
                -dplyr::starts_with("AGE"),
                -final_id)

init_dat <- init_dat[
  !duplicated(init_dat %>%
                dplyr::select(-Profile_ID, -lu_end, -lu_start)), ]

init_covars <- init_dat %>%
  dplyr::select(init_id, depth_avg, MAP, MAT,
                clay, nitrogen, irrigation,
                title_id, lu_start, lu_end) %>%
  dplyr::group_by(init_id) %>%
  dplyr::summarize(
    dplyr::across(dplyr::where(is.numeric),  mean),
    dplyr::across(dplyr::where(is.character), dplyr::first)
  )

message(sprintf("c_dat: %d rows | init_covars: %d rows (%d cropland, %d pasture)",
                nrow(c_dat), nrow(init_covars),
                sum(init_covars$lu_start == "C"),
                sum(init_covars$lu_start == "P")))

# =============================================================================
# Section 1: Helper functions
# =============================================================================

# -----------------------------------------------------------------------------
# 1a. Weighted quantile
# -----------------------------------------------------------------------------

wquantile <- function(x, w, probs) {
  ord <- order(x)
  x   <- x[ord]; w <- w[ord]
  cw  <- cumsum(w) / sum(w)
  sapply(probs, function(p) x[which(cw >= p)[1]])
}

# -----------------------------------------------------------------------------
# 1b. Weighted KDE on a fixed grid, returns a normalised probability vector
# -----------------------------------------------------------------------------

wkde <- function(x, w, from, to, n = 512) {
  kd <- density(x, weights = w / sum(w), from = from, to = to, n = n)
  kd$y / sum(kd$y)
}

# -----------------------------------------------------------------------------
# 1c. Joint (multivariate) coverage
#
# Evaluates what fraction of global area (weighted) falls within the
# field-data envelope, at several percentile thresholds:
#
#   Full range (0–100%):  has any field observation ever seen this condition?
#   99%  (0.5–99.5%):     full range minus the most extreme 0.5% tails
#   95%  (2.5–97.5%):     conventional central interval
#   90%  (5–95%):         more conservative central interval
#
# Coverage is reported:
#   - Per variable individually (univariate)
#   - Climate variables jointly (MAT + MAP)
#   - Soil variables jointly (clay + nitrogen)
#   - All four variables jointly (requires pixel-level joint data frames)
#
# df_rast / df_field use a split-row structure: each row carries either
# climate or soil values (with the other pair as NA). df_rast_joint /
# df_field_joint have all four variables on every row (stratum-weighted
# soil average joined with climate) and are required for the all-four metric.
# -----------------------------------------------------------------------------

joint_coverage <- function(df_rast, df_field, context_label,
                           df_rast_joint  = NULL,
                           df_field_joint = NULL,
                           probs = c(0.005, 0.995)) {
  
  var_groups <- list(climate = c("mat", "map"),
                     soil    = c("clay", "nitrogen"))
  
  # (i) Per-group joint coverage
  group_results <- map_dfr(names(var_groups), function(grp) {
    gvars <- var_groups[[grp]]
    
    rast_ok  <- rowSums(!is.na(df_rast[,  gvars, drop = FALSE])) == length(gvars)
    field_ok <- rowSums(!is.na(df_field[, gvars, drop = FALSE])) == length(gvars)
    df_r <- df_rast[rast_ok, ]; df_f <- df_field[field_ok, ]
    
    envelope <- map_dfr(gvars, function(v) {
      tibble(variable = v,
             lo = quantile(df_f[[v]], probs[1], na.rm = TRUE),
             hi = quantile(df_f[[v]], probs[2], na.rm = TRUE))
    })
    
    inside <- rep(TRUE, nrow(df_r))
    for (v in gvars) {
      lo_v   <- envelope$lo[envelope$variable == v]
      hi_v   <- envelope$hi[envelope$variable == v]
      inside <- inside & (df_r[[v]] >= lo_v) & (df_r[[v]] <= hi_v)
    }
    
    tibble(group     = grp,
           pct_joint = round(100 * sum(df_r$weight[inside], na.rm = TRUE) /
                               sum(df_r$weight, na.rm = TRUE), 1))
  })
  
  # (ii) All-four-variable joint coverage using the pixel-level joint df
  all_four <- if (!is.null(df_rast_joint) && !is.null(df_field_joint)) {
    envelope_all <- map_dfr(vars, function(v) {
      tibble(variable = v,
             lo = quantile(df_field_joint[[v]], probs[1], na.rm = TRUE),
             hi = quantile(df_field_joint[[v]], probs[2], na.rm = TRUE))
    })
    inside <- rep(TRUE, nrow(df_rast_joint))
    for (v in vars) {
      lo_v   <- envelope_all$lo[envelope_all$variable == v]
      hi_v   <- envelope_all$hi[envelope_all$variable == v]
      inside <- inside & !is.na(df_rast_joint[[v]]) &
        (df_rast_joint[[v]] >= lo_v) & (df_rast_joint[[v]] <= hi_v)
    }
    round(100 * sum(df_rast_joint$weight[inside], na.rm = TRUE) /
            sum(df_rast_joint$weight,             na.rm = TRUE), 1)
  } else {
    NA_real_
  }
  
  # (iii) Per-variable univariate coverage
  univariate <- map_dfr(vars, function(v) {
    grp      <- if (v %in% var_groups$climate) "climate" else "soil"
    field_ok <- rowSums(!is.na(df_field[, var_groups[[grp]], drop = FALSE])) ==
      length(var_groups[[grp]])
    df_f     <- df_field[field_ok, ]
    lo_v     <- quantile(df_f[[v]], probs[1], na.rm = TRUE)
    hi_v     <- quantile(df_f[[v]], probs[2], na.rm = TRUE)
    rast_ok  <- !is.na(df_rast[[v]])
    in_v     <- rast_ok & (df_rast[[v]] >= lo_v) & (df_rast[[v]] <= hi_v)
    tibble(variable   = var_labels[v],
           pct_in_env = round(100 * sum(df_rast$weight[in_v],    na.rm = TRUE) /
                                sum(df_rast$weight[rast_ok], na.rm = TRUE), 1))
  })
  
  list(context        = context_label,
       group_coverage = group_results,
       pct_all_four   = all_four,
       univariate     = univariate)
}

# =============================================================================
# Section 2: Variable definitions
# =============================================================================

vars <- c("mat", "map", "clay", "nitrogen")

var_labels <- c(
  mat      = "Mean annual temperature (°C)",
  map      = "Mean annual precipitation (mm)",
  clay     = "Clay content (g/kg)",
  nitrogen = "Soil nitrogen (mg/kg)"
)

# =============================================================================
# Section 3: Load and prepare raster data
#
# Raster loading, resampling, aggregation, and pixel extraction are the most
# time-consuming steps. If a cache file exists from a previous run, all
# extracted data frames are loaded directly. Delete the cache to force a
# full re-run (e.g. after changing thresholds or raster paths).
# =============================================================================

raster_cache <- "./data_inputs/env_raster_cache.RData"

if (file.exists(raster_cache)) {
  message("Loading raster data from cache: ", raster_cache)
  load(raster_cache)
  message("Cache loaded — skipping raster extraction.")
} else {
  message("No cache found — running full raster extraction (may take several minutes).")
  
  # ---------------------------------------------------------------------------
  # 3a. Cropland and pasture fraction rasters
  # Pixels below 5% cover excluded to reduce HYDE classification noise.
  # ---------------------------------------------------------------------------
  
  crop_frac    <- rast("./data_inputs/raster/cropland_frac_2010_HYDE.tif")
  pasture_frac <- rast("./data_inputs/raster/pasture_frac_2010_HYDE.tif")
  
  crop_frac[crop_frac       < 0.05] <- NA
  pasture_frac[pasture_frac < 0.05] <- NA
  
  # Preserve native resolution for resampling before aggregation
  crop_frac_native <- crop_frac
  
  # ---------------------------------------------------------------------------
  # 3b. Climate rasters — resample to native grid then aggregate
  # ---------------------------------------------------------------------------
  
  wc_mat <- rast("./data_inputs/worldclim/climate/wc2.1_10m/wc2.1_10m_bio_1.tif")
  wc_map <- rast("./data_inputs/worldclim/climate/wc2.1_10m/wc2.1_10m_bio_12.tif")
  
  wc_mat <- resample(wc_mat, crop_frac_native, method = "bilinear")
  wc_map <- resample(wc_map, crop_frac_native, method = "bilinear")
  
  agg_fact     <- 10
  crop_frac    <- aggregate(crop_frac,    fact = agg_fact, fun = "mean", na.rm = TRUE)
  pasture_frac <- aggregate(pasture_frac, fact = agg_fact, fun = "mean", na.rm = TRUE)
  wc_mat       <- aggregate(wc_mat,       fact = agg_fact, fun = "mean", na.rm = TRUE)
  wc_map       <- aggregate(wc_map,       fact = agg_fact, fun = "mean", na.rm = TRUE)
  
  env_stack        <- c(wc_mat, wc_map)
  names(env_stack) <- c("mat", "map")
  
  # ---------------------------------------------------------------------------
  # 3c. Depth-stratified soil rasters
  #
  # Clay and nitrogen are extracted from all six SoilGrids depth layers.
  # Each field observation already has its clay/nitrogen from the layer
  # matching its depth interval (assigned in 01_soc_data_preparation.R via
  # cut() with breaks c(0,5,15,30,60,100,500)). Keeping all strata and
  # weighting by field record depth proportions gives a global soil
  # distribution comparable to the actual depth distribution the model
  # was trained on.
  # ---------------------------------------------------------------------------
  
  depth_breaks <- c(0, 5, 15, 30, 60, 100, 500)
  depth_labels <- c("0_5cm", "5_15cm", "15_30cm", "30_60cm", "60_100cm", "100_200cm")
  
  assign_stratum <- function(depth_avg) {
    idx <- cut(depth_avg, breaks = depth_breaks,
               labels = FALSE, include.lowest = TRUE)
    depth_labels[idx]
  }
  
  load_soil_rast <- function(var, layer, native_rast, agg_fact) {
    path <- sprintf("./data_inputs/raster/via_gee/%s_%s_mean.tif", var, layer)
    r    <- rast(path)
    r    <- resample(r, native_rast, method = "bilinear")
    r    <- aggregate(r, fact = agg_fact, fun = "mean", na.rm = TRUE)
    r
  }
  
  soil_rasts <- list()
  for (lyr in depth_labels) {
    soil_rasts[[paste0("clay_",     lyr)]] <- load_soil_rast("clay",     lyr, crop_frac_native, agg_fact)
    soil_rasts[[paste0("nitrogen_", lyr)]] <- load_soil_rast("nitrogen", lyr, crop_frac_native, agg_fact)
  }
  
  # ---------------------------------------------------------------------------
  # 3d. Extract climate values weighted by fractional cover
  # ---------------------------------------------------------------------------
  
  extract_weighted <- function(env, frac_rast, label, n_max = 200000) {
    stk <- c(env, frac_rast)
    names(stk)[nlyr(stk)] <- "weight"
    df  <- as.data.frame(mask(stk, frac_rast), na.rm = TRUE)
    if (nrow(df) > n_max) df <- df[sample(nrow(df), n_max), ]
    df$source <- label
    df
  }
  
  df_crop_climate    <- extract_weighted(env_stack, crop_frac,    "Global cropland")
  df_pasture_climate <- extract_weighted(env_stack, pasture_frac, "Global pasture")
  
  # ---------------------------------------------------------------------------
  # 3e. Extract per-stratum soil values weighted by fractional cover
  #
  # Stratum weights from c_dat ensure deeper strata (rare in the field data)
  # contribute proportionally less to the density.
  # ---------------------------------------------------------------------------
  
  stratum_weights <- c_dat %>%
    dplyr::distinct(coord_id, depth_avg) %>%
    dplyr::mutate(stratum = assign_stratum(depth_avg)) %>%
    dplyr::filter(!is.na(stratum)) %>%
    dplyr::count(stratum) %>%
    dplyr::mutate(stratum_wt = n / sum(n))
  
  message("Stratum weights derived from c_dat field records:")
  print(stratum_weights)
  
  extract_soil_weighted <- function(frac_rast, label, n_max = 200000) {
    map_dfr(depth_labels, function(lyr) {
      sw <- stratum_weights$stratum_wt[stratum_weights$stratum == lyr]
      if (length(sw) == 0 || sw == 0) return(NULL)
      clay_r <- soil_rasts[[paste0("clay_",     lyr)]]
      nit_r  <- soil_rasts[[paste0("nitrogen_", lyr)]]
      stk    <- c(clay_r, nit_r, frac_rast)
      names(stk) <- c("clay", "nitrogen", "weight")
      df <- as.data.frame(mask(stk, frac_rast), na.rm = TRUE)
      if (nrow(df) > n_max) df <- df[sample(nrow(df), n_max), ]
      df$stratum <- lyr
      df$source  <- label
      df$weight  <- df$weight * sw
      df
    }) %>% filter(weight > 0)
  }
  
  df_crop_soil    <- extract_soil_weighted(crop_frac,    "Global cropland")
  df_pasture_soil <- extract_soil_weighted(pasture_frac, "Global pasture")
  
  # Combined split data frames (climate rows have NA soil; soil rows have NA climate)
  df_crop <- bind_rows(
    df_crop_climate %>% mutate(clay = NA_real_, nitrogen = NA_real_, stratum = NA_character_),
    df_crop_soil    %>% mutate(mat  = NA_real_, map      = NA_real_)
  )
  df_pasture <- bind_rows(
    df_pasture_climate %>% mutate(clay = NA_real_, nitrogen = NA_real_, stratum = NA_character_),
    df_pasture_soil    %>% mutate(mat  = NA_real_, map      = NA_real_)
  )
  
  # ---------------------------------------------------------------------------
  # 3f. Joint pixel-level data frames (all four variables per pixel)
  #
  # For the all-four-variable joint coverage metric, clay and nitrogen are
  # the stratum-weighted mean across depth layers, giving one representative
  # soil value per pixel. This is joined with climate values on shared pixel
  # coordinates to produce one row per pixel with all four variables valid
  # simultaneously.
  # ---------------------------------------------------------------------------
  
  extract_joint <- function(frac_rast, label) {
    soil_wide <- map(depth_labels, function(lyr) {
      sw <- stratum_weights$stratum_wt[stratum_weights$stratum == lyr]
      if (length(sw) == 0 || sw == 0) return(NULL)
      clay_r <- soil_rasts[[paste0("clay_",     lyr)]]
      nit_r  <- soil_rasts[[paste0("nitrogen_", lyr)]]
      stk    <- c(clay_r, nit_r, frac_rast)
      names(stk) <- c("clay", "nitrogen", "weight")
      df <- as.data.frame(mask(stk, frac_rast), xy = TRUE, na.rm = TRUE)
      df$clay     <- df$clay     * sw
      df$nitrogen <- df$nitrogen * sw
      df[, c("x", "y", "clay", "nitrogen")]
    }) %>%
      bind_rows() %>%
      group_by(x, y) %>%
      summarise(clay = sum(clay), nitrogen = sum(nitrogen), .groups = "drop")
    
    clim_xy <- as.data.frame(mask(c(env_stack, frac_rast), frac_rast),
                             xy = TRUE, na.rm = TRUE)
    names(clim_xy)[ncol(clim_xy)] <- "weight"
    
    inner_join(clim_xy, soil_wide, by = c("x", "y")) %>%
      mutate(source = label)
  }
  
  df_crop_joint    <- extract_joint(crop_frac,    "Global cropland")
  df_pasture_joint <- extract_joint(pasture_frac, "Global pasture")
  
  # Save cache
  save(df_crop_climate, df_pasture_climate,
       df_crop_soil,    df_pasture_soil,
       df_crop,         df_pasture,
       df_crop_joint,   df_pasture_joint,
       stratum_weights, env_stack,
       file = raster_cache)
  message("Cache saved to: ", raster_cache)
  
} # end raster extraction block

# =============================================================================
# Section 4: Prepare field data
# =============================================================================

# MAT/MAP: one record per title_id (study site) — climate does not vary with
# depth and title_id is the finest available site-level key in init_covars.
#
# Clay/nitrogen: all init_id rows are kept so the field soil distribution
# reflects the full range of depths the model was trained on. Each row's
# clay/nitrogen is already matched to the correct SoilGrids depth layer by
# 01_soc_data_preparation.R.

depth_breaks <- c(0, 5, 15, 30, 60, 100, 500)
depth_labels <- c("0_5cm", "5_15cm", "15_30cm", "30_60cm", "60_100cm", "100_200cm")
assign_stratum <- function(depth_avg) {
  idx <- cut(depth_avg, breaks = depth_breaks, labels = FALSE, include.lowest = TRUE)
  depth_labels[idx]
}

df_field_climate <- init_covars %>%
  rename(mat = MAT, map = MAP) %>%
  distinct(title_id, lu_start, .keep_all = TRUE) %>%
  dplyr::select(title_id, lu_start, mat, map) %>%
  drop_na() %>%
  mutate(weight = 1,
         clay = NA_real_, nitrogen = NA_real_, stratum = NA_character_)

df_field_soil <- init_covars %>%
  dplyr::select(init_id, depth_avg, lu_start, clay, nitrogen) %>%
  drop_na() %>%
  mutate(weight  = 1,
         stratum = assign_stratum(depth_avg),
         mat = NA_real_, map = NA_real_)

message(sprintf(
  "Climate records (one per title_id): %d cropland, %d pasture",
  sum(df_field_climate$lu_start == "C"),
  sum(df_field_climate$lu_start == "P")))
message(sprintf(
  "Soil records (all depth intervals): %d cropland, %d pasture",
  sum(df_field_soil$lu_start == "C"),
  sum(df_field_soil$lu_start == "P")))

df_field_crop <- bind_rows(
  df_field_climate %>% filter(lu_start == "C"),
  df_field_soil    %>% filter(lu_start == "C")
) %>% mutate(source = "Cropland field sites")

df_field_pasture <- bind_rows(
  df_field_climate %>% filter(lu_start == "P"),
  df_field_soil    %>% filter(lu_start == "P")
) %>% mutate(source = "Pasture field sites")

# Joint field data frames: one row per init_id with all four variables valid.
# MAT/MAP are present on every row of init_covars (they do not vary with
# depth), so they are extracted directly rather than via a separate join.
df_field_joint_base <- init_covars %>%
  rename(mat = MAT, map = MAP) %>%
  dplyr::select(init_id, depth_avg, lu_start, mat, map, clay, nitrogen) %>%
  drop_na() %>%
  mutate(weight = 1)

df_field_crop_joint <- df_field_joint_base %>%
  filter(lu_start == "C") %>%
  mutate(source = "Cropland field sites")

df_field_pasture_joint <- df_field_joint_base %>%
  filter(lu_start == "P") %>%
  mutate(source = "Pasture field sites")

message(sprintf(
  "Joint field records (all four vars): %d cropland, %d pasture",
  nrow(df_field_crop_joint),
  nrow(df_field_pasture_joint)))

# =============================================================================
# Section 5: Build density curves
# =============================================================================

make_density_df <- function(df_raster, df_field, raster_label, field_label) {
  items <- list(list(d = df_raster, lbl = raster_label),
                list(d = df_field,  lbl = field_label))
  map_dfr(vars, function(v) {
    map_dfr(items, function(item) {
      vals    <- item$d[[v]];  wts <- item$d$weight
      ok      <- !is.na(vals) & !is.na(wts)
      src_rng <- range(vals[ok])
      kd <- density(vals[ok],
                    weights = wts[ok] / sum(wts[ok]),
                    from = src_rng[1], to = src_rng[2],
                    cut = 0, n = 512)
      tibble(variable = var_labels[v],
             source   = item$lbl,
             x = kd$x, y = kd$y)
    })
  })
}

density_crop <- make_density_df(df_crop,    df_field_crop,
                                "Global cropland",    "Cropland field sites")
density_past <- make_density_df(df_pasture, df_field_pasture,
                                "Global pasture",     "Pasture field sites")

# =============================================================================
# Section 6: Combined density figure
#
# Cropland and pasture panels are stacked into a single figure for the
# Extended Data. Each row shows one land-use type; columns show the four
# predictor variables. Horizontal bars mark the 5–95% range of each
# distribution.
# =============================================================================

colours_crop <- c("Global cropland"      = "#d95f02",
                  "Cropland field sites" = "#1b9e77")
colours_past <- c("Global pasture"       = "#7570b3",
                  "Pasture field sites"  = "#1b9e77")

# -----------------------------------------------------------------------------
# 6a. 5–95% range bars
# -----------------------------------------------------------------------------

make_range_bars <- function(df_raster, raster_label, df_field, field_label) {
  items <- list(list(d = df_raster, lbl = raster_label, yoff = -0.02),
                list(d = df_field,  lbl = field_label,  yoff = -0.04))
  map_dfr(vars, function(v) {
    map_dfr(items, function(item) {
      vals <- item$d[[v]];  wts <- item$d$weight
      ok   <- !is.na(vals) & !is.na(wts)
      q    <- wquantile(vals[ok], wts[ok], c(0.05, 0.95))
      tibble(variable = factor(var_labels[v], levels = var_labels),
             source   = item$lbl,
             xmin = q[1], xmax = q[2],
             yoff = item$yoff)
    })
  })
}

# -----------------------------------------------------------------------------
# 6b. Single-panel plot function (one land-use type)
# -----------------------------------------------------------------------------

make_plot <- function(density_df, range_df, colour_scale, subtitle_txt) {
  density_df <- density_df %>%
    mutate(variable = factor(variable, levels = var_labels),
           source   = factor(source,   levels = names(colour_scale)))
  range_df <- range_df %>%
    mutate(source = factor(source, levels = names(colour_scale)))
  
  ymax_df <- density_df %>%
    group_by(variable) %>%
    summarise(ymax = max(y), .groups = "drop")
  
  range_df <- range_df %>%
    left_join(ymax_df, by = "variable") %>%
    mutate(y_bar = yoff * ymax)
  
  ggplot() +
    geom_area(data = density_df,
              aes(x = x, y = y, colour = source, fill = source),
              alpha = 0.20, position = "identity") +
    geom_line(data = density_df,
              aes(x = x, y = y, colour = source),
              linewidth = 0.8) +
    geom_segment(data = range_df,
                 aes(x = xmin, xend = xmax, y = y_bar, yend = y_bar,
                     colour = source),
                 linewidth = 1.5, lineend = "round") +
    geom_segment(data = range_df,
                 aes(x = xmin, xend = xmin,
                     y = y_bar - 0.004 * ymax,
                     yend = y_bar + 0.004 * ymax,
                     colour = source),
                 linewidth = 0.8) +
    geom_segment(data = range_df,
                 aes(x = xmax, xend = xmax,
                     y = y_bar - 0.004 * ymax,
                     yend = y_bar + 0.004 * ymax,
                     colour = source),
                 linewidth = 0.8) +
    facet_wrap(~ variable, scales = "free", ncol = 2) +
    scale_colour_manual(values = colour_scale, name = NULL) +
    scale_fill_manual(  values = colour_scale, name = NULL) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = "Density", subtitle = subtitle_txt) +
    theme_bw(base_size = 11) +
    theme(legend.position  = "bottom",
          strip.background = element_rect(fill = "grey92"),
          panel.grid.minor = element_blank())
}

# -----------------------------------------------------------------------------
# 6c. Build panels and combine into a single figure
# -----------------------------------------------------------------------------

ranges_crop <- make_range_bars(df_crop,    "Global cropland",
                               df_field_crop,    "Cropland field sites")
ranges_past <- make_range_bars(df_pasture, "Global pasture",
                               df_field_pasture, "Pasture field sites")

p_crop <- make_plot(density_crop, ranges_crop, colours_crop,
                    "Cropland: global raster (weighted by fractional cover) vs field sites")
p_past <- make_plot(density_past, ranges_past, colours_past,
                    "Pasture: global raster (weighted by fractional cover) vs field sites")

p_combined <- cowplot::plot_grid(p_crop, p_past, ncol = 1, labels = c("a", "b"))

ggsave("./figures/env_density_combined.pdf", p_combined, width = 8, height = 11)
ggsave("./figures/env_density_combined.png", p_combined, width = 8, height = 11, dpi = 300)
message("Combined density figure written to ./figures/env_density_combined.pdf/.png")

# =============================================================================
# Section 7: Joint (multivariate) coverage table
#
# Coverage is computed at four percentile thresholds:
#   Full range (0–100%):  has any field observation ever seen this condition?
#   99%  (0.5–99.5%):     full range minus the most extreme 0.5% tails
#   95%  (2.5–97.5%):     conventional central interval
#   90%  (5–95%):         more conservative central interval
#
# Results span: each variable individually, climate jointly (MAT + MAP),
# soil jointly (clay + nitrogen), and all four variables jointly.
# =============================================================================

coverage_intervals <- list(
  "Full range (0-100%)" = c(0,     1    ),
  "99% (0.5-99.5%)"     = c(0.005, 0.995),
  "95% (2.5-97.5%)"     = c(0.025, 0.975),
  "90% (5-95%)"         = c(0.05,  0.95 )
)

global_contexts <- list(
  list(label      = "Global cropland",
       df_rast    = df_crop,
       df_field   = df_field_crop,
       df_rast_j  = df_crop_joint,
       df_field_j = df_field_crop_joint),
  list(label      = "Global pasture",
       df_rast    = df_pasture,
       df_field   = df_field_pasture,
       df_rast_j  = df_pasture_joint,
       df_field_j = df_field_pasture_joint)
)

coverage_table <- map_dfr(names(coverage_intervals), function(interval_label) {
  probs <- coverage_intervals[[interval_label]]
  
  map_dfr(global_contexts, function(ctx) {
    cov <- joint_coverage(ctx$df_rast,  ctx$df_field,  ctx$label,
                          df_rast_joint  = ctx$df_rast_j,
                          df_field_joint = ctx$df_field_j,
                          probs = probs)
    
    # Univariate rows
    univ <- cov$univariate %>%
      rename(coverage_pct = pct_in_env) %>%
      mutate(coverage_type = "Univariate",
             metric        = variable)
    
    # Group rows
    grp <- cov$group_coverage %>%
      mutate(coverage_pct  = pct_joint,
             coverage_type = "Joint",
             metric        = case_when(
               group == "climate" ~ "Climate jointly (MAT + MAP)",
               group == "soil"    ~ "Soil jointly (clay + nitrogen)"
             )) %>%
      dplyr::select(coverage_type, metric, coverage_pct)
    
    # All-four row
    all4 <- tibble(
      coverage_type = "Joint",
      metric        = "All four variables jointly",
      coverage_pct  = cov$pct_all_four
    )
    
    bind_rows(univ %>% dplyr::select(coverage_type, metric, coverage_pct),
              grp, all4) %>%
      mutate(context        = ctx$label,
             interval_label = interval_label)
  })
}) %>%
  dplyr::select(context, interval_label, coverage_type, metric, coverage_pct) %>%
  arrange(context, coverage_type, metric, interval_label)

# Wide format: one row per context × metric, intervals as columns
coverage_wide <- coverage_table %>%
  pivot_wider(names_from  = interval_label,
              values_from = coverage_pct) %>%
  arrange(context, coverage_type, metric)

print(coverage_wide, n = Inf)

write.csv(coverage_wide,
          "./figures/env_joint_coverage_table.csv",
          row.names = FALSE)
message("Joint coverage table written to ./figures/env_joint_coverage_table.csv")