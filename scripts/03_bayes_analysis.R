# =============================================================================
# 03_bayes_analysis.R
#
# PURPOSE
# -------
# Fits the hierarchical Bayesian model of SOC accrual following cropland or
# pasture abandonment, evaluates convergence, runs three pre-planned
# sensitivity analyses, and produces all manuscript figures.
#
# All data preparation (cleaning, covariate extraction, BD gap-filling, unit
# harmonisation) is handled upstream by 01_soc_data_preparation.R.  This
# script reads the analysis-ready CSV produced there and proceeds directly
# to the JAGS model.
#
# SCRIPT DEPENDENCIES (must be run first)
# ----------------------------------------
#   01_soc_data_preparation.R  — produces c_dat_analysis_ready.csv and
#                                sc_df.csv
#   02_lc_priors.R             — produces no_lu_coefs.RData,
#                                cropland_coefs.RData, pasture_coefs.RData
#
# INPUTS
# ------
#   ./data_inputs/c_dat_analysis_ready.csv
#   ./data_inputs/sc_df.csv
#   ./data_output/no_lu_coefs.RData
#   ./data_output/cropland_coefs.RData
#   ./data_output/pasture_coefs.RData
#
# OUTPUTS
# -------
#   ./data_output/bayes_output.RData          — primary analysis MCMC samples
#   ./data_output/bayes_output_bd_sens.RData  — BD sensitivity MCMC samples
#   ./data_output/bayes_output_mgmt_sens.RData — MGMT sensitivity MCMC samples
#   ./data_output/bayes_output_som_sens.RData  — SOM sensitivity MCMC samples
#   ./data_output/bayes_output_weakprior.RData  — weak-prior sensitivity MCMC samples
#   ./data_output/bayes_output_origbd_sens.RData — original-BD sensitivity MCMC samples
#   ./figures/  (all figures described in section headers below)
#
# =============================================================================
# SENSITIVITY ANALYSIS SUBSETS
# =============================================================================
# Four sensitivity analyses evaluate the robustness of the primary results
# to specific data preparation decisions.  Each re-fits the identical JAGS
# model on a subset of c_dat defined by a boolean provenance column added
# during data preparation.  sc_df from the primary (full) dataset is used
# throughout, so all runs are on the same covariate scale and their
# posteriors are directly comparable.
#
# (a) BULK DENSITY GAP-FILLING
#     Primary:     c_dat (all rows)
#     Sensitivity: filter(c_dat, !bd_gap_filled)
#     Rationale:   Tests whether rows where BD was model-predicted rather than
#                  published materially affect the posterior effect sizes.
#
# (b) MANAGEMENT UNCERTAINTY
#     Primary:     c_dat (records with uncertain "?" management assigned to "N")
#     Sensitivity: filter(c_dat, !MGMT_is_uncertain)
#     Rationale:   Tests whether the conservative coding of ambiguous
#                  management records influences the recovery trajectories.
#
# (c) SOM-TO-SOC CONVERSION
#     Primary:     c_dat (SOM records converted via van Bemmelen factor 0.5)
#     Sensitivity: filter(c_dat, !som_converted)
#     Rationale:   Tests whether the fixed conversion factor introduces
#                  systematic bias in the recovered trajectories.
#
# (d) ORIGINAL BULK DENSITY
#     Primary:     c_dat (problematic BD replaced by model-predicted values)
#     Sensitivity: full c_dat, but c_init_obs / c_final_obs computed from
#                  Ctrl/Aban_c_dens_av_orig_bd (original reported BD, including
#                  repeated-value entries); rows where BD_raw is NA are dropped
#     Rationale:   Tests whether replacing implausible constant-BD profiles
#                  with model-predicted values drives the primary conclusions,
#                  versus simply trusting the published values as reported.
#
# Each sensitivity run produces a saved .RData file with the same structure
# as the primary output.  The sensitivity comparison figure at the end of
# the script overlays the effect size posteriors from all five runs using
# the same two-panel MCMCplot layout as the primary "effect sizes" figure.
# =============================================================================


# ── Packages ──────────────────────────────────────────────────────────────────
# ipak() installs and loads packages in one call
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(c("tidyverse", "rjags", "MCMCvis"))


# =============================================================================
# USER FLAG — REFIT MODELS
# =============================================================================
# Set refit_models <- TRUE  to force all JAGS models (primary + all
# sensitivity runs) to be re-estimated even if saved .RData files already
# exist on disk.
#
# Set refit_models <- FALSE (the default) to load cached .RData files
# whenever they are found, skipping the expensive jags.model / update /
# coda.samples calls.  If a cached file is NOT found for a given run,
# that run will be executed regardless of this flag.

refit_models <- FALSE


# =============================================================================
# SECTION 1 — READ PREPARED DATA
# =============================================================================
# All cleaning, covariate extraction, and BD gap-filling has been performed by
# 01_soc_data_preparation.R.  The CSV read here is fully analysis-ready.

c_dat <- read.csv("./data_inputs/c_dat_analysis_ready.csv",
                  stringsAsFactors = FALSE)

# sc_df provides centering/scaling parameters for all covariates.
# It was computed from the full unfiltered dataset and is used for all runs
# (primary and sensitivity) to keep posteriors on the same scale.
sc_df <- read.csv("./data_inputs/sc_df.csv",
                  row.names = 1, stringsAsFactors = FALSE)

message("c_dat loaded: ", nrow(c_dat), " rows")
message("Sensitivity subset sizes:")
message("  (a) !bd_gap_filled     : ",
        sum(!c_dat$bd_gap_filled,     na.rm = TRUE), " rows")
message("  (b) !MGMT_is_uncertain : ",
        sum(!c_dat$MGMT_is_uncertain, na.rm = TRUE), " rows")
message("  (c) !som_converted     : ",
        sum(!c_dat$som_converted,     na.rm = TRUE), " rows")



# =============================================================================
# SECTION 2 — SPATIAL DISTRIBUTION FIGURE
# =============================================================================
# Map of study sites and latitudinal density plots by land-use category.
# c_sf_short is built from c_dat (which already has LAT, LONG, lu_start,
# lu_end, MGMT from 01_soc_data_preparation.R) rather than from a spatial
# extraction object.
ipak(c("sf","ggplot2", "rnaturalearth", "rnaturalearthdata"))

projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

c_sf_short <- c_dat %>%
  dplyr::select(lu_start, lu_end, MGMT, LAT, LONG) %>%
  distinct() %>%
  st_as_sf(coords = c("LONG", "LAT"), crs = projcrs) %>%
  mutate(abs_lat = abs(st_coordinates(.)[, 2]))

world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(subregion != "Antarctica")

# Version of the map with overlapping plot locations
# p0 <- ggplot() +
#   geom_sf(data = world, fill = "gray90", color = NA) +
#   geom_sf(data = c_sf_short, stroke = F, size = 2, alpha = 0.4, color = "grey60") +
#   
#   # remove 'buffer'
#   coord_sf(expand = FALSE) +
#   
#   # minimal theme, but strip every background element
#   theme_minimal() +
#   theme(
#     # completely blank out background elements
#     panel.background = element_blank(),
#     plot.background  = element_blank(),
#     
#     # no grids or axis stuff
#     panel.grid       = element_blank(),
#     axis.title       = element_blank(),
#     axis.text        = element_blank(),
#     axis.ticks       = element_blank(),
#     
#     # zero out margins between panel and plot edge
#     plot.margin      = margin(0, 0, 0, 0, unit = "cm"),
#     
#     # if ever faceting, also zero spacing
#     panel.spacing    = unit(0.5, "cm")
#   )


# Version with squares that show the counts of cropland/pasture plots
grid   <- st_make_grid(world, cellsize = c(5,5), square = TRUE)
grid_sf <- st_sf(grid_id = seq_along(grid), geometry = grid)
pts_grid <- st_join(c_sf_short, grid_sf, join = st_within)
counts <- pts_grid %>%
  st_drop_geometry() %>%
  count(grid_id, lu_start) %>%
  pivot_wider(names_from = lu_start, values_from = n, values_fill = 0) %>%
  rename(C = C, P = P)
grid_counts <- grid_sf %>%
  left_join(counts, by = "grid_id") %>%
  replace_na(list(C = 0, P = 0)) %>%
  mutate(
    total  = C + P,
    prop_C = if_else(total > 0, C / total, NA_real_)
  )
centroids <- grid_counts %>%
  st_centroid() %>%
  mutate(
    lon = st_coordinates(.)[,1],
    lat = st_coordinates(.)[,2]
  )
p0 <- ggplot() +
  # base map
  geom_sf(data = world, fill = "gray90", color = NA) +
  
  # cropland squares
  geom_point(
    data    = centroids %>% filter(C > 0),
    aes(x = lon, y = lat, size = C/3),
    shape   = 15,
    color   = "#d95f02",
    alpha   = 0.5
  ) +
  # pasture squares
  geom_point(
    data    = centroids %>% filter(P > 0),
    aes(x = lon, y = lat, size = P/3),
    shape   = 15,
    color   = "#7570b3",
    alpha   = 0.5
  ) +
  
  # log‐scale sizing
  scale_size_continuous(
    name   = "Plots per cell",
    trans  = "log10",
    range  = c(0.7, 1.5) #
  ) +
  
  coord_sf(expand = FALSE) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

dens_cols <- c(
  "#d95f02",
  "#7570b3",
  "#1b9e77",
  "#e7298a",
  "#66a61e",
  "#e6ab02")

p1 <- ggplot(c_sf_short, aes(x = abs_lat, fill = lu_start, color = lu_start)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(breaks = c(0, 30, 60)) +
  scale_y_continuous(breaks = c(0, 0.025, 0.05), limits = c(0, 0.07)) +
  scale_fill_manual(
    name = "Agricultural\nland use",
    values = c("C" = dens_cols[1], "P" = dens_cols[2]),  # Adjust colors as needed
    labels = c("C" = "Cropland", "P" = "Pasture")
  ) +
  scale_color_manual(
    name = "Agricultural\nland use",
    values = c("C" = dens_cols[1], "P" = dens_cols[2]),  # Ensure colors match fill
    labels = c("C" = "Cropland", "P" = "Pasture")
  ) +
  labs(
    x = "",
    y = "Site\ndensity"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )



p2 <- ggplot(c_sf_short, aes(x = abs_lat, fill = lu_end, color = lu_end)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(breaks = c(0, 30, 60)) +
  scale_y_continuous(breaks = c(0, 0.025, 0.05), limits = c(0, 0.07)) +
  
  scale_fill_manual(
    name = "Restored\nlandcover", 
    values = c("F" = dens_cols[3], "O" = dens_cols[4]),  # Adjust colors as needed
    labels = c("F" = "Forest", "O" = "Open")
  ) +
  scale_color_manual(
    name = "Restored\nlandcover", 
    values = c("F" = dens_cols[3], "O" = dens_cols[4]),  # Adjust colors as needed
    labels = c("F" = "Forest", "O" = "Open")
  ) +
  labs(
    x = "",
    y = "Site\ndensity"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

p3 <- ggplot(c_sf_short, aes(x = abs_lat, fill = MGMT, color = MGMT)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(breaks = c(0, 30, 60)) +
  scale_y_continuous(breaks = c(0, 0.025, 0.05), limits = c(0, 0.07)) +
  scale_fill_manual(
    name = "Restoration\nstrategy", 
    values = c("N" = dens_cols[5], "A" = dens_cols[6]),  # Adjust colors as needed
    labels = c("N" = "Natural", "A" = "Active")
  ) +
  scale_color_manual(
    name = "Restoration\nstrategy",
    values = c("N" = dens_cols[5], "A" = dens_cols[6]),  # Adjust colors as needed
    labels = c("N" = "Natural", "A" = "Active")
  ) +
  labs(
    x = "Absolute Latitude",
    y = "Site\ndensity"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

# Now get a scatter plot of raw data
ipak(c("zoo", "boot")) # for rollapply and boot()

# window size and # of bootstrap replicates
window_size <- 75 # 25 and sd = 2 looks okay
nboot       <- 2000  

# helper: compute moving‐avg for one bootstrap draw of the data frame
movfun <- function(dat, idx, window_size) {
  d <- dat[idx, ]
  d$AGE <- d$AGE + rnorm(nrow(d), mean = 0, sd = 1)
  d <- d[order(d$AGE), ]
  resp <- d$Aban_c_dens_av / d$Ctrl_c_dens_av
  # returns a vector of length nrow(d), padded with NA at ends
  zoo::rollapply(resp, width = window_size, FUN = mean,
                 fill = NA, align = "center")
}

# split the data
crop_df    <- c_dat[c_dat$lu_start=="C", ]
pasture_df <- c_dat[c_dat$lu_start=="P", ]

# ages in sorted order (they’ll be identical positions to the mov‐avgs)
age_crop    <- sort(crop_df$AGE)
age_pasture <- sort(pasture_df$AGE)

set.seed(101)  # for reproducibility
boot_crop    <- boot(crop_df,    statistic = function(d,i) movfun(d, i, window_size), R = nboot)
boot_pasture <- boot(pasture_df, statistic = function(d,i) movfun(d, i, window_size), R = nboot)

# extract the bootstrap‐mean curves (and 95% CI if desired)
crop_mean    <- apply(boot_crop$t,    2, mean,   na.rm=TRUE)
crop_ci      <- apply(boot_crop$t,    2, quantile, probs = c(0.025, 0.975), na.rm=TRUE)

pasture_mean <- apply(boot_pasture$t, 2, mean,   na.rm=TRUE)
pasture_ci   <- apply(boot_pasture$t, 2, quantile, probs = c(0.025, 0.975), na.rm=TRUE)


p4 <- ggplot() +
  # POINTS: colored by lu_start, no legend
  geom_point(
    data = c_dat,
    aes(
      x = AGE,
      y = pmax(0.2, pmin(8, Aban_c_dens_av / Ctrl_c_dens_av)),
      color = lu_start
    ),
    alpha = 0.3,
    stroke = F,
    size = 2,
    show.legend = FALSE
  ) +
  
  # LINES: colored by group, with legend
  geom_line(
    data = data.frame(AGE = age_pasture, mean = pasture_mean, Group = "Pasture"),
    aes(x = AGE, y = mean, color = Group),
    size = 1
  ) +
  geom_line(
    data = data.frame(AGE = age_crop, mean = crop_mean, Group = "Cropland"),
    aes(x = AGE, y = mean, color = Group),
    size = 1
  ) +
  
  
  # COMBINED scale for both point colors (light tones) and line legend colors (darker)
  scale_color_manual(
    values = c(
      "C" = "#f4c29f",        # light cropland points
      "P" = "#c8c6e1",        # light pasture points
      "Cropland" = "#d95f02", # orange line
      "Pasture"  = "#7570b3"  # green line
    ),
    breaks = c("Cropland", "Pasture"),  # only show these in legend
    name = NULL
  ) +
  
  scale_y_log10(
    breaks = c(0.25, 0.5, 1, 2, 4, 8),
    labels = c("-75%", "-50%", "No change", "+100%", "+400%", ">+800%")
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  
  labs(
    x = "Years since restoration",
    y = "Percent change in SOC\n(Restored vs agriculture)"
  ) +
  theme_classic() +
  theme(
    legend.position = c(0.99, 0.01),
    legend.justification = c(1, 0),
    legend.direction = "vertical",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 10)
  )


ipak("patchwork")
# Some very detailed tweaking
p1 <- p1 + 
  theme(
    legend.key.size   = unit(0.5, "cm"),  # overall key size
    legend.key.height = unit(0.3, "cm"),  # optional fine control
    legend.key.width  = unit(0.3, "cm")
  )

p2 <- p2 + 
  theme(
    legend.key.size   = unit(0.5, "cm"),  # overall key size
    legend.key.height = unit(0.3, "cm"),  # optional fine control
    legend.key.width  = unit(0.3, "cm")
  )

p3 <- p3 + 
  theme(
    legend.key.size   = unit(0.5, "cm"),  # overall key size
    legend.key.height = unit(0.3, "cm"),  # optional fine control
    legend.key.width  = unit(0.3, "cm")
  )

p4a <- p4 +
  theme(
    plot.margin = margin(0, 0.5, 0, 0, "cm")
  )


p0  <- p0  + labs(tag = "A") + theme(plot.tag.position = c(0, 1))
p4a <- p4a + labs(tag = "B") + theme(plot.tag.position = c(-0.0375, 1))
p1  <- p1  + labs(tag = "C") + theme(plot.tag.position = c(0, 1))
p2  <- p2  + labs(tag = "D") + theme(plot.tag.position = c(0, 1))
p3  <- p3  + labs(tag = "E") + theme(plot.tag.position = c(0, 1))


# Misc additional tweaks
# 
# p4b <- p4a +
#   theme(
#     legend.position = c(0.5, 1.02),          # center-top, just above plot
#     legend.justification = c(0.5, 1),        # align bottom-center of legend box
#     legend.direction = "horizontal",
#     legend.background = element_blank(),
#     legend.box.background = element_blank(),
#     legend.key = element_blank(),
#     legend.margin = margin(0, 0, 0, 0),
#     plot.margin = margin(15, 0, 0, 0)        # extra top margin so legend fits
#   )




ipak("cowplot")

p0a <- p0 +
  annotate("point", x = -162, y = -20, size = 0.7, shape = 15, color = "grey") +
  annotate("text",  x = -156, y = -20, label = "1 site", hjust = 0, size = 3) +
  annotate("point", x = -162, y = -35, size = 1.51, shape = 15, color = "grey") +
  annotate("text",  x = -156, y = -35, label = "50 sites", hjust = 0, size = 3)

# # This is a histogram of chronosequence lengths, but it's really
# # busy here.
# chronoseq_lengths <- c_dat %>%
#   group_by(Profile_ID) %>%
#   summarise(length = max(AGE, na.rm = TRUE)) %>%
#   ungroup()
# ipak("cowplot")
# p_hist <- ggplot(chronoseq_lengths, aes(x = length)) +
#   geom_histogram(binwidth = 10, boundary = 0, fill = "gray", color = NA) +
#   scale_x_continuous(
#     breaks = seq(0, max(chronoseq_lengths$length), 50),
#     expand = c(0, 0)
#   ) +
#   scale_y_continuous(expand = c(0, 0)) +
#   coord_cartesian(clip = "off") +
#   labs(x = "Chronosequence length (y)", y = NULL) +
#   theme_minimal(base_size = 8) +
#   theme(
#     axis.title.y    = element_blank(),
#     axis.text.y     = element_blank(),
#     axis.ticks.y    = element_blank(),
#     panel.grid      = element_blank(),
#     panel.background = element_blank(),
#     plot.background  = element_rect(fill = NA, color = NA),  # <- fix here
#     axis.title.x    = element_text(size = 7, margin = margin(t = 2)),
#     axis.text.x     = element_text(size = 6, margin = margin(t = 0)),
#     plot.margin     = margin(0, 0, 0, 0)
#   )
# p0a <- ggdraw(p0a) +
#   draw_plot(p_hist, x = 0.45, y = 0, width = 0.3, height = 0.4)

# build the top row (just the map full width)
top_row <- p0a



# build the bottom row as [ blank | p4b ]
bottom_row <- (plot_spacer() | p4a) + 
  plot_layout(
    heights = c(1, 1),      # map twice as tall as the bottom row
    widths  = c(0, 1)      # bottom row: spacer is 0.1 of p4a’s width
  )

# now stack them, giving the map more height than the bottom row
left_col <- (top_row / bottom_row) + 
  plot_layout(
    heights = c(1, 1),      # map twice as tall as the bottom row
    widths  = c(1, 0)      # bottom row: spacer is 0.1 of p4a’s width
  )

# right column as before
right_col <- (p1 / p2 / p3)

# combine left & right
combined_plot <- (left_col | right_col) + 
  plot_layout(widths = c(6,1))


png("./figures/spatial_distribution.png",
    width  = 6.5, height = 4,
    unit   = "in", res = 440)
combined_plot
dev.off()

pdf("./figures/spatial_distribution.pdf",
    width  = 6.5, height = 4)
combined_plot
dev.off()




# =============================================================================
# SECTION 3 — PREPARE DATA FOR JAGS
# =============================================================================
# Reshapes c_dat into the init_dat / init_covars / final_dat / final_covars
# structure expected by the JAGS model.
#
# The core challenge this reshaping solves: the dataset was constructed so that
# each row represents a paired comparison (one reference site, one restored
# site).  When one reference site is compared to multiple restored sites, that
# reference site's measurements appear multiple times.  The model addresses
# this by treating reference and restored carbon as separate latent processes
# linked by shared IDs, so each measurement appears exactly once.

# Remove rows missing carbon density estimates or any spatial covariate.
# JAGS requires all data nodes to be fully defined; dropping here rather than
# imputing keeps the analysis honest about which observations are actually used.
# The warning messages from 01_soc_data_preparation.R identify which
# coordinates had extraction failures — those are the rows removed here.
n_before <- nrow(c_dat)
c_dat <- c_dat %>%
  filter(!is.na(Aban_c_dens_av), !is.na(Ctrl_c_dens_av),
         !is.na(clay), !is.na(nitrogen), !is.na(irrigation))
n_dropped <- n_before - nrow(c_dat)
if (n_dropped > 0) {
  message(n_dropped, " row(s) dropped for missing covariates or carbon density")
}



# Reshape dataset to accommodate analysis. 
# Essentially the issue is that the dataset was built while repeating measurements
# for comparison. We want a single measurement to appear in the dataset
# only once. As an example, if there was 1 crop site and 2 abandoned sites,
# currently the crop site's measurements are repeated twice for each of the two 
# rows that represent each crop-abandoned comparison.
# This data repetition problem is corrected by modeling each crop and abandoned
# site individually, and using the architecture of the model to enable the 
# direct comparison between the two.

# First, develop a comparison ID. This enables direct comparisons between
# the cropland and abandoned sites 
c_dat$init_id <- c_dat %>% 
  dplyr::select(coord_id, depth_avg, Authors, Year) %>% 
  apply(1, function(x) paste(x, collapse = "")) %>% 
  as.factor() %>% 
  as.numeric()

# Next, get an identifier representing the combination of coordinates, study,
# depth, age, and end land use. This is for when there are multiple observations of 
# abandoned carbon or multiple types of end land use.
c_dat$final_id <- c_dat %>% 
  dplyr::select(coord_id, depth_avg, Authors, Year, AGE, lu_end, MGMT) %>% 
  apply(1, function(x) paste(x, collapse = "")) %>% 
  as.factor() %>% 
  as.numeric()

# Third, get an identifier for the title
c_dat <- c_dat %>% 
  mutate(title_id = as.numeric(as.factor(Title)))

# To make the factor levels easier for jags, will sort by final_id
c_dat <- c_dat[order(c_dat$final_id),]

# Make a dataframe just for the cropland records
init_dat <- c_dat %>% 
  dplyr::select(-starts_with("Aban"),
                -contains("accrual"),
                -contains("Perc_change"),
                -contains("Cum_Aban"),
                -starts_with("SE"),
                -starts_with("SD"),
                -obs,
                -starts_with("AGE"),
                #-lu_end,
                -final_id)

# Only want the unique data values
init_dat <- init_dat[!duplicated(init_dat %>% dplyr::select(-Profile_ID, -lu_end, -lu_start)),]
dim(init_dat)

# Rename to init profile ID and re-factor them
init_dat <- init_dat %>% 
  mutate(init_profile_id = as.numeric(as.factor(Profile_ID)))


# Make a dataframe just for the abandoned records
final_dat <- c_dat %>% 
  dplyr::select(-starts_with("Ctrl"),
                -contains("accrual"),
                -contains("Perc_change"),
                -contains("Cum_Ctrl"),
                -starts_with("SE"),
                -starts_with("SD"),
                -obs
  )

final_dat <- final_dat[!duplicated(final_dat %>% dplyr::select(-Profile_ID)),]
dim(final_dat)

# Rename to final profile ID and re-factor the them
final_dat <- final_dat %>% 
  mutate(final_profile_id = as.numeric(as.factor(Profile_ID)))


# Need to make a tidy initial and final dataframe to store the covariate values.
# The number of rows is the number of unique values of init_id and
# final_id, respectively.

# In the very few cases where there are multiple MAP/MAT values for a site, we
# will take the average.
init_covars <- init_dat %>% 
  dplyr::select(init_id, depth_avg, MAP, MAT, 
                clay, 
                nitrogen, 
                irrigation,
                title_id, 
                lu_start, lu_end) %>% 
  group_by(init_id)  %>% 
  summarize(across(.cols = where(is.numeric),.fns = mean),
            across(.cols = where(is.character),.fns = first))

final_covars <- final_dat %>% 
  dplyr::rename(init_compare_id = init_id) %>% 
  dplyr::select(final_id, init_compare_id, depth_avg, MAP, MAT, 
                clay, 
                nitrogen, 
                irrigation,
                AGE, 
                title_id, coord_id,
                lu_start, lu_end, MGMT) %>% 
  group_by(final_id) %>% 
  summarize(across(.cols = where(is.numeric),.fns = mean),
            across(.cols = where(is.character),.fns = first))


# sc_df was computed from the full dataset and written by
# 01_soc_data_preparation.R.  Reading it here (rather than recomputing)
# ensures that sensitivity analysis subsets are centred and scaled using
# population-level parameters, keeping all runs on the same covariate scale
# and making their posteriors directly comparable.
#
# The sc_df object read at the top of this script is already available;
# this comment marks where the centering/scaling parameters are first used.

# Load prior information from previous estimates of natural soil carbon from
# Sanderman et al. 2017 PNAS and soil carbon in cropland and pasture from SoilGrids
load("./data_output/no_lu_coefs.RData")
load("./data_output/cropland_coefs.RData")
load("./data_output/pasture_coefs.RData")

jags_data <- list(
  # Data
  c_final_obs = final_dat$Aban_c_dens_av %>% log(),
  c_init_obs = init_dat$Ctrl_c_dens_av %>% log(),
  
  # Covariates for the initial data
  init_depth = (log(init_covars$depth_avg) - sc_df["depth", "center"]) / sc_df["depth", "scale"],
  init_age = rep(0, nrow(init_covars)), 
  
  init_map = (log(init_covars$MAP) - sc_df["map", "center"]) / sc_df["map", "scale"],
  init_mat = (init_covars$MAT - sc_df["mat", "center"]) / sc_df["mat", "scale"], # This is closer to normal (also contains negative numbers)
  
  init_clay = (init_covars$clay - sc_df["clay", "center"]) / sc_df["clay", "scale"], # Won't log scale as these are relatively normal already
  
  init_nitrogen = (log(init_covars$nitrogen) - sc_df["nitrogen", "center"]) / sc_df["nitrogen", "scale"],
  
  init_irrigation = (log(init_covars$irrigation) - sc_df["irrigation", "center"]) / sc_df["irrigation", "scale"],
  
  
  # Covariates for the final data
  final_depth = (log(final_covars$depth_avg) - sc_df["depth", "center"]) / sc_df["depth", "scale"], 
  final_age = final_covars$AGE,
  
  final_map = (log(final_covars$MAP) - sc_df["map", "center"]) / sc_df["map", "scale"],
  final_mat = (final_covars$MAT - sc_df["mat", "center"]) / sc_df["mat", "scale"], # This is closer to normal (also contains negative numbers)
  
  final_clay = (final_covars$clay - sc_df["clay", "center"]) / sc_df["clay", "scale"], # Won't log scale as these are relatively normal already
  
  final_nitrogen = (log(final_covars$nitrogen) - sc_df["nitrogen", "center"]) / sc_df["nitrogen", "scale"],
  
  
  # Land use (and management for the abandoned sites)
  # For the 'initial' data, this gives the start land use, being crop or pasture
  init_lu_start = init_covars$lu_start %>% as.factor() %>% as.numeric(),
  init_lu_start_is_cropland = ifelse(init_covars$lu_start == "C", 1, 0),
  init_lu_start_is_pasture = ifelse(init_covars$lu_start == "P", 1, 0),
  
  # For the 'final' data, this gives the start land use, being crop or pasture
  final_lu_start = final_covars$lu_start %>% as.factor() %>% as.numeric(),
  
  # Type: combinations of start and end land use and management, as shown in table below
  final_type = paste(final_covars$lu_start, final_covars$lu_end, final_covars$MGMT) %>% as.factor() %>% as.numeric(),
  
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
  
  
  # Study title for random effect
  init_title = init_covars$title_id,
  final_title = final_covars$title_id,
  
  
  # IDs
  init_id = init_dat$init_id,
  final_id = final_dat$final_id,
  
  init_lu_for_tau = init_dat$lu_start %>% as.factor() %>% as.numeric(),
  final_lu_for_tau = final_dat$lu_end %>% as.factor() %>% as.numeric(),
  
  init_compare_id = final_covars$init_compare_id,
  
  init_profile_id = init_dat$init_profile_id,
  final_profile_id = final_dat$final_profile_id,
  
  final_coord = as.numeric(as.factor(final_covars$coord_id)),
  
  # Indexing
  n_init_obs = nrow(init_dat),
  n_init_id = nrow(init_covars),
  
  n_final_obs = nrow(final_dat),
  n_final_id = nrow(final_covars),
  
  n_coord = length(unique(final_covars$coord_id)),
  
  n_lu_end = length(unique(paste(final_covars$lu_end, final_covars$MGMT))),
  n_lu_start = length(unique(final_dat$lu_start)),
  
  lu_start_vals = 1:length(unique(final_dat$lu_start)),
  type_vals = 1:length(unique(paste(final_covars$lu_start, final_covars$lu_end, final_covars$MGMT))),
  
  # Counts
  n_type = length(unique(paste(final_covars$lu_start, final_covars$lu_end, final_covars$MGMT))),
  
  n_title = length(unique(final_dat$title_id)),
  n_init_profile_ids = length(unique(init_dat$init_profile_id)),
  n_final_profile_ids = length(unique(final_dat$final_profile_id)),
  
  
  # Values to be used as moderately informative priors, with the tau values
  # converted to represent precision as implemented by JAGS.
  # Priors from Sanderman et al. 2017 PNAS 
  prior_max_intercept = no_lu_coefs$coefficients["(Intercept)", "Estimate"],
  prior_tau_max_intercept = 1/no_lu_coefs$coefficients["(Intercept)", "Std. Error"]^2/1e7,
  prior_max_depth = no_lu_coefs$coefficients["depth", "Estimate"],
  prior_tau_max_depth = 1/no_lu_coefs$coefficients["depth", "Std. Error"]^2/1e7,
  
  prior_max_clay = no_lu_coefs$coefficients["clay", "Estimate"],
  prior_tau_max_clay = 1/no_lu_coefs$coefficients["clay", "Std. Error"]^2/1e7,
  
  prior_max_nitrogen = no_lu_coefs$coefficients["nitrogen", "Estimate"],
  prior_tau_max_nitrogen = 1/no_lu_coefs$coefficients["nitrogen", "Std. Error"]^2/1e7,
  
  prior_max_map_mat = no_lu_coefs$coefficients["map:mat", "Estimate"],
  prior_tau_max_map_mat = 1/no_lu_coefs$coefficients["map:mat", "Std. Error"]^2/1e7,
  
  prior_max_map = no_lu_coefs$coefficients["map", "Estimate"],
  prior_tau_max_map = 1/no_lu_coefs$coefficients["map", "Std. Error"]^2/1e7,
  prior_max_mat = no_lu_coefs$coefficients["mat", "Estimate"],
  prior_tau_max_mat = 1/no_lu_coefs$coefficients["mat", "Std. Error"]^2/1e7,
  
  
  
  
  # Priors from SoilGrids for cropland sites
  prior_cropland_intercept = cropland_coefs$coefficients["(Intercept)", "Estimate"],
  prior_tau_cropland_intercept = 1/cropland_coefs$coefficients["(Intercept)", "Std. Error"]^2/1e7,
  prior_cropland_depth = cropland_coefs$coefficients["depth", "Estimate"],
  prior_tau_cropland_depth = 1/cropland_coefs$coefficients["depth", "Std. Error"]^2/1e7,
  
  prior_cropland_clay = cropland_coefs$coefficients["clay", "Estimate"],
  prior_tau_cropland_clay = 1/cropland_coefs$coefficients["clay", "Std. Error"]^2/1e7,
  
  prior_cropland_nitrogen = cropland_coefs$coefficients["nitrogen", "Estimate"],
  prior_tau_cropland_nitrogen = 1/cropland_coefs$coefficients["nitrogen", "Std. Error"]^2/1e7,
  
  prior_cropland_map_mat = cropland_coefs$coefficients["map:mat", "Estimate"],
  prior_tau_cropland_map_mat = 1/cropland_coefs$coefficients["map:mat", "Std. Error"]^2/1e7,
  
  prior_cropland_map = cropland_coefs$coefficients["map", "Estimate"],
  prior_tau_cropland_map = 1/cropland_coefs$coefficients["map", "Std. Error"]^2/1e7,
  prior_cropland_mat = cropland_coefs$coefficients["mat", "Estimate"],
  prior_tau_cropland_mat = 1/cropland_coefs$coefficients["mat", "Std. Error"]^2/1e7,
  
  
  # Priors from SoilGrids for pasture sites
  prior_pasture_intercept = pasture_coefs$coefficients["(Intercept)", "Estimate"],
  prior_tau_pasture_intercept = 1/pasture_coefs$coefficients["(Intercept)", "Std. Error"]^2/1e7,
  prior_pasture_depth = pasture_coefs$coefficients["depth", "Estimate"],
  prior_tau_pasture_depth = 1/pasture_coefs$coefficients["depth", "Std. Error"]^2/1e7,
  
  prior_pasture_clay = pasture_coefs$coefficients["clay", "Estimate"],
  prior_tau_pasture_clay = 1/pasture_coefs$coefficients["clay", "Std. Error"]^2/1e7,
  
  prior_pasture_nitrogen = pasture_coefs$coefficients["nitrogen", "Estimate"],
  prior_tau_pasture_nitrogen = 1/pasture_coefs$coefficients["nitrogen", "Std. Error"]^2/1e7,
  
  prior_pasture_map_mat = pasture_coefs$coefficients["map:mat", "Estimate"],
  prior_tau_pasture_map_mat = 1/pasture_coefs$coefficients["map:mat", "Std. Error"]^2/1e7,
  
  prior_pasture_map = pasture_coefs$coefficients["map", "Estimate"],
  prior_tau_pasture_map = 1/pasture_coefs$coefficients["map", "Std. Error"]^2/1e7,
  prior_pasture_mat = pasture_coefs$coefficients["mat", "Estimate"],
  prior_tau_pasture_mat = 1/pasture_coefs$coefficients["mat", "Std. Error"]^2/1e7
  
)


# =============================================================================
# SECTION 4 — JAGS MODEL AND PRIMARY ANALYSIS
# =============================================================================
# The model is written to a text file and then fit using rjags.
# See inline comments inside the model string for biological interpretation
# of each block.

# Write model
sink("./scripts/bayes_script.txt")
cat("
  model {

  ### Likelihood
  
  # ---------------------------------------------------------------------------
  # Data models
  # Each observed SOC value (on the log scale) is modeled as drawn from a normal
  # distribution around the true latent SOC value at that unique location/depth/
  # study combination. Observation-level precision is allowed to differ between
  # the reference (init) and restored (final) land use classes, reflecting
  # potentially different measurement variability in these two contexts.
  # ---------------------------------------------------------------------------
  
  for(i in 1:n_init_obs){
    c_init_obs[i] ~ dnorm(c_init_true[init_id[i]],
                          c_init_tau[init_lu_for_tau[i]])
  }
  
  for(i in 1:n_final_obs){
    c_final_obs[i] ~ dnorm(c_final_true[final_id[i]],
                           c_final_tau[final_lu_for_tau[i]])
  }
  
  
  # ---------------------------------------------------------------------------
  # Process model: reference (agricultural) carbon
  #
  # The true log SOC density at reference sites is modeled as a linear function
  # of environmental covariates (depth, climate, soil properties, irrigation),
  # with separate coefficients for cropland [1] and pasture [2].
  #
  # Note on study-level measurement bias: a study-level random effect on
  # c_init_true is not included.  The model cannot distinguish genuine
  # between-study differences in SOC from systematic analytical offsets
  # (e.g. loss-on-ignition vs. dry combustion, or bulk density estimation
  # method), making such an effect non-identifiable jointly with the recovery
  # trajectory.  The dominant form of inter-study measurement inconsistency
  # tends to be multiplicative on the raw SOC scale (additive on the log
  # scale).  Because the recovery trajectory models the ratio of restored to
  # reference carbon, multiplicative biases affect both the starting point
  # and the asymptote equally and largely cancel, so the model is reasonably
  # robust to this class of inconsistency without an explicit random effect.
  # ---------------------------------------------------------------------------
  
  for(i in 1:n_init_id){
    
    c_init_true[i] <- beta_init_intercept[init_lu_start[i]] +
      beta_init_depth[init_lu_start[i]] * init_depth[i] + 
      
      beta_init_clay[init_lu_start[i]] * init_clay[i] + 
      
      beta_init_nitrogen[init_lu_start[i]] * init_nitrogen[i] + 
      
      beta_init_map_mat[init_lu_start[i]] * init_map[i] * init_mat[i] +
      
      beta_init_irrigation[init_lu_start[i]] * init_irrigation[i] + 
      
      beta_init_map[init_lu_start[i]] * init_map[i] + 
      
      beta_init_mat[init_lu_start[i]] * init_mat[i]
    
  }
  
  
  # ---------------------------------------------------------------------------
  # Process model: restored (post-abandonment) carbon
  #
  # Recovery of SOC following agricultural abandonment is modeled using a
  # Chapman-Richards curve, which describes a flexible asymptotic trajectory:
  #
  #   C(t) = C_init + (C_max - C_init) * (1 - exp(-t / exp(r))) ^ exp(m)
  #
  # where:
  #   C_init = reference agricultural carbon (carried forward from above)
  #   C_max  = the potential maximum SOC under full recovery (environmental ceiling)
  #   r      = the recovery rate (time constant, on log scale)
  #   m      = the shape parameter controlling the inflection of the curve
  #   t      = years since abandonment
  #
  # All quantities are on the log scale; exp() transforms back to the raw scale
  # for the nonlinear arithmetic before log() returns to the log scale.
  #
  # r_val includes a coordinate-level random effect (rand_coord) to account
  # for the non-independence of multiple profiles measured at the same
  # geographic location.  A coordinate-level grouping is used rather than a
  # study-level grouping because spatial proximity is the primary source of
  # within-study clustering, and many studies contribute only one or two
  # coordinates, leaving a study-level variance component poorly identified.
  # ---------------------------------------------------------------------------
  
  for(i in 1:n_final_id){
    
    c_final_true[i] <- log(exp(c_init_true[init_compare_id[i]]) + 
                             (exp(c_max[i]) - exp(c_init_true[init_compare_id[i]])) * 
                             (1 - exp(-final_age[i] / exp(r_val[i]))) ^ exp(m_val[final_type[i]]))
    
    r_val[i] <- beta_r_intercept[final_type[i]] +
      beta_r_depth[final_lu_start[i]] * final_depth[i] +
      beta_r_clay[final_lu_start[i]] * final_clay[i] +
      beta_r_nitrogen[final_lu_start[i]] * final_nitrogen[i] +
      beta_r_map_mat[final_lu_start[i]] * final_map[i] * final_mat[i] +
      beta_r_map[final_lu_start[i]] * final_map[i] +
      beta_r_mat[final_lu_start[i]] * final_mat[i] +
      rand_coord[final_coord[i]]
    
    c_max[i] <- beta_max_intercept + 
      beta_max_depth * final_depth[i] +
      beta_max_clay * final_clay[i] +
      beta_max_nitrogen * final_nitrogen[i] +
      beta_max_map_mat * final_map[i] * final_mat[i] +
      beta_max_map * final_map[i] +
      beta_max_mat * final_mat[i]
    
  }
  
  
  ### Priors
  
  # ---------------------------------------------------------------------------
  # Reference carbon coefficients
  # Moderately informative priors derived from SoilGrids-based regressions for
  # cropland [1] and pasture [2] sites (see 01_lc_priors.R). The precision
  # values (tau) are divided by 1e7 to be weakly informative relative to the
  # data, using the external estimates only to regularize rather than dominate.
  # ---------------------------------------------------------------------------
  
  # Cropland [1]
  beta_init_intercept[1] ~ dnorm(prior_cropland_intercept, prior_tau_cropland_intercept)
  beta_init_depth[1]     ~ dnorm(prior_cropland_depth,     prior_tau_cropland_depth)
  beta_init_clay[1]      ~ dnorm(prior_cropland_clay,      prior_tau_cropland_clay)
  beta_init_nitrogen[1]  ~ dnorm(prior_cropland_nitrogen,  prior_tau_cropland_nitrogen)
  beta_init_map_mat[1]   ~ dnorm(prior_cropland_map_mat,   prior_tau_cropland_map_mat)
  beta_init_map[1]       ~ dnorm(prior_cropland_map,       prior_tau_cropland_map)
  beta_init_mat[1]       ~ dnorm(prior_cropland_mat,       prior_tau_cropland_mat)
  
  # Pasture [2]
  beta_init_intercept[2] ~ dnorm(prior_pasture_intercept, prior_tau_pasture_intercept)
  beta_init_depth[2]     ~ dnorm(prior_pasture_depth,     prior_tau_pasture_depth)
  beta_init_clay[2]      ~ dnorm(prior_pasture_clay,      prior_tau_pasture_clay)
  beta_init_nitrogen[2]  ~ dnorm(prior_pasture_nitrogen,  prior_tau_pasture_nitrogen)
  beta_init_map_mat[2]   ~ dnorm(prior_pasture_map_mat,   prior_tau_pasture_map_mat)
  beta_init_map[2]       ~ dnorm(prior_pasture_map,       prior_tau_pasture_map)
  beta_init_mat[2]       ~ dnorm(prior_pasture_mat,       prior_tau_pasture_mat)
  
  # Irrigation effect: hierarchical across land use classes
  for(i in lu_start_vals){
    beta_init_irrigation[i] ~ dnorm(beta_mean_init_irrigation, tau_mean_init_irrigation)
  }
  
  
  # ---------------------------------------------------------------------------
  # Recovery rate (r) covariate coefficients
  # Weakly informative independent priors (precision = 0.25, i.e. SD ~ 2) per
  # land use start class. A hypervariance structure across land use classes was
  # explored but found to be non-identifiable; independent priors are used instead.
  # ---------------------------------------------------------------------------
  
  for(i in 1:2){
    beta_r_depth[i]    ~ dnorm(0, 0.25)
    beta_r_clay[i]     ~ dnorm(0, 0.25)
    beta_r_nitrogen[i] ~ dnorm(0, 0.25)
    beta_r_map_mat[i]  ~ dnorm(0, 0.25)
    beta_r_map[i]      ~ dnorm(0, 0.25)
    beta_r_mat[i]      ~ dnorm(0, 0.25)
  }
  
  
  # ---------------------------------------------------------------------------
  # Recovery rate intercepts and shape parameters by transition type
  # There are 8 transition types defined by the combination of starting land use
  # (cropland/pasture), ending land cover (forest/open), and management strategy
  # (active/natural). Intercepts are drawn from a common hyperprior to enable
  # partial pooling across transition types.
  # Shape parameter m is given a tight prior around 1.1, reflecting approximately
  # sigmoidal recovery consistent with prior literature.
  # ---------------------------------------------------------------------------
  
  for(i in type_vals){
    beta_r_intercept[i] ~ dnorm(beta_mean_r_intercept, tau_mean_r_intercept)
    m_val[i] ~ dnorm(1.1, 5)
  }
  
  beta_mean_r_intercept ~ dnorm(3, 1)
  tau_mean_r_intercept <- pow(sd_mean_r_intercept, -2)
  sd_mean_r_intercept ~ dunif(0, 5)
  
  
  # ---------------------------------------------------------------------------
  # Maximum SOC (C_max) coefficients
  # Informative priors from Sanderman et al. 2017 PNAS, representing natural
  # (pre-disturbance) SOC as a function of environment. Sign constraints on MAP
  # (positive) and MAT (negative) are applied via truncation to enforce
  # ecologically expected directions for the maximum attainable SOC.
  # ---------------------------------------------------------------------------
  
  beta_max_intercept ~ dnorm(prior_max_intercept, prior_tau_max_intercept)
  beta_max_depth     ~ dnorm(prior_max_depth,     prior_tau_max_depth)
  beta_max_clay      ~ dnorm(prior_max_clay,      prior_tau_max_clay)
  beta_max_nitrogen  ~ dnorm(prior_max_nitrogen,  prior_tau_max_nitrogen)
  beta_max_map_mat   ~ dnorm(prior_max_map_mat,   prior_tau_max_map_mat)
  beta_max_map       ~ dnorm(prior_max_map,       prior_tau_max_map)T(0, )
  beta_max_mat       ~ dnorm(prior_max_mat,       prior_tau_max_mat)T(, 0)
  
  
  # ---------------------------------------------------------------------------
  # Irrigation hyperprior
  # ---------------------------------------------------------------------------
  
  beta_mean_init_irrigation ~ dnorm(0, 0.5)
  tau_mean_init_irrigation <- pow(sd_mean_init_irrigation, -2)
  sd_mean_init_irrigation ~ dunif(0, 5)
  
  
  # ---------------------------------------------------------------------------
  # Coordinate-level random effect on recovery rate (rand_coord)
  #
  # Captures non-independence among profiles measured at the same geographic
  # location.  Applied to r_val (the recovery rate) only; applying it to
  # c_max or c_init_true creates identifiability conflicts with those
  # parameters.  A sum-to-zero constraint is imposed for identifiability.
  # ---------------------------------------------------------------------------
  
  for(i in 1:(n_coord-1)){
    rand_coord[i] ~ dnorm(0, tau_coord)
  }
  rand_coord[n_coord] <- -sum(rand_coord[1:(n_coord-1)])
  tau_coord <- pow(coord_sd, -2)
  coord_sd ~ dunif(0, 5)
  
  
  # ---------------------------------------------------------------------------
  # Observation-level residual standard deviations
  # Separate SDs for reference (init) and restored (final) land use classes,
  # reflecting potentially different measurement variability in these contexts.
  # ---------------------------------------------------------------------------
  
  for(i in 1:2){
    c_init_tau[i] <- pow(c_init_sd[i], -2)
    c_init_sd[i]  ~ dunif(0, 5)
    c_final_tau[i] <- pow(c_final_sd[i], -2)
    c_final_sd[i]  ~ dunif(0, 5)
  }
  
  
  # ---------------------------------------------------------------------------
  # Derived quantities: random-effect-free predictions
  #
  # These quantities reproduce the recovery trajectory with the coordinate-level
  # random effect on r removed, yielding predictions that reflect only the
  # fixed-effect environmental and land-use structure. Used for inference on
  # population-average trajectories and for generating predictions at new
  # locations.
  # ---------------------------------------------------------------------------
  
  for(i in 1:n_final_id){
    
    c_final_norand[i] <- log(exp(c_init_true[init_compare_id[i]]) + 
                               (exp(c_max_norand[i]) - exp(c_init_true[init_compare_id[i]])) * 
                               (1 - exp(-final_age[i] / exp(r_val_norand[i]))) ^ exp(m_val[final_type[i]]))
    
    r_val_norand[i] <- beta_r_intercept[final_type[i]] +
      beta_r_depth[final_lu_start[i]] * final_depth[i] +
      beta_r_clay[final_lu_start[i]] * final_clay[i] +
      beta_r_nitrogen[final_lu_start[i]] * final_nitrogen[i] +
      beta_r_map_mat[final_lu_start[i]] * final_map[i] * final_mat[i] +
      beta_r_map[final_lu_start[i]] * final_map[i] +
      beta_r_mat[final_lu_start[i]] * final_mat[i]
    
    c_max_norand[i] <- beta_max_intercept +
      beta_max_depth * final_depth[i] +
      beta_max_clay * final_clay[i] +
      beta_max_nitrogen * final_nitrogen[i] +
      beta_max_map_mat * final_map[i] * final_mat[i] +
      beta_max_map * final_map[i] +
      beta_max_mat * final_mat[i]

  }

}

    ",fill=TRUE)
sink()

# =============================================================================
# MCMC SETTINGS
# =============================================================================
# All JAGS runs (primary analysis + all sensitivity runs) draw their adapt,
# burn-in, sampling, and thinning values from these four objects.  Change
# them here once to affect every run in the script.
#
# Recommended values:
#   Trial run  : n_adapt = 1000, n_burnin = 1000, n_iter = 10000,  n_thin = 10
#   Full run   : n_adapt = 2000, n_burnin = 2000, n_iter = 50000, n_thin = 50

n_chains <- 3
n_adapt  <- 2000    # adaptation steps
n_burnin <- 2000    # burn-in iterations (passed to update())
n_iter   <- 50000   # sampling iterations (passed to coda.samples())
n_thin   <- 50      # thinning interval


# =============================================================================
# INITIAL VALUES
# =============================================================================
# Providing starting values anchors all chains in the same plausible region.
#   beta_max_intercept = 1.5  (plausible log SOC density, above typical
#                              agricultural values to enforce c_max > c_init)
#   beta_r_intercept   = 3    (moderate recovery rate, log scale)
#   m_val              = 1.1  (approximately sigmoidal, consistent with prior)
#
# make_inits() is a function rather than a fixed list so that jags.model()
# calls it once per chain, producing slightly different starting values for
# each chain while still anchoring all chains near plausible values.
# Coverage is comprehensive: all major stochastic nodes are initialised,
# including covariate slopes, hyperparameters, and random effects.
# Random effects start near zero (the "no effect" prior expectation).
# SD parameters use abs() to prevent negative starts.
#
# n_type  — number of transition types in this run (may be <8 for BD subsets)
# n_coord — number of unique coordinates (rand_coord initialised near zero)

make_inits <- function(n_type, n_coord = 1) {
  list(
    # Maximum SOC coefficients
    beta_max_intercept = 1.5  + rnorm(1, 0, 0.10),
    beta_max_depth     = -0.5 + rnorm(1, 0, 0.05),
    beta_max_map       = 0.3  + rnorm(1, 0, 0.05),
    beta_max_mat       = -0.1 + rnorm(1, 0, 0.05),
    beta_max_clay      = 0.1  + rnorm(1, 0, 0.05),
    beta_max_nitrogen  = 0.1  + rnorm(1, 0, 0.05),
    beta_max_map_mat   = 0.1  + rnorm(1, 0, 0.05),
    
    # Recovery rate: intercepts by transition type and covariate slopes
    beta_r_intercept   = rep(3, n_type) + rnorm(n_type, 0, 0.2),
    beta_r_depth       = c(0, 0)        + rnorm(2, 0, 0.05),
    beta_r_clay        = c(0, 0)        + rnorm(2, 0, 0.05),
    beta_r_nitrogen    = c(0, 0)        + rnorm(2, 0, 0.05),
    beta_r_map         = c(0, 0)        + rnorm(2, 0, 0.05),
    beta_r_mat         = c(0, 0)        + rnorm(2, 0, 0.05),
    beta_r_map_mat     = c(0, 0)        + rnorm(2, 0, 0.05),
    
    # Shape parameter by transition type
    m_val              = rep(1.1, n_type) + rnorm(n_type, 0, 0.1),
    
    # Hyperparameters
    beta_mean_r_intercept    = 3   + rnorm(1, 0, 0.1),
    sd_mean_r_intercept      = abs(0.5 + rnorm(1, 0, 0.05)),
    beta_mean_init_irrigation = 0  + rnorm(1, 0, 0.05),
    sd_mean_init_irrigation  = abs(0.3 + rnorm(1, 0, 0.03)),
    
    # Reference carbon covariate slopes (one per land-use class)
    beta_init_intercept = c(0.1, 0.5)  + rnorm(2, 0, 0.05),
    beta_init_depth     = c(-0.3, -0.2) + rnorm(2, 0, 0.03),
    beta_init_map       = c(0, 0)      + rnorm(2, 0, 0.05),
    beta_init_mat       = c(0, 0)      + rnorm(2, 0, 0.05),
    beta_init_clay      = c(0, 0)      + rnorm(2, 0, 0.05),
    beta_init_nitrogen  = c(0, 0)      + rnorm(2, 0, 0.05),
    beta_init_map_mat   = c(0, 0)      + rnorm(2, 0, 0.05),
    beta_init_irrigation = c(0, 0)     + rnorm(2, 0, 0.05),
    
    # Coordinate-level random effects: start near zero
    rand_coord         = c(rep(0, n_coord - 1), NA),  # last is sum-to-zero deterministic
    
    # Observation-level SDs
    coord_sd           = abs(0.5  + rnorm(1, 0, 0.05)),
    c_init_sd          = abs(c(0.3, 0.3) + rnorm(2, 0, 0.02)),
    c_final_sd         = abs(c(0.3, 0.3) + rnorm(2, 0, 0.02))
  )
}


###
### JAGS analysis
###
ipak("rjags")

variables_to_sample <- c("beta_init_intercept"
                         ,"beta_init_depth"
                         ,"beta_init_clay"
                         ,"beta_init_nitrogen"
                         ,"beta_init_map_mat"
                         ,"beta_init_irrigation"
                         ,"beta_init_map"
                         ,"beta_init_mat"
                         
                         ,"beta_r_intercept"
                         ,"beta_r_depth"
                         ,"beta_r_clay"
                         ,"beta_r_nitrogen"
                         ,"beta_r_map_mat"
                         ,"beta_r_map"
                         ,"beta_r_mat"
                         
                         ,"m_val"
                         
                         ,"beta_max_intercept"
                         ,"beta_max_depth"
                         ,"beta_max_clay"
                         ,"beta_max_nitrogen"
                         ,"beta_max_map_mat"
                         ,"beta_max_map"
                         ,"beta_max_mat"
                         
                         ,"c_init_tau"
                         ,"c_init_sd"
                         ,"c_final_tau"
                         ,"c_final_sd"
                         
                         ,"c_init_true"
                         ,"c_final_true"
                         
                         ,"c_final_norand"
)

# Also sample the sd values
sd_variables <- gsub("^beta_mean_", "sd_mean_", variables_to_sample)

# Combine both sets of variables, ensuring no duplicates
variables_to_sample <- unique(c(variables_to_sample, sd_variables))


primary_outfile <- "./data_output/bayes_output.RData"

if (!refit_models && file.exists(primary_outfile)) {
  message("Loading cached primary analysis from ", primary_outfile)
  load(primary_outfile)
} else {
  if (refit_models) {
    message("refit_models = TRUE — fitting primary JAGS model from scratch.")
  } else {
    message("No cached primary output found — fitting primary JAGS model.")
  }
  set.seed(111)
  mod <- jags.model("./scripts/bayes_script.txt", data = jags_data,
                    inits    = make_inits(jags_data$n_type, jags_data$n_coord),
                    n.chains = n_chains, n.adapt = n_adapt)
  update(mod, n.iter = n_burnin)
  
  mod_samp <- coda.samples(mod,
                           variable.names = variables_to_sample,
                           n.iter = n_iter, thin = n_thin)
  
  # Evaluate convergence: trace plots for scalar parameters
  # (vectors such as c_init_true and c_final_norand are excluded for clarity)
  MCMCvis::MCMCtrace(mod_samp,
                     params = variables_to_sample[!(variables_to_sample %in% c(
                       "c_init_true", "c_final_true",
                       "c_init_norand", "c_final_norand"
                     ))] %>% sort())
  
  # R-hat convergence summary: values close to 1 indicate chain mixing
  mod_summary <- MCMCvis::MCMCsummary(mod_samp)
  hist(mod_summary$Rhat,
       main = "R-hat distribution — primary analysis",
       xlab = "R-hat")
  n_bad_rhat <- sum(mod_summary$Rhat > 1.1, na.rm = TRUE)
  if (n_bad_rhat > 0) {
    warning(n_bad_rhat, " parameters have R-hat > 1.1")
  }
  
  # Save primary analysis output
  save(file = primary_outfile,
       mod_samp, variables_to_sample, sc_df,
       init_dat, init_covars,
       final_dat, final_covars)
}

# Ensure mod_summary is available for the diagnostic block below, regardless
# of whether the model was just fitted or loaded from cache.
if (!exists("mod_summary")) {
  mod_summary <- MCMCvis::MCMCsummary(mod_samp)
}
if (!exists("n_bad_rhat")) {
  n_bad_rhat <- sum(mod_summary$Rhat > 1.1, na.rm = TRUE)
}


# -----------------------------------------------------------------------------
# MCMC DIAGNOSTIC SUMMARY — for reporting in manuscript
# -----------------------------------------------------------------------------
# Compile key convergence and mixing statistics from mod_summary, which was
# produced by MCMCvis::MCMCsummary() above and contains one row per parameter.

n_params    <- nrow(mod_summary)
max_rhat    <- round(max(mod_summary$Rhat,  na.rm = TRUE), 3)
n_bad_rhat_strict <- sum(mod_summary$Rhat > 1.05, na.rm = TRUE)  # stricter threshold
min_neff    <- round(min(mod_summary$n.eff, na.rm = TRUE))
median_neff <- round(median(mod_summary$n.eff, na.rm = TRUE))

message(
  "\n===== MCMC DIAGNOSTICS (primary analysis) =====",
  "\n  Total parameters monitored : ", n_params,
  "\n  Max R-hat                  : ", max_rhat,
  "\n  Parameters with R-hat > 1.1 : ", n_bad_rhat,   # from existing check above
  "\n  Parameters with R-hat > 1.05: ", n_bad_rhat_strict,
  "\n  Min effective sample size  : ", min_neff,
  "\n  Median effective sample size: ", median_neff,
  "\n===============================================\n"
)

# Parameters with the worst convergence (useful to inspect individually)
worst_rhat <- mod_summary[order(mod_summary$Rhat, decreasing = TRUE), ][1:10, ]
print(worst_rhat[, c("Rhat", "n.eff")])

# Parameters with lowest effective sample size
worst_neff <- mod_summary[order(mod_summary$n.eff), ][1:10, ]
print(worst_neff[, c("Rhat", "n.eff")])



# =============================================================================
# Two-panel figure combining the empirical variogram (left) and a global map
# of mean residuals per coordinate (right).  Both panels use fixed-effect-only
# residuals (c_final_norand and c_init_true) pooled across reference and
# restored sites.
#
# The variogram shows the scale at which residual spatial autocorrelation
# operates.  If semivariance reaches the marginal variance (sill) at short-to-
# moderate distances with no further rise at continental scales, unmodelled
# spatial dependence does not bias broad-scale covariate estimates.
#
# The map shows whether any region is systematically over- or under-predicted.
# A patchy, directionless mix of red and blue supports the claim that residual
# autocorrelation does not propagate broad-scale bias into the covariate
# estimates of interest.

ipak("geosphere")  # distGeo() — not loaded elsewhere

# ── Compute pooled residuals ───────────────────────────────────────────────────
norand_summary <- MCMCvis::MCMCsummary(mod_samp, params = "c_final_norand",
                                       probs = 0.5)
pred_final    <- norand_summary[, "50%"][jags_data$final_id]
resid_final_v <- jags_data$c_final_obs - pred_final

init_summary  <- MCMCvis::MCMCsummary(mod_samp, params = "c_init_true",
                                      probs = 0.5)
pred_init     <- init_summary[, "50%"][jags_data$init_id]
resid_init_v  <- jags_data$c_init_obs - pred_init

resid_pooled <- bind_rows(
  data.frame(resid = resid_final_v, LAT = final_dat$LAT, LONG = final_dat$LONG),
  data.frame(resid = resid_init_v,  LAT = init_dat$LAT,  LONG = init_dat$LONG)
) %>%
  filter(complete.cases(.))

# ── Variogram ─────────────────────────────────────────────────────────────────
compute_variogram <- function(df, max_dist_km = 12800,
                              breaks = c(0, 100, 200, 400, 800,
                                         1600, 3200, 6400, 12800)) {
  coords   <- as.matrix(df[, c("LONG", "LAT")])
  n        <- nrow(df)
  dist_mat <- matrix(NA_real_, n, n)
  for (i in seq_len(n - 1)) {
    dist_mat[i, (i + 1):n] <-
      geosphere::distGeo(coords[i, , drop = FALSE],
                         coords[(i + 1):n, , drop = FALSE]) / 1000
  }
  sv_mat   <- outer(df$resid, df$resid, FUN = function(a, b) 0.5 * (a - b)^2)
  idx      <- upper.tri(dist_mat)
  dist_vec <- dist_mat[idx]
  sv_vec   <- sv_mat[idx]
  keep     <- !is.na(dist_vec) & dist_vec <= max_dist_km
  dist_vec <- dist_vec[keep]
  sv_vec   <- sv_vec[keep]
  bin_labels <- cut(dist_vec, breaks = breaks, include.lowest = TRUE)
  data.frame(
    dist_mid = tapply(dist_vec, bin_labels, mean,   na.rm = TRUE),
    gamma    = tapply(sv_vec,   bin_labels, mean,   na.rm = TRUE),
    n_pairs  = tapply(sv_vec,   bin_labels, length)
  ) %>%
    tibble::rownames_to_column("bin_label") %>%
    filter(!is.na(dist_mid))
}

variog_df <- compute_variogram(resid_pooled)
sill_val  <- var(resid_pooled$resid, na.rm = TRUE)
n_pairs_total <- formatC(sum(variog_df$n_pairs, na.rm = TRUE),
                         format = "d", big.mark = ",")

p_variog <- ggplot(variog_df, aes(x = dist_mid, y = gamma)) +
  geom_line(colour = "grey30") +
  geom_point(colour = "grey30", size = 2) +
  geom_hline(yintercept = sill_val, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 12800, y = sill_val,
           label = "Sill", hjust = 1.1, vjust = -0.4,
           colour = "grey50", size = 3) +
  annotate("text", x = 12800, y = max(variog_df$gamma, na.rm = TRUE) * 1.05,
           label = paste0("n pairs: ", n_pairs_total),
           hjust = 1, size = 2.8, colour = "grey40") +
  scale_x_continuous(limits = c(0, 12800),
                     labels = scales::comma) +
  scale_y_continuous(limits = c(0, max(sill_val, variog_df$gamma,
                                       na.rm = TRUE) * 1.12)) +
  labs(x = "Distance (km)", y = "Semivariance",
       title = "Residual variogram") +
  theme_classic() +
  theme(
    plot.title   = element_text(size = 10, face = "bold"),
    axis.title   = element_text(size = 9),
    axis.text    = element_text(size = 8)
  )

# ── Map ────────────────────────────────────────────────────────────────────────
resid_coord_df <- resid_pooled %>%
  group_by(LAT, LONG) %>%
  summarise(resid_mean = mean(resid, na.rm = TRUE), .groups = "drop")

resid_sd  <- sd(resid_coord_df$resid_mean, na.rm = TRUE)
scale_lim <- 2 * resid_sd

world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(subregion != "Antarctica")

p_map <- ggplot() +
  geom_sf(data = world, fill = "grey90", colour = NA) +
  geom_point(
    data = resid_coord_df,
    aes(x = LONG, y = LAT, colour = resid_mean),
    size = 1.5, alpha = 0.85, stroke = 0
  ) +
  scale_colour_distiller(
    palette   = "RdBu",
    direction = 1,
    limits    = c(-scale_lim, scale_lim),
    oob       = scales::squish,
    name      = "Mean residual\n(log SOC density)",
    guide     = guide_colourbar(barwidth = 0.7, barheight = 5)
  ) +
  coord_sf(expand = FALSE) +
  labs(x = NULL, y = NULL, title = "Residual spatial pattern") +
  theme_minimal() +
  theme(
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    plot.title    = element_text(size = 10, face = "bold"),
    legend.title  = element_text(size = 8),
    legend.text   = element_text(size = 7)
  )

# ── Combine with patchwork ─────────────────────────────────────────────────────
combined <- p_variog + p_map +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(tag_levels = "A")

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png") {
    ggsave("./figures/residual_diagnostics.png", combined,
           width = 11, height = 4, dpi = 440)
  } else {
    ggsave("./figures/residual_diagnostics.pdf", combined,
           width = 11, height = 4)
  }
  message("Saved ./figures/residual_diagnostics.", dev_type)
}

rm(norand_summary, pred_final, resid_final_v,
   init_summary, pred_init, resid_init_v,
   resid_pooled, variog_df, sill_val, n_pairs_total,
   resid_coord_df, resid_sd, scale_lim, world,
   p_variog, p_map, combined)



# =============================================================================
# SECTION 5 — SENSITIVITY ANALYSES
# =============================================================================
# Each sensitivity analysis re-fits the model on a filtered subset of c_dat
# and is designed to vary only the data dimension under investigation, not
# the prior specification.  Accordingly, MGMT uncertainty and SOM conversion
# runs use the weak-prior model (bayes_script_weakprior.txt), and the
# weak-prior run in Section 6 tests prior sensitivity directly on the full
# dataset.
#
# The two BD sensitivity runs (a and d) use the informative model
# (bayes_script.txt) with the beta_max_* priors updated to be centred on
# the primary posterior means.  beta_max_* describes the maximum attainable
# SOC as a function of climate and soil properties — a background property
# of the landscape that is not under investigation in either BD sensitivity
# test and that is well identified by the full dataset.  Because the BD
# subsets retain substantially fewer unique geographic locations, anchoring
# these environmental scaling parameters to the full-data posterior is the
# appropriate prior choice: it reflects genuine prior knowledge about the
# SOC-environment relationship rather than an uninformative default, and
# does not constrain the recovery trajectory parameters (beta_r_*, m_val)
# that are the focus of the sensitivity comparison.
#
# Prior specification by run:
#   BD gap-filling (a):   bayes_script.txt, beta_max_* centred on primary posterior
#   Original-BD (d):      bayes_script.txt, beta_max_* centred on primary posterior
#   MGMT uncertainty (b): bayes_script_weakprior.txt
#   SOM conversion (c):   bayes_script_weakprior.txt
#   Weak-prior (Sec. 6):  bayes_script_weakprior.txt, full dataset
#
# jags_data is rebuilt from the filtered dataset using the same reshaping
# pipeline as Section 3; sc_df remains unchanged throughout.  A helper
# function encapsulates the reshaping + jags_data construction so that each
# sensitivity run is a one-liner filter + function call.

build_jags_data <- function(dat) {
  # ── Reshape ────────────────────────────────────────────────────────────────
  # Identical pipeline to Section 3, applied to the supplied (possibly filtered)
  # dataset.  See Section 3 for commentary on each step.
  
  dat$init_id <- dat %>%
    dplyr::select(coord_id, depth_avg, Authors, Year) %>%
    apply(1, function(x) paste(x, collapse = "")) %>%
    as.factor() %>% as.numeric()
  
  dat$final_id <- dat %>%
    dplyr::select(coord_id, depth_avg, Authors, Year, AGE, lu_end, MGMT) %>%
    apply(1, function(x) paste(x, collapse = "")) %>%
    as.factor() %>% as.numeric()
  
  dat <- dat %>% mutate(title_id = as.numeric(as.factor(Title)))
  dat <- dat[order(dat$final_id), ]
  
  i_dat <- dat %>%
    dplyr::select(-starts_with("Aban"), -contains("accrual"),
                  -contains("Perc_change"), -contains("Cum_Aban"),
                  -starts_with("SE"), -starts_with("SD"),
                  -obs, -starts_with("AGE"), -final_id)
  i_dat <- i_dat[!duplicated(i_dat %>%
                               dplyr::select(-Profile_ID, -lu_end, -lu_start)), ]
  i_dat <- i_dat %>% mutate(init_profile_id = as.numeric(as.factor(Profile_ID)))
  
  f_dat <- dat %>%
    dplyr::select(-starts_with("Ctrl"), -contains("accrual"),
                  -contains("Perc_change"), -contains("Cum_Ctrl"),
                  -starts_with("SE"), -starts_with("SD"), -obs)
  f_dat <- f_dat[!duplicated(f_dat %>% dplyr::select(-Profile_ID)), ]
  f_dat <- f_dat %>% mutate(final_profile_id = as.numeric(as.factor(Profile_ID)))
  
  i_covars <- i_dat %>%
    dplyr::select(init_id, depth_avg, MAP, MAT, clay, nitrogen, irrigation,
                  title_id, lu_start, lu_end) %>%
    group_by(init_id) %>%
    summarize(across(where(is.numeric), mean),
              across(where(is.character), first))
  
  f_covars <- f_dat %>%
    dplyr::rename(init_compare_id = init_id) %>%
    dplyr::select(final_id, init_compare_id, depth_avg, MAP, MAT,
                  clay, nitrogen, irrigation, AGE, title_id, coord_id,
                  lu_start, lu_end, MGMT) %>%
    group_by(final_id) %>%
    summarize(across(where(is.numeric), mean),
              across(where(is.character), first))
  
  # ── Assemble jags_data ──────────────────────────────────────────────────────
  # Uses the same sc_df from the primary (full) dataset throughout.
  list(
    c_final_obs = log(f_dat$Aban_c_dens_av),
    c_init_obs  = log(i_dat$Ctrl_c_dens_av),
    
    init_depth = (log(i_covars$depth_avg) - sc_df["depth","center"]) / sc_df["depth","scale"],
    init_age   = rep(0, nrow(i_covars)),
    init_map   = (log(i_covars$MAP)       - sc_df["map",  "center"]) / sc_df["map",  "scale"],
    init_mat   = (i_covars$MAT            - sc_df["mat",  "center"]) / sc_df["mat",  "scale"],
    init_clay  = (i_covars$clay           - sc_df["clay", "center"]) / sc_df["clay", "scale"],
    init_nitrogen  = (log(i_covars$nitrogen)   - sc_df["nitrogen",  "center"]) / sc_df["nitrogen",  "scale"],
    init_irrigation = (log(i_covars$irrigation) - sc_df["irrigation","center"]) / sc_df["irrigation","scale"],
    
    final_depth = (log(f_covars$depth_avg) - sc_df["depth","center"]) / sc_df["depth","scale"],
    final_age   = f_covars$AGE,
    final_map   = (log(f_covars$MAP)       - sc_df["map",  "center"]) / sc_df["map",  "scale"],
    final_mat   = (f_covars$MAT            - sc_df["mat",  "center"]) / sc_df["mat",  "scale"],
    final_clay  = (f_covars$clay           - sc_df["clay", "center"]) / sc_df["clay", "scale"],
    final_nitrogen = (log(f_covars$nitrogen) - sc_df["nitrogen","center"]) / sc_df["nitrogen","scale"],
    
    init_lu_start             = i_covars$lu_start %>% as.factor() %>% as.numeric(),
    init_lu_start_is_cropland = ifelse(i_covars$lu_start == "C", 1, 0),
    init_lu_start_is_pasture  = ifelse(i_covars$lu_start == "P", 1, 0),
    final_lu_start = f_covars$lu_start %>% as.factor() %>% as.numeric(),
    final_type     = paste(f_covars$lu_start, f_covars$lu_end, f_covars$MGMT) %>%
      as.factor() %>% as.numeric(),
    
    init_title  = i_covars$title_id,
    final_title = f_covars$title_id,
    
    init_id  = i_dat$init_id,
    final_id = f_dat$final_id,
    
    init_lu_for_tau  = i_dat$lu_start %>% as.factor() %>% as.numeric(),
    final_lu_for_tau = f_dat$lu_end   %>% as.factor() %>% as.numeric(),
    
    init_compare_id     = f_covars$init_compare_id,
    init_profile_id     = i_dat$init_profile_id,
    final_profile_id    = f_dat$final_profile_id,
    final_coord         = as.numeric(as.factor(f_covars$coord_id)),
    
    n_init_obs  = nrow(i_dat),
    n_init_id   = nrow(i_covars),
    n_final_obs = nrow(f_dat),
    n_final_id  = nrow(f_covars),
    n_coord     = length(unique(f_covars$coord_id)),
    n_lu_end    = length(unique(paste(f_covars$lu_end, f_covars$MGMT))),
    n_lu_start  = length(unique(f_dat$lu_start)),
    lu_start_vals = seq_len(length(unique(f_dat$lu_start))),
    type_vals     = seq_len(length(unique(paste(f_covars$lu_start,
                                                f_covars$lu_end, f_covars$MGMT)))),
    n_type    = length(unique(paste(f_covars$lu_start, f_covars$lu_end, f_covars$MGMT))),
    n_title   = length(unique(f_dat$title_id)),
    n_init_profile_ids  = length(unique(i_dat$init_profile_id)),
    n_final_profile_ids = length(unique(f_dat$final_profile_id)),
    
    # Mean log SOC across reference observations in this subset — used as the
    # prior mean for beta_init_intercept and beta_max_intercept in the
    # weak-prior model (bayes_script_weakprior.txt).  Computed from the
    # filtered data so it reflects the subset being analysed.
    obs_log_soc_mean = mean(log(i_dat$Ctrl_c_dens_av), na.rm = TRUE),
    
    # Prior parameters — same values used for primary run
    prior_max_intercept       = jags_data$prior_max_intercept,
    prior_tau_max_intercept   = jags_data$prior_tau_max_intercept,
    prior_max_depth           = jags_data$prior_max_depth,
    prior_tau_max_depth       = jags_data$prior_tau_max_depth,
    prior_max_clay            = jags_data$prior_max_clay,
    prior_tau_max_clay        = jags_data$prior_tau_max_clay,
    prior_max_nitrogen        = jags_data$prior_max_nitrogen,
    prior_tau_max_nitrogen    = jags_data$prior_tau_max_nitrogen,
    prior_max_map_mat         = jags_data$prior_max_map_mat,
    prior_tau_max_map_mat     = jags_data$prior_tau_max_map_mat,
    prior_max_map             = jags_data$prior_max_map,
    prior_tau_max_map         = jags_data$prior_tau_max_map,
    prior_max_mat             = jags_data$prior_max_mat,
    prior_tau_max_mat         = jags_data$prior_tau_max_mat,
    prior_cropland_intercept  = jags_data$prior_cropland_intercept,
    prior_tau_cropland_intercept = jags_data$prior_tau_cropland_intercept,
    prior_cropland_depth      = jags_data$prior_cropland_depth,
    prior_tau_cropland_depth  = jags_data$prior_tau_cropland_depth,
    prior_cropland_clay       = jags_data$prior_cropland_clay,
    prior_tau_cropland_clay   = jags_data$prior_tau_cropland_clay,
    prior_cropland_nitrogen   = jags_data$prior_cropland_nitrogen,
    prior_tau_cropland_nitrogen = jags_data$prior_tau_cropland_nitrogen,
    prior_cropland_map_mat    = jags_data$prior_cropland_map_mat,
    prior_tau_cropland_map_mat = jags_data$prior_tau_cropland_map_mat,
    prior_cropland_map        = jags_data$prior_cropland_map,
    prior_tau_cropland_map    = jags_data$prior_tau_cropland_map,
    prior_cropland_mat        = jags_data$prior_cropland_mat,
    prior_tau_cropland_mat    = jags_data$prior_tau_cropland_mat,
    prior_pasture_intercept   = jags_data$prior_pasture_intercept,
    prior_tau_pasture_intercept = jags_data$prior_tau_pasture_intercept,
    prior_pasture_depth       = jags_data$prior_pasture_depth,
    prior_tau_pasture_depth   = jags_data$prior_tau_pasture_depth,
    prior_pasture_clay        = jags_data$prior_pasture_clay,
    prior_tau_pasture_clay    = jags_data$prior_tau_pasture_clay,
    prior_pasture_nitrogen    = jags_data$prior_pasture_nitrogen,
    prior_tau_pasture_nitrogen = jags_data$prior_tau_pasture_nitrogen,
    prior_pasture_map_mat     = jags_data$prior_pasture_map_mat,
    prior_tau_pasture_map_mat = jags_data$prior_tau_pasture_map_mat,
    prior_pasture_map         = jags_data$prior_pasture_map,
    prior_tau_pasture_map     = jags_data$prior_tau_pasture_map,
    prior_pasture_mat         = jags_data$prior_pasture_mat,
    prior_tau_pasture_mat     = jags_data$prior_tau_pasture_mat
  )
}

run_sensitivity <- function(filtered_dat, label, outfile, obj_name,
                            chains             = n_chains,
                            adapt              = n_adapt,
                            burnin             = n_burnin,
                            iter               = n_iter,
                            thin               = n_thin,
                            model_file         = "./scripts/bayes_script_weakprior.txt",
                            jags_data_override = list()) {
  # Fits the same model on a filtered subset and saves the output.
  # Arguments:
  #   filtered_dat       - subset of c_dat produced by filter()
  #   label              - short string used in console messages
  #   outfile            - path for the saved .RData file
  #   obj_name           - character string naming the samples object inside
  #                        the saved .RData file (e.g. "mod_samp_bd"), so that
  #                        loading multiple sensitivity outputs into the same
  #                        session does not cause name conflicts
  #   chains/adapt/burnin/iter/thin - MCMC settings; default to the primary
  #                        analysis values but can be overridden if needed
  #   model_file         - path to the JAGS model text file; defaults to
  #                        bayes_script_weakprior.txt for all sensitivity runs
  #   jags_data_override - named list of jags_data elements to replace after
  #                        build_jags_data(); used by BD runs to substitute
  #                        primary-posterior-derived beta_max_* priors
  
  # If a cached file exists and the user has not requested a refit, load and return
  if (!refit_models && file.exists(outfile)) {
    message("\nLoading cached sensitivity output (", label, ") from ", outfile)
    load(outfile, envir = parent.env(environment()))
    return(invisible(get(obj_name, envir = parent.env(environment()))))
  }
  
  if (refit_models) {
    message("\n-- Sensitivity run: ", label,
            " (refit_models = TRUE)",
            " (", nrow(filtered_dat), " rows)",
            " | chains=", chains, " adapt=", adapt,
            " burnin=", burnin, " iter=", iter, " thin=", thin, " --")
  } else {
    message("\n-- Sensitivity run: ", label,
            " (no cached file found, fitting from scratch)",
            " (", nrow(filtered_dat), " rows)",
            " | chains=", chains, " adapt=", adapt,
            " burnin=", burnin, " iter=", iter, " thin=", thin, " --")
  }
  
  jd_sens <- build_jags_data(filtered_dat)
  
  # Apply any prior overrides before passing to JAGS
  if (length(jags_data_override) > 0) {
    override_names <- names(jags_data_override)
    message("  Overriding jags_data elements: ",
            paste(override_names, collapse = ", "))
    for (nm in override_names) {
      jd_sens[[nm]] <- jags_data_override[[nm]]
    }
  }
  
  set.seed(111)
  mod_sens <- jags.model(model_file,
                         data     = jd_sens,
                         inits    = make_inits(jd_sens$n_type, jd_sens$n_coord),
                         n.chains = chains,
                         n.adapt  = adapt)
  update(mod_sens, n.iter = burnin)
  
  mod_samp_out <- coda.samples(mod_sens,
                               variable.names = variables_to_sample,
                               n.iter = iter, thin = thin)
  
  mod_summary_out <- MCMCvis::MCMCsummary(mod_samp_out)
  n_bad <- sum(mod_summary_out$Rhat > 1.1, na.rm = TRUE)
  if (n_bad > 0) {
    warning(label, ": ", n_bad, " parameters with R-hat > 1.1")
  } else {
    message(label, ": convergence OK (all R-hat <= 1.1)")
  }
  
  # Save under the caller-specified name so loading multiple sensitivity
  # outputs into the same session does not cause object name conflicts
  assign(obj_name,              mod_samp_out)
  assign(paste0("mod_summary_", sub("mod_samp_", "", obj_name)), mod_summary_out)
  save(list = c(obj_name,
                paste0("mod_summary_", sub("mod_samp_", "", obj_name)),
                "variables_to_sample", "sc_df", "jd_sens"),
       file = outfile,
       envir = environment())
  
  invisible(mod_samp_out)
}

# drop_sparse_types() removes rows belonging to transition types with fewer
# than min_n observations before data are passed to build_jags_data().  This
# is necessary for the BD sensitivity subsets, where P->O (A) is absent or
# nearly so, leaving beta_r_intercept and m_val for that type unidentified.
#
# Crucially, the JAGS indexing remains internally consistent after dropping:
#   - final_type is re-factored from the types actually present, producing a
#     dense integer sequence (e.g. 1-7 when one type is dropped).
#   - type_vals and n_type are computed from the same reduced set, so the
#     model prior loop (for i in type_vals) covers exactly the types present.
#   - make_inits() uses the updated n_type and n_coord, so initial value
#     vectors (beta_r_intercept, m_val, rand_coord) match exactly.
#   - final_coord and n_coord are likewise recomputed from remaining rows.
# The dropped type simply does not appear in the sensitivity model or its
# output.  sens_summary() handles this correctly by filling NA for any
# transition type present in the primary model but absent in the sensitivity
# run, so the comparison figure shows a blank at the correct row position
# rather than a misaligned point.
drop_sparse_types <- function(dat, min_n = 20) {
  type_vec    <- paste(dat$lu_start, dat$lu_end, dat$MGMT)
  type_counts <- table(type_vec)
  sparse      <- names(type_counts[type_counts < min_n])
  if (length(sparse) == 0) return(dat)
  message("Dropping sparse transition type(s) (n < ", min_n, "): ",
          paste(sparse, collapse = ", "))
  dat[!type_vec %in% sparse, ]
}

# Compute beta_max_* priors for the BD sensitivity runs.  These are centred
# on the primary posterior means with precision matching the primary posterior,
# reflecting that the SOC-environment scaling captured by beta_max_* is well
# characterised by the full dataset and is not the quantity under investigation
# in either BD sensitivity test.  Computing these at runtime from mod_samp
# ensures they remain current if the primary model is rerun.
prim_max_summary <- MCMCvis::MCMCsummary(
  mod_samp,
  params = c("beta_max_intercept", "beta_max_map", "beta_max_map_mat",
             "beta_max_mat",       "beta_max_depth",
             "beta_max_clay",      "beta_max_nitrogen")
)
bd_tau <- function(param) 1 / prim_max_summary[param, "sd"]^2

bd_max_priors <- list(
  prior_max_intercept     = prim_max_summary["beta_max_intercept", "mean"],
  prior_tau_max_intercept = bd_tau("beta_max_intercept"),
  prior_max_map           = prim_max_summary["beta_max_map",       "mean"],
  prior_tau_max_map       = bd_tau("beta_max_map"),
  prior_max_map_mat       = prim_max_summary["beta_max_map_mat",   "mean"],
  prior_tau_max_map_mat   = bd_tau("beta_max_map_mat"),
  prior_max_mat           = prim_max_summary["beta_max_mat",       "mean"],
  prior_tau_max_mat       = bd_tau("beta_max_mat"),
  prior_max_depth         = prim_max_summary["beta_max_depth",     "mean"],
  prior_tau_max_depth     = bd_tau("beta_max_depth"),
  prior_max_clay          = prim_max_summary["beta_max_clay",      "mean"],
  prior_tau_max_clay      = bd_tau("beta_max_clay"),
  prior_max_nitrogen      = prim_max_summary["beta_max_nitrogen",  "mean"],
  prior_tau_max_nitrogen  = bd_tau("beta_max_nitrogen")
)
rm(prim_max_summary)

# ── (a) Bulk density gap-filling sensitivity ──────────────────────────────────
# Excludes rows where either Ctrl_BD or Aban_BD was model-predicted rather
# than published.  Tests whether gap-filled BD values affect the conclusions.
# drop_sparse_types() is applied because P->O (A) has only 9 observations
# in this subset — insufficient to identify beta_r_intercept and m_val for
# that type.  It is excluded from this sensitivity run and appears as a blank
# in the comparison figure (see sens_summary() and draw_sensitivity_overlay()).
mod_samp_bd <- run_sensitivity(
  filter(c_dat, !bd_gap_filled) %>% drop_sparse_types(),
  label              = "BD gap-filling",
  outfile            = "./data_output/bayes_output_bd_sens.RData",
  obj_name           = "mod_samp_bd",
  model_file         = "./scripts/bayes_script.txt",
  jags_data_override = bd_max_priors
)

# ── (b) Management uncertainty sensitivity ────────────────────────────────────
# Excludes rows where management was originally coded "?" (assigned to "N"
# in the primary analysis as a conservative assumption).
mod_samp_mgmt <- run_sensitivity(
  filter(c_dat, !MGMT_is_uncertain),
  label    = "MGMT uncertainty",
  outfile  = "./data_output/bayes_output_mgmt_sens.RData",
  obj_name = "mod_samp_mgmt"
)

# ── (c) SOM-to-SOC conversion sensitivity ─────────────────────────────────────
# Excludes rows where the SOM→SOC van Bemmelen conversion (factor 0.5) was
# applied, testing whether the fixed conversion factor introduces bias.
mod_samp_som <- run_sensitivity(
  filter(c_dat, !som_converted),
  label    = "SOM conversion",
  outfile  = "./data_output/bayes_output_som_sens.RData",
  obj_name = "mod_samp_som"
)

# ── (d) Original bulk density sensitivity ─────────────────────────────────────
# Uses the original reported BD values (Ctrl/Aban_BD_raw) for all rows,
# including profiles flagged as problematic (repeated constant values).
# Rows where the original study did not report BD at all (BD_raw is NA)
# are dropped since no alternative BD is available for stock computation.
# drop_sparse_types() is applied because P->O (A) is entirely absent from
# this subset — it would leave beta_r_intercept and m_val for that type
# completely unidentified.  It is excluded from this sensitivity run and
# appears as a blank in the comparison figure.
#
# The _orig_bd carbon density columns are renamed to match the column names
# that build_jags_data() reads (Ctrl/Aban_c_dens_av) in a new object
# (c_dat_origbd), leaving c_dat itself unmodified.  This ensures the
# deduplication and ordering inside build_jags_data() is handled identically
# to every other sensitivity run.
c_dat_origbd <- c_dat %>%
  filter(!is.na(Ctrl_c_dens_av_orig_bd), !is.na(Aban_c_dens_av_orig_bd)) %>%
  mutate(
    Ctrl_c_dens_av = Ctrl_c_dens_av_orig_bd,
    Aban_c_dens_av = Aban_c_dens_av_orig_bd
  ) %>%
  drop_sparse_types()

origbd_outfile <- "./data_output/bayes_output_origbd_sens.RData"

if (!refit_models && file.exists(origbd_outfile)) {
  message("\nLoading cached sensitivity output (original BD) from ", origbd_outfile)
  load(origbd_outfile)
} else {
  if (refit_models) {
    message("\n-- Sensitivity run: original BD (refit_models = TRUE) (",
            nrow(c_dat_origbd), " rows) ------")
  } else {
    message("\n-- Sensitivity run: original BD (no cached file found, fitting from scratch) (",
            nrow(c_dat_origbd), " rows) ------")
  }
  jd_origbd <- build_jags_data(c_dat_origbd)
  
  # Substitute the primary-posterior-centred beta_max_* priors before fitting,
  # consistent with the BD gap-filling run above.
  for (nm in names(bd_max_priors)) jd_origbd[[nm]] <- bd_max_priors[[nm]]
  
  set.seed(111)
  mod_origbd <- jags.model("./scripts/bayes_script.txt",
                           data     = jd_origbd,
                           inits    = make_inits(jd_origbd$n_type, jd_origbd$n_coord),
                           n.chains = n_chains,
                           n.adapt  = n_adapt)
  update(mod_origbd, n.iter = n_burnin)
  mod_samp_origbd <- coda.samples(mod_origbd,
                                  variable.names = variables_to_sample,
                                  n.iter = n_iter, thin = n_thin)
  mod_summary_origbd <- MCMCvis::MCMCsummary(mod_samp_origbd)
  n_bad_origbd <- sum(mod_summary_origbd$Rhat > 1.1, na.rm = TRUE)
  if (n_bad_origbd > 0) {
    warning("Original-BD run: ", n_bad_origbd, " parameters with R-hat > 1.1")
  } else {
    message("Original-BD run: convergence OK (all R-hat <= 1.1)")
  }
  save(file = origbd_outfile,
       mod_samp_origbd, variables_to_sample, sc_df,
       jd_origbd, mod_summary_origbd)
}




# =============================================================================
# SECTION 6 — WEAK-PRIOR SENSITIVITY ANALYSIS
# =============================================================================

# ── (d) Weak-prior sensitivity ────────────────────────────────────────────────
# Tests whether the informative priors from SoilGrids / Sanderman et al.
# (used for beta_init_* and beta_max_*) materially influence the posteriors.
# The weak-prior model replaces those gridded-product priors with:
#   - beta_init_* intercepts: dnorm(obs_log_soc_mean, 4) — anchored to the
#     observed mean log SOC density to keep intercepts in a plausible range;
#     precision = 4 (SD = 0.5) is intentionally broad
#   - All slope priors: dnorm(0, 1) — centred at zero on the standardised
#     covariate scale, weakly informative
#   - Sign constraints on beta_max_map T(0,) and beta_max_mat T(,0) are
#     RETAINED — these encode a basic physical expectation (more precipitation
#     → more SOC; more heat → less SOC) and are independent of any gridded
#     product
#
# All other model structure — the Chapman-Richards trajectory, coordinate-level
# random effects on r_val, hierarchical r intercepts — is identical to the
# primary model.  The full c_dat dataset is used (no rows filtered), so any
# difference in posteriors reflects prior sensitivity only.

# Write the weak-prior model text file
sink("./scripts/bayes_script_weakprior.txt")
cat("
  model {

  ### Likelihood - identical to primary model

  for(i in 1:n_init_obs){
    c_init_obs[i] ~ dnorm(c_init_true[init_id[i]],
                          c_init_tau[init_lu_for_tau[i]])
  }

  for(i in 1:n_final_obs){
    c_final_obs[i] ~ dnorm(c_final_true[final_id[i]],
                           c_final_tau[final_lu_for_tau[i]])
  }


  # Process model: reference carbon - identical structure to primary
  for(i in 1:n_init_id){
    c_init_true[i] <- beta_init_intercept[init_lu_start[i]] +
      beta_init_depth[init_lu_start[i]]      * init_depth[i] +
      beta_init_clay[init_lu_start[i]]       * init_clay[i] +
      beta_init_nitrogen[init_lu_start[i]]   * init_nitrogen[i] +
      beta_init_map_mat[init_lu_start[i]]    * init_map[i] * init_mat[i] +
      beta_init_irrigation[init_lu_start[i]] * init_irrigation[i] +
      beta_init_map[init_lu_start[i]]        * init_map[i] +
      beta_init_mat[init_lu_start[i]]        * init_mat[i]
  }


  # Process model: recovery trajectory - identical structure to primary,
  # including coordinate-level random effect on r_val
  for(i in 1:n_final_id){

    c_final_true[i] <- log(exp(c_init_true[init_compare_id[i]]) +
                             (exp(c_max[i]) - exp(c_init_true[init_compare_id[i]])) *
                             (1 - exp(-final_age[i] / exp(r_val[i]))) ^ exp(m_val[final_type[i]]))

    r_val[i] <- beta_r_intercept[final_type[i]] +
      beta_r_depth[final_lu_start[i]]    * final_depth[i] +
      beta_r_clay[final_lu_start[i]]     * final_clay[i] +
      beta_r_nitrogen[final_lu_start[i]] * final_nitrogen[i] +
      beta_r_map_mat[final_lu_start[i]]  * final_map[i] * final_mat[i] +
      beta_r_map[final_lu_start[i]]      * final_map[i] +
      beta_r_mat[final_lu_start[i]]      * final_mat[i] +
      rand_coord[final_coord[i]]

    c_max[i] <- beta_max_intercept +
      beta_max_depth    * final_depth[i] +
      beta_max_clay     * final_clay[i] +
      beta_max_nitrogen * final_nitrogen[i] +
      beta_max_map_mat  * final_map[i] * final_mat[i] +
      beta_max_map      * final_map[i] +
      beta_max_mat      * final_mat[i]
  }


  ### Priors - weakly informative, no gridded-product information

  # Reference carbon intercepts: anchored to observed mean log SOC to prevent
  # the observed mean log SOC density to keep intercepts in a plausible range;
  # precision = 4 (SD = 0.5) is intentionally broad. Slope priors:
  # dnorm(0, 1) - centred at zero, uninformative about direction.
  for(i in 1:2){
    beta_init_intercept[i] ~ dnorm(obs_log_soc_mean, 4)
    beta_init_depth[i]     ~ dnorm(0, 1)
    beta_init_clay[i]      ~ dnorm(0, 1)
    beta_init_nitrogen[i]  ~ dnorm(0, 1)
    beta_init_map_mat[i]   ~ dnorm(0, 1)
    beta_init_map[i]       ~ dnorm(0, 1)
    beta_init_mat[i]       ~ dnorm(0, 1)
  }

  for(i in lu_start_vals){
    beta_init_irrigation[i] ~ dnorm(beta_mean_init_irrigation, tau_mean_init_irrigation)
  }

  # r covariate priors - unchanged from primary (already weakly informative)
  for(i in 1:2){
    beta_r_depth[i]    ~ dnorm(0, 0.25)
    beta_r_clay[i]     ~ dnorm(0, 0.25)
    beta_r_nitrogen[i] ~ dnorm(0, 0.25)
    beta_r_map_mat[i]  ~ dnorm(0, 0.25)
    beta_r_map[i]      ~ dnorm(0, 0.25)
    beta_r_mat[i]      ~ dnorm(0, 0.25)
  }

  # r intercepts - hierarchical prior, unchanged
  for(i in type_vals){
    beta_r_intercept[i] ~ dnorm(beta_mean_r_intercept, tau_mean_r_intercept)
    m_val[i] ~ dnorm(1.1, 5)
  }

  beta_mean_r_intercept ~ dnorm(3, 1)
  tau_mean_r_intercept <- pow(sd_mean_r_intercept, -2)
  sd_mean_r_intercept ~ dunif(0, 5)

  # Max SOC: intercept anchored to observed mean log SOC.
  # Sign constraints on map and mat RETAINED - physical, not data-derived.
  beta_max_intercept ~ dnorm(obs_log_soc_mean, 4)
  beta_max_depth     ~ dnorm(0, 1)
  beta_max_clay      ~ dnorm(0, 1)
  beta_max_nitrogen  ~ dnorm(0, 1)
  beta_max_map_mat   ~ dnorm(0, 1)
  beta_max_map       ~ dnorm(0, 1)T(0, )
  beta_max_mat       ~ dnorm(0, 1)T(, 0)

  # Irrigation hyperprior - unchanged
  beta_mean_init_irrigation ~ dnorm(0, 0.5)
  tau_mean_init_irrigation <- pow(sd_mean_init_irrigation, -2)
  sd_mean_init_irrigation ~ dunif(0, 5)

  # Coordinate-level random effect on r - identical to primary model
  for(i in 1:(n_coord-1)){
    rand_coord[i] ~ dnorm(0, tau_coord)
  }
  rand_coord[n_coord] <- -sum(rand_coord[1:(n_coord-1)])
  tau_coord <- pow(coord_sd, -2)
  coord_sd ~ dunif(0, 5)

  # Observation-level SDs - unchanged
  for(i in 1:2){
    c_init_tau[i]  <- pow(c_init_sd[i],  -2)
    c_init_sd[i]   ~ dunif(0, 5)
    c_final_tau[i] <- pow(c_final_sd[i], -2)
    c_final_sd[i]  ~ dunif(0, 5)
  }

  # Derived quantities - identical to primary model
  for(i in 1:n_final_id){
    c_final_norand[i] <- log(exp(c_init_true[init_compare_id[i]]) +
                               (exp(c_max_norand[i]) - exp(c_init_true[init_compare_id[i]])) *
                               (1 - exp(-final_age[i] / exp(r_val_norand[i]))) ^ exp(m_val[final_type[i]]))

    r_val_norand[i] <- beta_r_intercept[final_type[i]] +
      beta_r_depth[final_lu_start[i]]    * final_depth[i] +
      beta_r_clay[final_lu_start[i]]     * final_clay[i] +
      beta_r_nitrogen[final_lu_start[i]] * final_nitrogen[i] +
      beta_r_map_mat[final_lu_start[i]]  * final_map[i] * final_mat[i] +
      beta_r_map[final_lu_start[i]]      * final_map[i] +
      beta_r_mat[final_lu_start[i]]      * final_mat[i]

    c_max_norand[i] <- beta_max_intercept +
      beta_max_depth    * final_depth[i] +
      beta_max_clay     * final_clay[i] +
      beta_max_nitrogen * final_nitrogen[i] +
      beta_max_map_mat  * final_map[i] * final_mat[i] +
      beta_max_map      * final_map[i] +
      beta_max_mat      * final_mat[i]
  }

}
    ", fill = TRUE)
sink()

# Build jags_data for the weak-prior model.
# Identical to jags_data but with all prior_* entries stripped (they are not
# referenced in the weak-prior model text) and obs_log_soc_mean added as the
# anchor for the intercept priors.
prior_keys     <- names(jags_data)[grepl("^prior_", names(jags_data))]
jags_data_wp   <- jags_data[setdiff(names(jags_data), prior_keys)]
jags_data_wp$obs_log_soc_mean <- mean(jags_data$c_init_obs)

weakprior_outfile <- "./data_output/bayes_output_weakprior.RData"

if (!refit_models && file.exists(weakprior_outfile)) {
  message("\nLoading cached sensitivity output (weak priors) from ", weakprior_outfile)
  load(weakprior_outfile)
} else {
  if (refit_models) {
    message("\n-- Sensitivity run: weak priors (refit_models = TRUE) (",
            nrow(c_dat), " rows, full dataset) ------")
  } else {
    message("\n-- Sensitivity run: weak priors (no cached file found, fitting from scratch) (",
            nrow(c_dat), " rows, full dataset) ------")
  }
  set.seed(111)
  mod_wp <- jags.model("./scripts/bayes_script_weakprior.txt",
                       data     = jags_data_wp,
                       inits    = make_inits(jags_data_wp$n_type, jags_data_wp$n_coord),
                       n.chains = n_chains,
                       n.adapt  = n_adapt)
  update(mod_wp, n.iter = n_burnin)
  
  mod_samp_wp <- coda.samples(mod_wp,
                              variable.names = variables_to_sample,
                              n.iter = n_iter, thin = n_thin)
  
  mod_summary_wp <- MCMCvis::MCMCsummary(mod_samp_wp)
  hist(mod_summary_wp$Rhat,
       main = "R-hat distribution -- weak-prior run",
       xlab = "R-hat")
  n_bad_wp <- sum(mod_summary_wp$Rhat > 1.1, na.rm = TRUE)
  if (n_bad_wp > 0) {
    warning("Weak-prior run: ", n_bad_wp, " parameters with R-hat > 1.1")
  } else {
    message("Weak-prior run: convergence OK (all R-hat <= 1.1)")
  }
  
  save(file = weakprior_outfile,
       mod_samp_wp, variables_to_sample, sc_df, jags_data_wp, mod_summary_wp)
}

# SECTION 7 — SENSITIVITY COMPARISON FIGURE
# =============================================================================
# Overlays effect size posteriors from the primary analysis and all five
# sensitivity runs using the same two-panel MCMCplot layout as the main
# "effect sizes" figure.  Filled circles with solid CIs = primary analysis;
# open circles with dashed CIs = sensitivity run.  Each sensitivity run is
# shown in a separate pair of panels so the figures remain readable.
#
# This section can be run independently after loading the four saved .RData
# files if the user wants to regenerate figures without re-running the MCMC.

# Helper: extract posterior summaries for a set of parameters from a
# coda.samples object, returned in the same row order as MCMCplot uses.
# sens_summary() extracts posterior summaries in the same row order that
# MCMCplot uses for the primary model.  It uses the primary model (mod_samp)
# to determine the full set of expected rows, then fills NA for any rows
# absent from the sensitivity run (e.g. a transition type not present in
# that subset).  This keeps the overlay correctly aligned with the primary
# MCMCplot even when n_type differs between runs.
sens_summary <- function(samp, params) {
  # Full row structure from the primary model
  sm_primary <- MCMCvis::MCMCsummary(mod_samp, params = params,
                                     probs = c(0.025, 0.5, 0.975))
  ordered_names <- unlist(lapply(params, function(p)
    grep(paste0("^", p, "(\\[|$)"), rownames(sm_primary), value = TRUE)))
  template <- sm_primary[rev(ordered_names), , drop = FALSE]
  template[, ] <- NA   # clear values; keep row names as alignment template
  
  # Fill in whatever the sensitivity run has
  sm_sens <- MCMCvis::MCMCsummary(samp, params = params,
                                  probs = c(0.025, 0.5, 0.975))
  common  <- intersect(rownames(template), rownames(sm_sens))
  template[common, ] <- sm_sens[common, ]
  template
}

draw_sensitivity_overlay <- function(sm_sens, n_rows, col_vec, offset = 0.18) {
  # Adds open-circle + dashed-CI points for one sensitivity run on top of
  # an existing MCMCplot (which must already be drawn on the active device).
  # Rows with NA (parameter absent from this sensitivity subset) are skipped.
  for (j in seq_len(n_rows)) {
    if (anyNA(sm_sens[j, c("2.5%", "50%", "97.5%")])) next
    y_s   <- j - offset
    col_j <- col_vec[n_rows - j + 1]
    lines(x = c(sm_sens[j, "2.5%"], sm_sens[j, "97.5%"]),
          y = c(y_s, y_s), col = col_j, lwd = 1, lty = 3)
    points(x = sm_sens[j, "50%"], y = y_s,
           pch = 21, bg = "white", col = col_j, cex = 1.1)
  }
}

params_to_plot1 <- c("beta_init_intercept", "beta_init_depth",
                     "beta_init_clay", "beta_init_nitrogen",
                     "beta_init_irrigation", "beta_init_map",
                     "beta_init_mat", "beta_init_map_mat",
                     "beta_max_intercept", "beta_max_depth",
                     "beta_max_clay", "beta_max_nitrogen",
                     "beta_max_map", "beta_max_mat", "beta_max_map_mat")

params_to_plot2 <- c("beta_r_intercept", "beta_r_depth",
                     "beta_r_clay", "beta_r_nitrogen",
                     "beta_r_map", "beta_r_mat", "beta_r_map_mat",
                     "m_val")

col_left  <- c(rep(c("#d95f02", "#7570b3"), times = 8), rep("black", 7))
col_right <- c(rep(c("#d95f02", "#7570b3"), each  = 4),
               rep(c("#d95f02", "#7570b3"), times = 6),
               rep(c("#d95f02", "#7570b3"), each  = 4))

sens_list <- list(
  list(samp  = mod_samp_bd,
       label = "BD gap-filling sensitivity",
       files = c("./figures/effect_sizes_bd_sens.png",
                 "./figures/effect_sizes_bd_sens.pdf")),
  list(samp  = mod_samp_mgmt,
       label = "Management uncertainty sensitivity",
       files = c("./figures/effect_sizes_mgmt_sens.png",
                 "./figures/effect_sizes_mgmt_sens.pdf")),
  list(samp  = mod_samp_som,
       label = "SOM conversion sensitivity",
       files = c("./figures/effect_sizes_som_sens.png",
                 "./figures/effect_sizes_som_sens.pdf")),
  list(samp  = mod_samp_wp,
       label = "Weak-prior sensitivity",
       files = c("./figures/effect_sizes_weakprior_sens.png",
                 "./figures/effect_sizes_weakprior_sens.pdf")),
  list(samp  = mod_samp_origbd,
       label = "Original-BD sensitivity",
       files = c("./figures/effect_sizes_origbd_sens.png",
                 "./figures/effect_sizes_origbd_sens.pdf"))
)

for (s in sens_list) {
  sm_left  <- sens_summary(s$samp, params_to_plot1)
  sm_right <- sens_summary(s$samp, params_to_plot2)
  
  for (dev_type in c("png", "pdf")) {
    if (dev_type == "png")
      png(s$files[1], width = 12, height = 8, unit = "in", res = 440)
    else
      pdf(s$files[2], width = 12, height = 8)
    
    par(mfrow = c(1, 2))
    
    # Left panel: beta_init + beta_max
    MCMCvis::MCMCplot(mod_samp,
                      params = params_to_plot1,
                      labels = c("Intercept","","Depth","","Clay","",
                                 "Nitrogen","","Irrigation","",
                                 "Precipitation","","Temperature","",
                                 "Precip*Temp","",
                                 "Intercept","Depth","Clay","Nitrogen",
                                 "Precipitation","Temperature","Precip*Temp"),
                      col = col_left,
                      mar = c(5.1, 6.6, 1.1, 2.1),
                      xlim = c(-2, 2))
    draw_sensitivity_overlay(sm_left, nrow(sm_left), col_left)
    segments(y0 = 7.5 + seq(0, 16, by = 2), y1 = 7.5 + seq(0, 16, by = 2),
             x0 = rep(-2, 10), x1 = rep(2, 10), lty = 2)
    text(-4.0, 7.5 / 2, "Maximum", srt = 90, xpd = TRUE, cex = 1.3)
    segments(x0 = -3.8, y0 = 0.5, y1 = 7, xpd = TRUE)
    text(-4.0, 15.5, "Agricultural land cover", srt = 90, xpd = TRUE, cex = 1.3)
    segments(x0 = -3.8, y0 = 8, y1 = 23.5, xpd = TRUE)
    text("Cropland", col = "#d95f02", cex = 1, x = 1.5, y = 22.95, xpd = TRUE)
    text("Pasture",  col = "#7570b3", cex = 1, x = 1.5, y = 22.05, xpd = TRUE)
    
    # Right panel: beta_r + m_val
    MCMCvis::MCMCplot(mod_samp,
                      params = params_to_plot2,
                      labels = c("C\u2192F (A)","C\u2192F (N)","C\u2192O (A)","C\u2192O (N)",
                                 "P\u2192F (A)","P\u2192F (N)","P\u2192O (A)","P\u2192O (N)",
                                 "Depth","","Clay","","Nitrogen","",
                                 "Precipitation","","Temperature","","Precip*Temp","",
                                 "C\u2192F (A)","C\u2192F (N)","C\u2192O (A)","C\u2192O (N)",
                                 "P\u2192F (A)","P\u2192F (N)","P\u2192O (A)","P\u2192O (N)"),
                      col  = col_right,
                      mar  = c(5.1, 8, 1.1, 2.1),
                      xlim = c(-5, 10))
    draw_sensitivity_overlay(sm_right, nrow(sm_right), col_right)
    abline(h = 0.5 + seq(0, 12, by = 2) + 8, lty = 2)
    text(-12.6, 4.25,  "m value",   font = 3, srt = 90, xpd = TRUE, cex = 1.3)
    segments(x0 = -12.0, y0 = 0.25, y1 = 8,     xpd = TRUE)
    text(-10.2, 4.25,  "Intercept",              srt = 90, xpd = TRUE, cex = 1.3)
    segments(x0 =  -9.6, y0 = 0.25, y1 = 8,     xpd = TRUE)
    text(-12.6, 18.5,  "r value",   font = 3, srt = 90, xpd = TRUE, cex = 1.3)
    segments(x0 = -12.0, y0 = 9,    y1 = 28.25, xpd = TRUE)
    text(-10.2, 24.25, "Intercept",              srt = 90, xpd = TRUE, cex = 1.3)
    segments(x0 =  -9.6, y0 = 20.75, y1 = 28.25, xpd = TRUE)
    
    # Legend distinguishing primary (filled) from sensitivity (open)
    legend(x = -5, y = 28.5, xpd = TRUE,
           legend = c("Primary analysis", s$label),
           pch    = c(19, 21),
           lty    = c(1, 3),
           pt.bg  = c(NA, "white"),
           col    = "grey30",
           bty    = "n", cex = 1, pt.cex = c(1.3, 1.1))
    
    dev.off()
    message("Saved ", s$files[which(c("png","pdf") == dev_type)])
  }
}

rm(s, sm_left, sm_right)

# =============================================================================
# COMBINED BD SENSITIVITY FIGURE
# =============================================================================
# Shows both BD sensitivity runs in a single figure using the same shapes
# as the omnibus figure for consistency:
# Primary is redrawn last (on top) at y = j; sensitivity runs sit below:
#   Primary analysis — filled circle (pch 19), solid,  y = j  (redrawn last)
#   BD gap-filling   — open circle   (pch 21), dotted (lty 3), offset -0.18
#   Original-BD      — open dntri   (pch 25), dashed (lty 2), offset -0.36
# This allows direct visual comparison of how the two different approaches
# to the BD sensitivity question affect the posteriors relative to one another
# and relative to the primary analysis.

draw_second_overlay <- function(sm_sens, n_rows, col_vec, offset = 0.36) {
  # Used for the Original-BD run in the combined BD figure.
  # Matches the omnibus figure: pch = 25 (down-triangle), lty = 2 (dashed).
  # Offset is negative (below primary) to match omnibus layout convention.
  for (j in seq_len(n_rows)) {
    if (anyNA(sm_sens[j, c("2.5%", "50%", "97.5%")])) next
    y_s   <- j - offset
    col_j <- col_vec[n_rows - j + 1]
    lines(x = c(sm_sens[j, "2.5%"], sm_sens[j, "97.5%"]),
          y = c(y_s, y_s), col = col_j, lwd = 1, lty = 2)
    points(x = sm_sens[j, "50%"], y = y_s,
           pch = 25, bg = "white", col = col_j, cex = 1.1)
  }
}

# Redraws the primary analysis points and CIs on top of all sensitivity
# overlays so they are not obscured.  Uses the same colour vector as MCMCplot.
draw_primary_on_top <- function(sm_primary, n_rows, col_vec, cex = 0.9) {
  # Redraws only the primary analysis point (not the CI line) on top of all
  # sensitivity overlays.  MCMCplot has already drawn the CI line correctly;
  # redrawing it would produce a double line at a different lwd.
  for (j in seq_len(n_rows)) {
    if (anyNA(sm_primary[j, c("2.5%", "50%", "97.5%")])) next
    col_j <- col_vec[n_rows - j + 1]
    points(x = sm_primary[j, "50%"], y = j,
           pch = 19, col = col_j, cex = cex)
  }
}

sm_bd_left     <- sens_summary(mod_samp_bd,     params_to_plot1)
sm_bd_right    <- sens_summary(mod_samp_bd,     params_to_plot2)
sm_origbd_left  <- sens_summary(mod_samp_origbd, params_to_plot1)
sm_origbd_right <- sens_summary(mod_samp_origbd, params_to_plot2)

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png") {
    png("./figures/effect_sizes_bd_combined.png",
        width = 12, height = 8, unit = "in", res = 440)
  } else {
    pdf("./figures/effect_sizes_bd_combined.pdf",
        width = 12, height = 8)
  }
  
  par(mfrow = c(1, 2))
  
  # Left panel
  MCMCvis::MCMCplot(mod_samp,
                    params = params_to_plot1,
                    labels = c("Intercept","","Depth","","Clay","",
                               "Nitrogen","","Irrigation","",
                               "Precipitation","","Temperature","",
                               "Precip*Temp","",
                               "Intercept","Depth","Clay","Nitrogen",
                               "Precipitation","Temperature","Precip*Temp"),
                    col = col_left,
                    sz_med = 0.75,
                    mar = c(5.1, 6.6, 1.1, 2.1),
                    xlim = c(-2, 2))
  draw_sensitivity_overlay(sm_bd_left,    nrow(sm_bd_left),    col_left)
  draw_second_overlay(     sm_origbd_left, nrow(sm_origbd_left), col_left)
  sm_prim_bd_left <- sens_summary(mod_samp, params_to_plot1)
  draw_primary_on_top(sm_prim_bd_left, nrow(sm_prim_bd_left), col_left)
  segments(y0 = 7.5 + seq(0, 16, by = 2), y1 = 7.5 + seq(0, 16, by = 2),
           x0 = rep(-2, 10), x1 = rep(2, 10), lty = 2)
  text(-4.0, 7.5 / 2, "Maximum", srt = 90, xpd = TRUE, cex = 1.3)
  segments(x0 = -3.8, y0 = 0.5, y1 = 7, xpd = TRUE)
  text(-4.0, 15.5, "Agricultural land cover", srt = 90, xpd = TRUE, cex = 1.3)
  segments(x0 = -3.8, y0 = 8, y1 = 23.5, xpd = TRUE)
  text("Cropland", col = "#d95f02", cex = 1, x = 1.5, y = 22.95, xpd = TRUE)
  text("Pasture",  col = "#7570b3", cex = 1, x = 1.5, y = 22.05, xpd = TRUE)
  
  # Right panel
  MCMCvis::MCMCplot(mod_samp,
                    params = params_to_plot2,
                    labels = c("C\u2192F (A)","C\u2192F (N)","C\u2192O (A)","C\u2192O (N)",
                               "P\u2192F (A)","P\u2192F (N)","P\u2192O (A)","P\u2192O (N)",
                               "Depth","","Clay","","Nitrogen","",
                               "Precipitation","","Temperature","","Precip*Temp","",
                               "C\u2192F (A)","C\u2192F (N)","C\u2192O (A)","C\u2192O (N)",
                               "P\u2192F (A)","P\u2192F (N)","P\u2192O (A)","P\u2192O (N)"),
                    col  = col_right,
                    sz_med = 0.75,
                    mar  = c(5.1, 8, 1.1, 2.1),
                    xlim = c(-5, 10))
  draw_sensitivity_overlay(sm_bd_right,    nrow(sm_bd_right),    col_right)
  draw_second_overlay(     sm_origbd_right, nrow(sm_origbd_right), col_right)
  sm_prim_bd_right <- sens_summary(mod_samp, params_to_plot2)
  draw_primary_on_top(sm_prim_bd_right, nrow(sm_prim_bd_right), col_right)
  abline(h = 0.5 + seq(0, 12, by = 2) + 8, lty = 2)
  text(-12.6, 4.25,  "m value",   font = 3, srt = 90, xpd = TRUE, cex = 1.3)
  segments(x0 = -12.0, y0 = 0.25, y1 = 8,     xpd = TRUE)
  text(-10.2, 4.25,  "Intercept",              srt = 90, xpd = TRUE, cex = 1.3)
  segments(x0 =  -9.6, y0 = 0.25, y1 = 8,     xpd = TRUE)
  text(-12.6, 18.5,  "r value",   font = 3, srt = 90, xpd = TRUE, cex = 1.3)
  segments(x0 = -12.0, y0 = 9,    y1 = 28.25, xpd = TRUE)
  text(-10.2, 24.25, "Intercept",              srt = 90, xpd = TRUE, cex = 1.3)
  segments(x0 =  -9.6, y0 = 20.75, y1 = 28.25, xpd = TRUE)
  
  legend(x = -5, y = 28.5, xpd = TRUE,
         legend = c("Primary analysis",
                    "BD gap-filling sensitivity",
                    "Original-BD sensitivity"),
         pch    = c(19, 21, 25),
         lty    = c(1, 3, 2),
         pt.bg  = c(NA, "white", "white"),
         col    = "grey30",
         bty    = "n", cex = 1, pt.cex = c(1.3, 1.1, 1.1))
  
  dev.off()
  message("Saved ./figures/effect_sizes_bd_combined.", dev_type)
}


# =============================================================================
# OMNIBUS SENSITIVITY FIGURE
# =============================================================================
# All five sensitivity runs overlaid on the primary analysis in a single figure.
# Each sensitivity run uses a distinct point shape and line type; colour follows
# the same cropland / pasture / maximum SOC scheme as all other figures.
# The figure is taller than the individual sensitivity figures to accommodate
# the larger legend.
#
# Sensitivity runs and their symbols:
# Sensitivity runs stacked below the primary (all offsets negative); the
# primary is redrawn last at y = j so it sits visually on top:
#   Primary analysis — filled circle (pch 19), solid line, y = j  (redrawn last)
#   Weak priors      — uptri   (pch 24), dashed line (lty 2), offset -0.10
#   BD gap-filling   — circle  (pch 21), dotted line (lty 3), offset -0.20
#   Original-BD      — dntri   (pch 25), dashed line (lty 2), offset -0.30
#   MGMT uncertainty — square  (pch 22), dotted line (lty 3), offset -0.40
#   SOM conversion   — diamond (pch 23), dotted line (lty 3), offset -0.50

# Offsets are all negative so all five sensitivity runs sit below the primary,
# which is redrawn last at y = j (integer) so it appears on top.
omnibus_specs <- list(
  list(samp = mod_samp_wp,     label = "Weak priors",       pch = 24, lty = 2, offset = -0.15),
  list(samp = mod_samp_bd,     label = "BD gap-filling",    pch = 21, lty = 3, offset = -0.2375),
  list(samp = mod_samp_origbd, label = "Original BD",       pch = 25, lty = 2, offset = -0.3250),
  list(samp = mod_samp_mgmt,   label = "MGMT uncertainty",  pch = 22, lty = 3, offset = -0.4125),
  list(samp = mod_samp_som,    label = "SOM conversion",    pch = 23, lty = 3, offset = -0.5)
)

draw_omnibus_overlay <- function(sm_sens, n_rows, col_vec,
                                 offset, pch_val, lty_val) {
  for (j in seq_len(n_rows)) {
    if (anyNA(sm_sens[j, c("2.5%", "50%", "97.5%")])) next
    y_s   <- j + offset
    col_j <- col_vec[n_rows - j + 1]
    lines(x = c(sm_sens[j, "2.5%"], sm_sens[j, "97.5%"]),
          y = c(y_s, y_s), col = col_j, lwd = 1, lty = lty_val)
    points(x = sm_sens[j, "50%"], y = y_s,
           pch = pch_val, bg = "white", col = col_j, cex = 0.75)
  }
}

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png") {
    png("./figures/effect_sizes_omnibus.png",
        width = 12, height = 14, unit = "in", res = 440)
  } else {
    pdf("./figures/effect_sizes_omnibus.pdf",
        width = 12, height = 14)
  }
  
  par(mfrow = c(1, 2))
  
  # Left panel
  MCMCvis::MCMCplot(mod_samp,
                    params = params_to_plot1,
                    labels = c("Intercept","","Depth","","Clay","",
                               "Nitrogen","","Irrigation","",
                               "Precipitation","","Temperature","",
                               "Precip*Temp","",
                               "Intercept","Depth","Clay","Nitrogen",
                               "Precipitation","Temperature","Precip*Temp"),
                    col = col_left,
                    sz_med = 0.75,
                    mar = c(5.1, 6.6, 1.1, 2.1),
                    xlim = c(-2, 2))
  for (spec in omnibus_specs) {
    sm <- sens_summary(spec$samp, params_to_plot1)
    draw_omnibus_overlay(sm, nrow(sm), col_left,
                         spec$offset, spec$pch, spec$lty)
  }
  # Redraw primary on top so it is not obscured by sensitivity overlays
  sm_prim_left <- sens_summary(mod_samp, params_to_plot1)
  draw_primary_on_top(sm_prim_left, nrow(sm_prim_left), col_left)
  segments(y0 = 7.25 + seq(0, 16, by = 2), y1 = 7.25 + seq(0, 16, by = 2),
           x0 = rep(-2, 10), x1 = rep(2, 10), lty = 2)
  text(-4.0, 7.25 / 2, "Maximum", srt = 90, xpd = TRUE, cex = 1.3)
  segments(x0 = -3.8, y0 = 0.5, y1 = 7, xpd = TRUE)
  text(-4.0, 15.5, "Agricultural land cover", srt = 90, xpd = TRUE, cex = 1.3)
  segments(x0 = -3.8, y0 = 8, y1 = 23.5, xpd = TRUE)
  text("Cropland", col = "#d95f02", cex = 1, x = 1.5, y = 22.95, xpd = TRUE)
  text("Pasture",  col = "#7570b3", cex = 1, x = 1.5, y = 22.05, xpd = TRUE)
  
  # Right panel
  MCMCvis::MCMCplot(mod_samp,
                    params = params_to_plot2,
                    labels = c("C\u2192F (A)","C\u2192F (N)","C\u2192O (A)","C\u2192O (N)",
                               "P\u2192F (A)","P\u2192F (N)","P\u2192O (A)","P\u2192O (N)",
                               "Depth","","Clay","","Nitrogen","",
                               "Precipitation","","Temperature","","Precip*Temp","",
                               "C\u2192F (A)","C\u2192F (N)","C\u2192O (A)","C\u2192O (N)",
                               "P\u2192F (A)","P\u2192F (N)","P\u2192O (A)","P\u2192O (N)"),
                    col  = col_right,
                    sz_med = 0.75,
                    mar  = c(5.1, 8, 1.1, 2.1),
                    xlim = c(-5, 10))
  for (spec in omnibus_specs) {
    sm <- sens_summary(spec$samp, params_to_plot2)
    draw_omnibus_overlay(sm, nrow(sm), col_right,
                         spec$offset, spec$pch, spec$lty)
  }
  # Redraw primary on top so it is not obscured by sensitivity overlays
  sm_prim_right <- sens_summary(mod_samp, params_to_plot2)
  draw_primary_on_top(sm_prim_right, nrow(sm_prim_right), col_right)
  abline(h = 0.5 + seq(0, 12, by = 2) + 7.75, lty = 2)
  text(-12.6, 4.25,  "m value",   font = 3, srt = 90, xpd = TRUE, cex = 1.3)
  segments(x0 = -12.0, y0 = 0.25, y1 = 8,     xpd = TRUE)
  text(-10.2, 4.25,  "Intercept",              srt = 90, xpd = TRUE, cex = 1.3)
  segments(x0 =  -9.6, y0 = 0.25, y1 = 8,     xpd = TRUE)
  text(-12.6, 18.5,  "r value",   font = 3, srt = 90, xpd = TRUE, cex = 1.3)
  segments(x0 = -12.0, y0 = 9,    y1 = 28.25, xpd = TRUE)
  text(-10.2, 24.25, "Intercept",              srt = 90, xpd = TRUE, cex = 1.3)
  segments(x0 =  -9.6, y0 = 20.75, y1 = 28.25, xpd = TRUE)
  
  legend(x = -5, y = 28.5, xpd = TRUE,
         legend = c("Primary analysis",
                    sapply(omnibus_specs, "[[", "label")),
         pch    = c(19, sapply(omnibus_specs, "[[", "pch")),
         lty    = c(1,  sapply(omnibus_specs, "[[", "lty")),
         pt.bg  = c(NA, rep("white", length(omnibus_specs))),
         col    = "grey30",
         bty    = "n", cex = 0.95, pt.cex = c(1.3, rep(1.0, length(omnibus_specs))))
  
  dev.off()
  message("Saved ./figures/effect_sizes_omnibus.", dev_type)
}
rm(sm_bd_left, sm_bd_right, sm_origbd_left, sm_origbd_right)



# =============================================================================
# SECTION 8 — PRIMARY ANALYSIS FIGURES
# =============================================================================
# The remaining figures use mod_samp from the primary analysis and are
# unchanged from the original 03_bayes_analysis.R.  Sensitivity comparison
# figures were written above (Section 7).

# Figures ----------------------------------------------------------------------
# These are figures showing the model results and various ways to show the
# effect sizes and relationships.

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png")
    png("./figures/effect sizes.png", width = 12, height = 8, unit = "in", res = 440)
  else
    pdf("./figures/effect sizes.pdf", width = 12, height = 8)
  
  
  par(mfrow=c(1,2))
  # params_to_plot1 and params_to_plot2 are defined in Section 7 and reused here
  
  MCMCvis::MCMCplot(mod_samp,
                    params = params_to_plot1,
                    labels = c("Intercept", "",
                               "Depth","",
                               "Clay","",
                               "Nitrogen","",
                               "Irrigation","",
                               "Precipitation","",
                               "Temperature","",
                               "Precip*Temp","",
                               "Intercept",
                               "Depth",
                               "Clay",
                               "Nitrogen",
                               "Precipitation",
                               "Temperature",
                               "Precip*Temp"),
                    col = c(rep(c("#d95f02","#7570b3"), 
                                times = 8), 
                            rep("black", 7)),
                    mar = c(5.1,6.6,1.1,2.1),
                    xlim = c(-2, 2))
  segments(y0 = 7.5 + seq(0, 16, by = 2),
           y1 = 7.5 + seq(0, 16, by = 2),
           x0 = rep(-2, 10),
           x1 = rep(2, 10),
           lty = 2)
  text(-4.0, 7.5/2, "Maximum", srt = 90, xpd = T, cex = 1.3)
  segments(x0 = -3.8,
           y0 = 0.5,
           y1 = 7,
           xpd = T)
  
  text(-4.0, 15.5, "Agricultural land cover", srt = 90, xpd = T, cex = 1.3)
  segments(x0 = -3.8,
           y0 = 8,
           y1 = 23.5, 
           xpd = T)
  
  text("Cropland",
       col = "#d95f02",
       cex = 1.1,
       x = 1.5,
       y = 22.95,
       xpd = T)
  text("Pasture",
       col = "#7570b3",
       cex = 1.1,
       x = 1.5,
       y = 22.05,
       xpd = T)
  
  
  MCMCvis::MCMCplot(mod_samp,
                    params = c("beta_r_intercept"
                               ,"beta_r_depth"
                               ,"beta_r_clay"
                               ,"beta_r_nitrogen"
                               ,"beta_r_map"
                               ,"beta_r_mat"
                               ,"beta_r_map_mat"
                               
                               ,"m_val"),
                    labels = c("C→F (A)",
                               "C→F (P)",
                               "C→O (A)",
                               "C→O (P)",
                               "P→F (A)",
                               "P→F (P)",
                               "P→O (A)",
                               "P→O (P)",
                               "Depth","",
                               "Clay","",
                               "Nitrogen","",
                               "Precipitation","",
                               "Temperature","",
                               "Precip*Temp","",
                               "C→F (A)",
                               "C→F (P)",
                               "C→O (A)",
                               "C→O (P)",
                               "P→F (A)",
                               "P→F (P)",
                               "P→O (A)",
                               "P→O (P)"),
                    col = c(rep(c("#d95f02","#7570b3"), 
                                each = 4),
                            rep(c("#d95f02","#7570b3"), 
                                times = 6), 
                            rep(c("#d95f02","#7570b3"), 
                                each = 4)
                    ),
                    mar = c(5.1,8,1.1,2.1),
                    xlim = c(-5, 10))
  abline(h = 0.5 + seq(0, 12, by = 2) + 8, lty = 2)
  
  text(-12.6, 4.25, "m value", font = 3, srt = 90, xpd = T, cex = 1.3)
  segments(x0 = -12.0,
           y0 = 0.25,
           y1 = 8, 
           xpd = T)
  
  text(-10.2, 4.25, "Intercept", srt = 90, xpd = T, cex = 1.3)
  segments(x0 = -9.6,
           y0 = 0.25,
           y1 = 8, 
           xpd = T)
  
  text(-12.6, 18.5, "r value", font = 3, srt = 90, xpd = T, cex = 1.3)
  segments(x0 = -12.0,
           y0 = 9,
           y1 = 28.25, 
           xpd = T)
  
  text(-10.2, 24.25, "Intercept", srt = 90, xpd = T, cex = 1.3)
  segments(x0 = -9.6,
           y0 = 20.75,
           y1 = 28.25, 
           xpd = T)
  
  dev.off()
  message("Saved ./figures/effect sizes.", dev_type)
}







# Prior-posterior overlap for slope coefficients
# Uses samples already in mod_samp — no rerunning needed

ipak("MCMCvis")  # loaded here for Section 8 marginal figures

# Parameters to show — slopes only, no intercepts
slope_params <- c(
  "beta_init_depth",   "beta_init_clay",  "beta_init_nitrogen",
  "beta_init_map",     "beta_init_mat",   "beta_init_map_mat",
  "beta_max_depth",    "beta_max_clay",   "beta_max_nitrogen",
  "beta_max_map",      "beta_max_mat",    "beta_max_map_mat"
)

# Extract posterior samples as a matrix
post_mat <- MCMCvis::MCMCchains(mod_samp, params = slope_params)

# Build a lookup of prior mean and SD for each parameter row name
# param_root -> list of [1] and [2] indices map to cropland/pasture
prior_lookup <- list(
  "beta_init_depth[1]"    = c(jags_data$prior_cropland_depth,    1/sqrt(jags_data$prior_tau_cropland_depth)),
  "beta_init_depth[2]"    = c(jags_data$prior_pasture_depth,     1/sqrt(jags_data$prior_tau_pasture_depth)),
  "beta_init_clay[1]"     = c(jags_data$prior_cropland_clay,     1/sqrt(jags_data$prior_tau_cropland_clay)),
  "beta_init_clay[2]"     = c(jags_data$prior_pasture_clay,      1/sqrt(jags_data$prior_tau_pasture_clay)),
  "beta_init_nitrogen[1]" = c(jags_data$prior_cropland_nitrogen, 1/sqrt(jags_data$prior_tau_cropland_nitrogen)),
  "beta_init_nitrogen[2]" = c(jags_data$prior_pasture_nitrogen,  1/sqrt(jags_data$prior_tau_pasture_nitrogen)),
  "beta_init_map[1]"      = c(jags_data$prior_cropland_map,      1/sqrt(jags_data$prior_tau_cropland_map)),
  "beta_init_map[2]"      = c(jags_data$prior_pasture_map,       1/sqrt(jags_data$prior_tau_pasture_map)),
  "beta_init_mat[1]"      = c(jags_data$prior_cropland_mat,      1/sqrt(jags_data$prior_tau_cropland_mat)),
  "beta_init_mat[2]"      = c(jags_data$prior_pasture_mat,       1/sqrt(jags_data$prior_tau_pasture_mat)),
  "beta_init_map_mat[1]"  = c(jags_data$prior_cropland_map_mat,  1/sqrt(jags_data$prior_tau_cropland_map_mat)),
  "beta_init_map_mat[2]"  = c(jags_data$prior_pasture_map_mat,   1/sqrt(jags_data$prior_tau_pasture_map_mat)),
  "beta_max_depth"        = c(jags_data$prior_max_depth,         1/sqrt(jags_data$prior_tau_max_depth)),
  "beta_max_clay"         = c(jags_data$prior_max_clay,          1/sqrt(jags_data$prior_tau_max_clay)),
  "beta_max_nitrogen"     = c(jags_data$prior_max_nitrogen,      1/sqrt(jags_data$prior_tau_max_nitrogen)),
  "beta_max_map"          = c(jags_data$prior_max_map,           1/sqrt(jags_data$prior_tau_max_map)),
  "beta_max_mat"          = c(jags_data$prior_max_mat,           1/sqrt(jags_data$prior_tau_max_mat)),
  "beta_max_map_mat"      = c(jags_data$prior_max_map_mat,       1/sqrt(jags_data$prior_tau_max_map_mat))
)

# Nice labels
param_labels <- c(
  "beta_init_depth[1]"    = "Depth — Cropland",
  "beta_init_depth[2]"    = "Depth — Pasture",
  "beta_init_clay[1]"     = "Clay — Cropland",
  "beta_init_clay[2]"     = "Clay — Pasture",
  "beta_init_nitrogen[1]" = "Nitrogen — Cropland",
  "beta_init_nitrogen[2]" = "Nitrogen — Pasture",
  "beta_init_map[1]"      = "Precip — Cropland",
  "beta_init_map[2]"      = "Precip — Pasture",
  "beta_init_mat[1]"      = "Temp — Cropland",
  "beta_init_mat[2]"      = "Temp — Pasture",
  "beta_init_map_mat[1]"  = "Precip×Temp — Cropland",
  "beta_init_map_mat[2]"  = "Precip×Temp — Pasture",
  "beta_max_depth"        = "Depth — Max SOC",
  "beta_max_clay"         = "Clay — Max SOC",
  "beta_max_nitrogen"     = "Nitrogen — Max SOC",
  "beta_max_map"          = "Precip — Max SOC",
  "beta_max_mat"          = "Temp — Max SOC",
  "beta_max_map_mat"      = "Precip×Temp — Max SOC"
)

# Section labels for facet strips
section <- ifelse(grepl("beta_max", names(param_labels)), "Maximum SOC", "Agricultural land cover")

# Build plotting data frame
plot_list <- lapply(names(prior_lookup), function(pname) {
  post_samples <- post_mat[, pname]
  pr           <- prior_lookup[[pname]]
  
  # x range: cover both prior and posterior
  x_lo <- min(pr[1] - 4*pr[2], quantile(post_samples, 0.001))
  x_hi <- max(pr[1] + 4*pr[2], quantile(post_samples, 0.999))
  x    <- seq(x_lo, x_hi, length.out = 500)
  
  rbind(
    data.frame(param = pname, x = x,
               y     = dnorm(x, pr[1], pr[2]),
               type  = "Prior",
               label = param_labels[pname],
               section = section[which(names(param_labels) == pname)]),
    data.frame(param = pname, x = density(post_samples)$x,
               y     = density(post_samples)$y,
               type  = "Posterior",
               label = param_labels[pname],
               section = section[which(names(param_labels) == pname)])
  )
})
plot_df <- do.call(rbind, plot_list)
plot_df$label   <- factor(plot_df$label,   levels = param_labels)
plot_df$section <- factor(plot_df$section, levels = c("Agricultural land cover", "Maximum SOC"))

p_overlap <- ggplot(plot_df, aes(x = x, y = y, color = type, linetype = type)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ label, scales = "free", ncol = 3) +
  scale_color_manual(values = c("Prior" = "grey60", "Posterior" = "#1b7837"),
                     name = NULL) +
  scale_linetype_manual(values = c("Prior" = "dashed", "Posterior" = "solid"),
                        name = NULL) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40", linewidth = 0.4) +
  labs(x = "Parameter value", y = "Density",
       title = "Prior vs. posterior: slope coefficients") +
  theme_minimal(base_size = 9) +
  theme(
    strip.text      = element_text(size = 7.5),
    legend.position = "bottom",
    panel.grid      = element_blank(),
    axis.text       = element_text(size = 6)
  )

png("./figures/prior_posterior_slopes.png",
    width = 7, height = 10, unit = "in", res = 400)
print(p_overlap)
dev.off()

pdf("./figures/prior_posterior_slopes.pdf",
    width = 7, height = 10)
print(p_overlap)
dev.off()














# Other ways to look at effect sizes
age_max <- 100
ylims <- c(-10, 300)

# Get posterior means and credible interval vals
mod_samp_df <- do.call(rbind, mod_samp)
bsm <- colMeans(mod_samp_df)
bs975 <-apply(mod_samp_df, 2, function(x) quantile(x, 0.975))
bs025 <-apply(mod_samp_df, 2, function(x) quantile(x, 0.025))
mod_res_df <- tibble(coef = names(bsm),
                     mean = bsm,
                     lo = bs025,
                     hi = bs975)

# Make a more streamlined posterior df to make some of the predictions easier
ms_df <- mod_samp_df %>% 
  as.data.frame() %>% 
  dplyr::select(-starts_with(c(
    "c_init", "c_final",
    "beta_mean",
    "rand_coord")))

# This helps deal with scaled/centered values of each variable
sc_fun <- function(var,
                   natural_val,
                   sd_val){
  if(!is.null(natural_val)){
    natural_val <- ifelse(sc_df[var, "log"],
                          log(natural_val),
                          natural_val)
    out_val <- (natural_val - sc_df[var, "center"]) / sc_df[var, "scale"]
    
  } else{
    out_val <- sd_val
  }
  return(out_val)
}

# Age Transformation function
age_fun <- function(x){
  x # No transformation, but could be modified
}



# Want to get values for multiple depths
depth_vals <- c(5,30,100)
depth_ints <- c(0,15,60,200)


lu_avg <- list()

for(i in 1:length(depth_vals)){
  lu_avg[[i]] <- final_covars %>% 
    filter(depth_avg >= depth_ints[i],
           depth_avg < depth_ints[i+1]) %>% 
    group_by(lu_end) %>% 
    summarise(clay = mean(clay, na.rm = TRUE),
              #sand = mean(sand, na.rm = TRUE),
              nitrogen = mean(nitrogen, na.rm = TRUE),
              #phh2o = mean(phh2o, na.rm = TRUE),
              map = exp(mean(log(MAP), na.rm = TRUE)), # These are more log-normally distributed so taking the mean of log-transformed values
              mat = mean(MAT, na.rm = TRUE),
              irrigation = exp(mean(log(irrigation), na.rm = TRUE)))
}

# Uncomment if the user wants to visualize more consistent env. covariates
# lu_avg <- list()
# 
# lu_avg[[1]] <- tibble(lu_end = c("F", "O"),
#                       clay = c(250, 250),
#                       sand = c(400, 400),
#                       nitrogen = c(4200, 3900),
#                       phh20 = c(60, 60),
#                       map = c(1000, 600),
#                       mat = c(13, 10),
#                       irrigation = c(0.7, 1.6))
# 
# lu_avg[[2]] <- tibble(lu_end = c("F", "O"),
#                       clay = c(250, 250),
#                       sand = c(400, 400),
#                       nitrogen = c(2300, 1700),
#                       phh20 = c(60, 60),
#                       map = c(1000, 600),
#                       mat = c(13, 10),
#                       irrigation = c(0.7, 1.6))
# 
# lu_avg[[3]] <- tibble(lu_end = c("F", "O"),
#                       clay = c(250, 250),
#                       sand = c(400, 400),
#                       nitrogen = c(1100, 800),
#                       phh20 = c(60, 60),
#                       map = c(1000, 600),
#                       mat = c(13, 10),
#                       irrigation = c(0.7, 1.6))




# Decide on what appears on the X axis of these plots
n_x <- 100
age_x <- seq(1,#sc_df["age", "min"],
             age_max,#100,#sc_df["age", "max"],
             length.out = n_x)

# init_x <- seq(sc_df["init", "min"],
#               sc_df["init", "max"],
#               length.out = n_x)

# For the initial carbon, let's make it relative between 50% of average to 100% above
init_x <- seq(0.5, 1.5, length.out = n_x)

depth_x <- seq(sc_df["depth", "min"],
               sc_df["depth", "max"],
               length.out = n_x)

# I think this is the most updated version of all this

# Age with different depth values ----------------------------------------------

age_out <- tibble(age = age_x) %>% 
  as.data.frame()



for(depth in depth_vals){
  for(age in age_x){
    for(lu in 1:4){
      
      if(depth == depth_vals[1]){
        col_name <- paste0("depth_neg_", lu)
      }
      if(depth == depth_vals[2]){
        col_name <- paste0("depth_avg_", lu)
      }
      if(depth == depth_vals[3]){
        col_name <- paste0("depth_pos_", lu)
      }
      
      # Helper function to find the correct column. Updates each iteration
      # for the appropriate land use
      col_fun <- function(x){
        paste0(x, "[", lu, "]")
      }
      
      init_vals <- ms_df[,"beta_init_intercept[1]"] + 
        ms_df[,"beta_init_depth[1]"] * sc_fun(var = "depth", natural_val = depth) + 
        
        ms_df[,"beta_init_clay[1]"] * sc_fun(var = "clay", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "clay"])) + 
        
        ms_df[,"beta_init_nitrogen[1]"] * sc_fun(var = "nitrogen", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "nitrogen"])) + 
        
        ms_df[,"beta_init_irrigation[1]"] * sc_fun(var = "irrigation", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "irrigation"])) + 
        
        ms_df[,"beta_init_map_mat[1]"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "map"])) * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "mat"])) +
        
        ms_df[,"beta_init_map[1]"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "map"])) + 
        
        ms_df[,"beta_init_mat[1]"] * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "mat"]))
      
      
      
      r_vals <- ms_df[,col_fun("beta_r_intercept")] +
        
        ms_df[,"beta_r_depth[1]"] * sc_fun(var = "depth", natural_val = depth) +
        
        ms_df[,"beta_r_clay[1]"] * sc_fun(var = "clay", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "clay"])) +
        
        ms_df[,"beta_r_nitrogen[1]"] * sc_fun(var = "nitrogen", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "nitrogen"])) +
        
        ms_df[,"beta_r_map_mat[1]"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "map"])) * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "mat"])) +
        
        ms_df[,"beta_r_map[1]"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "map"])) +
        
        ms_df[,"beta_r_mat[1]"] * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "mat"]))
      
      max_vals <- ms_df[,"beta_max_intercept"] +
        ms_df[,"beta_max_depth"] * sc_fun(var = "depth", natural_val = depth) + 
        
        ms_df[,"beta_max_clay"] * sc_fun(var = "clay", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "clay"])) + 
        
        ms_df[,"beta_max_nitrogen"] * sc_fun(var = "nitrogen", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "nitrogen"])) + 
        
        ms_df[,"beta_max_map_mat"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "map"])) * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "mat"])) +
        
        ms_df[,"beta_max_map"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "map"])) + 
        
        ms_df[,"beta_max_mat"] * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "mat"]))
      
      final_vals <- log(exp(init_vals) + (exp(max_vals) - exp(init_vals)) * (1 - exp(-age_fun(age) / exp(r_vals))) ^ exp(ms_df[,col_fun("m_val")]))
      
      perc_change_vals <- ((exp(final_vals) - exp(init_vals)) / exp(init_vals)) * 100
      
      ind <- which(age_out$age == age)
      age_out[ind, c(paste(col_name, "lo", sep = "_"),
                     paste(col_name, "med", sep = "_"),
                     paste(col_name, "hi", sep = "_"))] <- unlist(quantile(perc_change_vals,
                                                                           probs = c(0.025, 0.5, 0.975)))
      
    }
  }
}

for (dev_type in c("png", "pdf")) {
  if (dev_type == "png")
    png(filename = "./figures/age with varied depths.png", res = 440, unit = "in", width = 8, height = 5.5)
  else
    pdf(file = "./figures/age with varied depths.pdf", width = 8, height = 5.5)
  par(mfrow=c(2, 4))
  top_row_mar <- c(3.1, 4.1, 3.6, 0.3)
  bottom_row_mar <- c(5.1, 4.1, 1.6, 0.3)
  par(mar = top_row_mar)
  
  plot(age_out$depth_neg_1_med ~ age_out$age, type = "l",
       ylim = ylims,#c(-50, 250),
       frame = F,
       las = 1,
       ylab = "Percent change in C density",
       xlab = "")
  mtext(side = 3, line = -0.5, "Cropland → Forest\nActive", cex = 0.7)
  polygon(x = c(age_out$age, rev(age_out$age)),
          y = c(age_out$depth_neg_1_lo,
                rev(age_out$depth_neg_1_hi)),
          col = rgb(1, 0.647, 0, alpha = 0.2),
          border = NA)
  
  lines(age_out$depth_avg_1_med ~ age_out$age)
  polygon(x = c(age_out$age, rev(age_out$age)),
          y = c(age_out$depth_avg_1_lo,
                rev(age_out$depth_avg_1_hi)),
          col = rgb(0,0,0,0.2),
          border = NA)
  
  lines(age_out$depth_pos_1_med ~ age_out$age)
  polygon(x = c(age_out$age, rev(age_out$age)),
          y = c(age_out$depth_pos_1_lo,
                rev(age_out$depth_pos_1_hi)),
          col = rgb(0.5, 0, 0.5, alpha = 0.2),
          border = NA)
  abline(h = 0, lty = 2)
  
  legend("topleft",
         col = c(rgb(1, 0.647, 0, alpha = 0.4),
                 rgb(0,0,0, alpha = 0.4),
                 rgb(0.5, 0, 0.5, alpha = 0.4)),
         pt.cex = 2,
         bty = "n",
         pch = 15,
         legend = paste(c(depth_vals), "cm"))
  
  plot(age_out$depth_neg_2_med ~ age_out$age, type = "l",
       ylim = ylims,#c(-50, 250),
       frame = F,
       las = 1,
       ylab = "",
       xlab = "")
  mtext(side = 3, line = -0.5, "Cropland → Forest\nPassive", cex = 0.7)
  polygon(x = c(age_out$age, rev(age_out$age)),
          y = c(age_out$depth_neg_2_lo,
                rev(age_out$depth_neg_2_hi)),
          col = rgb(1, 0.647, 0, alpha = 0.2),
          border = NA)
  
  lines(age_out$depth_avg_2_med ~ age_out$age)
  polygon(x = c(age_out$age, rev(age_out$age)),
          y = c(age_out$depth_avg_2_lo,
                rev(age_out$depth_avg_2_hi)),
          col = rgb(0,0,0,0.2),
          border = NA)
  
  lines(age_out$depth_pos_2_med ~ age_out$age)
  polygon(x = c(age_out$age, rev(age_out$age)),
          y = c(age_out$depth_pos_2_lo,
                rev(age_out$depth_pos_2_hi)),
          col = rgb(0.5, 0, 0.5, alpha = 0.2),
          border = NA)
  abline(h = 0, lty = 2)
  
  
  
  plot(age_out$depth_neg_3_med ~ age_out$age, type = "l",
       ylim = ylims,#c(-50, 250),
       frame = F,
       las = 1,
       ylab = "",
       xlab = "")
  mtext(side = 3, line = -0.5, "Cropland → Open\nActive", cex = 0.7)
  polygon(x = c(age_out$age, rev(age_out$age)),
          y = c(age_out$depth_neg_3_lo,
                rev(age_out$depth_neg_3_hi)),
          col = rgb(1, 0.647, 0, alpha = 0.2),
          border = NA)
  
  lines(age_out$depth_avg_3_med ~ age_out$age)
  polygon(x = c(age_out$age, rev(age_out$age)),
          y = c(age_out$depth_avg_3_lo,
                rev(age_out$depth_avg_3_hi)),
          col = rgb(0,0,0,0.2),
          border = NA)
  
  lines(age_out$depth_pos_3_med ~ age_out$age)
  polygon(x = c(age_out$age, rev(age_out$age)),
          y = c(age_out$depth_pos_3_lo,
                rev(age_out$depth_pos_3_hi)),
          col = rgb(0.5, 0, 0.5, alpha = 0.2),
          border = NA)
  abline(h = 0, lty = 2)
  
  
  
  plot(age_out$depth_neg_4_med ~ age_out$age, type = "l",
       ylim = ylims,#c(-50, 250),
       frame = F,
       las = 1,
       ylab = "",
       xlab = "")
  mtext(side = 3, line = -0.5, "Cropland → Open\nPassive", cex = 0.7)
  polygon(x = c(age_out$age, rev(age_out$age)),
          y = c(age_out$depth_neg_4_lo,
                rev(age_out$depth_neg_4_hi)),
          col = rgb(1, 0.647, 0, alpha = 0.2),
          border = NA)
  
  lines(age_out$depth_avg_4_med ~ age_out$age)
  polygon(x = c(age_out$age, rev(age_out$age)),
          y = c(age_out$depth_avg_4_lo,
                rev(age_out$depth_avg_4_hi)),
          col = rgb(0,0,0,0.2),
          border = NA)
  
  lines(age_out$depth_pos_4_med ~ age_out$age)
  polygon(x = c(age_out$age, rev(age_out$age)),
          y = c(age_out$depth_pos_4_lo,
                rev(age_out$depth_pos_4_hi)),
          col = rgb(0.5, 0, 0.5, alpha = 0.2),
          border = NA)
  abline(h = 0, lty = 2)
  
  
  
  
  # Now types 5-8
  
  for(depth in depth_vals){
    for(age in age_x){
      for(lu in 5:8){
        
        if(depth == depth_vals[1]){
          col_name <- paste0("depth_neg_", lu)
        }
        if(depth == depth_vals[2]){
          col_name <- paste0("depth_avg_", lu)
        }
        if(depth == depth_vals[3]){
          col_name <- paste0("depth_pos_", lu)
        }
        
        
        # Helper function to find the correct column. Updates each iteration
        # for the appropriate land use
        col_fun <- function(x){
          paste0(x, "[", lu, "]")
        }
        
        init_vals <- ms_df[,"beta_init_intercept[2]"] + 
          ms_df[,"beta_init_depth[2]"] * sc_fun(var = "depth", natural_val = depth) + 
          
          ms_df[,"beta_init_clay[2]"] * sc_fun(var = "clay", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "clay"])) + 
          
          ms_df[,"beta_init_nitrogen[2]"] * sc_fun(var = "nitrogen", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "nitrogen"])) +
          
          ms_df[,"beta_init_irrigation[2]"] * sc_fun(var = "irrigation", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "irrigation"])) + 
          
          ms_df[,"beta_init_map_mat[2]"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "map"])) * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "mat"])) +
          
          ms_df[,"beta_init_map[2]"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "map"])) + 
          
          ms_df[,"beta_init_mat[2]"] * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "mat"]))
        
        
        
        r_vals <- ms_df[,col_fun("beta_r_intercept")] +
          
          ms_df[,"beta_r_depth[2]"] * sc_fun(var = "depth", natural_val = depth) +
          
          ms_df[,"beta_r_clay[2]"] * sc_fun(var = "clay", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "clay"])) +
          
          ms_df[,"beta_r_nitrogen[2]"] * sc_fun(var = "nitrogen", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "nitrogen"])) +
          
          ms_df[,"beta_r_map_mat[2]"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "map"])) * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "mat"])) +
          
          ms_df[,"beta_r_map[2]"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "map"])) +
          
          ms_df[,"beta_r_mat[2]"] * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "mat"]))
        
        
        max_vals <- ms_df[,"beta_max_intercept"] +
          ms_df[,"beta_max_depth"] * sc_fun(var = "depth", natural_val = depth) + 
          
          ms_df[,"beta_max_clay"] * sc_fun(var = "clay", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "clay"])) + 
          
          ms_df[,"beta_max_nitrogen"] * sc_fun(var = "nitrogen", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "nitrogen"])) + 
          
          ms_df[,"beta_max_map_mat"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "map"])) * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "mat"])) +
          
          ms_df[,"beta_max_map"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "map"])) + 
          
          ms_df[,"beta_max_mat"] * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "mat"]))
        
        final_vals <- log(exp(init_vals) + (exp(max_vals) - exp(init_vals)) * (1 - exp(-age_fun(age) / exp(r_vals))) ^ exp(ms_df[,col_fun("m_val")]))
        
        perc_change_vals <- ((exp(final_vals) - exp(init_vals)) / exp(init_vals)) * 100
        
        ind <- which(age_out$age == age)
        age_out[ind, c(paste(col_name, "lo", sep = "_"),
                       paste(col_name, "med", sep = "_"),
                       paste(col_name, "hi", sep = "_"))] <- unlist(quantile(perc_change_vals,
                                                                             probs = c(0.025, 0.5, 0.975)))
        
      }
    }
  }
  
  
  par(mar = bottom_row_mar)
  
  plot(age_out$depth_neg_5_med ~ age_out$age, type = "l",
       ylim = ylims,#c(-50, 250),
       frame = F,
       las = 1,
       ylab = "Percent change in C density",
       xlab = "Years since restoration")
  mtext(side = 3, line = -0.5, "Pasture → Forest\nActive", cex = 0.7)
  polygon(x = c(age_out$age, rev(age_out$age)),
          y = c(age_out$depth_neg_5_lo,
                rev(age_out$depth_neg_5_hi)),
          col = rgb(1, 0.647, 0, alpha = 0.2),
          border = NA)
  
  lines(age_out$depth_avg_5_med ~ age_out$age)
  polygon(x = c(age_out$age, rev(age_out$age)),
          y = c(age_out$depth_avg_5_lo,
                rev(age_out$depth_avg_5_hi)),
          col = rgb(0,0,0,0.2),
          border = NA)
  
  lines(age_out$depth_pos_5_med ~ age_out$age)
  polygon(x = c(age_out$age, rev(age_out$age)),
          y = c(age_out$depth_pos_5_lo,
                rev(age_out$depth_pos_5_hi)),
          col = rgb(0.5, 0, 0.5, alpha = 0.2),
          border = NA)
  abline(h = 0, lty = 2)
  
  legend("topleft",
         col = c(rgb(1, 0.647, 0, alpha = 0.4),
                 rgb(0,0,0, alpha = 0.4),
                 rgb(0.5, 0, 0.5, alpha = 0.4)),
         pt.cex = 2,
         bty = "n",
         pch = 15,
         legend = paste(c(depth_vals), "cm"))
  
  plot(age_out$depth_neg_6_med ~ age_out$age, type = "l",
       ylim = ylims,#c(-50, 250),
       frame = F,
       las = 1,
       ylab = "",
       xlab = "Years since restoration")
  mtext(side = 3, line = -0.5, "Pasture → Forest\nPassive", cex = 0.7)
  polygon(x = c(age_out$age, rev(age_out$age)),
          y = c(age_out$depth_neg_6_lo,
                rev(age_out$depth_neg_6_hi)),
          col = rgb(1, 0.647, 0, alpha = 0.2),
          border = NA)
  
  lines(age_out$depth_avg_6_med ~ age_out$age)
  polygon(x = c(age_out$age, rev(age_out$age)),
          y = c(age_out$depth_avg_6_lo,
                rev(age_out$depth_avg_6_hi)),
          col = rgb(0,0,0,0.2),
          border = NA)
  
  lines(age_out$depth_pos_6_med ~ age_out$age)
  polygon(x = c(age_out$age, rev(age_out$age)),
          y = c(age_out$depth_pos_6_lo,
                rev(age_out$depth_pos_6_hi)),
          col = rgb(0.5, 0, 0.5, alpha = 0.2),
          border = NA)
  abline(h = 0, lty = 2)
  
  
  
  plot(age_out$depth_neg_7_med ~ age_out$age, type = "l",
       ylim = ylims,#c(-50, 250),
       frame = F,
       las = 1,
       ylab = "",
       xlab = "Years since restoration")
  mtext(side = 3, line = -0.5, "Pasture → Open\nActive", cex = 0.7)
  
  polygon(x = c(age_out$age, rev(age_out$age)),
          y = c(age_out$depth_neg_7_lo,
                rev(age_out$depth_neg_7_hi)),
          col = rgb(1, 0.647, 0, alpha = 0.2),
          border = NA)
  
  lines(age_out$depth_avg_7_med ~ age_out$age)
  polygon(x = c(age_out$age, rev(age_out$age)),
          y = c(age_out$depth_avg_7_lo,
                rev(age_out$depth_avg_7_hi)),
          col = rgb(0,0,0,0.2),
          border = NA)
  
  lines(age_out$depth_pos_7_med ~ age_out$age)
  polygon(x = c(age_out$age, rev(age_out$age)),
          y = c(age_out$depth_pos_7_lo,
                rev(age_out$depth_pos_7_hi)),
          col = rgb(0.5, 0, 0.5, alpha = 0.2),
          border = NA)
  abline(h = 0, lty = 2)
  
  
  
  plot(age_out$depth_neg_8_med ~ age_out$age, type = "l",
       ylim = ylims,#c(-50, 250),
       frame = F,
       las = 1,
       ylab = "",
       xlab = "Years since restoration")
  mtext(side = 3, line = -0.5, "Pasture → Open\nPassive", cex = 0.7)
  
  polygon(x = c(age_out$age, rev(age_out$age)),
          y = c(age_out$depth_neg_8_lo,
                rev(age_out$depth_neg_8_hi)),
          col = rgb(1, 0.647, 0, alpha = 0.2),
          border = NA)
  
  lines(age_out$depth_avg_8_med ~ age_out$age)
  polygon(x = c(age_out$age, rev(age_out$age)),
          y = c(age_out$depth_avg_8_lo,
                rev(age_out$depth_avg_8_hi)),
          col = rgb(0,0,0,0.2),
          border = NA)
  
  lines(age_out$depth_pos_8_med ~ age_out$age)
  polygon(x = c(age_out$age, rev(age_out$age)),
          y = c(age_out$depth_pos_8_lo,
                rev(age_out$depth_pos_8_hi)),
          col = rgb(0.5, 0, 0.5, alpha = 0.2),
          border = NA)
  abline(h = 0, lty = 2)
  
  dev.off()
  message("Saved ./figures/age with varied depths.", dev_type)
}




















# Version for a single depth 
# This is a little clunky/overcomplicated but want to plot at 15 cm

depth_vals <- c(5,15,100)

lu_avg <- list()

lu_avg[[1]] <- tibble(lu_end = c("F", "O"),
                      clay = c(250, 250),
                      sand = c(400, 400),
                      nitrogen = c(4200, 3900),
                      phh20 = c(60, 60),
                      map = c(1000, 600),
                      mat = c(13, 10),
                      irrigation = c(0.7, 1.6))

lu_avg[[2]] <- tibble(lu_end = c("F", "O"),
                      clay = c(250, 250),
                      sand = c(400, 400),
                      nitrogen = c(2300, 1700),
                      phh20 = c(60, 60),
                      map = c(1000, 600),
                      mat = c(13, 10),
                      irrigation = c(0.7, 1.6))

lu_avg[[3]] <- tibble(lu_end = c("F", "O"),
                      clay = c(250, 250),
                      sand = c(400, 400),
                      nitrogen = c(1100, 800),
                      phh20 = c(60, 60),
                      map = c(1000, 600),
                      mat = c(13, 10),
                      irrigation = c(0.7, 1.6))




# Decide on what appears on the X axis of these plots
n_x <- 100
age_x <- seq(1,#sc_df["age", "min"],
             age_max,#100,#sc_df["age", "max"],
             length.out = n_x)

# init_x <- seq(sc_df["init", "min"],
#               sc_df["init", "max"],
#               length.out = n_x)

# For the initial carbon, let's make it relative between 50% of average to 100% above
init_x <- seq(0.5, 1.5, length.out = n_x)

depth_x <- seq(sc_df["depth", "min"],
               sc_df["depth", "max"],
               length.out = n_x)

# I think this is the most updated version of all this

# Age with different depth values ----------------------------------------------

age_out <- tibble(age = age_x) %>% 
  as.data.frame()



for(depth in depth_vals){
  for(age in age_x){
    for(lu in 1:4){
      
      if(depth == depth_vals[1]){
        col_name <- paste0("depth_neg_", lu)
      }
      if(depth == depth_vals[2]){
        col_name <- paste0("depth_avg_", lu)
      }
      if(depth == depth_vals[3]){
        col_name <- paste0("depth_pos_", lu)
      }
      
      # Helper function to find the correct column. Updates each iteration
      # for the appropriate land use
      col_fun <- function(x){
        paste0(x, "[", lu, "]")
      }
      
      init_vals <- ms_df[,"beta_init_intercept[1]"] + 
        ms_df[,"beta_init_depth[1]"] * sc_fun(var = "depth", natural_val = depth) + 
        
        ms_df[,"beta_init_clay[1]"] * sc_fun(var = "clay", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "clay"])) + 
        
        ms_df[,"beta_init_nitrogen[1]"] * sc_fun(var = "nitrogen", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "nitrogen"])) + 
        
        ms_df[,"beta_init_irrigation[1]"] * sc_fun(var = "irrigation", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "irrigation"])) + 
        
        ms_df[,"beta_init_map_mat[1]"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "map"])) * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "mat"])) +
        
        ms_df[,"beta_init_map[1]"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "map"])) + 
        
        ms_df[,"beta_init_mat[1]"] * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "mat"]))
      
      
      
      r_vals <- ms_df[,col_fun("beta_r_intercept")] +
        
        ms_df[,"beta_r_depth[1]"] * sc_fun(var = "depth", natural_val = depth) +
        
        ms_df[,"beta_r_clay[1]"] * sc_fun(var = "clay", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "clay"])) +
        
        ms_df[,"beta_r_nitrogen[1]"] * sc_fun(var = "nitrogen", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "nitrogen"])) +
        
        ms_df[,"beta_r_map_mat[1]"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "map"])) * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "mat"])) +
        
        ms_df[,"beta_r_map[1]"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "map"])) +
        
        ms_df[,"beta_r_mat[1]"] * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "mat"]))
      
      max_vals <- ms_df[,"beta_max_intercept"] +
        ms_df[,"beta_max_depth"] * sc_fun(var = "depth", natural_val = depth) + 
        
        ms_df[,"beta_max_clay"] * sc_fun(var = "clay", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "clay"])) + 
        
        ms_df[,"beta_max_nitrogen"] * sc_fun(var = "nitrogen", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "nitrogen"])) + 
        
        ms_df[,"beta_max_map_mat"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "map"])) * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "mat"])) +
        
        ms_df[,"beta_max_map"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "map"])) + 
        
        ms_df[,"beta_max_mat"] * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "mat"]))
      
      final_vals <- log(exp(init_vals) + (exp(max_vals) - exp(init_vals)) * (1 - exp(-age_fun(age) / exp(r_vals))) ^ exp(ms_df[,col_fun("m_val")]))
      
      perc_change_vals <- ((exp(final_vals) - exp(init_vals)) / exp(init_vals)) * 100
      
      ind <- which(age_out$age == age)
      age_out[ind, c(paste(col_name, "lo", sep = "_"),
                     paste(col_name, "med", sep = "_"),
                     paste(col_name, "hi", sep = "_"))] <- unlist(quantile(perc_change_vals,
                                                                           probs = c(0.025, 0.5, 0.975)))
      
    }
  }
}


# Now types 5-8

for(depth in depth_vals){
  for(age in age_x){
    for(lu in 5:8){
      
      if(depth == depth_vals[1]){
        col_name <- paste0("depth_neg_", lu)
      }
      if(depth == depth_vals[2]){
        col_name <- paste0("depth_avg_", lu)
      }
      if(depth == depth_vals[3]){
        col_name <- paste0("depth_pos_", lu)
      }
      
      
      # Helper function to find the correct column. Updates each iteration
      # for the appropriate land use
      col_fun <- function(x){
        paste0(x, "[", lu, "]")
      }
      
      init_vals <- ms_df[,"beta_init_intercept[2]"] + 
        ms_df[,"beta_init_depth[2]"] * sc_fun(var = "depth", natural_val = depth) + 
        
        ms_df[,"beta_init_clay[2]"] * sc_fun(var = "clay", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "clay"])) + 
        
        ms_df[,"beta_init_nitrogen[2]"] * sc_fun(var = "nitrogen", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "nitrogen"])) + 
        
        ms_df[,"beta_init_irrigation[2]"] * sc_fun(var = "irrigation", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "irrigation"])) + 
        
        ms_df[,"beta_init_map_mat[2]"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "map"])) * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "mat"])) +
        
        ms_df[,"beta_init_map[2]"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "map"])) + 
        
        ms_df[,"beta_init_mat[2]"] * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "mat"]))
      
      
      
      r_vals <- ms_df[,col_fun("beta_r_intercept")] +
        
        ms_df[,"beta_r_depth[2]"] * sc_fun(var = "depth", natural_val = depth) +
        
        ms_df[,"beta_r_clay[2]"] * sc_fun(var = "clay", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "clay"])) +
        
        ms_df[,"beta_r_nitrogen[2]"] * sc_fun(var = "nitrogen", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "nitrogen"])) +
        
        ms_df[,"beta_r_map_mat[2]"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "map"])) * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "mat"])) +
        
        ms_df[,"beta_r_map[2]"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "map"])) +
        
        ms_df[,"beta_r_mat[2]"] * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "mat"]))
      
      
      max_vals <- ms_df[,"beta_max_intercept"] +
        ms_df[,"beta_max_depth"] * sc_fun(var = "depth", natural_val = depth) + 
        
        ms_df[,"beta_max_clay"] * sc_fun(var = "clay", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "clay"])) + 
        
        ms_df[,"beta_max_nitrogen"] * sc_fun(var = "nitrogen", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "nitrogen"])) + 
        
        ms_df[,"beta_max_map_mat"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "map"])) * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "mat"])) +
        
        ms_df[,"beta_max_map"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "map"])) + 
        
        ms_df[,"beta_max_mat"] * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "mat"]))
      
      final_vals <- log(exp(init_vals) + (exp(max_vals) - exp(init_vals)) * (1 - exp(-age_fun(age) / exp(r_vals))) ^ exp(ms_df[,col_fun("m_val")]))
      
      perc_change_vals <- ((exp(final_vals) - exp(init_vals)) / exp(init_vals)) * 100
      
      ind <- which(age_out$age == age)
      age_out[ind, c(paste(col_name, "lo", sep = "_"),
                     paste(col_name, "med", sep = "_"),
                     paste(col_name, "hi", sep = "_"))] <- unlist(quantile(perc_change_vals,
                                                                           probs = c(0.025, 0.5, 0.975)))
      
    }
  }
}



# ylims <- c(-50, 300)

# 2×2 layout

ylims <- c(0, 300)
plot_panel <- function(active_i, natural_i, title, poly_col, line_col, trans_index){
  # Mapping from index to lu_start / lu_end
  trans_map <- data.frame(
    index    = 1:4,
    label    = c("Cropland → Forest", "Cropland → Open", "Pasture → Forest", "Pasture → Open"),
    lu_start = c("C", "C", "P", "P"),
    lu_end   = c("F", "O", "F", "O")
  )
  
  this_trans <- trans_map[trans_index, ]
  
  # Filter observations for this transition (active + passive)
  rug_df <- final_dat %>%
    filter(lu_start == this_trans$lu_start,
           lu_end   == this_trans$lu_end)
  
  rug_vals <- rug_df %>%
    group_by(final_profile_id) %>%
    summarise(max_age = max(AGE, na.rm = TRUE)) %>%
    pull(max_age)
  
  # Compute percent label
  trans_counts <- final_covars %>%
    mutate(group = paste(lu_start, lu_end)) %>%
    count(group)
  
  this_group <- paste(this_trans$lu_start, this_trans$lu_end)
  pct_val <- trans_counts %>%
    filter(group == this_group) %>%
    mutate(pct = round(100 * n / sum(trans_counts$n))) %>%
    pull(pct)
  
  pct_label <- if(length(pct_val) == 1) {
    paste0(pct_val, "% of records")
  } else {
    "NA% of records"
  }
  
  # plot setup
  plot(age_out$age, age_out[[paste0("depth_avg_", active_i, "_med")]],
       type = "n",
       ylim = ylims,
       xlab = "",
       xaxt = "n",
       yaxt = "n",
       ylab = ifelse(active_i %in% c(1,5), "Percent change in SOC density", ""),
       las = 1,
       frame = FALSE)
  mtext(side = 3, line = 0, title, cex = 0.8)
  
  if(active_i %in% c(1,5)){
    my_y_labels <- seq(0, ylims[2], by = 50)
  } else{
    my_y_labels <- rep("",length(seq(0, ylims[2], by = 50)))
  }
  
  axis(2, at = seq(0, ylims[2], by = 50), las = 1,
       labels = my_y_labels)
  
  if(active_i %in% c(5,7)){
    my_x_labels <- seq(0, 100, by = 20)
  } else{
    my_x_labels <- rep("",length(seq(0, 100, by = 20)))
  }
  
  axis(1, at = seq(0, 100, by = 20), las = 1,
       labels = my_x_labels)
  
  mtext(side = 1, line = 2.25, ifelse(active_i %in% c(5,7), "Years since restoration", ""), cex = 0.8)
  
  # Natural regrowth: shaded CI + solid median
  polygon(c(age_out$age, rev(age_out$age)),
          c(age_out[[paste0("depth_avg_", natural_i, "_lo")]],
            rev(age_out[[paste0("depth_avg_", natural_i, "_hi")]])),
          col   = poly_col,
          border= NA)
  lines(age_out$age,
        age_out[[paste0("depth_avg_", natural_i, "_med")]],
        lwd = 2,
        col = line_col)
  
  # Active restoration: dashed median + dashed CI
  lines(age_out$age,
        age_out[[paste0("depth_avg_", active_i, "_med")]],
        lty = 2, lwd = 1,
        col = line_col)
  lines(age_out$age,
        age_out[[paste0("depth_avg_", active_i, "_lo")]],
        lty = 2, lwd = 0.5,
        col = line_col)
  lines(age_out$age,
        age_out[[paste0("depth_avg_", active_i, "_hi")]],
        lty = 2, lwd = 0.5,
        col = line_col)
  
  # Rug: jitter, clip at 100, add alpha transparency
  if(length(rug_vals) > 0){
    rug((pmin(rug_vals, 100) + runif(length(rug_vals), -0.5,0.5)), ticksize = 0.02,
        col = adjustcolor(line_col, alpha.f = 0.4))
  }
  
  # Label: percent of all observations
  text(x = 100, y = 290, labels = pct_label, adj = 1, cex = 0.7)
}

add_label <- function(lab){
  u <- par("usr")
  text(
    x     = u[1] - 0.05*(u[2]-u[1]),
    y     = u[4] + 0.05*(u[4]-u[3]),
    labels = lab,
    font   = 1,
    cex    = 1.2,
    xpd    = NA
  )
}

# Make the plot
for (dev_type in c("png", "pdf")) {
  if (dev_type == "png")
    png(file = "./figures/SOC vs age four panel.png", width = 5, height = 5, units = "in", res = 440)
  else
    pdf(file = "./figures/SOC vs age four panel.pdf", width = 5, height = 5)
  par(mfrow = c(2,2))
  par(mar = c(2, 4, 2.5, 0.2) + 0.1)
  
  # Top row
  plot_panel(active_i  = 1, natural_i  = 2, title = "Cropland → Forest",
             line_col = "#d95f02",
             poly_col = "#f4c29f",
             trans_index = 1)
  add_label("A")
  legend(x = 0, y = ylims[2],
         legend = c("Passive", "Active"),
         lty    = c(1, 2),
         lwd    = c(2, 1),
         cex = 0.8,
         bty = "n")
  
  plot_panel(active_i  = 3, natural_i  = 4, title = "Cropland → Open",
             line_col = "#d95f02",
             poly_col = "#f4c29f",
             trans_index = 2)
  add_label("B")
  
  # Bottom row
  par(mar = c(4, 4, 0.5, 0.2) + 0.1)
  plot_panel(active_i  = 5, natural_i  = 6, title = "Pasture → Forest",
             line_col = "#7570b3",
             poly_col = "#c8c6e1",
             trans_index = 3)
  add_label("C")
  
  plot_panel(active_i  = 7, natural_i  = 8, title = "Pasture → Open",
             line_col = "#7570b3",
             poly_col = "#c8c6e1",
             trans_index = 4)
  add_label("D")
  
  dev.off()
  message("Saved ./figures/SOC vs age four panel.", dev_type)
}



# Summary stats
# 1) Grab the age = 100 row once
age100 <- age_out %>% filter(age == 100)

# 2) Define the eight transition labels and the matching column names
types    <- c(
  "C→F active",   "C→F passive",
  "C→O active",   "C→O passive",
  "P→F active",   "P→F passive",
  "P→O active",   "P→O passive"
)
med_cols <- paste0("depth_avg_", 1:8, "_med")
lo_cols  <- paste0("depth_avg_", 1:8, "_lo")
hi_cols  <- paste0("depth_avg_", 1:8, "_hi")

# 3) Build the summary table
summary_table <- tibble(
  type = types,
  lo   = map_dbl(lo_cols,  ~ age100[[.x]]),
  med  = map_dbl(med_cols, ~ age100[[.x]]),
  hi   = map_dbl(hi_cols, ~ age100[[.x]])
) %>%
  mutate(across(c(lo, med, hi), ~ signif(.x, 3)))

# 5) Preview
summary_table

# write it out
write.csv(summary_table,
          file = "./figures/table_age100_depth_avg_summary.csv",
          row.names = FALSE)





# Some additional posterior predictions. Specific contrasts (e.g., passive vs active)
post_pred_perc_change <- function(depth = 15, 
                                  age = 30,
                                  lu = 1){
  
  col_fun <- function(x){
    paste0(x, "[", lu, "]")
  }
  
  if(lu < 5){
    init_vals <- ms_df[,"beta_init_intercept[1]"] + # [1] is cropland, [1] is pasture
      ms_df[,"beta_init_depth[1]"] * sc_fun(var = "depth", natural_val = depth) + 
      
      ms_df[,"beta_init_clay[1]"] * sc_fun(var = "clay", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "clay"])) + 
      
      ms_df[,"beta_init_nitrogen[1]"] * sc_fun(var = "nitrogen", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "nitrogen"])) + 
      
      ms_df[,"beta_init_irrigation[1]"] * sc_fun(var = "irrigation", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "irrigation"])) + 
      
      ms_df[,"beta_init_map_mat[1]"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "map"])) * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "mat"])) +
      
      ms_df[,"beta_init_map[1]"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "map"])) + 
      
      ms_df[,"beta_init_mat[1]"] * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "mat"]))
    
    
    
    r_vals <- ms_df[,col_fun("beta_r_intercept")] +
      
      ms_df[,"beta_r_depth[1]"] * sc_fun(var = "depth", natural_val = depth) +
      
      ms_df[,"beta_r_clay[1]"] * sc_fun(var = "clay", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "clay"])) +
      
      ms_df[,"beta_r_nitrogen[1]"] * sc_fun(var = "nitrogen", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "nitrogen"])) +
      
      ms_df[,"beta_r_map_mat[1]"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "map"])) * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "mat"])) +
      
      ms_df[,"beta_r_map[1]"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "map"])) +
      
      ms_df[,"beta_r_mat[1]"] * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "mat"]))
    
    
    max_vals <- ms_df[,"beta_max_intercept"] +
      ms_df[,"beta_max_depth"] * sc_fun(var = "depth", natural_val = depth) + 
      
      ms_df[,"beta_max_clay"] * sc_fun(var = "clay", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "clay"])) + 
      
      ms_df[,"beta_max_nitrogen"] * sc_fun(var = "nitrogen", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "nitrogen"])) + 
      
      ms_df[,"beta_max_map_mat"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "map"])) * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "mat"])) +
      
      ms_df[,"beta_max_map"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "map"])) + 
      
      ms_df[,"beta_max_mat"] * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2), "mat"]))
    
    final_vals <- log(exp(init_vals) + (exp(max_vals) - exp(init_vals)) * (1 - exp(-age_fun(age) / exp(r_vals))) ^ exp(ms_df[,col_fun("m_val")]))
    
    perc_change_vals <- ((exp(final_vals) - exp(init_vals)) / exp(init_vals)) * 100
  } else{
    
    init_vals <- ms_df[,"beta_init_intercept[2]"] + # [1] is cropland, [2] is pasture
      ms_df[,"beta_init_depth[2]"] * sc_fun(var = "depth", natural_val = depth) + 
      
      ms_df[,"beta_init_clay[2]"] * sc_fun(var = "clay", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "clay"])) + 
      
      ms_df[,"beta_init_nitrogen[2]"] * sc_fun(var = "nitrogen", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "nitrogen"])) +
      
      ms_df[,"beta_init_irrigation[2]"] * sc_fun(var = "irrigation", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "irrigation"])) + 
      
      ms_df[,"beta_init_map_mat[2]"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "map"])) * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "mat"])) +
      
      ms_df[,"beta_init_map[2]"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "map"])) + 
      
      ms_df[,"beta_init_mat[2]"] * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "mat"]))
    
    
    
    r_vals <- ms_df[,col_fun("beta_r_intercept")] +
      
      ms_df[,"beta_r_depth[2]"] * sc_fun(var = "depth", natural_val = depth) +
      
      ms_df[,"beta_r_clay[2]"] * sc_fun(var = "clay", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "clay"])) +
      
      ms_df[,"beta_r_nitrogen[2]"] * sc_fun(var = "nitrogen", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "nitrogen"])) +
      
      ms_df[,"beta_r_map_mat[2]"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "map"])) * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "mat"])) +
      
      ms_df[,"beta_r_map[2]"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "map"])) +
      
      ms_df[,"beta_r_mat[2]"] * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "mat"]))
    
    
    max_vals <- ms_df[,"beta_max_intercept"] +
      ms_df[,"beta_max_depth"] * sc_fun(var = "depth", natural_val = depth) + 
      
      ms_df[,"beta_max_clay"] * sc_fun(var = "clay", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "clay"])) + 
      
      ms_df[,"beta_max_nitrogen"] * sc_fun(var = "nitrogen", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "nitrogen"])) + 
      
      ms_df[,"beta_max_map_mat"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "map"])) * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "mat"])) +
      
      ms_df[,"beta_max_map"] * sc_fun(var = "map", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "map"])) + 
      
      ms_df[,"beta_max_mat"] * sc_fun(var = "mat", natural_val = as.numeric(lu_avg[[which(depth_vals == depth)]][ceiling(lu/2)-2, "mat"]))
    
    final_vals <- log(exp(init_vals) + (exp(max_vals) - exp(init_vals)) * (1 - exp(-age_fun(age) / exp(r_vals))) ^ exp(ms_df[,col_fun("m_val")]))
    
    perc_change_vals <- ((exp(final_vals) - exp(init_vals)) / exp(init_vals)) * 100
    
  }
  return(perc_change_vals)
}


# 
prop_over_zero <- function(x) mean(x > 0)

# the four contrasts we care about:
contrasts <- tibble(
  name     = c("C→F", "C→O", "P→F", "P→O"),
  lu_act   = c(2,      4,      6,      8   ),
  lu_nat   = c(1,      3,      5,      7   )
)

# ages to evaluate
ages <- c(10, 100)

# build the table
results <- map_dfr(ages, function(age_val) {
  contrasts %>%
    mutate(age = age_val) %>%
    pmap_dfr(function(name, lu_act, lu_nat, age) {
      # get posterior draws of percent change
      post_act <- post_pred_perc_change(depth = 15, age = age, lu = lu_act)
      post_nat <- post_pred_perc_change(depth = 15, age = age, lu = lu_nat)
      diff     <- post_act - post_nat
      
      # summarise
      tibble(
        contrast = name,
        age       = age,
        lower     = quantile(diff, 0.025),
        median    = median(diff),
        upper     = quantile(diff, 0.975),
        prop_pos  = prop_over_zero(diff)
      )
    })
})

# reorder columns & rows nicely
results <- results %>% 
  dplyr::select(age, contrast, lower, median, upper, prop_pos) %>%
  arrange(age, contrast)

# write to CSV
write.csv(results,
          file = "./figures/table_active_vs_natural_restoration.csv",
          row.names = FALSE)










#—— 1. set up which env vars to do ——————————————————————————————
env_vars   <- c("mat", "map", "clay", "nitrogen")
env_labels <- c("Mean annual\ntemperature (°C)",
                "Mean annual\nprecipitation (mm)",
                "Clay content\n(g/kg)",
                "Soil nitrogen\n(mg/kg)")

#—— 2. helper to get a sequence from –2SD to +2SD in ORIGINAL units ———
make_grid <- function(var, n = 100) {
  # pull center & scale from sc_df
  ctr <- sc_df[var, "center"]
  scl <- sc_df[var, "scale"]
  # create a sequence of scaled values from -2 to +2
  scaled_vals <- seq(-1.5, 1.5, length.out = n)
  
  if (sc_df[var, "log"]) {
    # center/scale are on log(original), so
    # log_val = center + scaled * scale
    log_val    <- ctr + scaled_vals * scl
    original   <- exp(log_val)
  } else {
    # variables that weren’t logged: center/scale are on raw units
    original   <- ctr + scaled_vals * scl
  }
  
  tibble(
    scaled   = scaled_vals,
    original = original
  )
}

#—— 3. fixed depth and other covars at “mean” (i.e. scaled = 0) ————————
depth_s  <- sc_fun("depth", natural_val = 15)       # always 15 cm
sand_s   <- 0  # we drop sand, pH, etc. from the final model
phh2o_s  <- 0
irr_s    <- sc_fun("irrigation", natural_val = exp(mean(log(lu_avg[[2]]$irrigation)))) 

#—— 4. build a big data.frame of predictions ——————————————————————
predictions <-
  map2_dfr(env_vars, env_labels, function(var, lab) {
    
    grid_df <- make_grid(var, n = 100)
    
    
    map_dfr(c(2,6), function(transition_idx) {
      
      # for each land‐use (1=cropland, 2=pasture)
      
      lu_idx <- ifelse(transition_idx < 5, 1, 2)
      
      # run through each grid‐point
      map_dfr(seq_len(nrow(grid_df)), function(i) {
        x_s    <- grid_df$scaled[i]
        # set up all scaled covars
        clay_s   <- if (var=="clay")    x_s else sc_fun("clay",    natural_val = NULL, sd_val = 0)
        nit_s    <- if (var=="nitrogen") x_s else sc_fun("nitrogen", natural_val = NULL, sd_val = 0)
        mat_s    <- if (var=="mat")     x_s else sc_fun("mat",      natural_val = NULL, sd_val = 0)
        map_s    <- if (var=="map")     x_s else sc_fun("map",      natural_val = NULL, sd_val = 0)
        
        #—— 5. pull out all the posterior draws from ms_df ——————————————
        # beta_init_intercept[1] etc.
        init_vals <- 
          ms_df[, paste0("beta_init_intercept[", lu_idx, "]")] +
          ms_df[, paste0("beta_init_depth[",     lu_idx, "]")] * depth_s +
          ms_df[, paste0("beta_init_clay[",      lu_idx, "]")] * clay_s +
          ms_df[, paste0("beta_init_nitrogen[",  lu_idx, "]")] * nit_s +
          ms_df[, paste0("beta_init_map[",       lu_idx, "]")] * map_s +
          ms_df[, paste0("beta_init_mat[",       lu_idx, "]")] * mat_s +
          ms_df[, paste0("beta_init_map_mat[",   lu_idx, "]")] * map_s * mat_s +
          ms_df[, paste0("beta_init_irrigation[",lu_idx, "]")] * irr_s
        
        r_vals <- 
          ms_df[, paste0("beta_r_intercept[", transition_idx, "]")] +
          ms_df[, paste0("beta_r_depth[",     lu_idx, "]")] * depth_s +
          ms_df[, paste0("beta_r_clay[",      lu_idx, "]")] * clay_s +
          ms_df[, paste0("beta_r_nitrogen[",  lu_idx, "]")] * nit_s +
          ms_df[, paste0("beta_r_map[",       lu_idx, "]")] * map_s +
          ms_df[, paste0("beta_r_mat[",       lu_idx, "]")] * mat_s +
          ms_df[, paste0("beta_r_map_mat[",   lu_idx, "]")] * map_s * mat_s
        
        max_vals <- 
          ms_df[, "beta_max_intercept"] +
          ms_df[, "beta_max_depth"]   * depth_s +
          ms_df[, "beta_max_clay"]    * clay_s +
          ms_df[, "beta_max_nitrogen"]* nit_s +
          ms_df[, "beta_max_map"]     * map_s +
          ms_df[, "beta_max_mat"]     * mat_s +
          ms_df[, "beta_max_map_mat"] * map_s * mat_s
        
        m_vals <- ms_df[, paste0("m_val[", transition_idx, "]")]
        
        #—— 6. compute final C at age = 100, then percent deficit ————————
        age_val <- 100
        vals_100 <- log(
          exp(init_vals) +
            (exp(max_vals) - exp(init_vals)) *
            (1 - exp(-age_val / exp(r_vals)))^exp(m_vals)
        )
        
        pct_deficit <- -100 * (1 - exp(init_vals) / exp(vals_100))
        
        #—— 7. summarise posterior ————————————————————————————————
        tibble(
          var        = lab,
          original   = grid_df$original[i],
          lu         = if (lu_idx==1) "Cropland" else "Pasture",
          transition = transition_idx,
          median     = median(pct_deficit),
          lower      = quantile(pct_deficit, 0.025),
          upper      = quantile(pct_deficit, 0.975)
        )
      })
    })
  })

#—— now plot! ——————————————————————————————————————————————

# 1. panel definitions
env_list <- list(
  list(key = "Mean annual\ntemperature (°C)", var = "Mean annual\ntemperature (°C)"),
  list(key = "Mean annual\nprecipitation (mm)", var = "Mean annual\nprecipitation (mm)"),
  list(key = "Clay content\n(g/kg)",                var = "Clay content\n(g/kg)"),
  list(key = "Soil nitrogen\n(mg/kg)",         var = "Soil nitrogen\n(mg/kg)")
)

# 2. force y–limits 0–100
ylims <- c(-100, 0)

# 3. open device and do 2×2 panels, with generous bottom margin for x–labels
for (dev_type in c("png", "pdf")) {
  if (dev_type == "png")
    png(file = "./figures/deficits_marginal.png", width = 5, height = 5, units = "in", res = 440)
  else
    pdf(file = "./figures/deficits_marginal.pdf", width = 5, height = 5)
  par(mfrow = c(2,2),
      mar   = c(4, 4, 1.5, 0.3) + 0.1)
  
  for(i in seq_along(env_list)) {
    info <- env_list[[i]]
    dfv  <- subset(predictions, var == info$var)
    cro  <- subset(dfv, lu == "Cropland")
    pas  <- subset(dfv, lu == "Pasture")
    o1   <- order(cro$original); x1 <- cro$original[o1]
    o2   <- order(pas$original); x2 <- pas$original[o2]
    
    # empty plot with custom xlab/ylab but no axes yet
    plot(x1, cro$median[o1], type = "n", ylim = ylims,
         xlab = info$key,
         ylab = if(i %in% c(1,3)) "Percent deficit in SOC" else "",
         frame = FALSE, axes = FALSE)
    
    # draw both axes
    axis(1)        # bottom
    #axis(2, las = 1)
    
    # Cropland ribbon + line
    polygon(c(x1, rev(x1)),
            c(cro$lower[o1], rev(cro$upper[o1])),
            col   = rgb(217/255, 95/255,  2/255, 0.2),
            border= NA)
    lines(x1, cro$median[o1], col = "#d95f02", lwd = 2)
    
    # Pasture ribbon + line
    polygon(c(x2, rev(x2)),
            c(pas$lower[o2], rev(pas$upper[o2])),
            col   = rgb( 27/255,158/255,119/255, 0.2),
            border= NA)
    lines(x2, pas$median[o2], col = "#7570b3", lwd = 2)
    
    abline(h = 0, lty = 2)
    
    if(i %in% c(1,3)){
      axis(2, at  = seq(-100,0, by = 20), 
           labels = seq(100,0, by = -20),
           las = 1)
    }else{
      axis(2, at  = seq(-100,0, by = 20), 
           labels = rep("",6),
           las = 1)
    }
    
    # legend only in first panel
    if(i == 1) {
      legend("bottomright",
             legend = c("Cropland", "Pasture"),
             col    = c("#d95f02", "#7570b3"),
             lwd    = 2,
             bty    = "n")
    }
    add_label(c("A","B","C","D")[i])
  }
  
  dev.off()
  message("Saved ./figures/deficits_marginal.", dev_type)
}











# Alternate version where it's about restoration potential
predictions <-
  map2_dfr(env_vars, env_labels, function(var, lab) {
    
    grid_df <- make_grid(var, n = 100)
    
    
    map_dfr(c(2,6), function(transition_idx) {
      
      # for each land‐use (1=cropland, 2=pasture)
      
      lu_idx <- ifelse(transition_idx < 5, 1, 2)
      
      # run through each grid‐point
      map_dfr(seq_len(nrow(grid_df)), function(i) {
        x_s    <- grid_df$scaled[i]
        # set up all scaled covars
        clay_s   <- if (var=="clay")    x_s else sc_fun("clay",    natural_val = NULL, sd_val = 0)
        nit_s    <- if (var=="nitrogen") x_s else sc_fun("nitrogen", natural_val = NULL, sd_val = 0)
        mat_s    <- if (var=="mat")     x_s else sc_fun("mat",      natural_val = NULL, sd_val = 0)
        map_s    <- if (var=="map")     x_s else sc_fun("map",      natural_val = NULL, sd_val = 0)
        
        #—— 5. pull out all the posterior draws from ms_df ——————————————
        # beta_init_intercept[1] etc.
        init_vals <- 
          ms_df[, paste0("beta_init_intercept[", lu_idx, "]")] +
          ms_df[, paste0("beta_init_depth[",     lu_idx, "]")] * depth_s +
          ms_df[, paste0("beta_init_clay[",      lu_idx, "]")] * clay_s +
          ms_df[, paste0("beta_init_nitrogen[",  lu_idx, "]")] * nit_s +
          ms_df[, paste0("beta_init_map[",       lu_idx, "]")] * map_s +
          ms_df[, paste0("beta_init_mat[",       lu_idx, "]")] * mat_s +
          ms_df[, paste0("beta_init_map_mat[",   lu_idx, "]")] * map_s * mat_s +
          ms_df[, paste0("beta_init_irrigation[",lu_idx, "]")] * irr_s
        
        r_vals <- 
          ms_df[, paste0("beta_r_intercept[", transition_idx, "]")] +
          ms_df[, paste0("beta_r_depth[",     lu_idx, "]")] * depth_s +
          ms_df[, paste0("beta_r_clay[",      lu_idx, "]")] * clay_s +
          ms_df[, paste0("beta_r_nitrogen[",  lu_idx, "]")] * nit_s +
          ms_df[, paste0("beta_r_map[",       lu_idx, "]")] * map_s +
          ms_df[, paste0("beta_r_mat[",       lu_idx, "]")] * mat_s +
          ms_df[, paste0("beta_r_map_mat[",   lu_idx, "]")] * map_s * mat_s
        
        max_vals <- 
          ms_df[, "beta_max_intercept"] +
          ms_df[, "beta_max_depth"]   * depth_s +
          ms_df[, "beta_max_clay"]    * clay_s +
          ms_df[, "beta_max_nitrogen"]* nit_s +
          ms_df[, "beta_max_map"]     * map_s +
          ms_df[, "beta_max_mat"]     * mat_s +
          ms_df[, "beta_max_map_mat"] * map_s * mat_s
        
        m_vals <- ms_df[, paste0("m_val[", transition_idx, "]")]
        
        #—— 6. compute final C at age = 100, then percent deficit ————————
        age_val <- 100
        vals_100 <- log(
          exp(init_vals) +
            (exp(max_vals) - exp(init_vals)) *
            (1 - exp(-age_val / exp(r_vals)))^exp(m_vals)
        )
        #pct_deficit <- 100 * (1 - exp(vals_100) / exp(max_vals))
        
        pct_improvement <- 100 * (exp(vals_100) / exp(init_vals)) - 100
        
        #—— 7. summarise posterior ————————————————————————————————
        tibble(
          var        = lab,
          original   = grid_df$original[i],
          lu         = if (lu_idx==1) "Cropland" else "Pasture",
          transition = transition_idx,
          median     = median(pct_improvement),
          lower      = quantile(pct_improvement, 0.025),
          upper      = quantile(pct_improvement, 0.975)
        )
      })
    })
  })



for (dev_type in c("png", "pdf")) {
  if (dev_type == "png")
    png(file = "./figures/improvements_marginal.png", width = 5, height = 5, units = "in", res = 440)
  else
    pdf(file = "./figures/improvements_marginal.pdf", width = 5, height = 5)
  
  ylims <- c(0, 500)
  par(mfrow = c(2,2))
  
  
  for(i in seq_along(env_list)) {
    
    
    if(i %in% c(1,3)){
      par(mar   = c(4, 5, 1.5, 0.1) + 0.1)
    } else{
      par(mar   = c(4, 4, 1.5, 1.1) + 0.1)
    }
    
    
    info <- env_list[[i]]
    dfv  <- subset(predictions, var == info$var)
    cro  <- subset(dfv, lu == "Cropland")
    pas  <- subset(dfv, lu == "Pasture")
    o1   <- order(cro$original); x1 <- cro$original[o1]
    o2   <- order(pas$original); x2 <- pas$original[o2]
    
    # empty plot with custom xlab/ylab but no axes yet
    plot(x1, cro$median[o1], type = "n", ylim = ylims,
         xlab = info$key,
         ylab = if(i %in% c(1,3)) "SOC restoration potential\n(% change)" else "",
         frame = FALSE, axes = FALSE)
    
    # draw both axes
    axis(1)        # bottom
    #axis(2, las = 1)
    
    # Cropland ribbon + line
    polygon(c(x1, rev(x1)),
            c(cro$lower[o1], rev(cro$upper[o1])),
            col   = rgb(217/255, 95/255,  2/255, 0.2),
            border= NA)
    lines(x1, cro$median[o1], col = "#d95f02", lwd = 2)
    
    # Pasture ribbon + line
    polygon(c(x2, rev(x2)),
            c(pas$lower[o2], rev(pas$upper[o2])),
            col   = rgb( 27/255,158/255,119/255, 0.2),
            border= NA)
    lines(x2, pas$median[o2], col = "#7570b3", lwd = 2)
    
    abline(h = 0, lty = 2)
    
    y_axis_vals <- seq(0,500, by = 100)
    
    if(i %in% c(1,3)){
      axis(2, at  = y_axis_vals, 
           labels = y_axis_vals,
           las = 1)
    }else{
      axis(2, at  = y_axis_vals, 
           labels = rep("",length(y_axis_vals)),
           las = 1)
    }
    
    # legend only in first panel
    if(i == 1) {
      legend("topright",
             legend = c("Cropland", "Pasture"),
             col    = c("#d95f02", "#7570b3"),
             lwd    = 2,
             bty    = "n")
    }
    add_label(c("A","B","C","D")[i])
  }
  
  dev.off()
  message("Saved ./figures/improvements_marginal.", dev_type)
}


#—— 3. fixed depth and other covars at “mean” (i.e. scaled = 0) ————————
depth_s  <- sc_fun("depth", natural_val = 15)       # always 15 cm
sand_s   <- 0  # we drop sand, pH, etc. from the final model
phh2o_s  <- 0
irr_s    <- sc_fun("irrigation", natural_val = exp(mean(log(lu_avg[[2]]$irrigation)))) 

#—— 4. build a big data.frame of predictions ——————————————————————
predictions <-
  map2_dfr(env_vars, env_labels, function(var, lab) {
    
    grid_df <- make_grid(var, n = 100)
    
    
    map_dfr(c(2,6), function(transition_idx) {
      
      # for each land‐use (1=cropland, 2=pasture)
      
      lu_idx <- ifelse(transition_idx < 5, 1, 2)
      
      # run through each grid‐point
      map_dfr(seq_len(nrow(grid_df)), function(i) {
        x_s    <- grid_df$scaled[i]
        # set up all scaled covars
        clay_s   <- if (var=="clay")    x_s else sc_fun("clay",    natural_val = NULL, sd_val = 0)
        nit_s    <- if (var=="nitrogen") x_s else sc_fun("nitrogen", natural_val = NULL, sd_val = 0)
        mat_s    <- if (var=="mat")     x_s else sc_fun("mat",      natural_val = NULL, sd_val = 0)
        map_s    <- if (var=="map")     x_s else sc_fun("map",      natural_val = NULL, sd_val = 0)
        
        #—— 5. pull out all the posterior draws from ms_df ——————————————
        # beta_init_intercept[1] etc.
        init_vals <- 
          ms_df[, paste0("beta_init_intercept[", lu_idx, "]")] +
          ms_df[, paste0("beta_init_depth[",     lu_idx, "]")] * depth_s +
          ms_df[, paste0("beta_init_clay[",      lu_idx, "]")] * clay_s +
          ms_df[, paste0("beta_init_nitrogen[",  lu_idx, "]")] * nit_s +
          ms_df[, paste0("beta_init_map[",       lu_idx, "]")] * map_s +
          ms_df[, paste0("beta_init_mat[",       lu_idx, "]")] * mat_s +
          ms_df[, paste0("beta_init_map_mat[",   lu_idx, "]")] * map_s * mat_s +
          ms_df[, paste0("beta_init_irrigation[",lu_idx, "]")] * irr_s
        
        r_vals <- 
          ms_df[, paste0("beta_r_intercept[", transition_idx, "]")] +
          ms_df[, paste0("beta_r_depth[",     lu_idx, "]")] * depth_s +
          ms_df[, paste0("beta_r_clay[",      lu_idx, "]")] * clay_s +
          ms_df[, paste0("beta_r_nitrogen[",  lu_idx, "]")] * nit_s +
          ms_df[, paste0("beta_r_map[",       lu_idx, "]")] * map_s +
          ms_df[, paste0("beta_r_mat[",       lu_idx, "]")] * mat_s +
          ms_df[, paste0("beta_r_map_mat[",   lu_idx, "]")] * map_s * mat_s
        
        max_vals <- 
          ms_df[, "beta_max_intercept"] +
          ms_df[, "beta_max_depth"]   * depth_s +
          ms_df[, "beta_max_clay"]    * clay_s +
          ms_df[, "beta_max_nitrogen"]* nit_s +
          ms_df[, "beta_max_map"]     * map_s +
          ms_df[, "beta_max_mat"]     * mat_s +
          ms_df[, "beta_max_map_mat"] * map_s * mat_s
        
        m_vals <- ms_df[, paste0("m_val[", transition_idx, "]")]
        
        #—— 6. compute final C at age = 100, then percent deficit ————————
        age_val <- 100
        vals_100 <- log(
          exp(init_vals) +
            (exp(max_vals) - exp(init_vals)) *
            (1 - exp(-age_val / exp(r_vals)))^exp(m_vals)
        )
        #pct_deficit <- 100 * (1 - exp(vals_100) / exp(max_vals))
        
        pct_restoration <- 100 * (exp(vals_100) / exp(init_vals)) - 100
        
        #—— 7. summarise posterior ————————————————————————————————
        tibble(
          var        = lab,
          original   = grid_df$original[i],
          lu         = if (lu_idx==1) "Cropland" else "Pasture",
          transition = transition_idx,
          init_vals_median     = median(exp(init_vals)),
          init_vals_lower      = quantile(exp(init_vals), 0.025),
          init_vals_upper      = quantile(exp(init_vals), 0.975),
          
          vals_100_median     = median(exp(vals_100)),
          vals_100_lower      = quantile(exp(vals_100), 0.025),
          vals_100_upper      = quantile(exp(vals_100), 0.975),
          
          
          vals_restoration_median     = median(pct_restoration),
          vals_restoration_lower      = quantile(pct_restoration, 0.025),
          vals_restoration_upper      = quantile(pct_restoration, 0.975)
        )
      })
    })
  })

#—— now plot! ——————————————————————————————————————————————

# 1. panel definitions
env_list <- list(
  list(key = "Mean annual\ntemperature (°C)", var = "Mean annual\ntemperature (°C)"),
  list(key = "Mean annual\nprecipitation (mm)", var = "Mean annual\nprecipitation (mm)"),
  list(key = "Clay content\n(g/kg)",                var = "Clay content\n(g/kg)"),
  list(key = "Soil nitrogen\n(mg/kg)",         var = "Soil nitrogen\n(mg/kg)")
)

env_vars   <- c("mat","map","clay","nitrogen")
env_labs   <- c("Mean annual\ntemperature (°C)",
                "Mean annual\nprecipitation (mm)",
                "Clay content\n(g/kg)",
                "Soil nitrogen\n(mg/kg)")

# 2. force y–limits 0–100
ylims <- c(0, 10)

# 3. open device and do 3×4 panels, with generous bottom margin for x–labels
for (dev_type in c("png", "pdf")) {
  if (dev_type == "png")
    png(file = "./figures/c_dens_and_deficits_marginal.png", width = 7.5, height = 6.5, units = "in", res = 440)
  else
    pdf(file = "./figures/c_dens_and_deficits_marginal.pdf", width = 7.5, height = 6.5)
  par(mfrow = c(3,4))
  
  for(row in 1:3) {
    if(row == 1) par(mar = c(1, 3.5, 4, 0.5) + 0.1)
    
    if(row == 2) par(mar = c(3, 3.5, 2, 0.5) + 0.1)
    
    if(row == 3) par(mar = c(5, 3.5, 0, 0.5) + 0.1)
    
    for(col in 1:4) {
      key <- env_labs[col]
      
      # select data depending on row
      if (row == 1) {
        LU   <- "Cropland";     mycol <- "#d95f02"
        dfv  <- subset(predictions, var == key & lu == LU)
        dfv  <- dfv[order(dfv$original), ]
        y1   <- dfv$vals_100_median
        lo1  <- dfv$vals_100_lower
        hi1  <- dfv$vals_100_upper
        y2   <- dfv$init_vals_median
        lo2  <- dfv$init_vals_lower
        hi2  <- dfv$init_vals_upper
        x    <- dfv$original
        ylims <- range(c(lo1, hi1, lo2, hi2))
        
      } else if (row == 2) {
        LU   <- "Pasture";      mycol <- "#7570b3"
        dfv  <- subset(predictions, var == key & lu == LU)
        dfv  <- dfv[order(dfv$original), ]
        y1   <- dfv$vals_100_median
        lo1  <- dfv$vals_100_lower
        hi1  <- dfv$vals_100_upper
        y2   <- dfv$init_vals_median
        lo2  <- dfv$init_vals_lower
        hi2  <- dfv$init_vals_upper
        x    <- dfv$original
        ylims <- range(c(lo1, hi1, lo2, hi2))
        
      } else {
        # row 3: restoration potential
        dfC  <- subset(predictions, var == key & lu == "Cropland") %>% arrange(original)
        dfP  <- subset(predictions, var == key & lu == "Pasture")  %>% arrange(original)
        x    <- dfC$original
        yC   <- dfC$vals_restoration_median
        loC  <- dfC$vals_restoration_lower
        hiC  <- dfC$vals_restoration_upper
        yP   <- dfP$vals_restoration_median
        loP  <- dfP$vals_restoration_lower
        hiP  <- dfP$vals_restoration_upper
        ylims <- range(c(loC, hiC, loP, hiP))
      }
      
      # Instead, set ylims manually to be the same by row
      if(row == 1) ylims <- c(0, 5)
      if(row == 2) ylims <- c(0, 5)
      if(row == 3) ylims <- c(0, 350)
      
      # draw empty frame
      plot(x, x, type = "n", ylim = ylims,
           xlab = "",
           ylab = "",
           axes = FALSE, frame = FALSE)
      mtext(if(row == 3) key else "",
            line = 3.5,
            side = 1,
            cex = 0.8)
      
      mtext(if(col == 1) {
        if(row < 3) "Soil C density" else "Restoration potential (%)"
      } else "",
      line = 2.5,
      side = 2,
      cex = 0.8)
      
      ii <- ifelse(row == 1 & col == 1, 1, ii)
      add_label(LETTERS[ii])
      ii <- ii + 1
      
      
      # y‐axis ticks & labels only on col 1
      axis(2, las = 1)
      
      # x‐axis ticks & labels only on row 3
      axis(1)
      
      # plot data
      if (row < 3) {
        # 100‐yr ribbon + solid line
        polygon(c(x, rev(x)), c(lo1, rev(hi1)),
                col = adjustcolor(mycol, alpha.f = 0.2), border = NA, xpd = T)
        lines(x, y1, col = mycol, lwd = 2, lty = 1)
        # init ribbon + dashed line
        polygon(c(x, rev(x)), c(lo2, rev(hi2)),
                col = adjustcolor(mycol, alpha.f = 0.2), border = NA, xpd = T)
        lines(x, y2, col = mycol, lwd = 2, lty = 2)
        
        # legend in first panel
        if (row == 1 && col == 1) {
          legend("topright",
                 legend = c("Restored", "Cropland"),
                 lty    = c(1,2), lwd = 2,
                 col = "#d95f02",
                 bty    = "n", inset = 0.02)
        }
        
        # legend in first panel
        if (row == 2 && col == 1) {
          legend("topright",
                 legend = c("Restored", "Pasture"),
                 lty    = c(1,2), lwd = 2,
                 col = "#7570b3",
                 bty    = "n", inset = 0.02)
        }
        
      } else {
        # row 3: restoration potential CI + lines for both LU
        # Cropland
        polygon(c(x, rev(x)), c(loC, rev(hiC)),
                col = adjustcolor("#d95f02", alpha.f = 0.2), border = NA, xpd = T)
        lines(x, yC, col = "#d95f02", lwd = 2, lty = 1)
        # Pasture
        polygon(c(x, rev(x)), c(loP, rev(hiP)),
                col = adjustcolor("#7570b3", alpha.f = 0.2), border = NA, xpd = T)
        lines(x, yP, col = "#7570b3", lwd = 2, lty = 1)
        
        # legend in first bottom panel
        if (col == 1) {
          legend("topright",
                 legend = c("Cropland", "Pasture"),
                 col   = c("#d95f02", "#7570b3"),
                 lwd   = 2, bty = "n", inset = 0.02)
        }
      }
    }
  }
  
  dev.off()
  message("Saved ./figures/c_dens_and_deficits_marginal.", dev_type)
}