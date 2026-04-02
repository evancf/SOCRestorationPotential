# =============================================================================
# 01_soc_data_preparation.R
#
# PURPOSE
# -------
# This script reads the raw data spreadsheet, cleans and harmonises all
# measurements, extracts spatially-varying covariates (MAT, MAP, bulk density,
# clay, nitrogen, irrigation), applies gap-filling where observations are
# missing, and outputs two files that downstream scripts depend on:
#
#   ./data_inputs/c_dat_analysis_ready.csv   — analysis-ready observation-level
#                                              data (one row per depth interval
#                                              per chronosequence comparison)
#   ./data_inputs/sc_df.csv                  — centering / scaling parameters
#                                              used by 01_lc_priors.R and
#                                              02_bayes_analysis.R
#
# SCRIPT DEPENDENCIES (must be run first)
# ----------------------------------------
#   00_spatial_data.R  — downloads raster inputs into ./data_inputs/raster/
#
# INPUTS
# ------
#   ./data_inputs/Global_data_cleaned_region_260311.xlsx
#   ./data_inputs/missing_lu_end.csv
#   ./data_inputs/raster/cru_ts4.06.2011.2020.tmp.dat.nc
#   ./data_inputs/raster/cru_ts4.06.2011.2020.pre.dat.nc
#   ./data_inputs/raster/via_gee/bdod_*_mean.tif   (5 depth layers)
#   ./data_inputs/raster/via_gee/clay_*_mean.tif   (6 depth layers)
#   ./data_inputs/raster/via_gee/nitrogen_*_mean.tif (6 depth layers)
#   ./data_inputs/raster/gmia_plus_0.1.tif
#
# OUTPUTS
# -------
#   ./data_inputs/c_dat_analysis_ready.csv  — includes Ctrl/Aban_c_dens_av_orig_bd
#                                             columns for sensitivity run (d)
#   ./data_inputs/sc_df.csv
#   ./data_inputs/bd_validation_inputs.rda  — BD gap-filling model objects,
#                                             saved so that the validation
#                                             diagnostics in Section 11 can be
#                                             re-run from 03_bayes_analysis.R
#                                             without re-running this script
#
# =============================================================================
# SENSITIVITY ANALYSIS FILTER COLUMNS
# =============================================================================
# c_dat_analysis_ready.csv carries several boolean / categorical columns
# designed to make it straightforward to construct the subsetted datasets
# needed for the three pre-planned sensitivity analyses in 02_bayes_analysis.R.
# Filtering instructions are repeated at the top of 02_bayes_analysis.R.
#
# (a) BULK DENSITY GAP-FILLING SENSITIVITY
#     Tests whether gap-filled BD values materially affect conclusions.
#     Filter column: bd_gap_filled  (logical)
#       TRUE  — at least one of Ctrl_BD or Aban_BD was gap-filled (i.e. the
#               raw published value was absent or flagged as problematic and
#               the model-predicted value was substituted)
#       FALSE — both Ctrl_BD and Aban_BD come directly from the source
#               publication with no model involvement
#     Sensitivity run: filter(c_dat, !bd_gap_filled)
#     Supporting detail columns: Ctrl_BD_raw, Aban_BD_raw (published values,
#       NA where not reported), bd_SG (SoilGrids extraction),
#       Ctrl/Aban_BD_problematic (repeated-value flag),
#       Ctrl/Aban_BD_source ("original_study" or "model_soilgrids_climate")
#
# (b) MANAGEMENT UNCERTAINTY SENSITIVITY
#     Tests whether records with uncertain management ("?") bias results.
#     Filter column: MGMT_is_uncertain  (logical)
#       TRUE  — management was recorded as "?" in the source spreadsheet and
#               could not be resolved; the record has been assigned to "N"
#               (natural) for the primary analysis but its attribution is
#               genuinely ambiguous
#       FALSE — management is unambiguously "A" (assisted/active) or "N"
#               (natural/passive)
#     Primary analysis coding: MGMT column uses "A" / "N" throughout;
#       uncertain records are assigned "N" (conservative assumption)
#     Sensitivity run: filter(c_dat, !MGMT_is_uncertain)
#     Supporting detail column: MGMT_raw — the original code from the
#       spreadsheet before any recoding ("A", "N", "N/A", "?", etc.)
#
# (c) SOM-TO-SOC CONVERSION SENSITIVITY
#     Tests whether application of a fixed van Bemmelen factor (0.5)
#     to SOM records introduces systematic bias.
#     Filter column: som_converted  (logical)
#       TRUE  — the stock for this row was derived by applying the 0.5
#               SOM→SOC conversion factor (i.e. the original measurement was
#               organic matter, not organic carbon)
#       FALSE — the stock was measured directly as SOC or reported as a
#               carbon stock; no conversion factor was applied
#     Sensitivity run: filter(c_dat, !som_converted)
#
# (d) ORIGINAL BULK DENSITY SENSITIVITY
#     Tests whether the decision to replace problematic BD values with
#     model-predicted values affects conclusions, by instead using the
#     original reported BD throughout (including repeated-value entries).
#     Columns: Ctrl_c_dens_av_orig_bd / Aban_c_dens_av_orig_bd
#       Computed identically to Ctrl/Aban_c_dens_av but using Ctrl/Aban_BD_raw
#       (the published value, NA where not reported) instead of the final
#       gap-filled Ctrl/Aban_BD.  Rows where Ctrl_BD_raw or Aban_BD_raw is NA
#       will have NA in these columns and must be dropped before passing to
#       JAGS (see 03_bayes_analysis.R Section 5, sensitivity run (d)).
#     Note: rows flagged as Ctrl/Aban_BD_problematic retain their original
#       reported BD in these columns rather than the model-predicted value.
#
# =============================================================================
# OTHER "_raw" PROVENANCE COLUMNS
# =============================================================================
# Beyond the sensitivity filter columns above, one further column is modified
# during preparation and its pre-modification value is preserved:
#
#   MGMT_raw
#       Documented above under sensitivity (b).
#
#   Ctrl_BD_raw / Aban_BD_raw
#       Documented above under sensitivity (a) / BD column guide.
#
# The following transformations do NOT require separate _raw columns because
# the original values are already recoverable from columns that are retained
# unchanged in the output CSV:
#
#   lu_start / lu_end  — both derived from LUC, which is written as-is
#   Region             — single label fix; LUC encodes no regional grouping
#   MAT / MAP          — _source columns distinguish published vs. CRU values;
#                        raw values are NA (not overwritten) where CRU fills in
#
# =============================================================================
# BULK DENSITY COLUMN GUIDE (in c_dat_analysis_ready.csv)
# =============================================================================
#   Ctrl_BD_raw / Aban_BD_raw
#       The BD value exactly as reported in the source publication (after unit
#       correction to g cm⁻³).  NA where not reported.  Never modified after
#       ingestion.
#
#   bd_SG
#       SoilGrids 2.0 BD at the depth layer whose midpoint is closest to
#       depth_avg.  Extracted via GEE rasters; ocean/nodata pixels filled by
#       averaging over an expanding buffer (1 → 5 → 25 km).  Units: g cm⁻³.
#
#   Ctrl_BD_problematic / Aban_BD_problematic  (logical)
#       TRUE for multi-depth profiles where every interval shares the same
#       rounded BD — consistent with a single imputed value rather than
#       depth-resolved measurements.  Treated as missing for gap-filling.
#
#   Ctrl_BD / Aban_BD
#       Final BD used for stock conversion.  Equal to Ctrl/Aban_BD_raw when
#       available and not problematic; otherwise the fixed-effect prediction
#       from the chosen gap-filling model:
#         Ctrl_BD: SoilGrids BD + lu_start + MAT + MAP + study random effect
#         Aban_BD: SoilGrids BD + log(AGE) + study random effect
#       Fixed effects only are applied for prediction (no study random effect),
#       consistent with out-of-sample use.
#
#   Ctrl_BD_source / Aban_BD_source
#       "original_study" or "model_soilgrids_climate"
#
#   bd_gap_filled  (logical, derived from the two _source columns)
#       TRUE when either Ctrl_BD_source OR Aban_BD_source is
#       "model_soilgrids_climate".  The single column most convenient for
#       the BD sensitivity filter described above.
# =============================================================================


# ── Packages ──────────────────────────────────────────────────────────────────
# ipak() installs and loads packages in one call
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(c("tidyverse", "sf", "stars",   # spatial tools used for covariate extraction
       "openxlsx",                    # read raw Excel workbook
       "cruts",                       # convert CRU netCDF climate data to rasters
       "raster",                      # raster::extract() for BD and MAT/MAP
       "lme4"))                       # mixed models for BD gap-filling


# =============================================================================
# SECTION 1 — READ AND RENAME RAW DATA
# =============================================================================

# The raw workbook contains one observation per row, where each row represents
# a single depth interval from a paired (reference / restoration) comparison.
filename <- "./data_inputs/Global_data_cleaned_region_260311.xlsx"
df <- read.xlsx(filename, sheet = 1, startRow = 1) %>% as.data.frame()

# ── Step 1: rename the four duplicate SD/SE column pairs by position ──────────
# The spreadsheet contains four unlabelled SD and SE columns, one pair for each
# of: Ctrl carbon (SOC/SOM/stock), Ctrl BD, Aban carbon, Aban BD.  They appear
# in that order immediately after the corresponding value columns.  openxlsx
# preserves the duplicate names verbatim (all four are literally "SD" / "SE"),
# so dplyr::rename() cannot be used here — it refuses to operate on duplicate
# column names.  We assign by position instead, which is safe and unambiguous.
#
# Column positions confirmed from colnames(df) output (37 columns total):
#   13 / 14  → Ctrl_C_SD  / Ctrl_C_SE   (Ctrl carbon uncertainty)
#   18 / 19  → Ctrl_BD_SD / Ctrl_BD_SE  (Ctrl bulk density uncertainty)
#   23 / 24  → Aban_C_SD  / Aban_C_SE   (Aban carbon uncertainty)
#   32 / 33  → Aban_BD_SD / Aban_BD_SE  (Aban bulk density uncertainty)
colnames(df)[c(13, 14)] <- c("Ctrl_C_SD",  "Ctrl_C_SE")
colnames(df)[c(18, 19)] <- c("Ctrl_BD_SD", "Ctrl_BD_SE")
colnames(df)[c(23, 24)] <- c("Aban_C_SD",  "Aban_C_SE")
colnames(df)[c(32, 33)] <- c("Aban_BD_SD", "Aban_BD_SE")

# ── Step 2: rename the verbose header columns ─────────────────────────────────
# The exact source names (with backtick-quoted specials) reflect how openxlsx
# preserves the original Excel header text.
df <- df %>%
  dplyr::rename(
    MAT        = `MAT.(°C)`,
    MAP        = `MAP.(mm)`,
    LUC        = `LUC.(C=cropland,.P=pasture,.G=grassland,.S=shrubland,.F=forest)`,
    depth      = `Depth.(cm)`,
    MGMT       = `MGMT.(N=natural,.A=assisted,.N+A=natural+assisted,.G=grazed)`,
    Ctrl_stock = `Ctrl_stock.(Mg/ha)`,
    Aban_stock = `Aban_stock.(Mg/ha)`
  )

# ── Guard: confirm all expected columns exist after renaming ──────────────────
# If the Excel header text or column order changes in a future data update,
# this will error immediately with a clear message rather than failing silently
# much later in the script.
required_cols <- c("MAT", "MAP", "LUC", "depth", "MGMT",
                   "Ctrl_stock", "Aban_stock",
                   "Ctrl_SOC", "Ctrl_SOM", "Ctrl_BD",
                   "Aban_SOC",  "Aban_SOM",  "Aban_BD",
                   "AGE1", "LAT", "LONG", "Units",
                   "Ctrl_C_SD",  "Ctrl_C_SE",
                   "Ctrl_BD_SD", "Ctrl_BD_SE",
                   "Aban_C_SD",  "Aban_C_SE",
                   "Aban_BD_SD", "Aban_BD_SE",
                   "Title", "Authors", "Year", "Region")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("After renaming, the following expected columns are missing:\n  ",
       paste(missing_cols, collapse = ", "),
       "\nCheck that the Excel header text and column order match the ",
       "rename() calls above.\nRun glimpse(df) immediately after read.xlsx() ",
       "to inspect the raw names.")
}


# =============================================================================
# SECTION 2 — REMOVE MANUALLY IDENTIFIED INVALID STUDIES
# =============================================================================
# Each exclusion is documented with the reason supplied by the original data
# curators.

# Study excluded: SOC data are not properly matched between paired plots
df <- subset(df,
             Title != "Nitrogen cycle patterns during forest regrowth in an African Miombo woodland landscape. J. Geophys. Res. Biogeosciences 124, 1591–1603 (2019)." |
               is.na(Title))

# Study excluded: reference ("control") data is undisturbed prairie rather than
# an agricultural land use
df <- subset(df,
             Title != "Soil carbon sequestration across a chronosequence of tallgrass prairie restorations in the Ozark Highlands region of northwest Arkansas" |
               is.na(Title))

# Study excluded: abandoned plots were harvested once a year, which constitutes
# an ongoing agricultural practice rather than true abandonment
df <- subset(df,
             Title != "Restoration of ecosystem carbon stocks following exclosure establishment in communal grazing lands in Tigray, Ethiopia. Soil Sci. Soc. Am. J. 75, 246–256 (2011)." |
               is.na(Title))


# =============================================================================
# SECTION 3 — CORRECT WRONGLY REPORTED DEPTH INTERVALS
# =============================================================================
# A small number of studies reported cumulative depth stocks (e.g., 0–60 cm)
# rather than incremental interval stocks (e.g., 20–60 cm).  The corrections
# below convert cumulative values to incremental ones by subtraction.

# Study: "Land-use change and management effects on carbon sequestration in soils
# of Russia's South Taiga zone"
# Reported as 0–20 and 0–60; corrected to 0–20 and 20–60.
title1 <- "Land‐use change and management effects on carbon sequestration in soils of Russia's South Taiga zone"
tdf1 <- subset(df, Title == title1)
tdf1[c(2,4,6,8,10,12), "Ctrl_stock"] <-
  as.double(tdf1[c(2,4,6,8,10,12), "Ctrl_stock"]) - as.double(tdf1[c(1,3,5,7,9,11), "Ctrl_stock"])
tdf1[c(2,4,6,8,10,12), "Aban_stock"] <-
  as.double(tdf1[c(2,4,6,8,10,12), "Aban_stock"]) - as.double(tdf1[c(1,3,5,7,9,11), "Aban_stock"])
tdf1[c(2,4,6,8,10,12), "depth"] <- "20-60"
tdf1[, "Ctrl_C_SE"] <- NA
tdf1[, "Aban_C_SE"] <- NA
df[which(df$Title == title1), ] <- tdf1

# Study: "Changes in Physical Properties and Carbon Stocks of Gray Forest Soils..."
# Reported as 0–10, 10–20, 0–60; corrected to 0–10, 10–20, 20–60.
title2 <- "Changes in Physical Properties and Carbon Stocks of Gray Forest Soils in the Southern Part of Moscow Region during Postagrogenic Evolution"
tdf2 <- subset(df, Title == title2)
tdf2[c(3,6,9,12), "Ctrl_stock"] <-
  as.double(tdf2[c(3,6,9,12), "Ctrl_stock"]) -
  as.double(tdf2[c(1,4,7,10), "Ctrl_stock"]) -
  as.double(tdf2[c(2,5,8,11), "Ctrl_stock"])
tdf2[c(3,6,9,12), "Aban_stock"] <-
  as.double(tdf2[c(3,6,9,12), "Aban_stock"]) -
  as.double(tdf2[c(1,4,7,10), "Aban_stock"]) -
  as.double(tdf2[c(2,5,8,11), "Aban_stock"])
tdf2[c(2,5,8,11), "depth"] <- "10-20"
tdf2[c(3,6,9,12), "depth"] <- "20-60"
df[which(df$Title == title2), ] <- tdf2

# Study: "Self-restoration of post-agrogenic chernozems of Russia..."
# Rows 1–4 were reported as 0–50; rows 5–8 as 0–20.  Corrected rows 1–4 to
# 20–50 by subtracting the 0–20 values.
title3 <- "Self-restoration of post-agrogenic chernozems of Russia: Soil development, carbon stocks, and dynamics of carbon pools"
tdf3 <- subset(df, Title == title3)
tdf3[c(1,2,3,4), "Ctrl_stock"] <-
  as.double(tdf3[c(1,2,3,4), "Ctrl_stock"]) - as.double(tdf3[c(5,6,7,8), "Ctrl_stock"])
tdf3[c(1,2,3,4), "Aban_stock"] <-
  as.double(tdf3[c(1,2,3,4), "Aban_stock"]) - as.double(tdf3[c(5,6,7,8), "Aban_stock"])
tdf3[c(1,2,3,4), "depth"] <- "20-50"
df[which(df$Title == title3), ] <- tdf3

# Assign missing units for one study
title5 <- "Influence of Abandoning Agricultural Land Use on Hydrophysical Properties of Sandy Soil"
df[which(df$Title == title5), ]$Units <- "%"

# Fix a decimal-place error in stock value (1831 → 18.31 Mg/ha)
title6 <- "Soil Organic Carbon Beneath Croplands and Re-established Grasslands in the North Dakota Prairie Pothole Region"
df[which(df$Title == title6 & df$Ctrl_stock == 1831), ]$Ctrl_stock <- 18.31


# =============================================================================
# SECTION 4 — CLEAN NUMERIC COLUMNS
# =============================================================================

# Some SOC/SOM/BD cells contain "value ± SD" strings.  The helper below
# separates the central value and moves the uncertainty into its dedicated
# column.
format_sdse <- function(df, colname, sdsename) {
  col      <- df[, colname]
  col_sdse <- df[, sdsename]
  fobs     <- grep("±", col)
  temp     <- str_split(col, "±")
  col_sdse[fobs] <- unlist(sapply(temp[fobs], "[", 2))
  col      <- unlist(sapply(temp, "[", 1))
  col_sdse <- gsub("sd", "", gsub("se", "", col_sdse))
  df[, colname]  <- col
  df[, sdsename] <- col_sdse
  return(df)
}

df <- format_sdse(df, "Ctrl_SOM", "Ctrl_C_SD")
df <- format_sdse(df, "Ctrl_BD",  "Ctrl_BD_SD")
df <- format_sdse(df, "Aban_SOM", "Aban_C_SD")
df <- format_sdse(df, "Aban_BD",  "Aban_BD_SD")

# Strip comparison operators from age and depth fields
df$AGE1  <- gsub("[><= ]", "", df$AGE1)
df$AGE2  <- gsub("[><= ]", "", df$AGE2)
df$depth <- gsub("[><= ]", "", df$depth)

# Strip trailing letter suffixes from publication years (e.g., "2003a" → "2003")
df$Year <- gsub("[ab]", "", df$Year)

# Handle a handful of coordinates stored in exponential notation in the
# source spreadsheet, which would be silently corrupted by as.numeric()
df$LAT <- gsub("7.0000000000000007E-2",  "0.07",      df$LAT)
df$LAT <- gsub("2.4167000000000001E-2",  "0.024167",  df$LAT)
df$LAT <- gsub("9.4167000000000001E-2",  "0.094167",  df$LAT)

if (length(df$LAT[grepl("E",  df$LAT)])  > 0 |
    length(df$LONG[grepl("E", df$LONG)]) > 0) {
  stop("Unexpected exponential notation in LAT or LONG — check raw data.")
}

# Some coordinates and ages were reported as ranges (e.g., "34.5-36.2").
# This helper parses them into separate v1 / v2 columns; the midpoint is
# used downstream.
split_range_num <- function(x) {
  x <- str_squish(str_replace_all(x, "\u2013", "-"))
  x[is.na(x)] <- NA_character_
  parts <- str_split(x, "(?<=\\d)\\s*-\\s*(?=-?\\d)", n = 2)
  get1 <- function(p) if (length(p) >= 1 && nzchar(p[1])) p[1] else NA_character_
  get2 <- function(p) if (length(p) >= 2 && nzchar(p[2])) p[2] else NA_character_
  tibble(
    v1 = suppressWarnings(as.numeric(map_chr(parts, get1))),
    v2 = suppressWarnings(as.numeric(map_chr(parts, get2)))
  )
}

df <- df %>%
  mutate(
    LONG_clean = str_squish(str_replace_all(LONG, "\u2013", "-")),
    LAT_clean  = str_squish(str_replace_all(LAT,  "\u2013", "-"))
  ) %>%
  bind_cols(
    split_range_num(.$LONG_clean) %>% rename(LONG1 = v1, LONG2 = v2),
    split_range_num(.$LAT_clean)  %>% rename(LAT1  = v1, LAT2  = v2)
  )

# Convert all measurement columns to numeric (non-numeric strings → NA)
cols_num <- c("Year", "Ctrl_SOM", "Ctrl_SOC", "Ctrl_stock",
              "Aban_SOM", "Aban_SOC", "Aban_stock", "Ctrl_BD", "Aban_BD",
              "MAT", "MAP",
              "Ctrl_C_SD",  "Ctrl_C_SE",
              "Ctrl_BD_SD", "Ctrl_BD_SE",
              "Aban_C_SD",  "Aban_C_SE",
              "Aban_BD_SD", "Aban_BD_SE",
              "AGE1", "AGE2", "LAT1", "LAT2", "LONG1", "LONG2")
df[cols_num] <- sapply(df[cols_num], as.double)


# =============================================================================
# SECTION 5 — FORMAT LAND-USE CODES AND DEPTH INTERVALS
# =============================================================================

# Remove parenthetical notes from LUC strings (e.g., "C-F (plantation)")
df$LUC <- gsub("\\s*\\([^\\)]+\\)", "", df$LUC)

# Split LUC into source (Ctrl_LUC) and destination (Aban_LUC) codes
df <- df %>%
  mutate(LUC_sep  = str_split(LUC, "-"),
         Ctrl_LUC = sapply(LUC_sep, "[", 1),
         Aban_LUC = sapply(LUC_sep, "[", 2))

# A small number of depth values were stored as Excel serial dates in some
# exports.  Values > 1000 are converted via date arithmetic to recover the
# intended numeric depth.
single_id <- which(
  as.double(df[(!grepl("-", df$depth) & !grepl("~", df$depth)), ]$depth) > 1000)
depth_as_date <- as.double(
  df[(!grepl("-", df$depth) & !grepl("~", df$depth)), ]$depth[single_id])
depth_as_date <- format(as.Date(depth_as_date, origin = "1899-12-30", "%Y-%m-%d"), "%m-%d")
df[(!grepl("-", df$depth) & !grepl("~", df$depth)), ]$depth[single_id] <- depth_as_date

# Separate depth strings into upper (depth1) and lower (depth2) bounds.
# Single values (no separator) are treated as the lower bound with 0 as upper.
df$depth <- str_split(df$depth, "-|~")
df <- df %>%
  mutate(
    depth1 = ifelse(lengths(depth) == 1, 0,
                    sapply(depth, "[", 1)),
    depth2 = ifelse(lengths(depth) == 1,
                    sapply(depth, "[", 1),
                    sapply(depth, "[", 2))
  )

# Any unparseable depth values become NA and are dropped with a diagnostic note
idx_bad <- which(
  (!is.na(df$depth1) & is.na(suppressWarnings(as.double(as.character(df$depth1))))) |
    (!is.na(df$depth2) & is.na(suppressWarnings(as.double(as.character(df$depth2)))))
)
if (length(idx_bad) > 0) {
  message(length(idx_bad),
          " rows have unparseable depth values and will be treated as NA.")
}
df[, c("depth1", "depth2")] <- sapply(df[, c("depth1", "depth2")], as.double)

df <- df %>%
  mutate(
    D         = depth2 - depth1,            # interval thickness (cm)
    depth_avg = (depth2 + depth1) / 2       # midpoint depth (cm)
  )


# =============================================================================
# SECTION 6 — FILTER BY AGE RANGE AND AVERAGE RANGES
# =============================================================================

# Where age or coordinate ranges span more than a factor of 3 AND more than
# 10 years, the pair is too heterogeneous to represent a single time-point;
# these rows are removed.
df <- df %>%
  filter(!((AGE2 / AGE1 > 3 & AGE2 - AGE1 > 10)) | is.na(AGE2))

# Represent ranged values by their midpoints
df <- df %>%
  mutate(
    AGE  = ifelse(is.na(AGE2),  AGE1,              (AGE1  + AGE2)  / 2),
    LAT  = ifelse(is.na(LAT2),  LAT1,              (LAT1  + LAT2)  / 2),
    LONG = ifelse(is.na(LONG2), LONG1,             (LONG1 + LONG2) / 2)
  )

# Drop rows lacking spatial coordinates, depth, or age — these cannot be
# linked to spatial covariates or placed on a chronosequence
df <- df %>%
  filter(!is.na(LAT), !is.na(LONG), !is.na(depth_avg), !is.na(AGE))

# Assign sequential observation IDs (used to link back to spatial extractions)
df <- df %>% mutate(obs = row_number())


# =============================================================================
# SECTION 7 — ASSIGN SOIL PROFILE IDs
# =============================================================================
# A profile is defined by a unique combination of study, year, coordinates,
# land-use transition, and sampling age.  Multiple depth intervals from the
# same physical core share the same Profile_ID.

df <- df %>% mutate(Profile_ID = NA_integer_)
max_ID <- 0L

for (i in seq_len(nrow(distinct(df, Year, Title,
                                round(LAT, 9), round(LONG, 9),
                                LUC, AGE1)))) {
  temp_vals <- distinct(df, Year, Title,
                        round(LAT, 9), round(LONG, 9),
                        LUC, AGE1)[i, ]
  temp_df  <- subset(df,
                     Year == temp_vals$Year &
                       Title == temp_vals$Title &
                       round(LAT, 9) == temp_vals[, 3] &
                       round(LONG, 9) == temp_vals[, 4] &
                       LUC  == temp_vals$LUC &
                       AGE1 == temp_vals$AGE1)
  temp_obs <- temp_df$obs
  if (nrow(temp_df) == 0) next
  
  if (nrow(temp_df) == 1) {
    seq_ID <- max_ID + 1L
    df$Profile_ID[temp_obs] <- seq_ID
    max_ID <- max(seq_ID)
    
  } else if (nrow(temp_df) %% length(unique(temp_df$depth)) == 0 &&
             identical(temp_df$depth[1], temp_df$depth[2])) {
    nprofile  <- nrow(temp_df) / length(unique(temp_df$depth))
    ninterval <- length(unique(temp_df$depth))
    seq_ID    <- max_ID + rep(seq_len(nprofile), times = ninterval)
    df$Profile_ID[temp_obs] <- seq_ID
    max_ID <- max(seq_ID)
    
  } else if (!identical(temp_df$depth[1], temp_df$depth[2]) && nrow(temp_df) > 1) {
    temp_ID <- max_ID + 1L
    seq_ID  <- temp_ID
    for (j in 2:nrow(temp_df)) {
      if (temp_df$depth1[j] == 0 || temp_df$depth1[j] == temp_df$depth1[1]) {
        temp_ID <- temp_ID + 1L
      }
      seq_ID <- c(seq_ID, temp_ID)
    }
    df$Profile_ID[temp_obs] <- seq_ID
    max_ID <- max(seq_ID)
    
  } else {
    message("Profile ID assignment edge case — rows: ", paste(temp_obs, collapse = ", "))
    next
  }
}

# Manual override for one study whose profile structure is not captured by the
# general algorithm above
title_manual <- "Tree species and time since afforestation drive soil C and N mineralization on former cropland"
temp_df_m  <- subset(df, Title == title_manual)
temp_obs_m <- temp_df_m$obs
seq_ID_m   <- max_ID + c(rep(seq_len(8), times = 2), 8L + rep(seq_len(7), times = 2))
df$Profile_ID[temp_obs_m] <- seq_ID_m
max_ID <- max(seq_ID_m)


# =============================================================================
# SECTION 8 — FORMAT MEASUREMENT UNITS
# =============================================================================

# Harmonise BD unit labels
df[(df$Ctrl_BD <= 10  & df$Ctrl_BD >= 0.1 & !is.na(df$Ctrl_BD)), ]$Ctrl_BD.units <- "g/cm3"
df[(df$Aban_BD <= 10  & df$Aban_BD >= 0.1 & !is.na(df$Aban_BD)), ]$Aban_BD.units <- "g/cm3"
df[(df$Ctrl_BD >= 200 & !is.na(df$Ctrl_BD)), ]$Ctrl_BD.units <- "mg/cm3"
df[(df$Aban_BD >= 200 & !is.na(df$Aban_BD)), ]$Aban_BD.units <- "mg/cm3"
df[is.na(df$Ctrl_BD), ]$Ctrl_BD.units <- NA
df[is.na(df$Aban_BD), ]$Aban_BD.units <- NA

# Harmonise SOC/SOM/stock unit labels
df[df$Units %in% c("g/kg","mg/g","g kg\u20131","g kg-1","mg g-1"), ]$Units <- "g/kg"
df[df$Units %in% c("Mg ha-1","Mg ha\u20131","Mg/ha","t/ha","tons/ha","t ha-1"), ]$Units <- "Mg/ha"
df[df$Units %in% c("kg/m²","kg/m2","kg/m3","mg/cm3","kg/(m²*15cm)","kg m\u20132","kg C m-2","kg m-2","kg m-3"), ]$Units <- "kg/m2"
df[df$Units %in% c("g/m²","g C m-2","g m-2"), ]$Units <- "g/m2"

# Stocks reported in the SOM column (wrong_unit profiles) are moved to the
# stock column after applying the SOM→SOC conversion factor
converter      <- 0.5
wrong_unit_obs <- which(
  df$Profile_ID %in% {
    wrong_unit <- c()
    for (k in seq_len(nrow(distinct(df, Profile_ID)))) {
      kprofile   <- distinct(df, Profile_ID)[k, ]
      profile_df <- subset(df, Profile_ID == kprofile) %>%
        .[order(.$depth2), ]
      measurements <- cbind(!is.na(profile_df[, c("Ctrl_SOC","Ctrl_SOM","Ctrl_stock",
                                                  "Aban_SOC","Aban_SOM","Aban_stock")]))
      if (TRUE %in% measurements[, c("Ctrl_SOC","Ctrl_SOM","Aban_SOC","Aban_SOM")] &&
          TRUE %in% !(profile_df$Units %in% c("g/kg","mg/kg","%"))) {
        wrong_unit <- c(wrong_unit, kprofile)
      }
      if (TRUE %in% measurements[, c("Ctrl_stock","Aban_stock")] &&
          TRUE %in% !(profile_df$Units %in% c("Mg/ha","kg/m2","kg/ha","g/m2"))) {
        wrong_unit <- c(wrong_unit, kprofile)
      }
    }
    wrong_unit
  } &
    rowSums(!is.na(df[, c("Ctrl_SOM","Aban_SOM")])) >= 1
)
mmsmt_list <- df[wrong_unit_obs, c("Ctrl_SOM","Aban_SOM")] * converter
df[wrong_unit_obs, c("Ctrl_stock","Aban_stock")] <- mmsmt_list
df[wrong_unit_obs, c("Ctrl_SOM","Aban_SOM")]     <- NA

# Final removal of rows lacking coordinates, depth, or age
df <- df %>%
  filter(!is.na(LAT), !is.na(LONG), !is.na(depth_avg), !is.na(AGE))
df <- df %>% mutate(obs = row_number())   # refresh observation IDs


# =============================================================================
# SECTION 9 — EXTRACT CLIMATE COVARIATES (MAT, MAP) FROM CRU
# =============================================================================
# Original study-reported values are preferred; raster values fill gaps.
# A _source column records the origin of each value.

# Build a SpatialPointsDataFrame for raster::extract()
points <- dplyr::select(df, obs, LAT, LONG)
coordinates(points) <- c("LONG", "LAT")

# MAT — CRU TS 4.06 mean monthly temperature, averaged 2012–2020
temp_nc <- cruts2raster("./data_inputs/raster/cru_ts4.06.2011.2020.tmp.dat.nc",
                        timeRange = c("2012-01-01", "2020-12-31"))
mat_raster <- mean(temp_nc)
mat_extracted <- data.frame(obs = points$obs,
                            mat = raster::extract(mat_raster, points))
df$MAT_source <- ifelse(is.na(df$MAT), "cru_ts", "original_study")
df$MAT        <- ifelse(is.na(df$MAT), mat_extracted$mat, df$MAT)
rm(temp_nc, mat_raster, mat_extracted)

# MAP — CRU TS 4.06 total monthly precipitation, summed then averaged to
# annual (divide by 9 for the 9-year 2012–2020 window)
prec_nc  <- cruts2raster("./data_inputs/raster/cru_ts4.06.2011.2020.pre.dat.nc",
                         timeRange = c("2012-01-01", "2020-12-31"))
map_raster <- sum(prec_nc)
map_extracted <- data.frame(obs = points$obs,
                            map = raster::extract(map_raster, points))
df$MAP_source <- ifelse(is.na(df$MAP), "cru_ts", "original_study")
df$MAP        <- ifelse(is.na(df$MAP), map_extracted$map / 9, df$MAP)
rm(prec_nc, map_raster, map_extracted)

message("MAT missing after CRU fill: ", sum(is.na(df$MAT)))
message("MAP missing after CRU fill: ", sum(is.na(df$MAP)))


# =============================================================================
# SECTION 10 — EXTRACT ALL SPATIAL COVARIATES
# =============================================================================
# lu_start is parsed first because it is needed as a predictor in the BD
# gap-filling models.  All other spatially-varying covariates (BD, clay,
# nitrogen, irrigation) are then extracted so that the complete covariate set
# is assembled on df before any modelling takes place.


# ── 10a. Parse initial land use ───────────────────────────────────────────────
# LUC_sep, Ctrl_LUC, and Aban_LUC were already parsed in Section 5.
# lu_start is the additional column needed here as a BD gap-filling predictor.
df <- df %>% mutate(lu_start = substr(LUC, 1, 1))


# ── 10b. Correct BD units and preserve raw published values ───────────────────
# A small number of published BD values were recorded in mg/cm³ (× 1000 too
# large).  Correct to g/cm³ before snapshotting the raw values.
df$Aban_BD[which(df$Aban_BD.units == "mg/cm3")] <-
  df$Aban_BD[which(df$Aban_BD.units == "mg/cm3")] / 1000
df$Ctrl_BD[which(df$Ctrl_BD.units == "mg/cm3")] <-
  df$Ctrl_BD[which(df$Ctrl_BD.units == "mg/cm3")] / 1000
df$Aban_BD.units[df$Aban_BD.units == "mg/cm3"] <- "g/cm3"
df$Ctrl_BD.units[df$Ctrl_BD.units == "mg/cm3"] <- "g/cm3"

# Snapshot unmodified published BD values; these are never altered after this
# point and serve as the validation targets in the gap-filling CV below.
df$Ctrl_BD_raw <- df$Ctrl_BD
df$Aban_BD_raw <- df$Aban_BD


# ── 10c. Extract SoilGrids BD at each observation's depth ────────────────────
# SoilGrids provides BD (bdod) at five standard depth intervals.  The layer
# whose midpoint is closest to depth_avg is selected; at exact boundaries the
# average of the two flanking layers is used.  Where the exact pixel is nodata
# (ocean, urban), the mean within an expanding buffer (1 → 5 → 25 km) is used.

extract_bd_buffered <- function(raster_obj, points_sp, buf_km = c(1, 5, 25)) {
  vals        <- raster::extract(raster_obj, points_sp)
  missing_idx <- which(is.na(vals))
  if (length(missing_idx) == 0) return(vals)
  pts_sf     <- st_as_sf(as.data.frame(points_sp),
                         coords = c("LONG", "LAT"),
                         crs    = "+proj=longlat +datum=WGS84")
  rast_stars <- st_as_stars(raster_obj)
  for (buf in buf_km) {
    still_missing <- missing_idx[is.na(vals[missing_idx])]
    if (length(still_missing) == 0) break
    for (i in still_missing) {
      buf_geom <- st_buffer(st_geometry(pts_sf[i, ]), dist = buf * 1000)
      vals[i]  <- tryCatch(
        st_extract(rast_stars, buf_geom, FUN = mean, na.rm = TRUE) %>% pull(1),
        error = function(e) NA_real_)
    }
  }
  if (any(is.na(vals)))
    warning(sum(is.na(vals)), " point(s) still have NA in BD extraction after ",
            max(buf_km), " km buffer — these will remain NA and be gap-filled ",
            "by the model in Section 11.")
  return(vals)
}

bd_rasters <- list(
  "0_5"    = raster("./data_inputs/raster/via_gee/bdod_0_5cm_mean.tif"),
  "5_15"   = raster("./data_inputs/raster/via_gee/bdod_5_15cm_mean.tif"),
  "15_30"  = raster("./data_inputs/raster/via_gee/bdod_15_30cm_mean.tif"),
  "30_60"  = raster("./data_inputs/raster/via_gee/bdod_30_60cm_mean.tif"),
  "60_100" = raster("./data_inputs/raster/via_gee/bdod_60_100cm_mean.tif")
)

d       <- df$depth_avg
bd_vals <- rep(NA_real_, nrow(df))

# Interior intervals
bd_vals[d <  5]          <- extract_bd_buffered(bd_rasters[["0_5"]],   points[d <  5,  ])
bd_vals[d >  5 & d < 15] <- extract_bd_buffered(bd_rasters[["5_15"]],  points[d >  5 & d < 15, ])
bd_vals[d > 15 & d < 30] <- extract_bd_buffered(bd_rasters[["15_30"]], points[d > 15 & d < 30, ])
bd_vals[d > 30 & d < 60] <- extract_bd_buffered(bd_rasters[["30_60"]], points[d > 30 & d < 60, ])
bd_vals[d >= 60]          <- extract_bd_buffered(bd_rasters[["60_100"]],points[d >= 60, ])

# Boundary depths: interpolate between flanking layers
if (any(d == 5))
  bd_vals[d == 5]  <- (extract_bd_buffered(bd_rasters[["0_5"]],   points[d == 5,  ]) +
                         extract_bd_buffered(bd_rasters[["5_15"]],  points[d == 5,  ])) / 2
if (any(d == 15))
  bd_vals[d == 15] <- (extract_bd_buffered(bd_rasters[["5_15"]],  points[d == 15, ]) +
                         extract_bd_buffered(bd_rasters[["15_30"]], points[d == 15, ])) / 2
if (any(d == 30))
  bd_vals[d == 30] <- (extract_bd_buffered(bd_rasters[["15_30"]], points[d == 30, ]) +
                         extract_bd_buffered(bd_rasters[["30_60"]], points[d == 30, ])) / 2
if (any(d == 60))
  bd_vals[d == 60] <- (extract_bd_buffered(bd_rasters[["30_60"]],  points[d == 60, ]) +
                         extract_bd_buffered(bd_rasters[["60_100"]], points[d == 60, ])) / 2

# SoilGrids bdod is stored as cg/dm³ (= kg/m³); divide by 100 to get g/cm³
bd.df <- data.frame(obs = df$obs, depth_avg = df$depth_avg,
                    bd_SG = bd_vals / 100)
rm(bd_rasters, bd_vals, d)

# Scale check: bd_SG should be ~0.8–1.8 g/cm³
message("bd_SG range after unit conversion (expect ~0.8–1.8 g/cm³): ",
        round(min(bd.df$bd_SG, na.rm = TRUE), 3), " – ",
        round(max(bd.df$bd_SG, na.rm = TRUE), 3))
message("Ctrl_BD_raw range (expect ~0.8–1.8 g/cm³): ",
        round(min(df$Ctrl_BD_raw, na.rm = TRUE), 3), " – ",
        round(max(df$Ctrl_BD_raw, na.rm = TRUE), 3))

# Attach bd_SG to df now so it is available as a covariate throughout
df <- df %>% left_join(bd.df %>% dplyr::select(obs, bd_SG), by = "obs")


# ── 10d. Extract clay and nitrogen from SoilGrids ────────────────────────────
# Six depth layers each for clay (g/kg) and nitrogen (mg/kg); the layer whose
# midpoint brackets depth_avg is selected, with buffer-based gap-filling for
# nodata pixels (urban locations, coastal points).
# Sand and pH are omitted due to collinearity with clay in preliminary analyses.

projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
df_sf <- df %>% st_as_sf(coords = c("LONG", "LAT"), crs = projcrs)

extract_stars_buffered <- function(rast_path, points_sf, buf_km = c(1, 5, 25)) {
  rast  <- read_stars(rast_path)
  vals  <- st_extract(rast, points_sf) %>% pull(1) %>% st_drop_geometry()
  missing_idx <- which(is.na(vals))
  if (length(missing_idx) == 0) return(vals)
  for (buf in buf_km) {
    still_missing <- missing_idx[is.na(vals[missing_idx])]
    if (length(still_missing) == 0) break
    for (i in still_missing) {
      buf_geom <- st_buffer(st_geometry(points_sf[i, ]), dist = buf * 1000)
      vals[i]  <- tryCatch(
        st_extract(rast, buf_geom, FUN = mean, na.rm = TRUE) %>% pull(1),
        error = function(e) NA_real_)
    }
  }
  if (any(is.na(vals)))
    warning(sum(is.na(vals)), " point(s) still have NA after ", max(buf_km),
            " km buffer — these will be NA in the output.  Check coordinates ",
            "for these observations or increase buf_km.")
  return(vals)
}

depth_layers <- c("0_5cm", "5_15cm", "15_30cm", "30_60cm", "60_100cm", "100_200cm")

for (layer in depth_layers) {
  df_sf[[paste0("clay_",     layer, "_mean")]] <- extract_stars_buffered(
    paste0("./data_inputs/raster/via_gee/clay_",     layer, "_mean.tif"), df_sf)
  df_sf[[paste0("nitrogen_", layer, "_mean")]] <- extract_stars_buffered(
    paste0("./data_inputs/raster/via_gee/nitrogen_", layer, "_mean.tif"), df_sf)
}
gc()

clay_cols     <- paste0("clay_",     depth_layers, "_mean")
nitrogen_cols <- paste0("nitrogen_", depth_layers, "_mean")

interval_id <- cut(df_sf$depth_avg,
                   breaks = c(0, 5, 15, 30, 60, 100, 500),
                   labels = FALSE, include.lowest = TRUE)

# Drop geometry once upfront so row/column indexing in the loop is fast
df_sf_tbl <- df_sf %>% st_drop_geometry()

df_sf$clay <- as.numeric(sapply(seq_len(nrow(df_sf_tbl)), function(x)
  df_sf_tbl[x, clay_cols[interval_id[x]]] %>% as.numeric()))

df_sf$nitrogen <- as.numeric(sapply(seq_len(nrow(df_sf_tbl)), function(x)
  df_sf_tbl[x, nitrogen_cols[interval_id[x]]] %>% as.numeric()))

# ── 10e. Extract irrigation (GMIA) ───────────────────────────────────────────
# The GMIA raster has been pre-augmented by 0.1 to allow log-transformation.
df_sf$irrigation <- extract_stars_buffered(
  "./data_inputs/raster/gmia_plus_0.1.tif", df_sf)

# Add all spatial covariates to df, explicitly coercing to atomic numeric to
# prevent list-column propagation from sf subsetting operations above.
df <- df %>%
  bind_cols(
    df_sf %>%
      st_drop_geometry() %>%
      dplyr::select(clay, nitrogen, irrigation) %>%
      mutate(across(everything(), as.numeric))
  )
df$clay_source       <- "soilgrids_via_gee"
df$nitrogen_source   <- "soilgrids_via_gee"
df$irrigation_source <- "gmia_augmented_0.1"

rm(df_sf, df_sf_tbl)


# =============================================================================
# SECTION 11 — BULK DENSITY GAP-FILLING: MODEL SELECTION AND VALIDATION
# =============================================================================
# All covariates (MAT, MAP, AGE, depth, clay, nitrogen, irrigation, lu_start)
# are now available on df, so the full candidate set of gap-filling models can
# be evaluated before a final model is chosen and applied.
#
# Workflow:
#   11a. Flag implausible repeated BD values
#   11b. Build CV datasets (observation-level folds, fixed-effects-only
#        prediction — see note below on why)
#   11c. Evaluate all candidate specifications for ctrl and aban BD
#   11d. Repeatability analysis (upper-bound R²)
#   11e. Print results — review these to choose the final gap-filling model
#   11f. Fit the chosen final models and apply gap-filling
#
# Note on CV design: predictions use fixed effects only (re.form = NA), because
# the study random effect cannot be estimated for genuinely gap-filled rows
# (those where no BD was published).  This gives a conservative but realistic
# estimate of gap-filling accuracy.  The Pearson r between SoilGrids BD and
# field BD (~0.26–0.31) sets the practical ceiling for fixed-effects R².

set.seed(42)
K <- 10
my_transform <- function(x) log(x)   # log(AGE) used throughout

# ── 11a. Flag implausible repeated BD values ──────────────────────────────────
# Multi-depth profiles where every interval shares the same rounded BD are
# flagged as likely imputed by the original authors; treated as missing.
flag_repeated_bd <- function(bd_vec, profile_id_vec) {
  tibble(bd = bd_vec, Profile_ID = profile_id_vec) %>%
    mutate(rounded_bd = round(bd, 1)) %>%
    group_by(Profile_ID) %>%
    mutate(
      n_samples      = sum(!is.na(rounded_bd)),
      n_distinct_bd  = n_distinct(rounded_bd, na.rm = TRUE),
      is_problematic = (n_samples > 1 & n_distinct_bd == 1 & !is.na(rounded_bd))
    ) %>%
    ungroup() %>%
    pull(is_problematic)
}

ctrl_bd_problematic <- flag_repeated_bd(df$Ctrl_BD, df$Profile_ID)
aban_bd_problematic <- flag_repeated_bd(df$Aban_BD, df$Profile_ID)

df$Ctrl_BD_problematic <- ctrl_bd_problematic
df$Aban_BD_problematic <- aban_bd_problematic

message("Ctrl_BD rows flagged as problematic (repeated value): ",
        sum(ctrl_bd_problematic, na.rm = TRUE))
message("Aban_BD rows flagged as problematic (repeated value): ",
        sum(aban_bd_problematic, na.rm = TRUE))

# Correlation check: confirms bd_SG is on the same scale as published BD
ctrl_cor <- cor(df$Ctrl_BD_raw[!ctrl_bd_problematic & !is.na(df$Ctrl_BD_raw)],
                df$bd_SG[!ctrl_bd_problematic & !is.na(df$Ctrl_BD_raw)],
                use = "complete.obs")
aban_cor <- cor(df$Aban_BD_raw[!aban_bd_problematic & !is.na(df$Aban_BD_raw)],
                df$bd_SG[!aban_bd_problematic & !is.na(df$Aban_BD_raw)],
                use = "complete.obs")
message("Pearson r(Ctrl_BD_raw, bd_SG): ", round(ctrl_cor, 3))
message("Pearson r(Aban_BD_raw, bd_SG): ", round(aban_cor, 3))


# ── 11b. Build CV datasets ────────────────────────────────────────────────────
# Folds are assigned at the observation level, stratified within studies.
# Only rows with source == "original_study" are used as validation targets.

assign_folds <- function(titles, K) {
  titles_fac <- as.factor(titles)
  folds <- integer(length(titles_fac))
  for (lvl in levels(titles_fac)) {
    idx <- which(titles_fac == lvl)
    folds[idx] <- sample(rep(seq_len(K), length.out = length(idx)))
  }
  folds
}

ctrl_keep_cv <- !replace(ctrl_bd_problematic, is.na(ctrl_bd_problematic), FALSE)
aban_keep_cv <- !replace(aban_bd_problematic, is.na(aban_bd_problematic), FALSE)

ctrl_cv_df <- df[ctrl_keep_cv, ] %>%
  filter(!is.na(Ctrl_BD_raw)) %>%
  transmute(
    row_id    = obs,
    Ctrl_BD   = as.numeric(Ctrl_BD_raw),
    bd_SG     = as.numeric(bd_SG),
    AGE_trans = as.numeric(my_transform(AGE)),
    lu_start  = as.character(lu_start),
    MAT       = as.numeric(MAT),
    MAP       = as.numeric(MAP),
    depth_avg = as.numeric(depth_avg),
    clay      = as.numeric(clay),
    nitrogen  = as.numeric(nitrogen),
    irrigation = as.numeric(irrigation),
    Title     = as.character(Title)
  ) %>%
  filter(complete.cases(row_id, Ctrl_BD, bd_SG, lu_start, Title)) %>%
  mutate(Title = as.factor(Title),
         fold  = assign_folds(Title, K))

aban_cv_df <- df[aban_keep_cv, ] %>%
  filter(!is.na(Aban_BD_raw)) %>%
  transmute(
    row_id    = obs,
    Aban_BD   = as.numeric(Aban_BD_raw),
    bd_SG     = as.numeric(bd_SG),
    AGE_trans = as.numeric(my_transform(AGE)),
    lu_start  = as.character(lu_start),
    MAT       = as.numeric(MAT),
    MAP       = as.numeric(MAP),
    depth_avg = as.numeric(depth_avg),
    clay      = as.numeric(clay),
    nitrogen  = as.numeric(nitrogen),
    irrigation = as.numeric(irrigation),
    Title     = as.character(Title)
  ) %>%
  filter(complete.cases(row_id, Aban_BD, bd_SG, AGE_trans, lu_start, Title)) %>%
  mutate(Title = as.factor(Title),
         fold  = assign_folds(Title, K))

message("ctrl_cv_df rows: ", nrow(ctrl_cv_df),
        "  Pearson r(Ctrl_BD, bd_SG): ",
        round(cor(ctrl_cv_df$Ctrl_BD, ctrl_cv_df$bd_SG, use = "complete.obs"), 3))
message("aban_cv_df rows: ", nrow(aban_cv_df),
        "  Pearson r(Aban_BD, bd_SG): ",
        round(cor(aban_cv_df$Aban_BD, aban_cv_df$bd_SG, use = "complete.obs"), 3))


# ── 11c. Evaluate all candidate models ───────────────────────────────────────
# Candidates span from SoilGrids-only up to the full covariate set.
# lu_start is included as a candidate predictor for both ctrl and aban BD
# because cropland and pasture may have systematically different BD.

calc_metrics <- function(obs, pred) {
  if (length(obs) == 0 || length(pred) == 0)
    return(tibble(n = 0L, R2 = NA_real_, RMSE = NA_real_, bias = NA_real_))
  cc <- complete.cases(obs, pred)
  obs <- obs[cc]; pred <- pred[cc]; n <- length(obs)
  if (n < 3) return(tibble(n = n, R2 = NA, RMSE = NA, bias = NA))
  residuals <- obs - pred
  tibble(n = n,
         R2   = 1 - sum(residuals^2) / sum((obs - mean(obs))^2),
         RMSE = sqrt(mean(residuals^2)),
         bias = mean(residuals))
}

run_cv <- function(cv_df, formula, response_var) {
  imap_dfr(split(cv_df, cv_df$fold), function(test, k) {
    train <- cv_df %>% filter(fold != as.integer(k))
    vars_needed <- intersect(all.vars(formula), colnames(train))
    # Coerce to a plain data frame of atomic columns before complete.cases()
    # to avoid errors from any residual list-type columns
    train_sub <- as.data.frame(lapply(train[, vars_needed, drop = FALSE],
                                      function(x) if (is.list(x)) NA else x))
    test_sub  <- as.data.frame(lapply(test[,  vars_needed, drop = FALSE],
                                      function(x) if (is.list(x)) NA else x))
    train <- train[complete.cases(train_sub), ]
    test  <- test[ complete.cases(test_sub),  ]
    if (nrow(train) < 5 || nrow(test) == 0) return(tibble())
    mod_k <- tryCatch(
      suppressWarnings(lmer(formula, data = train, REML = FALSE)),
      error = function(e) NULL)
    if (is.null(mod_k)) return(tibble())
    tibble(obs  = test[[response_var]],
           pred = predict(mod_k, newdata = test, re.form = NA,
                          allow.new.levels = TRUE))
  })
}

candidate_models <- list(
  # ── Control BD ──────────────────────────────────────────────────────────────
  "M1: SoilGrids only [ctrl]" = list(
    type = "ctrl", response = "Ctrl_BD",
    formula = Ctrl_BD ~ bd_SG + (1 | Title)),
  "M2: + lu_start [ctrl]" = list(
    type = "ctrl", response = "Ctrl_BD",
    formula = Ctrl_BD ~ bd_SG + lu_start + (1 | Title)),
  "M3: + climate [ctrl]" = list(
    type = "ctrl", response = "Ctrl_BD",
    formula = Ctrl_BD ~ bd_SG + lu_start + MAT + MAP + (1 | Title)),
  "M4: + depth [ctrl]" = list(
    type = "ctrl", response = "Ctrl_BD",
    formula = Ctrl_BD ~ bd_SG + lu_start + MAT + MAP + depth_avg + (1 | Title)),
  "M5: + clay [ctrl]" = list(
    type = "ctrl", response = "Ctrl_BD",
    formula = Ctrl_BD ~ bd_SG + lu_start + MAT + MAP + depth_avg + clay + (1 | Title)),
  "M6: + nitrogen + irrigation [ctrl]" = list(
    type = "ctrl", response = "Ctrl_BD",
    formula = Ctrl_BD ~ bd_SG + lu_start + MAT + MAP + depth_avg + clay +
      nitrogen + irrigation + (1 | Title)),
  # ── Abandoned BD ────────────────────────────────────────────────────────────
  "M7: SoilGrids + log(AGE) [aban]" = list(
    type = "aban", response = "Aban_BD",
    formula = Aban_BD ~ bd_SG + AGE_trans + (1 | Title)),
  "M8: + lu_start [aban]" = list(
    type = "aban", response = "Aban_BD",
    formula = Aban_BD ~ bd_SG + AGE_trans + lu_start + (1 | Title)),
  "M9: + log(AGE)^2 [aban]" = list(
    type = "aban", response = "Aban_BD",
    formula = Aban_BD ~ bd_SG + AGE_trans + I(AGE_trans^2) + lu_start + (1 | Title)),
  "M10: + climate [aban]" = list(
    type = "aban", response = "Aban_BD",
    formula = Aban_BD ~ bd_SG + AGE_trans + lu_start + MAT + MAP + (1 | Title)),
  "M11: + depth [aban]" = list(
    type = "aban", response = "Aban_BD",
    formula = Aban_BD ~ bd_SG + AGE_trans + lu_start + MAT + MAP + depth_avg + (1 | Title)),
  "M12: + clay [aban]" = list(
    type = "aban", response = "Aban_BD",
    formula = Aban_BD ~ bd_SG + AGE_trans + lu_start + MAT + MAP + depth_avg +
      clay + (1 | Title)),
  "M13: + nitrogen + irrigation [aban]" = list(
    type = "aban", response = "Aban_BD",
    formula = Aban_BD ~ bd_SG + AGE_trans + lu_start + MAT + MAP + depth_avg +
      clay + nitrogen + irrigation + (1 | Title))
)

model_comparison_bd <- imap_dfr(candidate_models, function(spec, name) {
  cv_df  <- if (spec$type == "ctrl") ctrl_cv_df else aban_cv_df
  cv_out <- run_cv(cv_df, spec$formula, spec$response)
  bind_cols(tibble(model = name, BD_type = spec$type),
            calc_metrics(cv_out$obs, cv_out$pred))
})

message("\n── BD model comparison (Control BD) ─────────────────────")
print(model_comparison_bd %>% filter(BD_type == "ctrl") %>%
        dplyr::select(model, n, R2, RMSE, bias), n = Inf)
message("\n── BD model comparison (Abandoned BD) ───────────────────")
print(model_comparison_bd %>% filter(BD_type == "aban") %>%
        dplyr::select(model, n, R2, RMSE, bias), n = Inf)

# Optionally save:
# write.csv(model_comparison_bd, "./figures/bd_model_comparison.csv",
#           row.names = FALSE)


# ── 11d. Repeatability — upper-bound R² ──────────────────────────────────────
depth_bins   <- c(0, 10, 30, 60, Inf)
depth_labels <- c("0-10 cm", "10-30 cm", "30-60 cm", ">60 cm")

repeat_bd <- df %>%
  mutate(site_key  = paste(Title, round(LAT, 3), round(LONG, 3), sep = "_"),
         depth_bin = cut(depth_avg, breaks = depth_bins, labels = depth_labels,
                         include.lowest = TRUE))

icc_r2 <- function(repeat_df, bd_col) {
  d <- repeat_df %>%
    filter(!is.na(.data[[bd_col]])) %>%
    group_by(site_key, depth_bin) %>% filter(n() >= 2) %>% ungroup()
  if (nrow(d) == 0) return(NA_real_)
  mod <- lmer(reformulate(c("1", "(1|site_key)", "(1|depth_bin)"),
                          response = bd_col), data = d)
  vc <- as.data.frame(VarCorr(mod))
  (vc$vcov[vc$grp == "site_key"] + vc$vcov[vc$grp == "depth_bin"]) / sum(vc$vcov)
}

ctrl_repeat_R2 <- icc_r2(
  repeat_bd[ctrl_keep_cv & !is.na(repeat_bd$Ctrl_BD_raw), ], "Ctrl_BD_raw")
aban_repeat_R2 <- icc_r2(
  repeat_bd[aban_keep_cv & !is.na(repeat_bd$Aban_BD_raw), ], "Aban_BD_raw")

message("Repeatability R² — Control BD (upper bound): ",  round(ctrl_repeat_R2, 3))
message("Repeatability R² — Abandoned BD (upper bound): ", round(aban_repeat_R2, 3))

rm(ctrl_cv_df, aban_cv_df, repeat_bd, model_comparison_bd, run_cv, calc_metrics)


# ── 11e. Fit chosen gap-filling models and apply ──────────────────────────────
# Chosen specifications based on the CV comparison above:
#
#   ctrl: M3 — SoilGrids BD + lu_start + MAT + MAP
#     Adding climate (MAT + MAP) raises CV R² from ~0.03 to ~0.18.
#     Depth and soil texture predictors (M4–M6) add nothing further.
#
#   aban: M7 — SoilGrids BD + log(AGE)
#     The baseline AGE model performs well (R² = 0.174). lu_start and
#     climate predictors do not improve it, and in most cases hurt slightly.
#     M13 (+ nitrogen + irrigation) gives a marginal gain (0.191) at the
#     cost of two extra predictors that may be collinear; M7 is preferred
#     for parsimony.
#
# To change the specification, update the lmer() formula in each block below
# and the predict() call will automatically use the new fixed effects.

ctrl_mod_df <- df[ctrl_keep_cv, ] %>%
  filter(!is.na(Ctrl_BD_raw), !is.na(bd_SG), !is.na(lu_start),
         !is.na(MAT), !is.na(MAP)) %>%
  mutate(Title = as.factor(Title))

ctrl_bd_mod <- lmer(Ctrl_BD_raw ~ bd_SG + lu_start + MAT + MAP + (1 | Title),
                    data = ctrl_mod_df, REML = FALSE)

aban_mod_df <- df[aban_keep_cv, ] %>%
  filter(!is.na(Aban_BD_raw), !is.na(bd_SG), !is.na(AGE)) %>%
  mutate(AGE_trans = my_transform(AGE),
         Title     = as.factor(Title))

aban_bd_mod <- lmer(Aban_BD_raw ~ bd_SG + AGE_trans + (1 | Title),
                    data = aban_mod_df, REML = FALSE)

# Apply gap-filling using fixed effects only (no study random effect).
# Using predict() with re.form = NA is simpler and safer than manually
# extracting individual fixed-effect coefficients, especially if lu_start
# has more than two levels or the reference level varies.
needs_ctrl_fill <- is.na(df$Ctrl_BD) | ctrl_bd_problematic
needs_aban_fill <- is.na(df$Aban_BD) | aban_bd_problematic

df$Ctrl_BD_source <- ifelse(needs_ctrl_fill, "model_soilgrids_climate", "original_study")
ctrl_pred_df <- df %>%
  mutate(Title = as.factor(Title))
df$Ctrl_BD <- ifelse(
  needs_ctrl_fill,
  predict(ctrl_bd_mod, newdata = ctrl_pred_df, re.form = NA, allow.new.levels = TRUE),
  df$Ctrl_BD
)

df$Aban_BD_source <- ifelse(needs_aban_fill, "model_soilgrids_climate", "original_study")
aban_pred_df <- df %>%
  mutate(AGE_trans = my_transform(AGE),
         Title     = as.factor(Title))
df$Aban_BD <- ifelse(
  needs_aban_fill,
  predict(aban_bd_mod, newdata = aban_pred_df, re.form = NA, allow.new.levels = TRUE),
  df$Aban_BD
)

message("Ctrl_BD NAs remaining after gap-fill: ", sum(is.na(df$Ctrl_BD)))
message("Aban_BD NAs remaining after gap-fill: ", sum(is.na(df$Aban_BD)))

# Save model objects so diagnostics can be re-run without rerunning this script
save(bd.df, ctrl_bd_mod, aban_bd_mod,
     ctrl_bd_problematic, aban_bd_problematic, my_transform,
     file = "./data_inputs/bd_validation_inputs.rda")


# =============================================================================
# SECTION 12 — COMPUTE SOC STOCKS
# =============================================================================
# All covariates and final BD values are now on df.  Unit harmonisation and
# concentration → stock conversion are performed here.
# A som_converted flag records which rows used the SOM → SOC conversion factor.

var_columns <- c("Ctrl_C_SD", "Ctrl_C_SE", "Aban_C_SD", "Aban_C_SE")

# Convert stock units to Mg/ha
df[df$Units == "kg/m2", c("Ctrl_stock","Aban_stock", var_columns)] <-
  10    * df[df$Units == "kg/m2",  c("Ctrl_stock","Aban_stock", var_columns)]
df[df$Units == "g/m2",  c("Ctrl_stock","Aban_stock", var_columns)] <-
  0.01  * df[df$Units == "g/m2",   c("Ctrl_stock","Aban_stock", var_columns)]
df[df$Units == "kg/ha", c("Ctrl_stock","Aban_stock", var_columns)] <-
  0.001 * df[df$Units == "kg/ha",  c("Ctrl_stock","Aban_stock", var_columns)]

# Convert SOC/SOM concentrations to g/kg
df[df$Units == "%",     c("Ctrl_SOC","Aban_SOC","Ctrl_SOM","Aban_SOM", var_columns)] <-
  10    * df[df$Units == "%",     c("Ctrl_SOC","Aban_SOC","Ctrl_SOM","Aban_SOM", var_columns)]
df[df$Units == "mg/kg", c("Ctrl_SOC","Aban_SOC","Ctrl_SOM","Aban_SOM", var_columns)] <-
  0.001 * df[df$Units == "mg/kg", c("Ctrl_SOC","Aban_SOC","Ctrl_SOM","Aban_SOM", var_columns)]

# Separate by measurement type and flag SOM rows before conversion
df_SOC   <- subset(df, df$Ctrl_SOC  != "na" | df$Aban_SOC  != "na")
df_SOM   <- subset(df, df$Ctrl_SOM  != "na" | df$Aban_SOM  != "na")
df_stock <- subset(df, df$Ctrl_stock != "na" | df$Aban_stock != "na")

df_SOC$som_converted   <- FALSE
df_SOM$som_converted   <- TRUE
df_stock$som_converted <- FALSE

# SOM → SOC (van Bemmelen factor 0.5)
df_SOM[, c("Ctrl_SOC","Aban_SOC", var_columns)] <-
  df_SOM[, c("Ctrl_SOM","Aban_SOM", var_columns)] * converter
df_SOC2 <- rbind(df_SOC, df_SOM)

# SOC concentration × BD × depth → stock (Mg/ha)
# Formula: SOC (g/kg) × BD (g/cm³) × D (cm) / 10 = Mg C ha⁻¹
df_SOC2 <- df_SOC2 %>%
  mutate(
    Ctrl_stock_computed    = Ctrl_SOC * Ctrl_BD * D / 10,
    Aban_stock_computed    = Aban_SOC * Aban_BD * D / 10,
    Ctrl_stock_computed_SD = Ctrl_C_SD * Ctrl_BD * D / 10,
    Ctrl_stock_computed_SE = Ctrl_C_SE * Ctrl_BD * D / 10,
    Aban_stock_computed_SD = Aban_C_SD * Aban_BD * D / 10,
    Aban_stock_computed_SE = Aban_C_SE * Aban_BD * D / 10
  )

# Parallel stock calculation using original reported BD (before gap-filling).
# BD_raw retains the published value for all rows, including those flagged as
# problematic (repeated values).  Where the original study did not report BD
# at all, BD_raw is NA and these columns will also be NA — those rows are
# dropped in 03_bayes_analysis.R when the sensitivity run is assembled.
df_SOC2 <- df_SOC2 %>%
  mutate(
    Ctrl_stock_computed_orig_bd = Ctrl_SOC * Ctrl_BD_raw * D / 10,
    Aban_stock_computed_orig_bd = Aban_SOC * Aban_BD_raw * D / 10
  )

# Records already in stock units pass through directly.
# _orig_bd columns are set to NA for stock-unit rows: no BD multiplication
# was performed, so there is no alternative BD to substitute.
df_stock <- df_stock %>%
  mutate(
    Ctrl_stock_computed    = Ctrl_stock,
    Aban_stock_computed    = Aban_stock,
    Ctrl_stock_computed_orig_bd = NA_real_,
    Aban_stock_computed_orig_bd = NA_real_,
    Ctrl_stock_computed_SD = Ctrl_C_SD,
    Ctrl_stock_computed_SE = Ctrl_C_SE,
    Aban_stock_computed_SD = Aban_C_SD,
    Aban_stock_computed_SE = Aban_C_SE
  )

df_computed_stock <- rbind(df_SOC2, df_stock) %>%
  mutate(
    Stock_accrual     = Aban_stock_computed - Ctrl_stock_computed,
    Perc_change_stock = Stock_accrual / Ctrl_stock_computed,
    Accrual_rate      = Stock_accrual / AGE
  )


# =============================================================================
# SECTION 13 — FINALISE c_dat
# =============================================================================
# Derive analysis columns, apply recoding, add sensitivity filter flags,
# and drop rows that cannot be used in the Bayesian model.

c_dat <- df_computed_stock

# Coordinate ID (integer key for unique lat/long pairs)
c_dat <- c_dat %>%
  mutate(coord_id = paste(LAT, LONG) %>% as.factor() %>% as.numeric())

# Carbon density (Mg ha⁻¹ cm⁻¹): stock divided by interval thickness
c_dat <- c_dat %>%
  mutate(Aban_c_dens_av = Aban_stock_computed / D,
         Ctrl_c_dens_av = Ctrl_stock_computed / D,
         # Parallel density columns using original reported BD (sensitivity (d)).
         # NA where BD_raw was not reported or where the row uses stock units.
         Aban_c_dens_av_orig_bd = Aban_stock_computed_orig_bd / D,
         Ctrl_c_dens_av_orig_bd = Ctrl_stock_computed_orig_bd / D)

# ── BD gap-filling sensitivity flag ──────────────────────────────────────────
c_dat <- c_dat %>%
  mutate(bd_gap_filled = (Ctrl_BD_source == "model_soilgrids_climate") |
           (Aban_BD_source == "model_soilgrids_climate"))
message("Rows with at least one gap-filled BD value: ",
        sum(c_dat$bd_gap_filled, na.rm = TRUE), " of ", nrow(c_dat))

# Region label correction (Lithuania is not itself a region)
c_dat <- c_dat %>%
  mutate(Region = recode_factor(Region,
                                "Lithuania" = "Europe (west+east but not Russia)"))

# lu_start was already parsed in Section 10a.  Parse lu_end here and apply
# lookup-table gap-fill and binary collapse.
c_dat <- c_dat %>%
  mutate(lu_end = word(LUC, 2, sep = "-"))

missing_lu <- read.csv("./data_inputs/missing_lu_end.csv") %>%
  filter(!is.na(lu_end))
for (i in which(c_dat$lu_end == "?")) {
  c_dat$lu_end[i] <- ifelse(
    c_dat$Title[i] %in% missing_lu$Title,
    filter(missing_lu, Title == c_dat$Title[i]) %>% pull(lu_end),
    NA)
}

c_dat <- c_dat %>%
  filter(!is.na(lu_end)) %>%
  mutate(lu_end = ifelse(lu_end == "F", "F", "O"))

# ── Management sensitivity flags ──────────────────────────────────────────────
c_dat <- c_dat %>%
  mutate(MGMT_raw          = MGMT,
         MGMT_is_uncertain = (MGMT == "?"))
message("Rows with uncertain MGMT ('?'): ",
        sum(c_dat$MGMT_is_uncertain, na.rm = TRUE), " of ", nrow(c_dat))

# Remove grazed sites; recode remaining MGMT to "A" (assisted) / "N" (natural)
c_dat <- c_dat %>%
  filter(MGMT != "G") %>%
  mutate(MGMT = ifelse(grepl("A", MGMT), "A", "N"))

# Drop rows with missing carbon density
c_dat <- c_dat %>%
  filter(!is.na(Aban_c_dens_av), !is.na(Ctrl_c_dens_av))


# =============================================================================
# SECTION 14 — CENTERING / SCALING PARAMETERS
# =============================================================================
# Used by 01_lc_priors.R and 02_bayes_analysis.R.  Variables marked log = TRUE
# are log-transformed before centering and scaling.

scaled_vars <- c("depth_avg", "MAP", "MAT", "clay", "nitrogen", "irrigation")
# Coerce to plain numeric vectors — clay/nitrogen/irrigation can be list-typed
# if they were extracted via sf subsetting without explicit as.numeric()
all_df <- c_dat %>%
  dplyr::select(all_of(scaled_vars)) %>%
  mutate(across(everything(), as.numeric))

meanlog <- function(x) mean(log(x), na.rm = TRUE)
sdlog   <- function(x) sd(log(x),   na.rm = TRUE)

make_sc_row <- function(var, x, log_tf) {
  x <- as.numeric(x)   # guard against list columns
  tibble(var = var, log = log_tf,
         center = if (log_tf) meanlog(x) else mean(x, na.rm = TRUE),
         scale  = if (log_tf) sdlog(x)   else sd(x,   na.rm = TRUE),
         min  = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE),
         lo99 = quantile(x, 0.005, na.rm = TRUE),
         hi99 = quantile(x, 0.995, na.rm = TRUE))
}

sc_df <- bind_rows(
  tibble(var = "init", log = TRUE, center = NA, scale = NA,
         min  = min(c_dat$Ctrl_c_dens_av, na.rm = TRUE),
         max  = max(c_dat$Ctrl_c_dens_av, na.rm = TRUE),
         lo99 = quantile(c_dat$Ctrl_c_dens_av, 0.005, na.rm = TRUE),
         hi99 = quantile(c_dat$Ctrl_c_dens_av, 0.995, na.rm = TRUE)),
  make_sc_row("depth",      all_df$depth_avg,   TRUE),
  make_sc_row("map",        all_df$MAP,          TRUE),
  make_sc_row("mat",        all_df$MAT,          FALSE),
  make_sc_row("clay",       all_df$clay,         FALSE),
  make_sc_row("nitrogen",   all_df$nitrogen,     TRUE),
  make_sc_row("irrigation", all_df$irrigation,   TRUE),
  make_sc_row("age",        c_dat$AGE,           TRUE)
) %>% as.data.frame()
rownames(sc_df) <- sc_df$var


# =============================================================================
# SECTION 15 — WRITE OUTPUTS
# =============================================================================
# Sensitivity filters for 02_bayes_analysis.R:
#   (a) BD gap-filling:       filter(c_dat, !bd_gap_filled)
#   (b) Management certainty: filter(c_dat, !MGMT_is_uncertain)
#   (c) SOM conversion:       filter(c_dat, !som_converted)

c_dat_csv <- c_dat %>% dplyr::select(-any_of(c("depth", "LUC_sep")))

write.csv(c_dat_csv, "./data_inputs/c_dat_analysis_ready.csv", row.names = FALSE)
write.csv(sc_df,     "./data_inputs/sc_df.csv",                row.names = TRUE)

message("Wrote ./data_inputs/c_dat_analysis_ready.csv  (",
        nrow(c_dat_csv), " rows, ", ncol(c_dat_csv), " columns)")
message("Wrote ./data_inputs/sc_df.csv")

message("\n── Sensitivity subset sizes (of ", nrow(c_dat_csv), " total rows) ────")
message("  (a) !bd_gap_filled       : ", sum(!c_dat_csv$bd_gap_filled,     na.rm = TRUE), " rows")
message("  (b) !MGMT_is_uncertain   : ", sum(!c_dat_csv$MGMT_is_uncertain, na.rm = TRUE), " rows")
message("  (c) !som_converted       : ", sum(!c_dat_csv$som_converted,     na.rm = TRUE), " rows")
message("  (d) orig BD available    : ",
        sum(!is.na(c_dat_csv$Ctrl_c_dens_av_orig_bd) &
              !is.na(c_dat_csv$Aban_c_dens_av_orig_bd), na.rm = TRUE),
        " rows (both Ctrl and Aban BD_raw non-NA)")