### SOCRestorationPotential

This repository contains the data and code used to produce the analyses and figures for the study **"Soils match biomass carbon gains across restored farmlands"** (Terrer et al.).

The analysis quantifies how soil organic carbon (SOC) accrues following the abandonment or planned retirement of cropland and pasture globally. A hierarchical Bayesian model is fitted to a compiled chronosequence dataset and then projected onto global raster layers to produce spatially explicit SOC stock maps. These maps are combined with estimates of above- and belowground biomass carbon and four land-use restoration scenarios to compute global carbon sequestration potentials, including under a projected future climate (SSP3-7.0, 2041–2060).

------------------------------------------------------------------------

### Important notes for users

**Spatial data downloads are required before the analysis can run.** The first script (`00_spatial_data.R`) downloads rasters from Google Earth Engine, WorldClim/CMIP6, SoilGrids, and several other sources. This step requires a working Google Earth Engine account, an `rgee` installation, and substantial download time — expect several hours depending on server load, with each GEE raster taking roughly 10 minutes. Several datasets cannot be downloaded programmatically and must be obtained manually; the script comments point to the relevant sources.

**The full pipeline is computationally intensive.** The Bayesian model (`03_bayes_analysis.R`) is fit via JAGS MCMC and, depending on hardware, can take many hours. The uncertainty propagation script (`05_bayes_uncertainty.R`) runs the full posterior over global rasters using parallel workers and requires substantial RAM. Scripts 04–06 apply model outputs to global rasters and also benefit from high-memory machines. The analysis was originally executed on an Apple MacBook Pro M2 Max with 32 GB RAM running R 4.2.2, and runtimes will vary considerably on other hardware.

Scripts include `refit_models` flags (default `FALSE`) so that cached intermediate outputs can be reused without re-running expensive steps.

------------------------------------------------------------------------

### Explanation of scripts

Scripts are numbered and should be run in order.

-   **`00_spatial_data.R`** — Downloads and preprocesses all spatial inputs. Acquires contemporary climate (WorldClim V1 via GEE), future climate (CMIP6 SSP3-7.0, 2041–2060, averaged over 11 GCMs), SoilGrids soil properties at multiple depth intervals, land cover and land use layers, tree cover, irrigation, and biome classifications. All rasters are resampled to match a 1 km reference grid. Several datasets must be downloaded manually; the script provides the relevant URLs and instructions.

-   **`01_soc_data_preparation.R`** — Reads the raw chronosequence spreadsheet, harmonises units, and prepares the analysis-ready dataset. Key steps include: extracting spatially-varying covariates (MAT, MAP, bulk density, clay, nitrogen, irrigation) from the rasters produced in `00`; gap-filling missing bulk density values using mixed-effects models informed by SoilGrids; converting SOM measurements to SOC; and adding provenance flags for four pre-planned sensitivity analyses (bulk density gap-filling, management uncertainty, SOM conversion, original reported BD). Outputs two files used by all downstream scripts: `c_dat_analysis_ready.csv` and `sc_df.csv` (covariate centering/scaling parameters).

-   **`02_lc_priors.R`** — Derives moderately informative prior distributions for the Bayesian model by regressing log SOC density against environmental covariates on three globally gridded datasets: a no-land-use potential SOC scenario (Sanderman et al. 2017), current cropland areas, and current pasture areas. The resulting regression coefficients inform the `beta_init` and `beta_max` parameter priors in the JAGS model. Outputs are cached as `.RData` files.

-   **`03_bayes_analysis.R`** — Fits the hierarchical Bayesian SOC accrual model via JAGS (MCMC). The model describes how SOC density at a given depth recovers from an initial (agricultural) state toward a potential maximum capacity, following a saturating trajectory governed by environmental covariates (MAP, MAT, clay, nitrogen, irrigation, depth) and land-use covariates (starting land use, restoration pathway, management type). The script evaluates convergence, produces posterior visualisations, and reruns the model under four sensitivity subsets. All six MCMC outputs are saved as `.RData` files for use by downstream scripts.

-   **`03a_env_representativeness.R`** — A standalone supplementary analysis that compares the distribution of environmental predictor variables between the field chronosequence sites and global cropland/pasture areas (weighted by fractional cover). Produces density plots and a multivariate coverage table.

-   **`04_bayes_results.R`** — Applies posterior median coefficients to global raster layers to produce spatially explicit SOC stock maps. For each combination of land-use transition type (8 combinations of starting land use × restoration pathway × management), depth interval, and time horizon, the model's linear predictor is evaluated cell-by-cell using the covariate rasters from `00`. Depth-integrated stocks are corrected for coarse fragment volume to yield final 0–30 cm and 0–100 cm maps. Can be run for both contemporary and SSP3-7.0 2041–2060 climate inputs; results for both scenarios are required by `06`. An `agg_fact` parameter controls spatial resolution (default 10×, i.e. \~10 km); reducing this substantially increases compute time.

-   **`05_bayes_uncertainty.R`** — Propagates posterior uncertainty through the global mapping. Iterates over individual MCMC posterior samples and applies each to the raster stack, producing a set of per-sample global stock maps. Parallelised via `mclapply`. Must be run for both climate scenarios before executing `06`. Configure `n_workers` and `agg_fact` to match available RAM.

-   **`06_scenarios_estimates_maps.R`** — The final synthesis script. Combines SOC stock maps with aboveground and belowground biomass carbon estimates and four land-use restoration scenario rasters (current abandoned cropland, cropland sparing, cropland relocation, pasture sparing) to compute global carbon sequestration estimates in Pg C over time. Credible intervals are constructed from the per-sample maps produced by `05`. Also loads SSP3-7.0 outputs from `04` and `05` to compare SOC sequestration potential under contemporary versus projected future climates. Produces all manuscript and supplementary scenario figures and tables.

------------------------------------------------------------------------

### Additional software setup

#### Google Earth Engine

`00_spatial_data.R` downloads data from the Earth Engine catalog using the `rgee` package (v. 1.1.5). Users must have a Google Earth Engine account and a configured Google Cloud project. Detailed setup instructions: [Earth Engine access](https://developers.google.com/earth-engine/guides/access) and [rgee](https://github.com/r-spatial/rgee).

#### Spatial R packages

The scripts use `sf` (v. 1.0-21) and `terra` (v. 1.8-70). If you have not used these before, follow the installation instructions at [r-spatial.github.io/sf](https://r-spatial.github.io/sf/).

#### JAGS

The Bayesian model requires JAGS (Just Another Gibbs Sampler, v. 4.3.1) and the `rjags` R package. Download JAGS from [mcmc-jags.sourceforge.io](https://mcmc-jags.sourceforge.io/).

------------------------------------------------------------------------

### Data notes

All datasets used in the analyses are publicly available. Where possible, the scripts download data directly from source. In cases where programmatic download is not possible (data access requires a request to the data owner, or data are not shared in a machine-downloadable format), the data are either provided within this repository or the script comments direct users to the source.

------------------------------------------------------------------------

### Computing notes

The full analysis was executed on an Apple MacBook Pro M2 Max with 32 GB RAM running R 4.5.2. Runtimes commented in the scripts reflect this hardware; times will vary on other systems. Data download speeds also depend on server load.
