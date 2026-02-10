# Calanda JSDM - Workflow

Joint Species Distribution Model (JSDM) analysis of alpine vegetation in the Calanda region, Switzerland. The workflow processes vegetation survey data, remote-sensing environmental variables, and functional traits to fit a spatial JSDM and analyze variance components.

Authoritative model results come from `results_from_Max/` (GPU-fitted). Analysis and visualization scripts load from there.

---

## Directory Structure

```
R/
├── 00_setup/
│   ├── 00_workflow.R              # Coordinator: sources all scripts in order
│   └── functions_calanda.R        # Shared helper functions
│
├── 01_data_prep/
│   ├── 01_prepare_data.R          # Environmental data pipeline
│   ├── 02_prepare_trait_data.R    # Trait processing (isotopes, LMA, LDMC, PCA)
│   └── 03_assess_trait_coverage.R # QA: trait coverage + bias assessment
│
├── 02_model/
│   └── 04_jsdm.R                  # Fits sjSDM model (GPU required)
│
├── 03_analysis/
│   ├── 05_variance_partitioning.R # Venn diagram + ternary plots
│   ├── 06_community_postJSDM.R    # Community-level environmental regressions
│   └── 07_species_postJSDM.R      # Species-level coefficient analysis
│
├── 04_visualization/
│   ├── 08_map.R                   # Study area maps
│   └── 09_map_rgb_results.R       # 3D DEM with RGB variance coloring
│
└── archive/                       # Deprecated scripts
```

---

## Execution Order

Run `00_workflow.R` to execute the full pipeline, or source scripts individually:

```r
# Step 1: Data preparation
source("R/01_data_prep/01_prepare_data.R")
source("R/01_data_prep/02_prepare_trait_data.R")
source("R/01_data_prep/03_assess_trait_coverage.R")

# Step 2: Model fitting (GPU required) -- or use results_from_Max/
source("R/02_model/04_jsdm.R")

# Step 3: Post-model analysis
source("R/03_analysis/05_variance_partitioning.R")
source("R/03_analysis/06_community_postJSDM.R")
source("R/03_analysis/07_species_postJSDM.R")

# Step 4: Visualization
source("R/04_visualization/08_map.R")
source("R/04_visualization/09_map_rgb_results.R")
```

---

## Scripts

### `00_setup/functions_calanda.R`

Shared helper functions sourced by other scripts. Not run directly.

| Function | What it does | Key library |
|----------|-------------|-------------|
| `process_gfsc_data()` | Extracts Copernicus snow cover from ZIP archives, applies BISE correction, spline interpolation, and Savitzky-Golay smoothing | `terra::rast()`, `terra::extract()`, `zoo::rollapply()`, `zoo::na.spline()`, `signal::sgolayfilt()` |
| `calculate_snow_metrics()` | Computes snow disappearance date and snow-cover days from smoothed time series | — |
| `calculate_temp_metrics()` | Computes summer temperature and spring warming date from ECOSTRESS LST | — |
| `impute_environmental_data()` | Random-forest imputation of missing environmental variables | `miceRanger::miceRanger()` |
| `impute_functional_traits()` | Random-forest imputation of missing trait values using nearest-neighbour distances | `miceRanger::miceRanger()`, `sf::st_distance()` |
| `extract_mowing_events()` | Extracts grassland-use intensity from raster time series | `terra::rast()`, `terra::writeRaster()` |
| `calculate_cwm()` | Computes community-weighted mean traits | — |
| `plot_anova_custom()` | Draws Venn diagram of variance partitioning | — |
| `plot_tern_species()` / `plot_tern_sites()` | Ternary plots of variance components colored by altitude | `ggrepel::geom_text_repel()` |
| `label_env_var()` | Formats environmental variable names for plot labels | — |

---

### `01_data_prep/01_prepare_data.R`

Builds the environmental matrix (X) and species matrix (Y) for the JSDM.

**What it does:** Loads raw vegetation surveys, processes remote-sensing variables (Copernicus snow, ECOSTRESS LST, MODIS ET, SwissALTI3D topography), computes derived variables (freezing degree days, snow metrics, temperature metrics), imputes missing values, creates a land-use proxy via PCA on community-weighted mean traits, scales and assembles the final X and Y matrices.

| Inputs | |
|--------|--|
| `data/vegetation/2024_CAPHE_SpeDis_CleanData_20240214.csv` | Vegetation survey data |
| `data/vegetation/veg.coord.csv` | Plot coordinates |
| `data/mask/study_region_2024.shp`, `data/mask/study_region.shp` | Study area boundaries |
| `data/ecostress/*.csv` | ECOSTRESS land surface temperature |
| `data/modis/*.csv` | MODIS evapotranspiration |
| `data/wekeo/` | Copernicus GFSC snow cover ZIPs |
| `data/traits/try.quantitative_traits_*.csv` | TRY trait data |
| `data/traits/indicators_cleaned_calanda_*.csv` | Ecological indicator values |
| `data/traits/dispersal_cleaned_calanda_*.csv` | Dispersal traits |

| Outputs | |
|---------|--|
| `output/veg_clim.csv` | Merged vegetation + environmental data |
| `output/traits.csv` | Processed trait data |
| `output/starter_data_25.04.25.RData` | All R objects for downstream scripts |
| `output/calanda_mask.shp` | Merged study area polygon |
| `output/snow_metrics.csv` | Snow disappearance dates and cover days |
| `plot/land_use_variable.pdf` | PCA biplot of land-use proxy |
| `plot/climate_calanda_surveys.pdf` | Climate variable overview |

| Key library usage | |
|-------------------|--|
| `sf::st_read()`, `sf::st_union()`, `sf::st_write()` | Read, merge, and write study area shapefiles |
| `terra` (via helper functions) | Raster I/O for snow cover, topography |
| `stringi::stri_trans_general()` | Transliterate species names from Latin to ASCII |
| `miceRanger::miceRanger()` (via helpers) | Impute missing environmental and trait values |
| `FactoMineR::PCA()` | PCA on topography, land-use proxy, and final environmental variables |
| `factoextra::fviz_pca_biplot()` | PCA biplot visualization |
| `corrplot::corrplot()` | Correlation matrices for variable selection |
| `gt::gt()` | Render imputation performance tables |

---

### `01_data_prep/02_prepare_trait_data.R`

Processes raw functional trait measurements into species-level summaries.

**What it does:** Parses leaf N/C isotope data from Excel files, maps CN analysis IDs to plant IDs, merges with leaf area measurements, calculates derived traits (LMA, LDMC), validates data quality (weight errors, nutrient outliers), computes species-level means and coefficients of variation (kCV), runs PCA on mean traits and kCV values.

| Inputs | |
|--------|--|
| `data/traits/N+C_Samples*.xls` | Isotope analysis results (multiple files) |
| `data/traits/Trial*.xls` | Additional isotope runs |
| `data/traits/2025_CAPHE_traits_sample_id.csv` | Sample ID mapping |
| `data/traits/2025_TRAITS_clean_20251023.csv` | Cleaned trait measurements |
| `data/traits/leaf_area.csv` | Scanned leaf areas |

| Outputs | |
|---------|--|
| `output/species_trait_summary.csv` | Species-level trait means and kCV |
| `output/final_traits_clean.csv` | Individual-level cleaned traits |
| `output/all_traits_merged.csv` | All traits merged before cleaning |
| `output/pca_means_eigenvalues.csv`, `output/pca_kcv_eigenvalues.csv` | PCA eigenvalues |
| `output/pca_means_var_contrib.csv`, `output/pca_kcv_var_contrib.csv` | PCA variable contributions |
| `plot/corrplot_all_traits.pdf` | Trait correlation matrix |
| `plot/pca_trait_means_biplot.pdf` | PCA biplot of trait means |
| `plot/pca_kcv_biplot.pdf` | PCA biplot of trait variability |

| Key library usage | |
|-------------------|--|
| `readxl::read_excel()` | Parse multi-section Excel isotope data sheets |
| `FactoMineR::PCA()` | PCA on species-level mean traits and kCV values |
| `corrplot::corrplot()` | Correlation matrix of all traits |
| `ggrepel::geom_text_repel()` | Non-overlapping species labels on PCA biplots |

---

### `01_data_prep/03_assess_trait_coverage.R`

Quality assurance: tests whether trait-sampled communities are representative.

**What it does:** Computes per-plot trait coverage (proportion of species with trait data weighted by abundance). Tests whether plots with high trait coverage differ systematically in JSDM variance components from the full dataset using t-tests, Wilcoxon tests, and Cramer-von Mises tests at both community and species levels.

| Inputs | |
|--------|--|
| `output/species_trait_summary.csv` | Species-level trait summary |
| `output/starter_data_25.04.25.RData` | Abundance matrix and JSDM results |

| Outputs | |
|---------|--|
| `output/community_trait_coverage.csv` | Per-plot trait coverage |
| `output/coverage_bias_summary.csv` | Bias test results (community level) |
| `output/species_bias_summary.csv` | Bias test results (species level) |
| `plot/trait_coverage_assessment.pdf` | Coverage diagnostic plots |
| `plot/trait_coverage_bias.pdf` | Community-level bias plots |
| `plot/species_trait_coverage_bias.pdf` | Species-level bias plots |

| Key library usage | |
|-------------------|--|
| `dgof::cvm.test()` | Cramer-von Mises goodness-of-fit test on variance component distributions |
| `patchwork` | Multi-panel plot composition |

---

### `02_model/04_jsdm.R`

Fits the spatial joint species distribution model. **Requires GPU.**

**What it does:** Defines a linear environmental predictor, a deep neural network for the spatial component (coordinates as input, 2 hidden layers of 30 units with SELU activation), and a biotic covariance structure. Fits the model with RMSprop optimizer and learning rate scheduling. Runs variance partitioning (R-squared, ANOVA, internal structure).

| Inputs | |
|--------|--|
| `output/data_calanda_jsdm.rds` | Scaled X (environment + coordinates) and Y (species) matrices |

| Outputs | |
|---------|--|
| `output/model_sjsdm_calanda.rds` | Fitted sjSDM model object |
| `output/R2_sjsdm_calanda.rds` | Model R-squared |
| `output/an_sjsdm_calanda.rds` | ANOVA variance partitioning |
| `output/res_sjsdm_calanda.rds` | Site-level and species-level variance fractions |

| Key library usage | |
|-------------------|--|
| `sjSDM::sjSDM()` | Fit the joint species distribution model |
| `sjSDM::linear()` | Define environmental predictor with elastic-net regularization |
| `sjSDM::DNN()` | Define deep neural network for spatial component |
| `sjSDM::bioticStruct()` | Define species covariance structure |
| `sjSDM::Rsquared()` | Compute model R-squared |
| `sjSDM::anova()` | Partition variance across environment, space, and biotic components |
| `sjSDM::internalStructure()` | Extract site-level and species-level variance fractions |

---

### `03_analysis/05_variance_partitioning.R`

Visualizes the overall variance partitioning results.

**What it does:** Creates a Venn diagram showing shared and unique variance explained by environment, species associations, and space. Creates ternary plots for both sites and species, colored by altitude.

| Inputs | |
|--------|--|
| `output/starter_data_25.04.25.RData` | Plot metadata (altitude, coordinates) |
| `results_from_Max/an_sjsdm_calanda.rds` | ANOVA results |
| `results_from_Max/model_sjsdm_calanda.rds` | Fitted model |
| `results_from_Max/res_sjsdm_calanda.rds` | Variance fractions |

| Outputs | |
|---------|--|
| `plot/venn.pdf` | Venn diagram of variance partitioning |
| `plot/ternary_sites.pdf` | Ternary plot of site-level variance |
| `plot/ternary_species.pdf` | Ternary plot of species-level variance |

| Key library usage | |
|-------------------|--|
| `sjSDM` | Loaded for class compatibility when reading model objects |
| Custom `plot_anova_custom()` | Venn diagram from ANOVA results |
| Custom `plot_tern_sites()` / `plot_tern_species()` | Ternary plots |

---

### `03_analysis/06_community_postJSDM.R`

Tests what drives community-level variance components along environmental gradients.

**What it does:** Regresses site-level variance proportions (environment, species associations) against environmental gradients (summer temperature, freezing degree days, ET, land use) using quadratic models with backward stepwise AIC selection. Separately regresses spatial variance against topographic variables (slope, TPI, roughness).

| Inputs | |
|--------|--|
| `output/starter_data_25.04.25.RData` | Environmental data per site |
| `results_from_Max/res_sjsdm_calanda.rds` | Site-level variance fractions |

| Outputs | |
|---------|--|
| `plot/community_regressions.pdf` | Variance vs. environment scatter plots |
| `plot/spatial_topography_regressions.pdf` | Spatial variance vs. topography |

| Key library usage | |
|-------------------|--|
| `stats::step()` | Backward stepwise AIC model selection |
| `ggrepel::geom_text_repel()` | Annotate effect sizes on plots |

---

### `03_analysis/07_species_postJSDM.R`

Tests what drives species-level co-distribution variance.

**What it does:** Extracts per-species environmental regression coefficients (betas) from the fitted sjSDM model. Regresses species-level co-distribution variance against those betas to test whether species with stronger environmental responses have higher or lower species-association variance. Uses stepwise AIC-selected quadratic models.

| Inputs | |
|--------|--|
| `results_from_Max/model_sjsdm_calanda.rds` | Fitted model (for species betas) |
| `results_from_Max/res_sjsdm_calanda.rds` | Species-level variance fractions |

| Outputs | |
|---------|--|
| `plot/species_regressions.pdf` | Species variance vs. beta scatter plots |

| Key library usage | |
|-------------------|--|
| `sjSDM::summary()` | Extract species-level environmental coefficients |
| `stats::step()` | Backward stepwise AIC model selection |
| `ggrepel::geom_text_repel()` | Annotate effect sizes on plots |

---

### `04_visualization/08_map.R`

Study area maps.

**What it does:** Creates an overview map of Switzerland highlighting the Calanda study region, and a zoomed map of vegetation survey points colored by altitude.

| Inputs | |
|--------|--|
| `output/starter_data_25.04.25.RData` | Plot coordinates and altitude |
| `output/calanda_mask.shp` | Study area polygon |

| Outputs | |
|---------|--|
| `plot/map_calanda.pdf` | Overview + zoomed study area maps |

| Key library usage | |
|-------------------|--|
| `rnaturalearth::ne_countries()` | Download Switzerland country boundary |
| `sf::st_as_sf()` | Convert coordinates to spatial points |
| `sf::st_bbox()` | Compute bounding box for zoom extent |

---

### `04_visualization/09_map_rgb_results.R`

3D and 2D RGB visualizations of variance components on the landscape.

**What it does:** Maps variance component proportions (environment, space, species associations) to RGB color channels (blue, red, green). Creates a ternary plot, a 3D rayshader rendering over the DEM with elevated RGB points and vertical stalks, a 2D map, an RGB color legend, and a DEM contour plot with overlaid RGB points.

| Inputs | |
|--------|--|
| `output/starter_data_25.04.25.RData` | Plot coordinates and metadata |
| `results_from_Max/res_sjsdm_calanda.rds` | Site-level variance fractions |
| `output/calanda_mask.shp` | Study area polygon |
| `output/metrics/dem.tif` | High-resolution DEM |

| Outputs | |
|---------|--|
| `plot/ternary_plot_variance_components.pdf` | Ternary plot with habitat type |
| `plot/calanda_rgb_3d.png` | 3D rayshader rendering |
| `plot/calanda_rgb_2d_map.pdf` | 2D map with RGB points |
| `plot/rgb_legend.pdf` | RGB color legend |
| `plot/calanda_dem_contour_rgb.pdf` | DEM contour plot with RGB points |

| Key library usage | |
|-------------------|--|
| `ggtern::ggtern()` | Ternary plot of variance proportions |
| `elevatr::get_elev_raster()` | Download elevation raster tiles |
| `rayshader::sphere_shade()`, `ray_shade()`, `ambient_shade()` | Compute surface textures and shadows |
| `rayshader::plot_3d()` | Render 3D elevation surface |
| `rayshader::render_points()` | Add RGB points above the terrain |
| `rayshader::render_path()` | Draw vertical stalks from terrain to points |
| `rayshader::render_snapshot()` | Save 3D scene as PNG |
| `terra::rast()`, `terra::crop()`, `terra::aggregate()` | Load and process DEM raster |
| `tidyterra::geom_spatraster()`, `geom_spatraster_contour()` | Plot raster and contours in ggplot |

---

## Data Flow

```
Raw data (data/)
    │
    ├── Vegetation surveys ──┐
    ├── Topography (SwissALTI3D) ──┤
    ├── Snow (Copernicus GFSC) ──┤
    ├── LST (ECOSTRESS) ──┤        01_prepare_data.R
    ├── ET (MODIS) ──┤                    │
    └── Traits (TRY + field) ──┘          │
                                          ▼
                              output/starter_data_25.04.25.RData
                              output/veg_clim.csv
                                          │
    Trait samples (data/traits/) ──► 02_prepare_trait_data.R
                                          │
                                          ▼
                              output/species_trait_summary.csv
                                          │
                                          ▼
                              03_assess_trait_coverage.R (QA)
                                          │
                                          ▼
                              04_jsdm.R / results_from_Max/
                                          │
                              ┌───────────┼───────────┐
                              ▼           ▼           ▼
                      05_variance   06_community  07_species
                      _partitioning _postJSDM     _postJSDM
                              │           │           │
                              └───────────┼───────────┘
                                          ▼
                                      plot/ PDFs
                                          │
                              ┌───────────┴───────────┐
                              ▼                       ▼
                        08_map.R           09_map_rgb_results.R
                              │                       │
                              ▼                       ▼
                      plot/map_calanda.pdf    plot/calanda_rgb_*.pdf/png
```

---

## Coding Conventions

| Rule | Convention |
|------|-----------|
| **Naming** | `snake_case` everywhere (objects, functions, files, columns). No dots in names. |
| **Assignment** | `=` (not `<-`) |
| **Booleans** | `TRUE` / `FALSE` (not `T` / `F`) |
| **CSV I/O** | `read_csv()` / `write_csv()` (tidyverse) |
| **Paths** | `here::here()` (no `setwd()`) |
| **Plots** | `theme_bw()` as default |
| **Libraries** | Each script loads its own (file-based, self-contained) |
| **Conflicts** | `library(conflicted)` + `conflict_prefer()` in scripts that load conflicting packages |
| **No** | `View()`, `setwd()`, `library(ggplot2)` when tidyverse is loaded |

---

*Last updated: 2026-02-02*
