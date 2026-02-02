# Calanda JSDM - R Scripts Workflow

## Overview

Joint Species Distribution Model (JSDM) analysis of alpine vegetation in the Calanda region, Switzerland. The workflow processes vegetation survey data, environmental variables, and functional traits to fit spatial JSDMs and analyze variance components.

Model outputs used for analysis come from `results_from_Max/` (authoritative model results).

---

## Coding Conventions

### Data passing: File-based

Every script loads its own inputs from saved files and saves its outputs. No script depends on another having been run in the same R session. This makes scripts self-contained and reproducible.

### Naming conventions

| Element | Convention | Example |
|---------|-----------|---------|
| **Objects / variables** | `snake_case` with `_` | `veg_clim`, `snow_metrics`, `species_betas` |
| **Functions** | `snake_case` with `_` | `process_snow_data()`, `calculate_temp_metrics()` |
| **File names** | `snake_case` with `_`, numbered prefix | `01_prepare_data.R`, `03_variance_partitioning.R` |
| **Column names** | `snake_case` with `_` | `plot_id_releve`, `summer_temp`, `et_annual` |
| **Constants / colors** | `snake_case` with `_` | `color_env`, `color_spa` |

**No dots in names.** Dots in R names (e.g., `veg.clim`) can be confused with S3 method dispatch. All existing dotted names will be renamed to underscores during cleanup:
- `veg.clim` → `veg_clim`
- `veg.env` → `veg_env`
- `veg.PA` → `veg_pa`
- `veg.abund` → `veg_abund`
- `veg.comm` → `veg_comm`
- `veg.rare` → `veg_rare`
- `imp.clim` → `imp_clim`
- `imp.traits` → `imp_traits`
- `et.annual` → `et_annual` (column name)
- `res.pca` → `pca_result`

### Other style rules

- Assignment: `=` (not `<-`)
- Booleans: `TRUE` / `FALSE` (not `T` / `F`)
- CSV I/O: `read_csv()` / `write_csv()` (tidyverse, not base R)
- Paths: `here::here()` (not `setwd()`)
- Plots: `theme_bw()` as default theme
- Libraries: each script loads its own (file-based approach)
- No `library(ggplot2)` when `library(tidyverse)` is already loaded
- No `View()` calls (fails non-interactively)

---

## Directory Structure

```
R/
├── 00_setup/
│   ├── 00_workflow.R              # Coordinator: loads libraries, sets paths
│   └── functions_calanda.R        # Helper functions (sourced by other scripts)
│
├── 01_data_prep/
│   ├── 01_prepare_data.R          # Environmental data pipeline (topo, snow, LST, ET)
│   ├── 01_prepare_trait_data.R    # Trait processing (N/C isotopes, LMA, LDMC, PCA)
│   └── assess_trait_coverage.R    # QA: trait coverage per community + bias assessment
│
├── 02_model/
│   └── 02_jsdm.R                  # Fits sjSDM model (GPU-accelerated)
│
├── 03_analysis/
│   ├── 03_variance_partitioning.R # Venn diagram + ternary plots
│   ├── 05_community_postJSDM.R    # Environmental regressions on variance components
│   └── 06_species_postJSDM.R      # Species-level coefficient analysis
│
├── 04_visualization/
│   ├── map.R                      # Geographic maps of study area
│   └── map_rgb_results.R          # 3D DEM with RGB variance coloring
│
└── archive/                       # Deprecated / external scripts
```

---

## Execution Order

```r
# 1. Data preparation (each script is self-contained)
source("R/01_data_prep/01_prepare_data.R")
source("R/01_data_prep/01_prepare_trait_data.R")
source("R/01_data_prep/assess_trait_coverage.R")

# 2. Model fitting (GPU required) -- or use results_from_Max/
source("R/02_model/02_jsdm.R")

# 3. Post-model analysis (loads from results_from_Max/ and output/)
source("R/03_analysis/03_variance_partitioning.R")
source("R/03_analysis/05_community_postJSDM.R")
source("R/03_analysis/06_species_postJSDM.R")

# 4. Visualization
source("R/04_visualization/map.R")
source("R/04_visualization/map_rgb_results.R")
```

---

## Key Data Flow

```
Raw data (data/)
    │
    ├── Vegetation surveys ──► 01_prepare_data.R ──► output/starter_data.RData
    ├── Topography (SwissALTI3D)  ──┘                     │
    ├── Snow (Copernicus)         ──┘                     │
    ├── LST (ECOSTRESS)           ──┘                     │
    └── ET (MODIS)                ──┘                     │
                                                          ▼
    Trait data (data/traits/) ──► 01_prepare_trait_data.R ──► output/species_trait_summary.csv
                                                          │
                                                          ▼
                                      02_jsdm.R / results_from_Max/
                                                          │
                                                          ▼
                                      03_analysis/ scripts ──► plot/ PDFs
```

---

## Cleanup Plan: Step-by-Step

After screening every active script, here is a comprehensive list of issues and tasks to address, organized script by script. We will work through these one at a time.

### Required R packages (full audit)

Each script will load only the packages it needs. Here is the full inventory:

| Package | Used by | Notes |
|---------|---------|-------|
| `tidyverse` | all scripts | Core framework |
| `sjSDM` | 02_jsdm, 03_variance_partitioning, 05/06_postJSDM | JSDM fitting & analysis |
| `here` | all scripts | Project-relative paths |
| `terra` | 01_prepare_data, map_rgb_results | Raster processing |
| `sf` | 01_prepare_data, map, map_rgb_results | Spatial vectors |
| `conflicted` | 00_workflow (if kept) | Namespace conflict resolution |
| `stringi` | 01_prepare_data | String transliteration |
| `miceRanger` | functions_calanda (imputation) | Missing data imputation |
| `corrplot` | 01_prepare_data, 01_prepare_trait_data | Correlation plots |
| `FactoMineR` | 01_prepare_data, 01_prepare_trait_data | PCA |
| `factoextra` | 01_prepare_data, 01_prepare_trait_data | PCA visualization |
| `gt` | 01_prepare_data | Tables |
| `ggrepel` | 03_variance_partitioning, 05/06_postJSDM, 01_prepare_trait_data | Label repulsion |
| `patchwork` | assess_trait_coverage | Plot composition |
| `scales` | 03_variance_partitioning | Rescaling |
| `readxl` | 01_prepare_trait_data | Excel file reading |
| `dgof` | assess_trait_coverage | Cramer-von Mises test |
| `rnaturalearth` / `rnaturalearthdata` | map | Country outlines |
| `rayshader` | map_rgb_results | 3D rendering |
| `raster` | map_rgb_results | raster_to_matrix (rayshader needs it) |
| `elevatr` | map_rgb_results | Elevation downloads |
| `tidyterra` | map_rgb_results | ggplot + terra integration |
| `ggtern` | map_rgb_results | Ternary plots (evaluate if needed) |

**Packages to REMOVE** (no longer needed after archiving):
`qgam`, `mvabund`, `randomForest`, `caret`, `ggpubr`, `gridExtra`, `grid`, `ggforce`, `TNRS`, `spatialEco`

---

### TASK 1: Overhaul `00_setup/00_workflow.R`

**Current problems:**
- Loads ~30 libraries with many duplicates (`tidyverse` x2, `sf` x2, `terra` x2, `patchwork` x2, `factoextra` x2, `FactoMineR` x2, `ggrepel` x2, `mvabund` x2)
- Loads libraries no longer needed (archived scripts' dependencies)
- Has `setwd("Calanda_JSDM/")` hardcoded
- Sources `functions_calanda.R` twice (lines 50 and 88)
- Contains inline data processing (veg_coord extraction, ecostress coord export) that belongs in `01_prepare_data.R`
- Loads `TransPlantNetwork` data (archived)
- Has commented-out reticulate/Python config code
- Path references broken after reorganization

**Decision: This script may become unnecessary.** With the file-based approach, each script loads its own libraries and data. `00_workflow.R` could be reduced to just sourcing `functions_calanda.R`, or eliminated entirely if functions are sourced directly by scripts that need them.

**Tasks:**
- [ ] Decide: keep as minimal coordinator or eliminate entirely
- [ ] If kept: strip down to only `here` setup + `source(functions_calanda.R)`
- [ ] Move data processing code to `01_prepare_data.R`
- [ ] Remove all unused library calls, duplicates, commented-out code
- [ ] Remove TransPlant data loading

---

### TASK 2: Clean `00_setup/functions_calanda.R`

**Current problems:**
- Starts with `bioclim_data` + `match_clim_var()` (lines 1-59) -- never called by any active script
- ~2100 lines mixing data processing, imputation, spatial processing, and custom plotting
- Contains base-R ternary plot functions -- inconsistent with ggplot2 elsewhere (but still called by `03_variance_partitioning.R`)
- Some functions may no longer be called after archiving

**Tasks:**
- [ ] Remove `bioclim_data` and `match_clim_var()` (unused)
- [ ] Audit every function: which are still called by active scripts?
- [ ] Remove functions only used by archived scripts
- [ ] Rename all dotted object names inside functions (`veg.clim` → `veg_clim` etc.)
- [ ] Add roxygen-style documentation headers to remaining functions
- [ ] Consider splitting into `functions_data.R` and `functions_plotting.R`

---

### TASK 3: Clean `01_data_prep/01_prepare_data.R`

**Current problems:**
- No library loading -- relies on `00_workflow.R` workspace
- Massive commented-out code blocks (lines 212-276): old land-use, climate extraction
- Inconsistent `dplyr::select()` vs `select()`
- Uses `View()` (line 395)
- Uses deprecated `size` in ggplot2 (should be `linewidth`)
- `if(file.exists(...))` caching makes logic flow hard to follow
- Hardcoded output path with date: `"output/starter_data_25.04.25.RData"`
- Trait loading/imputation duplicated with archived `03_explain_variation.R`
- `save.image()` commented out -- RData not actually saved
- Mixed `write.csv` / `write_csv`
- Dotted object names throughout

**Tasks:**
- [ ] Add library loading: `tidyverse`, `sf`, `terra`, `stringi`, `corrplot`, `FactoMineR`, `factoextra`, `gt`, `here`
- [ ] Add `source(here(..., "functions_calanda.R"))` for custom functions
- [ ] Remove all commented-out code blocks
- [ ] Remove `View()` call
- [ ] Fix `size` → `linewidth`
- [ ] Standardize to `write_csv()` everywhere
- [ ] Rename all dotted objects: `veg.clim` → `veg_clim`, etc.
- [ ] Rename `et.annual` column → `et_annual`
- [ ] Properly save the RData output
- [ ] Add clear section headers
- [ ] Remove `select(-...1)` hacks -- fix by using `write_csv()` upstream

---

### TASK 4: Clean `01_data_prep/01_prepare_trait_data.R`

**Current problems:**
- Loads its own libraries (good)
- Best-structured script in the repo
- Uses `<-` assignment in some places, `=` in others
- Contains `clean_nc_data()` function inline -- could go in `functions_calanda.R`

**Tasks:**
- [ ] Standardize `<-` → `=`
- [ ] Move `clean_nc_data()` to `functions_calanda.R`
- [ ] Add `source(here(..., "functions_calanda.R"))` call
- [ ] Source `assess_trait_coverage.R` at the end (or keep separate)
- [ ] Use `here()` for all file paths

---

### TASK 5: Clean `01_data_prep/assess_trait_coverage.R`

**Current problems:**
- Standalone script that loads its own data (good for file-based)
- Species name remapping (lines 13-22) should be centralized
- Loads JSDM `res` for bias assessment -- tight coupling with analysis layer

**Tasks:**
- [ ] Centralize species name remapping in `functions_calanda.R`
- [ ] Split into two parts: (1) pure trait coverage (data_prep), (2) JSDM bias assessment (move to 03_analysis/)
- [ ] Use `here()` for paths
- [ ] Rename dotted objects

---

### TASK 6: Clean `02_model/02_jsdm.R`

**Current problems:**
- Clean and focused
- Hardcoded output filenames with dates (`260425`)
- No library loading
- `se=T` should be `se = TRUE`

**Tasks:**
- [ ] Add library loading: `sjSDM`, `tidyverse`, `here`
- [ ] Add script header documentation
- [ ] `T` → `TRUE`
- [ ] Use `here()` for paths
- [ ] Parameterize output filenames (or remove date suffix)

---

### TASK 7: Clean `03_analysis/03_variance_partitioning.R` (NEW)

**Status:** Just created from the kept code of archived `03_explain_variation.R`.

**Tasks:**
- [ ] Already uses `here()` and loads its own data -- good
- [ ] Rename dotted objects when loaded (`veg.clim` → `veg_clim`)
- [ ] Verify `plot_tern_sites()`, `plot_tern_species()`, `plot.anova.custom()` still work
- [ ] Rename `plot.anova.custom()` → `plot_anova_custom()` in functions_calanda.R

---

### TASK 8: Clean `03_analysis/05_community_postJSDM.R`

**Current problems:**
- No library loading
- References `res`, `veg.env` as globals
- Duplicated helpers (`label_env_var`, `predict_model`, `extract_coefs`) with `06_species_postJSDM.R`

**Tasks:**
- [ ] Add library loading and file-based data loading
- [ ] Extract shared helpers to `functions_calanda.R`
- [ ] Rename dotted objects
- [ ] Use `here()` for paths

---

### TASK 9: Clean `03_analysis/06_species_postJSDM.R`

**Current problems:**
- No library loading
- References `model_jsdm` (line 7) but model loaded as `model` elsewhere -- naming inconsistency
- Duplicated helpers with `05_community_postJSDM.R`

**Tasks:**
- [ ] Add library loading and file-based data loading
- [ ] Standardize model variable name
- [ ] Extract shared helpers to `functions_calanda.R`
- [ ] Rename dotted objects
- [ ] Use `here()` for paths

---

### TASK 10: Clean `04_visualization/map.R`

**Current problems:**
- References `veg.clim` and `veg_bbox` as globals
- `veg_bbox` is NEVER DEFINED anywhere -- script is broken
- Short script but incomplete

**Tasks:**
- [ ] Add file-based data loading
- [ ] Define/compute `veg_bbox` from data
- [ ] Rename dotted objects
- [ ] Use `here()` for paths
- [ ] Make fully self-contained

---

### TASK 11: Clean `04_visualization/map_rgb_results.R`

**Current problems:**
- Uses `setwd(here("Calanda_JSDM"))` -- fragile
- Uses `raster` (superseded by `terra`) -- but needed for `raster_to_matrix`
- Uses `ggtern` which had installation issues
- References `veg`, `veg.env` as globals
- `raster::extract()` conflicts with `terra::extract()`

**Tasks:**
- [ ] Replace `setwd()` with `here()` paths
- [ ] File-based data loading
- [ ] Evaluate `ggtern` -- replace with custom ternary if possible
- [ ] Rename dotted objects
- [ ] Document `raster` dependency (required by rayshader)

---

### TASK 12: Cross-cutting standardization (final pass)

After all individual scripts are cleaned:

- [ ] Rename all dotted objects across entire codebase (see naming table above)
- [ ] Rename `et.annual` → `et_annual` in all CSVs and code
- [ ] Rename `plot.anova.custom` → `plot_anova_custom`
- [ ] Verify all `source()` paths work with new folder structure
- [ ] Verify all `here()` paths resolve correctly
- [ ] Run each script independently to confirm file-based approach works
- [ ] Remove `00_workflow.R` if no longer needed
- [ ] Final check: no `T`/`F`, no `write.csv`, no `View()`, no `setwd()`

---

### Suggested execution order for cleanup

| Order | Task | Script | Priority |
|-------|------|--------|----------|
| 1 | TASK 2 | `functions_calanda.R` (audit + clean) | High |
| 2 | TASK 1 | `00_workflow.R` (strip down or remove) | High |
| 3 | TASK 3 | `01_prepare_data.R` | High |
| 4 | TASK 4 | `01_prepare_trait_data.R` | Medium |
| 5 | TASK 5 | `assess_trait_coverage.R` | Medium |
| 6 | TASK 6 | `02_jsdm.R` | Low |
| 7 | TASK 7 | `03_variance_partitioning.R` | Medium |
| 8 | TASK 8 | `05_community_postJSDM.R` | Medium |
| 9 | TASK 9 | `06_species_postJSDM.R` | Medium |
| 10 | TASK 10 | `map.R` | Low |
| 11 | TASK 11 | `map_rgb_results.R` | Low |
| 12 | TASK 12 | Cross-cutting standardization | Final pass |

---

*Last updated: 2026-02-02*
