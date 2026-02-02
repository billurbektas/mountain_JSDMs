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

### TASK 1: Overhaul `00_setup/00_workflow.R` ✅ DONE

**Completed.** Rewritten as a minimal coordinator script. Only loads `library(here)`, then sources all pipeline scripts in execution order (Steps 1-4). Model fitting (Step 2) is gated behind `run_model = FALSE` flag since authoritative results come from `results_from_Max/`. Removed all library loading (~30 duplicated/unused), `setwd()`, inline data processing, TransPlant data, commented-out reticulate code, and duplicate `source()` calls.

---

### TASK 2: Clean `00_setup/functions_calanda.R` ✅ DONE

**Completed.** Removed 8 unused functions/objects, renamed `plot.anova.custom` → `plot_anova_custom`, `addA` → `add_alpha`, renamed parameters (`veg.clim` → `veg_clim`, `num.trees` → `num_trees`, `et.annual` → `et_annual`), added roxygen docs, organized into 2 sections. Reduced from ~2120 to ~1693 lines.

---

### TASK 3: Clean `01_data_prep/01_prepare_data.R` ✅ DONE

**Completed.** Added self-contained library loading, script header, `here()` paths throughout. Renamed all dotted objects (`veg.clim` → `veg_clim`, `veg.comm` → `veg_comm`, `veg.PA` → `veg_pa`, `veg.abund` → `veg_abund`, `veg.rare` → `veg_rare`, `veg.env` → `veg_env`, `imp.clim` → `imp_clim`, `imp.traits` → `imp_traits`, `et.annual` → `et_annual`, `res.pca` → `res_pca_*`). Removed `View()`, commented-out blocks, `select(-...1)` hacks. `write.csv` → `write_csv`. `size` → `linewidth`. Output CSV renamed `veg.clim.csv` → `veg_clim.csv`. `try` → `try_data` to avoid shadowing base R. `theme_bw()` applied.

**Note:** Output CSV renamed from `veg.clim.csv` → `veg_clim.csv`. No other active scripts load this CSV directly (they use RData), so no breakage.

---

### TASK 4: Clean `01_data_prep/01_prepare_trait_data.R` ✅ DONE

**Completed.** Standardized `<-` → `=` throughout. All paths now use `here()`. Removed redundant `library(ggplot2)` (already in tidyverse), added `library(here)`. `theme_minimal()` → `theme_bw()`. Kept `clean_nc_data()` and `kcv()` inline (single-use functions, not worth abstracting to functions_calanda.R).

---

### TASK 5: Clean `01_data_prep/assess_trait_coverage.R` ✅ DONE

**Completed.** Added script header with inputs/outputs. Added `library(here)`, `library(patchwork)` at top, removed mid-script `library()` calls and redundant `library(ggplot2)`. All paths now use `here()`. Replaced `veg.abund` → `veg_abund` (with rename-on-load from RData). Made data loading unconditional (file-based, no `exists()` checks). Kept as single script (coverage + bias assessment are tightly coupled analytically).

---

### TASK 6: Clean `02_model/02_jsdm.R` ✅ DONE

**Completed.** Added script header with inputs/outputs. Added `library(here)`. All paths now use `here()`. `se=T` → `se = TRUE`. Removed date suffixes from output filenames. Renamed dotted hyperparameters (`lambda.env` → `lambda_env`, etc.). `et.annual` → `et_annual` in model formula. Removed `library(conflicted)` / `conflict_prefer` calls (using `dplyr::select` explicitly instead). Removed commented-out code.

---

### TASK 7: Clean `03_analysis/03_variance_partitioning.R` ✅ DONE

**Completed.** Updated function calls to match renamed functions (`plot.anova.custom` → `plot_anova_custom`), renamed `veg.clim` → `veg_clim` (with rename-on-load from RData), updated parameter names in ternary function calls.

---

### TASK 8: Clean `03_analysis/05_community_postJSDM.R` ✅ DONE

**Completed.** Added script header, self-contained library loading (`tidyverse`, `ggrepel`, `here`), file-based data loading with `here()`. Sources `functions_calanda.R` for shared `label_env_var()`. `veg.env` → `veg_env` (rename-on-load). `et.annual` → `et_annual` in formulas, variable references, and helper functions. All paths use `here()`. Removed inline `label_env_var` (moved to `functions_calanda.R`). Kept `predict_model` and `extract_coefs` inline (closure over script-local `models` list).

---

### TASK 9: Clean `03_analysis/06_species_postJSDM.R` ✅ DONE

**Completed.** Added script header, self-contained library loading (`tidyverse`, `sjSDM`, `ggrepel`, `here`), file-based data loading from `results_from_Max/`. Sources `functions_calanda.R` for shared `label_env_var()`. `et.annual` → `et_annual` in variable lists and `extract_coefs`. Removed inline `label_env_var` (using shared version). All paths use `here()`.

---

### TASK 10: Clean `04_visualization/map.R` ✅ DONE

**Completed.** Fixed broken `veg_bbox` — added `veg_bbox = st_bbox(veg_points)`. Added script header, self-contained library loading. `veg.clim` → `veg_clim` (rename-on-load). `<-` → `=`. `theme_minimal()` → `theme_bw()`. Removed redundant `library(ggplot2)`. All paths use `here()`. Renamed `my_shapefile` → `calanda_mask`.

---

### TASK 11: Clean `04_visualization/map_rgb_results.R` ✅ DONE

**Completed.** Removed `setwd()`. Added script header with inputs/outputs. All paths now use `here()`. `veg.env` → `veg_env` (rename-on-load). `theme_minimal()` → `theme_bw()`. Kept `raster` dependency (required by rayshader's `raster_to_matrix`). Kept `ggtern` dependency (used for ternary plot in this script).

---

### TASK 12: Cross-cutting standardization (final pass) ✅ DONE

All items verified:

- [x] Rename `plot.anova.custom` → `plot_anova_custom` (functions + 03_variance_partitioning)
- [x] Rename `et.annual` → `et_annual` (all active scripts: 01_prepare_data, 02_jsdm, 05_community, 06_species)
- [x] All active scripts use `veg_env` (not `veg.env`) — 05_community, map_rgb use rename-on-load
- [x] All active scripts use `veg_clim` (not `veg.clim`) — 03_variance, map use rename-on-load
- [x] All active scripts use `veg_abund` (not `veg.abund`) — assess_trait_coverage uses rename-on-load
- [x] All `source()` paths use `here()` with new folder structure
- [x] All file I/O paths use `here()`
- [x] `00_workflow.R` kept as coordinator script (sources all scripts in order)
- [x] Final check passed: no `T`/`F`, no `write.csv`, no `View()`, no `setwd()` in active scripts
- [x] Shared `label_env_var()` moved to `functions_calanda.R`, used by 05_community and 06_species

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

### TASK 13: Style audit and cosmetic fixes ✅ DONE

Post-cleanup style audit across all 11 active scripts. Fixes applied:

- [x] `01_prepare_trait_data.R` line 228: added spaces around operators (`LCC>70|LCC<20` → `LCC > 70 | LCC < 20`)
- [x] `functions_calanda.R`: replaced `require()` → `library()` for consistency
- [x] `01_prepare_data.R`: added `library(conflicted)` + `conflict_prefer()` for `select`, `filter`, `extract` (tidyverse vs terra). Removed `dplyr::select` qualifications.
- [x] `map_rgb_results.R`: added `library(conflicted)` + `conflict_prefer()` for `select`, `filter` (tidyverse vs terra)
- [x] `02_jsdm.R`: added `library(conflicted)` + `conflict_prefer()` for `select`, `filter` (tidyverse vs sjSDM). Removed `dplyr::select` qualification.
- [x] `functions_calanda.R`: kept explicit `dplyr::` and `terra::` qualifications (appropriate for a sourced helper file that doesn't load its own libraries)

---

*Last updated: 2026-02-02*
