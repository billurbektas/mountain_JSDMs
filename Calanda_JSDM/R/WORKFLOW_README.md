# Calanda JSDM - R Scripts Workflow

## Overview

This repository contains scripts for a Joint Species Distribution Model (JSDM) analysis of alpine vegetation in the Calanda region, Switzerland. The workflow processes vegetation survey data, environmental variables, and functional traits to fit spatial JSDMs and analyze variance components.

---

## Workflow Diagram

```
┌─────────────────────────────────────────────────────────────────────┐
│                     DATA PREPARATION LAYER                          │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  00_workflow.R ──────► Coordinator (loads libraries, sets paths)    │
│       │                                                             │
│       └──► functions_calanda.R (helper functions)                   │
│                                                                     │
│  01_prepare_data.R ──────► Main environmental data pipeline         │
│       │                    (topography, snow, LST, ET)              │
│       │                                                             │
│       └──► OUTPUT: starter_data_25.04.25.RData                      │
│                                                                     │
│  01_prepare_trait_data.R ──────► Trait data processing (parallel)   │
│       │                                                             │
│       └──► OUTPUT: species_trait_summary.csv                        │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────────┐
│                       MODEL FITTING LAYER                           │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  02_jsdm.R ──────► Fits sjSDM model (GPU-accelerated)               │
│       │                                                             │
│       └──► OUTPUTS: model_sjsdm_calanda_*.RDS                       │
│                     an_sjsdm_calanda_*.RDS (ANOVA)                  │
│                     R2_sjsdm_calanda_*.RDS                          │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────────┐
│                    POST-MODEL ANALYSIS LAYER                        │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  03_explain_variation.R ──────► Community-level variance analysis   │
│       └──► Ternary plots, CWM traits, PCA                           │
│                                                                     │
│  05_community_postJSDM.R ──────► Environmental effects on variance  │
│       └──► Community regression plots                               │
│                                                                     │
│  06_species_postJSDM.R ──────► Species-specific responses           │
│       └──► Species regression plots                                 │
│                                                                     │
│  explain_variation_species.R ──────► Detailed species perspective   │
│       └──► PCA biplots with RGB variance coloring                   │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────────┐
│                    SPECIALIZED ANALYSES                             │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  04_nmds.R ──────► NMDS ordination of species communities           │
│                                                                     │
│  04_predict_transplant_jsdm.R ──────► Predictions for TransPlant    │
│       └──► 8 scenarios: high/low elevation, marginal/conditional    │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────────┐
│                       VISUALIZATION LAYER                           │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  map.R ──────► Geographic maps of study area                        │
│                                                                     │
│  map_rgb_results.R ──────► 3D DEM with RGB variance coloring        │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Script Inventory

### ACTIVELY USED (Core Pipeline)

| Script | Purpose | Key Outputs |
|--------|---------|-------------|
| `00_workflow.R` | Main coordinator, loads libraries and paths | Coordinates CSVs |
| `01_prepare_data.R` | Environmental data integration (topography, snow, LST, ET, imputation) | `starter_data_25.04.25.RData` |
| `01_prepare_trait_data.R` | Trait data processing (N/C isotopes, LMA, LDMC, etc.) | `species_trait_summary.csv` |
| `02_jsdm.R` | Fits spatial JSDM with sjSDM package | Model RDS files |
| `03_explain_variation.R` | Community-level variance decomposition | Ternary plots, PDFs |
| `04_nmds.R` | NMDS ordination | `nmds.pdf` |
| `04_predict_transplant_jsdm.R` | Predictions for transplant experiments | Multiple scenario outputs |
| `05_community_postJSDM.R` | Environmental regressions on variance | `community_regressions.pdf` |
| `06_species_postJSDM.R` | Species-level coefficient analysis | `species_regressions.pdf` |
| `explain_variation_species.R` | Detailed species variance perspective | PCA biplots |
| `map.R` | Geographic mapping | `map_calanda.pdf` |
| `map_rgb_results.R` | 3D RGB variance visualization | 3D elevation maps |
| `functions_calanda.R` | Helper functions library | (sourced by other scripts) |

### POTENTIALLY OBSOLETE / NEEDS REVIEW

| Script | Status | Notes |
|--------|--------|-------|
| `01_prepare_vegetation_data.R` | **REVIEW** | May be older version of `01_prepare_data.R` - check if superseded |
| `helper_functions.R` | **REVIEW** | Contains `match_clim_var()` and `aucplot()` - verify if still used |
| `taken_out.R` | **UNCLEAR** | Contains ternary plot code snippets - may be extracted from other scripts |
| `get_data_chelsa.R` | **INACTIVE** | CHELSA climate projections - not part of current workflow |
| `assess_trait_coverage.R` | **QA/SUPPORT** | Trait coverage assessment - run occasionally, not core pipeline |

### EXTERNAL / DIFFERENT PROJECT

| Script | Status | Notes |
|--------|--------|-------|
| `Cai_2024.R` | **EXTERNAL** | sjSDM on aquatic (FHT) data - different project entirely |

### ARCHIVED

| Script | Location | Notes |
|--------|----------|-------|
| `03_explain_variation.R` | `archieve/` | Old version, superseded by main `03_explain_variation.R` |
| `03_explain_variation_mvt.R` | `archieve/` | Old multivariate version |

---

## Key Data Dependencies

### Input Data Required
- `data/2024_CAPHE_SpeDis_CleanData_20240214.csv` - Raw vegetation surveys
- `data/ch.swisstopo.swissalti3d*.csv` - Topographic data
- Copernicus snow data (various GeoTIFFs)
- ECOSTRESS LST data
- MODIS ET data
- Trait Excel files (`N+C_Samples*.xls`, `Trial*.xls`)
- `data/2025_TRAITS_clean_20251023.csv`

### Intermediate Outputs
- `output/veg.clim.csv` - Environmental variables with imputation
- `output/starter_data_25.04.25.RData` - Complete prepared data for JSDM
- `output/species_trait_summary.csv` - Species-level trait means and variability

### Model Outputs
- `output/model_sjsdm_calanda_*.RDS` - Fitted JSDM model
- `output/an_sjsdm_calanda_*.RDS` - ANOVA results
- `output/R2_sjsdm_calanda_*.RDS` - Model performance
- `results_from_Max/` - External model results (from collaborator)

---

## Suggested Reorganization

### Option A: Keep Current Structure (Clean Up)
1. Move `Cai_2024.R` out of R/ (different project)
2. Move `taken_out.R` to archive
3. Delete or confirm `01_prepare_vegetation_data.R` if superseded
4. Consolidate `helper_functions.R` into `functions_calanda.R`

### Option B: Reorganize by Stage
```
R/
├── 00_setup/
│   ├── 00_workflow.R
│   └── functions_calanda.R
├── 01_data_prep/
│   ├── 01_prepare_data.R
│   └── 01_prepare_trait_data.R
├── 02_model/
│   └── 02_jsdm.R
├── 03_analysis/
│   ├── 03_explain_variation.R
│   ├── 04_nmds.R
│   ├── 04_predict_transplant_jsdm.R
│   ├── 05_community_postJSDM.R
│   ├── 06_species_postJSDM.R
│   └── explain_variation_species.R
├── 04_visualization/
│   ├── map.R
│   └── map_rgb_results.R
└── archive/
    └── (obsolete scripts)
```

### Option C: Prefix-Based Naming (Minimal Change)
Rename scripts for clearer ordering:
- `03a_explain_variation_community.R`
- `03b_explain_variation_species.R`
- `04a_nmds.R`
- `04b_predict_transplant.R`

---

## Questions to Resolve

1. **01_prepare_vegetation_data.R vs 01_prepare_data.R**: Are both needed? Which is current?
2. **helper_functions.R**: Is this still sourced anywhere? Can it be merged?
3. **results_from_Max/**: Are these the authoritative model results, or should 02_jsdm.R be re-run?
4. **get_data_chelsa.R**: Is future climate projection work planned, or should this be archived?
5. **assess_trait_coverage.R**: Keep as QA tool or integrate into main pipeline?

---

## Execution Order (Current Workflow)

```bash
# 1. Setup and data preparation
source("R/00_workflow.R")
source("R/01_prepare_data.R")
source("R/01_prepare_trait_data.R")

# 2. Model fitting (requires GPU)
source("R/02_jsdm.R")

# 3. Post-model analysis
source("R/03_explain_variation.R")
source("R/05_community_postJSDM.R")
source("R/06_species_postJSDM.R")
source("R/explain_variation_species.R")

# 4. Specialized analyses
source("R/04_nmds.R")
source("R/04_predict_transplant_jsdm.R")

# 5. Visualization
source("R/map.R")
source("R/map_rgb_results.R")
```

---

*Last updated: 2026-01-27*
