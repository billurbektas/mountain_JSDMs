# Calanda JSDM - R Scripts Workflow

## Overview

Joint Species Distribution Model (JSDM) analysis of alpine vegetation in the Calanda region, Switzerland. The workflow processes vegetation survey data, environmental variables, and functional traits to fit spatial JSDMs and analyze variance components.

Model outputs used for analysis come from `results_from_Max/` (authoritative model results).

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
│   ├── 03_explain_variation.R     # Community-level variance decomposition
│   ├── 05_community_postJSDM.R    # Environmental regressions on variance components
│   └── 06_species_postJSDM.R      # Species-level coefficient analysis
│
├── 04_visualization/
│   ├── map.R                      # Geographic maps of study area
│   └── map_rgb_results.R          # 3D DEM with RGB variance coloring
│
└── archive/                       # Deprecated / external scripts
    ├── 01_prepare_vegetation_data.R   # Superseded by 01_prepare_data.R
    ├── 03_explain_variation_OLD.R     # Old version
    ├── 03_explain_variation_mvt.R     # Old multivariate version
    ├── 04_nmds.R                      # NMDS ordination (archived)
    ├── 04_predict_transplant_jsdm.R   # TransPlant predictions (archived)
    ├── explain_variation_species.R    # Species variance perspective (archived)
    ├── helper_functions.R             # Unused (content already in functions_calanda.R)
    ├── get_data_chelsa.R              # CHELSA climate projections (not used)
    ├── Cai_2024.R                     # External project (aquatic species)
    ├── taken_out.R                    # Code snippets
    └── 2025_TRAITS_Cleaning_20251023.R # Upstream trait cleaning
```

---

## Execution Order

```r
# 1. Setup
source("R/00_setup/00_workflow.R")           # Libraries, paths, functions

# 2. Data preparation
source("R/01_data_prep/01_prepare_data.R")          # -> output/starter_data_25.04.25.RData
source("R/01_data_prep/01_prepare_trait_data.R")     # -> output/species_trait_summary.csv
source("R/01_data_prep/assess_trait_coverage.R")     # -> output/community_trait_coverage.csv (QA)

# 3. Model fitting (GPU required) -- or use results_from_Max/
source("R/02_model/02_jsdm.R")                      # -> output/model_sjsdm_calanda_*.RDS

# 4. Post-model analysis (uses results_from_Max/)
source("R/03_analysis/03_explain_variation.R")       # -> plot/calanda_jsdm_outputs.pdf
source("R/03_analysis/05_community_postJSDM.R")      # -> plot/community_regressions.pdf
source("R/03_analysis/06_species_postJSDM.R")        # -> plot/species_regressions.pdf

# 5. Visualization
source("R/04_visualization/map.R")                   # -> plot/map_calanda.pdf
source("R/04_visualization/map_rgb_results.R")       # -> 3D RGB elevation maps
```

---

## Key Data Flow

```
Raw data (data/)
    │
    ├── Vegetation surveys ──► 01_prepare_data.R ──► starter_data_25.04.25.RData
    ├── Topography (SwissALTI3D)  ──┘                         │
    ├── Snow (Copernicus)         ──┘                         │
    ├── LST (ECOSTRESS)           ──┘                         │
    └── ET (MODIS)                ──┘                         │
                                                              ▼
    Trait data (data/traits/) ──► 01_prepare_trait_data.R ──► species_trait_summary.csv
                                                              │
                                                              ▼
                                          02_jsdm.R / results_from_Max/
                                                              │
                                                              ▼
                                          03_analysis/ scripts ──► plot/ PDFs
```

---

*Last updated: 2026-02-02*
