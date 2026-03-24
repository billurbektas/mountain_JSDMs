# Calanda JSDM — Analysis Code

Code for the manuscript:

> **Simultaneous disentangling processes behind species ranges and community assembly along an elevational gradient**
>
> Billur Bektas, Maximilian Pichler, Florian Hartig, Mikko Tiusanen, Camille Brioschi, Jake M. Alexander, Janneke Hille Ris Lambers

## Model framework

Species distributions and community assembly are modeled using the **sjSDM** package (Pichler & Hartig, 2021), a GPU-accelerated joint species distribution model that decomposes variance into environment, spatial, and biotic (co-distribution) components.

All model fitting and experiments are run on **Google Colab** (GPU: A100, default RAM) via the notebook:
- [`R/02_model/06_sjsdm_experiments_colab.ipynb`](https://colab.research.google.com/github/billurbektas/TP_JSDM/blob/main/Calanda_JSDM/R/02_model/06_sjsdm_experiments_colab.ipynb)

The Colab notebook runs 6 experiments: (1) full grid search over spatial form, alpha, and lambda; (2) decoupled lambda sensitivity; (3) anova sampling saturation; (4) model fit sampling saturation; (5) drop-one environmental variable; (6) 10-fold cross-validation with species AUC and site log-loss evaluation. Model fitting uses RMSprop optimizer with learning rate scheduling, early stopping, and 5000 sampling iterations. Anova decomposition uses 80,000 sampling iterations.

### Anova bug fix

The sjSDM anova function (v1.0.6) has a bug on **line 346 of `anova.R`**: `inherits(object, "spatial ")` has a trailing space, causing spatial models to never receive `spatial_formula = ~0` in the null model. This is patched at runtime in the Colab notebook (cell 13) by replacing `"spatial "` with `"spatial"` via `assignInNamespace`. See: https://github.com/TheoreticalEcology/s-jSDM/issues/172

## Data dimensions

After all filtering steps in the data preparation pipeline:

| Stage | Species | Sites |
|-------|---------|-------|
| Raw vegetation data | 723 binomial | 624 plots |
| After 5% prevalence + 70% coverage + env matching (Y matrix) | 168 | 529 |
| After AUC >= 0.7 filter | 158 | — |
| After logloss <= log(2) + ci_width < 0.1 filter | — | 503 |
| Environmental gradient analysis (script 09) | 158 species, 503 sites | |
| Functional trait analysis (script 10) | 115 species (woody + missing traits removed), 374 sites (woody-dominated plots excluded) | |

Full species-by-species and site-by-site attrition details are in `output/species_attrition.csv`, `output/site_attrition.csv`, and `output/attrition_summary.csv`.

---

## Directory structure

```
R/
├── 00_setup/
│   ├── 00_workflow.R              # Coordinator: sources all scripts in order
│   └── functions_calanda.R        # Shared helper functions
│
├── 01_data_prep/
│   ├── 01_prepare_TRY_traits.R
│   ├── 02_prepare_field_trait_data.R
│   ├── 03_merge_and_assess_traits.R
│   └── 04_prepare_climate_vegetation_data.R
│
├── 02_model/
│   ├── 05_test_cpu_jsdm.R                    # Local CPU test of model structure
│   └── 06_sjsdm_experiments_colab.ipynb       # Full model fitting + experiments (Colab GPU)
│
├── 03_analysis/
│   ├── 07.1–07.7_post_exp*.R                 # Model experiment post-analysis
│   ├── 07_modeling_experiments.R              # Consolidated experiment figures
│   ├── 08_variance_partitioning.R
│   ├── 09_environmental_gradient_analysis.R
│   └── 10_functional_post_analysis.R
│
├── 04_visualization/
│   ├── 11_map.R
│   └── 12_map_rgb_results.R
│
├── 05_validation/
│   ├── validation_checks.R        # Pipeline integrity checks
│   └── pipeline_attrition.R       # Species/site attrition tracking
│
└── archive/                       # Deprecated scripts
```

---

## Pipeline steps

### Step 1: Data preparation

| Script | What it does | Key inputs | Key outputs |
|--------|-------------|------------|-------------|
| `01_prepare_TRY_traits.R` | Fetches trait data from TRY database and ecological indicators/dispersal from floraveg.eu | `data/vegetation/*.csv`, `data/try/*.txt` | `output/try_traits_individual.csv`, `output/indicators.csv`, `output/dispersal.csv` |
| `02_prepare_field_trait_data.R` | Processes field-collected trait measurements (leaf area, LMA, LDMC, CN) with quality validation | `data/traits/*.csv` | `output/field_traits_clean.csv` |
| `03_merge_and_assess_traits.R` | Merges TRY + field traits, computes species medians and kCV, imputes missing values via random forest | `output/try_traits_individual.csv`, `output/field_traits_clean.csv`, `output/indicators.csv` | `output/traits_medians_imputed.csv`, `output/traits_kcv_imputed.csv`, `output/all_traits_individual.csv` |
| `04_prepare_climate_vegetation_data.R` | Builds environmental (X) and species (Y) matrices from remote sensing (snow, LST, ET, topography) and vegetation surveys; computes community traits and functional distinctiveness | `data/vegetation/*.csv`, `data/ecostress/*.csv`, `data/modis/*.csv`, `data/wekeo/`, `output/traits_medians_imputed.csv` | `output/data_calanda_jsdm_<date>.rds`, `output/veg_clim.csv`, `output/community_traits_*.csv`, `output/species_functional_distinctiveness.csv` |

### Step 2: Model fitting

| Script | What it does | Key inputs | Key outputs |
|--------|-------------|------------|-------------|
| `05_test_cpu_jsdm.R` | Local CPU test of the sjSDM model structure (low iterations, for debugging) | `output/data_calanda_jsdm_<date>.rds` | — |
| `06_sjsdm_experiments_colab.ipynb` | Full model fitting and 6 experiments on Google Colab GPU (A100); includes runtime anova bug fix | `output/data_calanda_jsdm_<date>.rds` (uploaded to Colab) | `output/results/cv/fold_*.rds`, `output/results/exp6_*_cv_*.csv`, `output/final_model_se_*.rds` |

### Step 3: Model experiment post-analysis

| Script | What it does | Key inputs | Key outputs |
|--------|-------------|------------|-------------|
| `07.1`–`07.7` | Individual post-analysis for each of the 6 experiments (grid search, decoupled, anova saturation, species stability, fit saturation, drop-one, CV) | `output/results/runs/`, `output/results/decoupled/`, etc. | `output/results/exp*_*.csv` |
| `07_modeling_experiments.R` | Consolidated publication-ready figures from all experiments | `output/results/exp*_*.csv` | `plot/exp_*.pdf` |

### Step 4: Post-model analysis

| Script | What it does | Key inputs | Key outputs |
|--------|-------------|------------|-------------|
| `08_variance_partitioning.R` | Aggregates species and site VP across 10-fold CV, computes CIs, produces Venn diagram and scatter plots | `output/results/cv/fold_*.rds`, `output/results/exp6_*_cv_*.csv` | `output/results/vp_species_summary_*.csv`, `output/results/vp_sites_summary_*.csv`, `plot/vp_combined.pdf` |
| `09_environmental_gradient_analysis.R` | Stepwise AIC regressions of VP components against 12 environmental predictors (sites, unweighted) and species betas (weighted) | `output/results/vp_*_summary_*.csv`, `output/data_calanda_jsdm_*.rds`, `output/final_model_se_*.rds` | `plot/env_gradient_climate.pdf`, `plot/env_gradient_other.pdf` |
| `10_functional_post_analysis.R` | Stepwise AIC regressions of VP against species trait medians, kCV, and distinctiveness (weighted); community trait means and variances (unweighted) | `output/results/vp_*_summary_*.csv`, `output/traits_*.csv`, `output/community_traits_unweighted_imputed.csv`, `output/species_functional_distinctiveness.csv` | `plot/functional_gradient.pdf` |

### Step 5: Visualization

| Script | What it does | Key inputs | Key outputs |
|--------|-------------|------------|-------------|
| `11_map.R` | Study area map of Switzerland with Calanda survey points colored by altitude | `output/starter_data_*.RData`, `output/calanda_mask.shp` | `plot/map_calanda.pdf` |
| `12_map_rgb_results.R` | 3D rayshader and 2D maps with variance components mapped to RGB channels on the DEM | `output/starter_data_*.RData`, `output/calanda_mask.shp`, `output/metrics/dem.tif` | `plot/calanda_rgb_3d.png`, `plot/calanda_rgb_2d_map.pdf` |

### Step 6: Validation (optional)

| Script | What it does | Key inputs | Key outputs |
|--------|-------------|------------|-------------|
| `validation_checks.R` | Pipeline integrity checks: data alignment, join correctness, VIF, residual normality, heteroscedasticity | All output files | Console diagnostics |
| `pipeline_attrition.R` | Tracks species and site attrition from raw data through all filtering stages with reasons | `data/vegetation/*.csv`, all output files | `output/species_attrition.csv`, `output/site_attrition.csv`, `output/attrition_summary.csv` |

---

## Execution

Run `00_workflow.R` to execute the full pipeline. Toggle steps with `TRUE`/`FALSE` at the top of the script:

```r
run_data_prep     = TRUE   # Step 1
run_model         = FALSE  # Step 2 (GPU, use Colab instead)
run_experiments   = FALSE  # Step 3
run_analysis      = TRUE   # Step 4
run_visualization = TRUE   # Step 5
run_validation    = FALSE  # Step 6
```

Or source individual scripts directly.

---

## Coding conventions

| Rule | Convention |
|------|-----------|
| Assignment | `=` (not `<-`) |
| Naming | `snake_case` |
| Paths | `here::here("Calanda_JSDM", ...)` |
| CSV I/O | `read_csv()` / `write_csv()` |
| Plots | `theme_bw()` default; colors: env=`#81caf3`, codist=`#00bd89`, spa=`#d00000` |
| Libraries | Each script loads its own |

---

*Last updated: 2026-03-17*
