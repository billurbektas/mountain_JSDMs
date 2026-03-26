# Mountain JSDMs

Joint Species Distribution Models for alpine vegetation along elevational gradients.

## Calanda JSDM

Code for the manuscript:

> **Simultaneous disentangling processes behind species ranges and community assembly along an elevational gradient**
>
> Billur Bektas, Maximilian Pichler, Florian Hartig, Mikko Tiusanen, Camille Brioschi, Jake M. Alexander, Janneke Hille Ris Lambers

### Model framework

Species distributions and community assembly are modeled using the **sjSDM** package (Pichler & Hartig, 2021), a GPU-accelerated joint species distribution model that decomposes variance into environment, spatial, and species association components.

All model fitting and experiments are run on **Google Colab** (GPU: A100, default RAM) via the notebook:
- [`Calanda_JSDM/R/02_model/06_sjsdm_experiments_colab.ipynb`](https://colab.research.google.com/github/billurbektas/TP_JSDM/blob/main/Calanda_JSDM/R/02_model/06_sjsdm_experiments_colab.ipynb)

The Colab notebook runs 6 experiments: (1) full grid search over spatial form, alpha, and lambda; (2) decoupled lambda sensitivity; (3) anova sampling saturation; (4) model fit sampling saturation; (5) drop-one environmental variable; (6) 10-fold cross-validation with species AUC and site log-loss evaluation.

#### Anova bug fix

The sjSDM anova function (v1.0.6) has a bug on **line 346 of `anova.R`**: `inherits(object, "spatial ")` has a trailing space, causing spatial models to never receive `spatial_formula = ~0` in the null model. This is patched at runtime in the Colab notebook (cell 13). See: https://github.com/TheoreticalEcology/s-jSDM/issues/172

### Data dimensions

| Stage | Species | Sites |
|-------|---------|-------|
| Raw vegetation data | 723 binomial | 624 plots |
| After 5% prevalence + 70% coverage + env matching (Y matrix) | 168 | 529 |
| After AUC >= 0.7 filter | 158 | -- |
| After logloss <= log(2) + ci_width < 0.1 filter | -- | 503 |
| Environmental gradient analysis (script 09) | 158 species, 503 sites | |
| Functional trait analysis (script 10) | 115 species, 374 sites | |

---

### Directory structure

```
Calanda_JSDM/
├── R/
│   ├── 00_setup/
│   │   ├── 00_workflow.R              # Coordinator with TRUE/FALSE toggles
│   │   └── functions_calanda.R        # Shared helper functions
│   │
│   ├── 01_data_prep/
│   │   ├── 01_prepare_TRY_traits.R
│   │   ├── 02_prepare_field_trait_data.R
│   │   ├── 03_merge_and_assess_traits.R
│   │   └── 04_prepare_climate_vegetation_data.R
│   │
│   ├── 02_model/
│   │   ├── 05_test_cpu_jsdm.R                    # Local CPU test
│   │   └── 06_sjsdm_experiments_colab.ipynb       # Colab GPU fitting + experiments
│   │
│   ├── 03_analysis/
│   │   ├── 07.1–07.7_post_exp*.R                 # Model experiment post-analysis
│   │   ├── 07_modeling_experiments.R              # Consolidated experiment figures
│   │   ├── 08_variance_partitioning.R             # Venn diagram, VP violins, scatter plots
│   │   ├── 09_environmental_gradient_analysis.R   # VP vs env predictors and species betas
│   │   └── 10_functional_post_analysis.R          # VP vs species traits and community traits
│   │
│   ├── 04_visualization/
│   │   └── 11_map.R                               # Study area map with RGB climate triangle
│   │
│   ├── 05_validation/
│   │   ├── validation_checks.R        # Pipeline integrity checks
│   │   ├── pipeline_attrition.R       # Species/site attrition tracking
│   │   └── vp_diagnostic_figures.R    # Discard vs proportional VP diagnostics
│   │
│   └── archive/                       # Deprecated scripts
│
├── plot/                              # All figures (PDFs)
├── output/                            # Data outputs (gitignored)
├── data/                              # Raw data (gitignored)
└── docs/                              # Documentation (gitignored)
```

---

### Pipeline steps

#### Step 1: Data preparation

| Script | What it does |
|--------|-------------|
| `01_prepare_TRY_traits.R` | Fetches trait data from TRY database and ecological indicators/dispersal from floraveg.eu |
| `02_prepare_field_trait_data.R` | Processes field-collected trait measurements (leaf area, LMA, LDMC, CN) |
| `03_merge_and_assess_traits.R` | Merges TRY + field traits, computes species medians and kCV, imputes missing values |
| `04_prepare_climate_vegetation_data.R` | Builds X and Y matrices from remote sensing and vegetation surveys; computes community traits and functional distinctiveness |

#### Step 2: Model fitting

| Script | What it does |
|--------|-------------|
| `05_test_cpu_jsdm.R` | Local CPU test of model structure |
| `06_sjsdm_experiments_colab.ipynb` | Full model fitting + 6 experiments on Colab GPU (A100) |

#### Step 3: Model experiment post-analysis

| Script | What it does |
|--------|-------------|
| `07.1`--`07.7` | Individual post-analysis for each experiment |
| `07_modeling_experiments.R` | Consolidated publication-ready figures |

#### Step 4: Post-model analysis

| Script | What it does |
|--------|-------------|
| `08_variance_partitioning.R` | Aggregates VP across 10-fold CV, Venn diagram with mean R² values, species/site violins and scatter plots |
| `09_environmental_gradient_analysis.R` | AIC stepwise regressions: community VP vs 12 env predictors (unweighted); species VP vs species betas (weighted) |
| `10_functional_post_analysis.R` | AIC stepwise regressions: species VP vs trait medians, kCV, distinctiveness (weighted); community VP vs trait means and variances (unweighted) |

#### Step 5: Visualization

| Script | What it does |
|--------|-------------|
| `11_map.R` | Study area map: greyscale DEM, RGB climate triangle (summer temp / ET / FDD), open/closed habitat shapes, species richness sizes |

#### Step 6: Validation (optional)

| Script | What it does |
|--------|-------------|
| `validation_checks.R` | Data alignment, join correctness, VIF, residual normality, heteroscedasticity checks |
| `pipeline_attrition.R` | Species and site attrition from raw data through all filters, with reasons |
| `vp_diagnostic_figures.R` | Discard vs proportional allocation comparison; raw shared fraction distributions |

---

### Execution

Run `00_workflow.R` to execute the full pipeline. Toggle steps at the top:

```r
run_data_prep     = TRUE   # Step 1
run_model         = FALSE  # Step 2 (GPU, use Colab)
run_experiments   = FALSE  # Step 3
run_analysis      = TRUE   # Step 4
run_visualization = TRUE   # Step 5
run_validation    = FALSE  # Step 6
```

---

*Last updated: 2026-03-24*
