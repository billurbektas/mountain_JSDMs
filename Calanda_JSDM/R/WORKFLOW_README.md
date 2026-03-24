# Calanda JSDM — Workflow Details

See the [main README](../../README.md) for an overview. This document contains additional workflow details.

## Pipeline outputs

### From script 08 (variance partitioning):
- `output/results/vp_species_summary_*.csv` — per-species VP with name-based identification
- `output/results/vp_sites_summary_*.csv` — per-site VP with name-based identification
- `output/results/vp_anova_summary_*.csv` — landscape-level anova decomposition
- `plot/vp_combined.pdf` — Venn diagram + violins + scatter plots

### From script 09 (environmental gradient):
- `output/env_gradient_coefficients_*.csv` — all model coefficients, SE, t, p
- `output/env_gradient_aic_steps_*.csv` — AIC step history
- `output/env_gradient_model_summary_*.csv` — R², F-stat, df, AIC per model
- `plot/env_gradient_climate.pdf` — summer temp, FDD, ET
- `plot/env_gradient_cover.pdf` — rocks, trees, shrubs, soil depth
- `plot/env_gradient_topography.pdf` — slope, TPI, flow direction
- `plot/env_gradient_indices.pdf` — nutrient, disturbance

### From script 10 (functional traits):
- `output/functional_coefficients_*.csv` — all model coefficients
- `output/functional_aic_steps_*.csv` — AIC step history
- `output/functional_model_summary_*.csv` — model-level stats
- `plot/functional_gradient.pdf` — vegetative height, LNC, seed mass
- `plot/functional_gradient_other.pdf` — LDMC

### From script 11 (map):
- `plot/map_calanda.pdf` — study area with RGB climate triangle, DEM relief

### From validation scripts:
- `output/species_attrition.csv` — per-species removal stage and reason
- `output/site_attrition.csv` — per-site removal stage
- `output/attrition_summary.csv` — combined stage table
- `plot/vp_diagnostic_discard_vs_proportional*.pdf` — allocation method comparison
- `plot/vp_diagnostic_raw_shared_fractions*.pdf` — raw fraction distributions

## Key methodological notes

### Variance partitioning levels
The Venn diagram (landscape) and violins (per-unit) use different aggregation. See `docs/variance_partitioning_notes.md` for a detailed explanation of why they differ.

### Site-level regressions
Community-level regressions in scripts 09 and 10 are **unweighted** (sites with ci_width >= 0.1 are excluded instead). Species-level regressions remain weighted by 1/ci_width².

### Name-based joins
All species and site matching uses name-based joins (not positional indices) to prevent silent data misalignment across folds.

---

*Last updated: 2026-03-24*
