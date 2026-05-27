# Salt × CNT priming in Barley and Cowpea

[![Live report](https://img.shields.io/badge/Live-report-blue?logo=github)](https://uta666xyz.github.io/plant-salt-cnt-stress/)

Portfolio piece — a self-contained demonstration of an analysis pipeline
for a factorial plant stress-physiology experiment, built end-to-end on a
fully synthetic dataset.

The example scenario is a two-factor design — a NaCl salinity gradient
crossed with a seed-priming treatment (labelled **CNT**) — measured in two
crop species, **barley** and **cowpea**. The data are generated from
scratch (see `R/01_simulate_data.R`) with response patterns set to
qualitatively mirror the published seed-priming literature, so the repo
showcases the full workflow — simulation → multivariate analysis →
reporting — without relying on any real experimental data.

## What's in the report

The report is organised into five analysis blocks:

1. **Percentage germination** — daily cumulative germination % over
   Day 0–4, NaCl 0–300 mM, both species, with and without CNT priming.
2. **Univariate seedling growth** — violin + boxplot composite for
   shoot length, root length and fresh weight, faceted by species ×
   NaCl × CNT.
3. **Three-way ANOVA across all responses** — `species × NaCl ×
   CNT` ANOVA table covering all 10 measured variables (4 growth + 6
   stress-response biochemistry).
4. **Growth × biochemistry correlation heatmap** — Pearson correlations
   between the four growth indices and the six stress markers,
   confirming the priming-hypothesis pattern.
5. **Per-species PCA (barley and cowpea)** — two parallel PCA
   layouts:
   * **growth-index PCA** on shoot length, root length, shoot-root axis
     weight (correlation circle + treatment biplot per species);
   * **biochemistry PCA** on SOD, CAT, APX, MDA, proline, H₂O₂
     (correlation circle + treatment biplot per species).
   Treatment biplots use 68 % confidence ellipses; CNT+ ellipses are
   drawn with dashed outlines so the priming offset is visible at a
   glance.

## What's in the repo

* `R/01_simulate_data.R` — generates `data/synthetic.csv` (2 species ×
  4 NaCl × 2 CNT × 6 reps = 96 plants; 10 response variables: 4 growth
  + SOD, CAT, APX, MDA, proline, H₂O₂).
* `R/02_simulate_germination.R` — generates `data/germination.csv`
  (2 species × 6 NaCl × 2 CNT × 4 dishes × 25 seeds, daily counts
  Day 0–4).
* `analysis.Rmd` — knits to `analysis.html` (self-contained).
* `figures/` — PNG exports of every figure.
* `index.html` — redirects the GitHub Pages root to `analysis.html`.

## Reproducing

```r
# from the repo root
install.packages(c("ggplot2", "dplyr", "tidyr", "forcats",
                   "FactoMineR", "factoextra", "gridExtra",
                   "broom", "knitr", "RColorBrewer", "rmarkdown"))
source("R/01_simulate_data.R")            # data/synthetic.csv
source("R/02_simulate_germination.R")     # data/germination.csv
rmarkdown::render("analysis.Rmd")         # analysis.html + figures/*.png
```

## Disclaimer

The numbers are **simulated**. The simulation is not calibrated against
any specific paper and should not be cited as evidence for any
biological claim. The repo demonstrates **the analysis workflow**, not
a finding.
