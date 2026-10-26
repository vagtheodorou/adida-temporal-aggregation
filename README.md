# On the benefits of temporal aggregation on forecasting and stock control

[![DOI](https://img.shields.io/badge/DOI-10.1080%2F00207543.2026.2747126-blue)](https://doi.org/10.1080/00207543.2026.2747126)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

R code accompanying the paper:

> Theodorou, E., Bougioukos, V., & Nikolopoulos, K. (2026). On the benefits of temporal aggregation on forecasting and stock control: Empirical evidence from statistical and machine learning methods. *International Journal of Production Research*. https://doi.org/10.1080/00207543.2026.2747126

## Authors

- **Evangelos Theodorou** — Forecasting and Strategy Unit, School of Electrical and Computer Engineering, National Technical University of Athens, Greece
- **Vasileios Bougioukos** — Richmond Business School, Richmond American University London, United Kingdom
- **Konstantinos Nikolopoulos** (corresponding author) — IHRR Forecasting Lab, Durham University Business School, Durham University, United Kingdom

## Overview

The study treats ADIDA (Aggregate–Disaggregate Intermittent Demand Approach) as a general aggregate–forecast–disaggregate framework. Each series is aggregated into non-overlapping buckets of size *k*, forecast at the aggregated level by a base method, and disaggregated back with equal weights. Seven statistical methods (Naive, Naive with drift, ARIMA, ES, CES, Theta, AutoTheta) and five global machine learning methods (linear model, GBM, XGBoost, LightGBM, Random Forest) are evaluated on the M3, M4 and M5 competition data.

## Repository structure

| File | Description |
|---|---|
| `00_adida_core.R` | Shared engine: ADIDA function, base forecasters, combination, error metrics (sMAPE, MASE, scaled ME/MSE, RMSSE), trend and seasonality tests |
| `01_m3_replication.R` | M3 monthly experiments (D / TD / NTD subsets; sMAPE, scaled ME, scaled MSE) |
| `02_m4_extension.R` | M4 monthly and quarterly experiments (sMAPE, MASE) |
| `03_m5_extension.R` | M5 Level 12 experiments (WRMSSE) |
| `bibliometrics/scholar_citation_analysis.R` | Bibliometric analysis of works citing the original ADIDA paper (Section 2) |

## Requirements

R (≥ 4.1) with the following packages:

```r
install.packages(c("forecast", "smooth", "zoo", "tsfeatures", "dplyr", "tidyr",
                   "data.table", "gbm", "xgboost", "lightgbm", "ranger", "Mcomp"))
```

The bibliometric script installs its own dependencies (`rvest`, `httr`, `tm`, `textstem`, `wordcloud`, among others).

## Data

No data are included in this repository. All datasets are publicly available and are downloaded by the scripts:

- **M3** — from the `Mcomp` R package
- **M4** — from the `M4comp2018` R package (release tarball)
- **M5** — from the `m5` R package

## Usage

Run each script from the repository root; every experiment script sources `00_adida_core.R`:

```r
source("01_m3_replication.R")   # M3
source("02_m4_extension.R")     # M4
source("03_m5_extension.R")     # M5
```

## Citation

If you wish to cite the paper, please use:

```bibtex
@article{theodorou2026benefits,
  title   = {On the benefits of temporal aggregation on forecasting and stock control: Empirical evidence from statistical and machine learning methods},
  author  = {Theodorou, Evangelos and Bougioukos, Vasileios and Nikolopoulos, Konstantinos},
  journal = {International Journal of Production Research},
  year    = {2026},
  doi     = {10.1080/00207543.2026.2747126}
}
```

## License

This code is released under the [MIT License](LICENSE).
