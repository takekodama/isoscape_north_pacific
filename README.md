# isoscape_north_pacific
The code and data for creating the isoscape of the North Pacific.

# README: North Pacific Isoscape & Trophic Transfer Pipeline

This repository contains the main R analysis script (`Isoscape_North_Pacific_final.R`) for modeling particulate organic matter (POM) nitrogen isotope ratios ($\delta^{15}\text{N}_{\text{POM}}$) and small pelagic fish muscle isotope ratios ($\delta^{15}\text{N}_{\text{fish}}$) across the North Pacific. The pipeline utilizes a two-stage hierarchical Bayesian framework to quantify oceanographic drivers of baseline isotopic variation and direct/indirect environmental effects on trophic transfer.

---

## 1. Pipeline Overview

The analysis workflow consists of the following sequential components:

* **Data Processing & Integration**: Cleans primary POM isotope database records (`SI_POM_Database_230511.xlsx`, `SI_POM_saury.xlsx`, `SI_POM_Database_211210.csv`), fish stable isotope data (`fish_SI.csv`), environmental datasets (`temp`, `sal`, `CHL`, `no3`), and published nitrate isotope records[cite: 1].
* **Outlier Removal & Regional Classification**: Detects and removes $\delta^{13}\text{C}$ outliers via Grubbs-like t-test criteria (`remove_outliers`) and categorizes stations into distinct oceanic regions (e.g., Kuroshio, Kuroshio Extension, Sea of Japan, Sea of Okhotsk)[cite: 1].
* **Stage 1: Bayesian GAM for POM Baseline ($\delta^{15}\text{N}_{\text{POM}}$)**: Models surface POM $\delta^{15}\text{N}$ using `brms` with a spatial smoother `s(long, lat)` and oceanographic covariates (`temp`, `sal`, `nitrate`, `log10_chla`) via the CmdStan backend[cite: 1]. Model comparison is performed using WAIC / LOO-CV[cite: 1].
* **Stage 2: Hierarchical Bayesian Model for Fish ($\delta^{15}\text{N}_{\text{fish}}$)**: Connects predicted latent baseline estimates (`d15NPOM_model1`) to fish muscle $\delta^{15}\text{N}$ alongside fish body length (BL), oceanographic variables, and genus-level random intercepts via Stan (`model_stage3_final.stan`)[cite: 1].
* **Spatial Variance Decomposition**: Calculates posterior proportions of site-level spatial variance attributable to baseline isotopic variation versus direct environmental effects (chlorophyll-$a$, nitrate, temperature, salinity)[cite: 1].
* **Publication-Ready Visualization**: Generates manuscript figures (Figure 01, Figure 02) and supporting figures (Figure S01–S04) incorporating IDW spatial interpolation, degree-based buffer masking (2.5°), and land polygon masking[cite: 1].

---

## 2. Directory Structure

* `Isoscape_North_Pacific_3.R`: Main R execution script[cite: 1].
* `data/`: Input Excel (`.xlsx`), CSV (`.csv`), and NetCDF climatology (`.nc`) datasets[cite: 1].
* `models/`: Stan model specification files (`model_stage3_final.stan`)[cite: 1].
* `output/`: Generated publication-quality PDF figures (Figure 01a–g, Figure 02, Figures S01–S04)[cite: 1].

---

## 3. Required R Packages

* **Data Manipulation & Core**: `tidyverse`, `readxl`, `lubridate`, `magrittr`, `data.table`[cite: 1]
* **Spatial & Oceanographic Processing**: `sf`, `terra`, `akima`, `mapdata`, `tidync`, `metR`, `gstat`, `stars`, `rnaturalearth`[cite: 1]
* **Bayesian Modeling & Statistics**: `brms`, `cmdstanr`, `posterior`, `mgcv`, `mgcViz`, `MuMIn`[cite: 1]
* **Visualization & Layout**: `ggplot2`, `patchwork`, `ggpubr`, `ggrastr`, `ggsci`, `MetBrewer`[cite: 1]

---

## 4. Role of Gemini (AI Collaborator)

Gemini acted as an AI collaborator throughout script development and pipeline optimization:

* **Transition to Bayesian GAM**: Assisted in migrating Stage 1 modeling from frequentist `mgcv::gam` to `brms::brm` to enable full posterior propagation and WAIC/LOO-CV comparison[cite: 1].
* **Spatial Geometry & Masking Logic**: Designed degree-based spatial buffering (`st_buffer`) and 180°-shifted land masking (`rnaturalearth`) to prevent extrapolation artifacts in unobserved ocean regions and landmasses[cite: 1].
* **Two-Stage Pipeline Coupling**: Structured the seamless integration linking Stage 1 `brms` latent baseline predictions (`fitted(..., re_formula = NA)`) to Stage 2 CmdStan hierarchical sampling[cite: 1].
* **Academic Documentation**: Standardized code formatting, unified English comments, and assisted in drafting manuscript-level Methods and figure captions[cite: 1].
