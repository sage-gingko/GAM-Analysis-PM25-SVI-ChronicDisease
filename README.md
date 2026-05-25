# GAM-Analysis-PM25-SVI-ChronicDisease

R code for Generalized Additive Model (GAM) analysis of $PM_{2.5}$ exposure, social vulnerability, and chronic disease prevalence across 498 U.S. cities (2016–2022).

## Data Availability

The analytical dataset used in this study is too large to host directly on GitHub. It is permanently archived and publicly accessible on Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20383589.svg)](https://doi.org/10.5281/zenodo.20383589)

### Data Sources Included in the Merged Dataset:
* **$PM2.5 Estimates:** Atmospheric Composition Analysis Group (Washington University in St. Louis)
* **Social Vulnerability Index (SVI):** CDC/ATSDR
* **Chronic Disease Prevalence:** CDC PLACES Dataset (Asthma, Diabetes, COPD, Coronary Heart Disease, and Cancer)

---

## Repository Structure

```text
├── data/
│   └── [Place Zenodo CSV file here]
├── .gitignore
├── GAMAnalysis_PM25_SVI_ChronicDisease.R
├── LICENSE
└── README.md
