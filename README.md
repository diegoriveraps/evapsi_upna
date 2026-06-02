# EVAPSI Project

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue.svg)](https://www.r-project.org/)

## Overview

This repository contains the R analysis code for the **EVAPSI** project, a psychometric research program developing and validating instruments related to academic stress and psychological well-being in university students.

The project spans multiple manuscripts, each in its own self-contained folder with its own documentation, data requirements, and script(s).

---

## Repository Structure

```
evapsi/
├── data/                          # Shared or raw data files (not tracked by git)
├── manuscript_1_stress/           # Psychometric validation of the EVAPSI stress scale
│   ├── evapsi_stress_final.R
│   └── README.md
├── manuscript_2/                  # [Next manuscript — title TBD]
│   ├── script.R
│   └── README.md
├── manuscript_3/                  # [Next manuscript — title TBD]
│   ├── script.R
│   └── README.md
├── .gitignore
├── LICENSE
└── README.md                      ← you are here
```

---

## Manuscripts

| # | Folder | Topic | Status |
|---|--------|-------|--------|
| 1 | [`manuscript_1_stress/`](manuscript_1_stress/) | Psychometric validation & measurement invariance of the EVAPSI stress scale | ✅ Complete |
| 2 | [`manuscript_2/`](manuscript_2/) | TBD | 🔄 In progress |
| 3 | [`manuscript_3/`](manuscript_3/) | TBD | 🔄 In progress |

---

## General Requirements

All scripts are written in **R** and rely on the following packages:

| Package | Purpose |
|---------|---------|
| `lavaan` | Structural equation modeling (CFA) |
| `semTools` | Measurement invariance, reliability |
| `EGAnet` | Exploratory Graph Analysis |
| `psych` | EFA, polychoric correlations |
| `MVN` | Multivariate normality testing |
| `dplyr` | Data wrangling |
| `table1` | Descriptive statistics tables |
| `knitr` | Report generation |

Install all at once:

```r
install.packages(c("lavaan", "semTools", "EGAnet", "psych",
                   "MVN", "dplyr", "table1", "knitr", "readxl"))
```

---

## Data

Raw data files are **not included** in this repository to protect participant privacy.  
Data access may be requested by contacting the corresponding author.

Place any required data files in the `data/` folder or in the relevant manuscript subfolder before running the scripts. Each manuscript README specifies the expected file names.

---

## Citation

If you use this code, please cite the corresponding manuscript(s). Citation details will be added upon publication.

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## Contact

For questions about the code or data, please open an issue or contact the project maintainers.
