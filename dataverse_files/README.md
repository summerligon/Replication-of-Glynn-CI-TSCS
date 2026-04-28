# Replication: Blackwell & Glynn (2018) APSR

Replication code and data for:

> Blackwell, M. & Glynn, A. N. (2018). How to Make Causal Inferences with Time-Series Cross-Sectional Data under Selection on Observables. *American Political Science Review*, 112(4), 1067–1082.

---

## Repository Structure

```
.
├── .devcontainer/
│   ├── devcontainer.json       # GitHub Codespaces environment config
│   └── install_packages.R      # Auto-installs all required R packages
├── code/
│   ├── panel-utils.R           # Helper functions (IPTW, panel lags, balance)
│   ├── burgoon.R               # Empirical reanalysis: Burgoon (2006)
│   ├── swank-steinmo.R         # Appendix replication: Swank & Steinmo
│   └── causal-tscs-sim.R       # Monte Carlo simulations + figures
├── data/
│   ├── burgoon.csv             # Data for burgoon.R
│   └── swank.dta               # Stata data for swank-steinmo.R
├── output/                     # Tables and figures written here
├── run_all.R                   # Master script to run all analyses
└── README.md
```

---

## How to Run in GitHub Codespaces

1. **Open the repo** in GitHub and click **Code → Codespaces → Create codespace on main**
2. Wait for the container to build (~3–5 min first time). It installs R 4.3.2 + RStudio Server + all packages automatically.
3. Codespaces will open a browser tab with **RStudio Server** on port 8787.
4. In RStudio, open `run_all.R` and source it — or run each script in `code/` individually.
5. All outputs (figures as PDFs/PNGs, tables as `.tex`) are written to `output/`.

> **Tip:** You can also open a Terminal in Codespaces and run `Rscript run_all.R` from the repo root.

---

## Scripts & Outputs

| Script | What it produces |
|---|---|
| `burgoon.R` | Tables 2–3 and Figures 3–5 (Burgoon welfare/terrorism reanalysis) |
| `swank-steinmo.R` | Appendix Table A2 and Figures A1–A2 (Swank & Steinmo replication) |
| `causal-tscs-sim.R` | Figures 1–2 and Appendix simulation figures |

---

## R Environment

| Item | Value |
|---|---|
| R version | 4.3.2 (container) — original paper used 3.4.3 |
| Key packages | `plm`, `geepack`, `mgcv`, `sandwich`, `lmtest`, `stargazer`, `ggplot2` |
| OS (original) | macOS 13.4 |
| OS (Codespaces) | Ubuntu 22.04 |

> **Note on R version:** The original analyses used R 3.4.3. Minor numerical differences in output are expected when using a newer R version, but all substantive results should replicate.

---

## Citation

```
@article{blackwell2018causal,
  title={How to Make Causal Inferences with Time-Series Cross-Sectional Data under Selection on Observables},
  author={Blackwell, Matthew and Glynn, Adam N.},
  journal={American Political Science Review},
  volume={112},
  number={4},
  pages={1067--1082},
  year={2018}
}
```
