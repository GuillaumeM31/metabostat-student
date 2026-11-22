# MetaboStat: an MSNet plug-in (Student edition)

A Shiny application for statistical analysis of metabolomics data, combining univariate and multivariate methods with optional feature metadata support. Developed at Metatoul-AgromiX, Toulouse University. Based on mixOmics and ggpubr packages. Initial script from G. Marti, shiny using vibe coding.
Directly compatible with MSNet outputs.

**Use it online (nothing to install):** https://guillaumem31.github.io/metabostat-student/

- R runs inside your web browser (Shinylive / webR, including `mixOmics`): no R installation needed.
- The first load takes 30-90 s (R and the packages are downloaded into the browser), later visits are faster.
- **Your data never leaves your computer**: files are processed locally in the browser.
- Click *Load example data* in the sidebar to try it with the Cannabis chemotypes dataset (`Dataset/Cannabis/`).
- Heavy computations (tuning, repeated cross-validation) are slower in the browser than in desktop R; the default number of repeats is reduced accordingly. For large datasets, run the app in desktop R (see below).

## Features

- **Flexible data import**: CSV/TXT with auto-detected separator (comma, semicolon, tab)
- **Optional feature metadata**: pathway, class, m/z, RT, confidence, etc. — used to enrich heatmaps and volcano plots
- **Preprocessing**: missing value imputation (none, Min/10, NIPALS), normalization (TSS), transformation (Log10, CLR), and scaling (centering, UV, Pareto)
- **Group selection**: include/exclude experimental groups for pairwise or subset comparisons
- **Univariate analysis**: ANOVA/Kruskal-Wallis testing, p-value correction, results table, heatmap, and volcano plots
- **Multivariate analysis (categorical groups)**: PCA, PLS-DA, and sPLS-DA (auto or manual tuning), with classification error rate/BER quality metrics and biomarker heatmaps
- **Multivariate analysis (continuous response)**: PLS Regression and sparse PLS Regression against a continuous variable (e.g. a phenotype, a dose, a viral load) — scores colored by the response, Arrow and Correlation Circle plots, VIP scores, and cross-validated MSEP/R²/Q² quality metrics
- **Model quality metrics everywhere supervised**: PLS-DA, sPLS-DA, PLS-R and sPLS-R all report cross-validated performance (error rate/BER, or MSEP/R²/Q²) for the actual model you end up using, not just the tuning curve
- **Bidirectional feature filtering**: send significant univariate features into sPLS-DA, or send sPLS-DA/sPLS-R-selected features back into the univariate view, to cross-validate biomarkers
- **Exports**: plots as PDF, tables and processed data as CSV

## Requirements & Installation (Mac & Windows)

- **R**: Version ≥ 4.0 recommended (fully tested on the latest R 4.4 and 4.3).
- **RStudio**: Strongly recommended for an interactive student experience.

### Option A: One-Line Automatic Installation (Recommended)

In R or RStudio (with working directory set to this project folder):

```r
source("install_dependencies.R")
```

This script will automatically detect your operating system (Mac or Windows), install any missing CRAN packages, install `BiocManager` and `mixOmics` from Bioconductor, and verify that everything is working properly.

> **Note on Pareto Scaling**: Unlike older versions that depended on the archived `MetabolAnalyze` package, this student edition incorporates native base-R Pareto scaling directly, avoiding installation roadblocks on modern R versions.

### Option B: Manual Package Installation

If you prefer to install packages manually in R:

```r
# 1. Install CRAN packages
install.packages(c(
  "shiny", "shinyjs", "shinythemes", "dplyr", "purrr", "tibble",
  "tidyr", "ggplot2", "data.table", "vegan", "cluster", "ggrepel",
  "ggpubr", "multcompView", "DT", "RColorBrewer", "viridis", "scales"
))

# 2. Install Bioconductor package (mixOmics)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("mixOmics", update = FALSE, ask = FALSE)
```

## Running the App

### From R / RStudio:
Open the project in RStudio, then:
- Simply open `app.R` and click the green **"Run App"** button at the top right, or
- Run in the R console:
  ```r
  shiny::runApp()
  ```
  or
  ```r
  source("run_app.R")
  ```

The application itself is `app/app.R` (the same file is used by the desktop and by the browser version).

### From System File Explorer:
- **macOS**: Double-click `run_mac.command`.
- **Windows**: Double-click `run_windows.bat`.

## Student Tutorial & Manual

A comprehensive step-by-step tutorial is available in **[TUTORIAL.md](TUTORIAL.md)**. It guides students through:
- Data import and feature metadata matching (`Alignment ID`)
- Preprocessing (Imputation, Normalization, Transformation, Pareto scaling)
- Univariate statistics (ANOVA / Kruskal-Wallis, FDR correction, Volcano plots, Tukey post-hoc tests)
- Unsupervised Multivariate exploration (PCA, Scree plot, PERMANOVA, Silhouette score)
- Supervised Multivariate modeling (PLS-DA, Cross-validation, Classification error rates, AUROC)
- Biomarker discovery (sPLS-DA tuning, Loadings, Clustered Image Maps / Heatmaps)
- Biological interpretation using the included **Cannabis chemotypes dataset** (`Dataset/Cannabis/`).

## Input data format

**Statistical data CSV** — rows are samples, columns are a sample ID, a grouping column, and one column per feature:

```
Sample, Group, Feature1, Feature2, ...
```

**Feature metadata CSV (optional)** — one row per feature, with a feature ID column matching the feature column names in the statistical data file. Additional columns (pathway, class, m/z, RT, annotation confidence, etc.) can be used to annotate heatmaps and size volcano plot points.

## Typical workflow

1. Upload the statistical data file (and optional feature metadata), select the grouping column, and validate.
2. Preprocess: choose imputation, normalization, transformation, and scaling.
3. Select which groups to include in the analysis.
4. Run univariate analysis (ANOVA/Kruskal-Wallis) to identify significant features.
5. Run multivariate analysis to explore structure and select biomarkers:
   - Categorical groups: PCA → PLS-DA → sPLS-DA
   - Continuous response (e.g. a phenotype, dose, or viral load): choose the response variable in sidebar section 7, then PLS-R → sPLS-R
6. Check the Quality metrics / Performance panel for the model you plan to interpret (error rate/BER for PLS-DA/sPLS-DA, MSEP/R²/Q² for PLS-R/sPLS-R) before trusting its scores or selected features.
7. Use bidirectional filtering to cross-check univariate and multivariate results.
8. Export plots and tables for reporting.

See the in-app **Guide** tab for detailed step-by-step instructions.

## Citation

If you use this software in your research, please cite:

> G. Marti. (2025). Metabolomics Statistical Analysis Platform - AgromiX. Zenodo. https://doi.org/10.5281/zenodo.17319833

## Author

Guillaume Marti (Toulouse University)

Contributors: Claude AI (Shiny application and metadata enhancement support)

## Disclaimer

> **For educational purposes only.** This application was developed for teaching
> metabolomics statistics to Master students. It is provided "as is", without
> warranty of any kind, express or implied. The authors do not guarantee the accuracy,
> completeness or fitness for a particular purpose of the results and cannot be held
> liable for any consequence of their use. It is not a validated tool (GLP/GMP, clinical
> or regulatory use): always verify results independently.

## Changes in v1.5.0 (browser version, review of the statistics)

- **Browser (standalone) version** with Shinylive, example datasets and a *Load example data* button.
- **PCA is now always centered** (it was not centered with the default scaling "None", so PC1 mostly captured the mean).
- **NIPALS imputation works** (it called an argument that no longer exists and silently fell back to min/10); only missing cells are replaced.
- **TSS keeps the intensity scale** (rescaled to the median total signal): `log10(x+1)` is no longer almost linear and CLR pseudo-counts no longer dominate. CLR uses a data-driven pseudo-count.
- **Scaling is computed on the selected groups only**, and sPLS-DA / sPLS-R / tuning no longer apply an extra UV scaling on top of the method you chose.
- The **grouping column, PLS-R response and numeric sample-index column are never treated as features** (they used to be normalized and analysed like metabolites).
- **Leave-one-out cross-validation works** (PLS-DA / PLS-R), folds are limited to the smallest group size, and a **random seed** makes cross-validation, tuning and PERMANOVA reproducible.
- **PERMANOVA uses the Euclidean distance** (the default Bray-Curtis is meaningless on scores with negative values).
- **"sPLS-DA -> Univariate" no longer filters the multivariate models** (and vice versa); the univariate p-values stay adjusted over all tested features; cached models/tuning/performance are cleared when the data change.
- Sortable p-value columns, correct sample names in exported tables, same figure on screen and in PDF for boxplots/volcano, Tukey letters for group names containing "-", decimal-comma files detected correctly.
- **Feature labels**: when metadata are loaded, any metadata column (metabolite name, class, m/z...) can replace the feature IDs in the volcano plot, heatmaps, loadings, VIP plot, boxplots and tables (truncated to 10 characters by default; CSV exports keep the IDs).
- Caveats are now shown in the app: PLS-DA PERMANOVA and AUROC are computed on training data (optimistic); pre-selecting features before cross-validation biases the error rate.

## License

This project's own code is released under the CeCILL-B License (see [LICENSE](LICENSE) in English and [LICENCE-fr.txt](LICENCE-fr.txt) in French), Copyright (c) 2025-2026 Guillaume Marti, University of Toulouse.

It depends on third-party R packages distributed under their own licenses, including mixOmics and ggpubr (GPL (≥ 2)). These are used as installed dependencies, not redistributed or modified, so they remain separately licensed and do not affect the CeCILL-B license of this code.


## For maintainers

- `app/app.R` is the application (plus the example data in `app/`).
- `docs/` is the generated website served by GitHub Pages; do not edit it by hand.
- After editing `app/app.R`: `Rscript build.R`, check locally with
  `Rscript -e "httpuv::runStaticServer('docs')"`, then commit and push.
