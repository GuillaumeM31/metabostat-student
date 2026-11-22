#==============================================================================
# MetaboStat / MSnet-Stat - Student Edition
# Automated Dependency Installer and Environment Checker
#
# Works seamlessly on macOS (Apple Silicon / Intel) and Windows (10 / 11)
# Compatible with modern R (R >= 4.0, including R 4.3 and R 4.4+)
#==============================================================================

cat("\n===================================================================\n")
cat("  MetaboStat / MSnet-Stat - Automated Dependency Installer\n")
cat("===================================================================\n\n")

# 1. System and R Version Information
os_name <- Sys.info()[["sysname"]]
r_ver <- paste(R.version$major, R.version$minor, sep = ".")
cat(sprintf("Operating System : %s (%s)\n", os_name, Sys.info()[["machine"]]))
cat(sprintf("R Version        : %s\n\n", r_ver))

if (getRversion() < "4.0.0") {
  warning("Your R version is older than 4.0.0. An upgrade to the latest R version is strongly recommended.")
}

# 2. Package Lists
cran_packages <- c(
  "shiny",
  "shinyjs",
  "shinythemes",
  "dplyr",
  "purrr",
  "tibble",
  "tidyr",
  "ggplot2",
  "data.table",
  "vegan",
  "cluster",
  "ggrepel",
  "ggpubr",
  "multcompView",
  "DT",
  "RColorBrewer",
  "viridis",
  "scales"
)

bioc_packages <- c(
  "mixOmics"
)

# Note: 'MetabolAnalyze' was previously listed in older tutorials, but has been
# archived from CRAN. In this student version, Pareto scaling is natively
# implemented in base R, eliminating that problematic dependency entirely.

# 3. Check and Install CRAN Packages
installed <- rownames(installed.packages())
missing_cran <- setdiff(cran_packages, installed)

if (length(missing_cran) > 0) {
  cat(sprintf("--> Installing %d missing CRAN package(s): %s\n", 
              length(missing_cran), paste(missing_cran, collapse = ", ")))
  
  # Set a reliable CRAN mirror if none is set
  repos <- getOption("repos")
  if (is.null(repos) || repos["CRAN"] == "@CRAN@") {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }
  
  install.packages(missing_cran, dependencies = TRUE)
} else {
  cat("--> All required CRAN packages are already installed.\n")
}

# 4. Check and Install Bioconductor & mixOmics
missing_bioc <- setdiff(bioc_packages, rownames(installed.packages()))

if (length(missing_bioc) > 0) {
  cat("\n--> Checking Bioconductor manager...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("    Installing BiocManager from CRAN...\n")
    install.packages("BiocManager")
  }
  
  cat(sprintf("--> Installing %d Bioconductor package(s) via BiocManager: %s\n",
              length(missing_bioc), paste(missing_bioc, collapse = ", ")))
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
} else {
  cat("--> All required Bioconductor packages (mixOmics) are already installed.\n")
}

# 5. Diagnostic Verification Test
cat("\n===================================================================\n")
cat("  Verifying Package Installations...\n")
cat("===================================================================\n")

all_pkgs <- c(cran_packages, bioc_packages)
all_ok <- TRUE

for (pkg in all_pkgs) {
  pkg_ok <- suppressWarnings(suppressPackageStartupMessages(
    requireNamespace(pkg, quietly = TRUE)
  ))
  if (pkg_ok) {
    ver <- as.character(packageVersion(pkg))
    cat(sprintf("  [OK] %-16s (v%s)\n", pkg, ver))
  } else {
    cat(sprintf("  [MISSING] %-16s FAILED TO LOAD\n", pkg))
    all_ok <- FALSE
  }
}

cat("-------------------------------------------------------------------\n")
if (all_ok) {
  cat("SUCCESS: All dependencies are installed and verified!\n")
  cat("You can now launch the application by running:\n\n")
  cat("   shiny::runApp()\n")
  cat("   or\n")
  cat("   source(\"run_app.R\")\n\n")
} else {
  cat("WARNING: Some packages could not be installed automatically.\n")
  cat("Please check the error messages above or refer to README.md.\n\n")
}
cat("===================================================================\n\n")
