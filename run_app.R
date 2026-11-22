#==============================================================================
# MetaboStat / MSnet-Stat - Student Edition
# One-Click Application Launcher (run_app.R)
#
# Usage in R / RStudio:
#   source("run_app.R")
#
# Usage from Terminal:
#   Rscript run_app.R
#==============================================================================

# Determine directory of this script
script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg)))
  } else if (!is.null(sys.frame(1)$ofile)) {
    dirname(normalizePath(sys.frame(1)$ofile))
  } else {
    getwd()
  }
}, error = function(e) getwd())

if (dir.exists(script_dir)) {
  setwd(script_dir)
}

# Check essential packages
essential_pkgs <- c("shiny", "mixOmics", "DT", "ggplot2", "dplyr", "data.table")
installed <- rownames(installed.packages())
missing <- setdiff(essential_pkgs, installed)

if (length(missing) > 0) {
  message("Some required packages are missing: ", paste(missing, collapse = ", "))
  message("Running automatic dependency installer...")
  source("install_dependencies.R")
}

message("\n--------------------------------------------------------------")
message("Launching MetaboStat Shiny Application...")
message("Opening your default web browser...")
message("To stop the application, press Ctrl+C or Esc in R/RStudio.")
message("--------------------------------------------------------------\n")

shiny::runApp(appDir = ".", launch.browser = TRUE)
