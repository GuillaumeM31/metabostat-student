#==============================================================================
# MetaboStat / MSnet-Stat - Student Edition
# Standard Shiny Application Entry Point (app.R)
#
# Allows running with:
#   shiny::runApp()
# or clicking "Run App" in RStudio
#==============================================================================

app_file <- "MSnet-Stat_1.260911_app.r"

if (!file.exists(app_file)) {
  # If run from another working directory, locate relative to script
  this_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
  app_file <- file.path(this_dir, "MSnet-Stat_1.260911_app.r")
}

if (!file.exists(app_file)) {
  stop("Could not locate ", app_file, ". Please set working directory to the project folder.")
}

source(app_file, local = FALSE)
shiny::shinyApp(ui = ui, server = server)
