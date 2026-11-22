#==============================================================================
# MetaboStat / MSnet-Stat - Student Edition
# Metabolomics Statistical Analysis Platform (Metatoul-AgromiX)
# Unified Shiny Application for Statistical Analysis of Metabolomics Data
#
# Version: 1.5.0 - runs both in desktop R and in the browser (Shinylive / webR)
# Date: November 2026
#
# Author: Guillaume Marti (Toulouse University)
#
# Contributors:
#   - Claude AI for shiny application and metadata enhancement
#
# License: CeCILL-B (BSD-like, compatible with French law) - see LICENSE
# Copyright (c) 2025-2026 Guillaume Marti, University of Toulouse
#
# DISCLAIMER: for educational purposes only. Provided "as is", without any
# warranty; the authors do not guarantee the accuracy of the results and accept
# no liability for their use. Not a validated tool (GLP/GMP, clinical or
# regulatory use).
#
# Citation:
# If you use this software in your research, please cite:
# G.Marti. (2025). Metabolomics Statistical Analysis Platform - AgromiX.
# Zenodo repository: 10.5281/zenodo.17319833
#
#==============================================================================

options(shiny.maxRequestSize = 200 * 1024^2)  # 200 Mo

# TRUE when the app runs inside the browser (Shinylive / webR): slower
# computations and no multi-core, so heavy defaults are reduced.
IS_WEBR <- identical(R.version$os, "emscripten")

# --- Load libraries ---
suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)
  library(shinythemes)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(data.table)
  library(mixOmics)
  library(vegan)
  library(cluster)
  library(ggrepel)
  library(ggpubr)
  library(multcompView)
  library(DT)
  library(RColorBrewer)
  library(viridis)
  library(scales)
})

#==============================================================================
# USER INTERFACE
#==============================================================================

# Downloads: fetch the file inside the page and save it as a blob. In Shinylive
# (webR in the browser), a plain link click opens a new tab that can miss the
# in-browser R server and saves an HTML page instead of the PDF/CSV.
download_js <- "
document.addEventListener('click', async function (e) {
  var a = e.target.closest('a.shiny-download-link');
  if (!a) return;
  e.preventDefault();
  e.stopPropagation();
  var href = a.getAttribute('href');
  if (!href || a.classList.contains('disabled')) {
    alert('The file is not ready yet: load and validate your data first.');
    return;
  }
  try {
    var r = await fetch(href);
    var type = r.headers.get('content-type') || '';
    if (!r.ok || type.indexOf('text/html') === 0) throw new Error('the app could not create the file (check that the analysis has been run)');
    var cd = r.headers.get('content-disposition') || '';
    var m = cd.match(/filename\\*?=(?:UTF-8'')?\"?([^\";]+)\"?/i);
    var name = m ? decodeURIComponent(m[1]) : 'download';
    var url = URL.createObjectURL(await r.blob());
    var tmp = document.createElement('a');
    tmp.href = url; tmp.download = name; tmp.style.display = 'none';
    document.body.appendChild(tmp); tmp.click(); tmp.remove();
    setTimeout(function () { URL.revokeObjectURL(url); }, 30000);
  } catch (err) {
    alert('Download failed: ' + err.message);
  }
}, true);
"

# Example datasets shipped next to the app (relative to the app folder).
EXAMPLES <- list(
  "Cannabis chemotypes (LC-MS, 18 samples, with feature metadata)" =
    list(data = "StatTableNEGPOS.csv", meta = "NodesNEGPOS.csv", group = "Class")
)

ui <- fluidPage(
  theme = shinytheme("flatly"),
  useShinyjs(),
  tags$head(tags$script(HTML(download_js))),
  # Load the DataTables JavaScript with the page: otherwise, in the browser
  # version, tables rendered at the same time as the library is first fetched
  # stay empty ("$table.DataTable is not a function").
  htmltools::attachDependencies(
    tags$div(style = "display: none;"),
    htmltools::findDependencies(DT::datatable(data.frame(x = 1)))
  ),
  
  titlePanel(
    div(
      h2("Metabolomics Statistical Analysis Platform"),
      h5("Integrated Univariate and Multivariate Analysis with Metadata Support")
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      h4("1. Data Input"),
      uiOutput("example_ui"),
      fileInput("file1", "Statistical Data CSV", accept = c(".csv", ".txt")),
      radioButtons("sep", "Separator:", 
                   choices = c("Comma" = ",", "Semicolon" = ";", "Tab" = "\t"),
                   selected = ",", inline = TRUE),
      selectInput("group_col", "Grouping Column", choices = NULL),
      tags$hr(),
      
      h4("2. Feature Metadata (Optional)"),
      fileInput("metadata_file", "Feature Metadata CSV", accept = c(".csv", ".txt")),
      radioButtons("meta_sep", "Metadata Separator:", 
                   choices = c("Comma" = ",", "Semicolon" = ";", "Tab" = "\t"),
                   selected = ",", inline = TRUE),
      selectInput("feature_id_col", "Feature ID Column", choices = NULL),
      checkboxInput("check_metadata_match", "Verify Matching", value = FALSE),
      conditionalPanel(
        condition = "input.check_metadata_match",
        verbatimTextOutput("metadata_match_info")
      ),
      actionButton("validate_data", "Validate All Data", class = "btn-primary"),
      tags$hr(),
      
      h4("3. Preprocessing"),
      selectInput("impute_method", "Missing value imputation:",
                  choices = c("None (remaining NAs set to 0)" = "none",
                              "Min/10 (detection limit)" = "min10",
                              "NIPALS (mixOmics)" = "nipals")),
      sliderInput("max_na_prop", "Remove variables with more than this % missing:",
                  min = 10, max = 80, value = 50, step = 5,
                  post = "%"),
      selectInput("norm_method", "Normalization:",
                  choices = c("None" = "none",
                              "TSS (rescaled to median total signal)" = "tss")),
      selectInput("trans_method", "Transformation:",
                  choices = c("None" = "none",
                              "Log10(x+1)" = "log10",
                              "CLR (Centered Log Ratio)" = "clr")),
      selectInput("scale_method", "Scaling:",
                  choices = c("None (Raw data)" = "none",
                              "Centering only" = "center", 
                              "UV" = "uv", 
                              "Pareto" = "pareto")),
      tags$hr(),
      
      h4("4. Group Selection"),
      uiOutput("group_selector_ui"),
      fluidRow(
        column(6, actionButton("reset_groups", "Reset", class = "btn-warning", style = "width: 100%;")),
        column(6, actionButton("apply_groups", "Apply", class = "btn-success", style = "width: 100%;"))
      ),
      tags$hr(),
      
      conditionalPanel(
        condition = "output.metadata_loaded",
        h4("5. Metadata Annotation"),
        selectInput("heatmap_annotation_col", "Heatmap Annotation (used in Univariate & Multivariate heatmaps):", choices = NULL),
        selectInput("feature_label_col", "Feature labels in plots and tables:", choices = NULL),
        helpText("Replace the feature IDs by a metadata column (e.g. the metabolite name) in the volcano plot, heatmaps, loadings, VIP and boxplots. Features without a name keep their ID; CSV exports keep the IDs."),
        numericInput("feature_label_max", "Max label length (characters):", value = 10, min = 5, max = 200),
        tags$hr()
      ),

      conditionalPanel(
        condition = "input.main_tabs == 'Multivariate'",
        h4("6. Multivariate Parameters"),
        sliderInput("ncomp_mv", "Number of components:", value = 5, min = 2, max = 15),
        
        radioButtons("splsda_mode", "sPLS-DA:",
                     choices = c("Auto" = "auto", "Manual" = "manual"),
                     selected = "auto"),
        
        conditionalPanel(
          condition = "input.splsda_mode == 'auto'",
          radioButtons("splsda_strategy", "Strategy:",
                       choices = c("Sparse" = "sparse", "Robust" = "robust"),
                       selected = "sparse"),
          actionButton("run_splsda_tuning", "Run Tuning", class = "btn-warning")
        ),
        
        conditionalPanel(
          condition = "input.splsda_mode == 'manual'",
          numericInput("splsda_ncomp_manual", "Nb comp:", value = 2, min = 1, max = 10),
          textInput("splsda_keepx_manual", "keepX:", value = "25,20")
        ),

        tags$hr(),
        h4("7. PLS Regression Parameters"),
        selectInput("plsr_y_col", "Response variable (Y):", choices = c("None" = "")),
        helpText("Numeric column to regress against (e.g. a continuous phenotype). ",
                 "Automatically excluded from the predictor matrix used for PLS-R/sPLS-R."),

        radioButtons("splsr_mode", "sPLS-R:",
                     choices = c("Auto" = "auto", "Manual" = "manual"),
                     selected = "auto"),

        conditionalPanel(
          condition = "input.splsr_mode == 'auto'",
          radioButtons("splsr_strategy", "Strategy:",
                       choices = c("Sparse" = "sparse", "Robust" = "robust"),
                       selected = "sparse"),
          selectInput("splsr_tune_measure", "Tuning criterion:",
                      choices = c("MSE" = "MSE", "R2" = "R2")),
          actionButton("run_splsr_tuning", "Run Tuning", class = "btn-warning")
        ),

        conditionalPanel(
          condition = "input.splsr_mode == 'manual'",
          numericInput("splsr_ncomp_manual", "Nb comp:", value = 2, min = 1, max = 10),
          textInput("splsr_keepx_manual", "keepX:", value = "10,10")
        )
      ),

      tags$hr(),
      numericInput("seed", "Random seed (cross-validation, tuning, PERMANOVA):",
                   value = 123, min = 1, step = 1),
      numericInput("pdf_width", "PDF width:", value = 10, min = 5, max = 20),
      numericInput("pdf_height", "PDF height:", value = 8, min = 5, max = 20)
    ),
    
    mainPanel(
      width = 9,
      
      tabsetPanel(
        id = "main_tabs",
        
        tabPanel(
          "Guide",
          div(class = "well", style = "background-color: #e8f4f8; border-left: 5px solid #2980b9;",
              h3("🔬 Statistical Analysis Platform for Metabolomics"),
              p(style = "font-size: 16px;", "Welcome! This platform provides comprehensive statistical analysis tools for metabolomics data with support for feature metadata and bidirectional analysis workflows.")
          ),
          
          tags$hr(),
          
          h3("📋 Step-by-Step Workflow"),
          
          div(class = "panel panel-default",
              div(class = "panel-heading", style = "background-color: #3498db; color: white;",
                  h4(icon("upload"), " Step 1: Data Upload")),
              div(class = "panel-body",
                  tags$ol(
                    tags$li(strong("Try it first:"), " pick an example dataset and click 'Load example data' (Cannabis chemotypes)"),
                    tags$li(strong("Main Data File:"), " Upload your metabolomics data matrix as CSV/TXT"),
                    tags$ul(
                      tags$li("Rows = samples, Columns = features + one grouping column"),
                      tags$li("Numeric feature values with optional missing data"),
                      tags$li("Example: Sample1, Sample2... | Group | Feature1, Feature2..."),
                      tags$li("The grouping column, the PLS-R response and a numeric sample-index column are never used as features")
                    ),
                    tags$li(strong("Select Separator:"), " The app auto-detects, but verify comma/semicolon/tab"),
                    tags$li(strong("Choose Grouping Column:"), " Select which column contains your experimental groups"),
                    tags$li(strong("Optional - Feature Metadata:"), " Upload additional feature information"),
                    tags$ul(
                      tags$li("One row per feature with Feature ID matching data column names"),
                      tags$li("Include: Pathway, Class, m/z, RT, Confidence, etc."),
                      tags$li("Used for enhanced visualizations (color-coded heatmaps, sized volcano points)")
                    ),
                    tags$li(strong("Validate:"), " Click 'Validate All Data' to proceed")
                  )
              )
          ),
          
          div(class = "panel panel-default",
              div(class = "panel-heading", style = "background-color: #9b59b6; color: white;",
                  h4(icon("cogs"), " Step 2: Preprocessing")),
              div(class = "panel-body",
                  tags$ol(
                    tags$li(strong("Missing Value Imputation:")),
                    tags$ul(
                      tags$li("None: Remove features with >X% missing"),
                      tags$li("Min/10: Replace with detection limit estimate"),
                      tags$li("NIPALS: Multivariate imputation (recommended for < 50% missing)")
                    ),
                    tags$li(strong("Normalization:"), " TSS (Total Sum Scaling, rescaled to the median total signal so the data stay in intensity units) for sample-to-sample variation"),
                    tags$li(strong("Transformation:"), " Log10 or CLR for variance stabilization"),
                    tags$li(strong("Scaling:"), " UV (Unit Variance) or Pareto for feature weighting. It is computed on the selected groups only, so excluding a group re-scales the data"),
                    tags$li(tags$em("Tip:"), " Typical untargeted metabolomics workflow: Min/10 → Log10 → Pareto (or UV)")
                  )
              )
          ),
          
          div(class = "panel panel-default",
              div(class = "panel-heading", style = "background-color: #e67e22; color: white;",
                  h4(icon("filter"), " Step 3: Group Selection")),
              div(class = "panel-body",
                  tags$ul(
                    tags$li("Select 2+ groups to include in your analysis"),
                    tags$li("Useful for pairwise comparisons or excluding control/outlier groups"),
                    tags$li("Analyses follow your selection; cached models / tuning / performance results are cleared whenever the data change"),
                    tags$li("Use 'Reset' to return to all groups")
                  )
              )
          ),
          
          div(class = "panel panel-default",
              div(class = "panel-heading", style = "background-color: #27ae60; color: white;",
                  h4(icon("chart-bar"), " Step 4: Univariate Analysis")),
              div(class = "panel-body",
                  tags$ol(
                    tags$li(strong("Statistical Tests:"), " ANOVA/Kruskal-Wallis for multi-group comparison"),
                    tags$li(strong("P-value Threshold:"), " Filter significant features (default p < 0.05)"),
                    tags$li(strong("Results Table:"), " Click rows to visualize individual features"),
                    tags$li(strong("Heatmap:"), " Shows top N significant features with group patterns"),
                    tags$li(strong("Volcano Plot:"), " Visualize fold-change vs. significance (2 groups)"),
                    tags$li(strong(icon("exchange-alt"), " Feature Export:"),
                            " Send significant features to sPLS-DA for focused biomarker discovery",
                            div(class = "alert alert-info", style = "margin-top: 10px;",
                                "Click 'Use Significant Features in sPLS-DA' to pre-filter multivariate analysis"))
                  )
              )
          ),
          
          div(class = "panel panel-default",
              div(class = "panel-heading", style = "background-color: #2980b9; color: white;",
                  h4(icon("project-diagram"), " Step 5: Multivariate Analysis")),
              div(class = "panel-body",
                  tags$ol(
                    tags$li(strong("PCA (Unsupervised):"), " Explore overall data structure without group labels"),
                    tags$ul(
                      tags$li("Check for outliers, batch effects, and natural clustering"),
                      tags$li("PERMANOVA tests group separation significance")
                    ),
                    tags$li(strong("PLS-DA (Supervised, categorical groups):"), " Maximize separation between known groups"),
                    tags$ul(
                      tags$li("Run Performance validation to assess model quality"),
                      tags$li("Check the 'Quality metrics' panel for classification error rate / BER per component"),
                      tags$li("Check AUC curves for discrimination ability")
                    ),
                    tags$li(strong("sPLS-DA (Sparse PLS-DA):"), " Select discriminant biomarkers"),
                    tags$ul(
                      tags$li(strong("Auto mode:"), " Algorithm optimizes number of features per component"),
                      tags$li(strong("Manual mode:"), " Specify ncomp and keepX values"),
                      tags$li("View loadings plots to identify key features"),
                      tags$li("Click 'Evaluate Final Model' for the error rate / BER of the model you actually end up using (the tuning curve only reports the criterion used to pick keepX)"),
                      tags$li("Heatmap shows selected biomarkers, columns colored by group and rows by metadata (if provided)")
                    ),
                    tags$li(strong("PLS-R (Supervised, continuous response):"), " Regress features against a continuous variable (e.g. a phenotype, a dose, a viral load)"),
                    tags$ul(
                      tags$li("Choose the response variable (Y) in sidebar section 7 - it is automatically excluded from the predictor matrix"),
                      tags$li("Scores (colored by Y), Arrow plot, and Correlation Circle show how features relate to the response"),
                      tags$li("Run Performance for cross-validated MSEP / R2 / Q2 - a negative Q2 means the model does not generalize to new samples"),
                      tags$li("VIP scores above 1 highlight the most influential predictors")
                    ),
                    tags$li(strong("sPLS-R (Sparse PLS Regression):"), " Select the predictors most associated with the response"),
                    tags$ul(
                      tags$li(strong("Auto mode:"), " Tunes ncomp and keepX by cross-validated MSE or R2"),
                      tags$li(strong("Manual mode:"), " Specify ncomp and keepX values"),
                      tags$li("Click 'Evaluate Final Model' for the Q2 / R2 / MSEP of the actual sparse model"),
                      tags$li("Heatmap shows selected features, columns colored by the response variable (blue = low, red = high)")
                    ),
                    tags$li(strong(icon("exchange-alt"), " Feature Import:"),
                            " Send sPLS-DA or sPLS-R selected features back to Univariate tab",
                            div(class = "alert alert-info", style = "margin-top: 10px;",
                                "Click 'Show sPLS-DA/sPLS-R Features in Univariate' to examine individual feature statistics"))
                  )
              )
          ),
          
          div(class = "panel panel-success",
              div(class = "panel-heading",
                  h4(icon("lightbulb"), " Best Practices & Tips")),
              div(class = "panel-body",
                  tags$ul(
                    tags$li(strong("Data Quality:"), " Always check the 'Data' tab to verify import success"),
                    tags$li(strong("Reproducibility:"), " cross-validation, tuning and PERMANOVA use the random seed set in the sidebar"),
                    tags$li(strong("Sample Size:"), " Minimum 4-5 samples per group recommended"),
                    tags$li(strong("Feature Selection:"), " Start with univariate → filter → multivariate workflow"),
                    tags$li(strong("Validation:"), " Check the Quality metrics / Performance panels (PLS-DA, sPLS-DA, PLS-R, sPLS-R) before interpreting any supervised model - a poor error rate/BER or a negative Q2 means the model is not trustworthy, however good the scores plot looks"),
                    tags$li(strong("Choosing DA vs. Regression:"), " Use PLS-DA/sPLS-DA when comparing categorical groups, and PLS-R/sPLS-R when relating features to a continuous variable (dose, time, phenotype score, viral load, etc.)"),
                    tags$li(strong("Metadata:"), " Add pathway/class info for more interpretable visualizations"),
                    tags$li(strong("Bidirectional Analysis:"), " Combine univariate and multivariate for comprehensive insights:"),
                    tags$ul(
                      tags$li("Univariate → sPLS-DA: Pre-filter noise before biomarker selection"),
                      tags$li("sPLS-DA / sPLS-R → Univariate: Validate multivariate biomarkers individually")
                    ),
                    tags$li(strong("Downloads:"), " Export all plots as PDF and tables as CSV for publications")
                  )
              )
          ),
          
          tags$hr(),
          
          h3("🔄 Bidirectional Feature Filtering"),
          
          div(class = "panel panel-warning",
              div(class = "panel-heading", h4("Cross-Analysis Workflow")),
              div(class = "panel-body",
                  fluidRow(
                    column(6,
                           div(style = "border: 2px solid #f39c12; padding: 15px; border-radius: 5px; background-color: #fef5e7;",
                               h5(icon("arrow-right"), strong(" Univariate → sPLS-DA")),
                               p("Pre-filter features before multivariate analysis"),
                               tags$ol(
                                 tags$li("Run ANOVA to identify significant features"),
                                 tags$li("Click 'Use Significant Features in sPLS-DA'"),
                                 tags$li("sPLS-DA will focus on pre-filtered biomarkers"),
                                 tags$li("Reduces noise and improves model performance")
                               ),
                               p(strong("Use when:"), " You have many features and want to reduce complexity")
                           )
                    ),
                    column(6,
                           div(style = "border: 2px solid #3498db; padding: 15px; border-radius: 5px; background-color: #ebf5fb;",
                               h5(icon("arrow-left"), strong(" sPLS-DA / sPLS-R → Univariate")),
                               p("Examine individual biomarker statistics"),
                               tags$ol(
                                 tags$li("Run sPLS-DA (categorical groups) or sPLS-R (continuous response) to select features"),
                                 tags$li("Click 'Show sPLS-DA Features in Univariate' or 'Show sPLS-R Features in Univariate'"),
                                 tags$li("View detailed statistics for selected biomarkers"),
                                 tags$li("Create focused visualizations and comparisons")
                               ),
                               p(strong("Use when:"), " You want to validate multivariate biomarkers individually")
                           )
                    )
                  ),
                  tags$hr(),
                  div(class = "alert alert-success",
                      icon("check-circle"),
                      strong(" Both directions can be used iteratively to refine your biomarker discovery!"))
              )
          ),
          
          tags$hr(),
          
          div(class = "well",
              h4(icon("question-circle"), " Need Help?"),
              p("For questions or issues, please refer to the documentation or contact the development team."),
              p(strong("Citation:"), " G. Marti. (2025). Metabolomics Statistical Analysis Platform - AgromiX. Zenodo: 10.5281/zenodo.17319833")
          )
        ),
        
        tabPanel(
          "Data",
          h4("Data Summary"),
          verbatimTextOutput("data_summary"),
          tags$hr(),
          h5("Feature Metadata Summary"),
          verbatimTextOutput("metadata_summary"),
          tags$hr(),
          h5("Raw data preview"),
          DT::dataTableOutput("data_preview_raw"),
          tags$hr(),
          h5("Metadata preview"),
          DT::dataTableOutput("metadata_preview"),
          tags$hr(),
          h5("Preprocessed data preview"),
          helpText(textOutput("processed_note")),
          DT::dataTableOutput("data_preview_processed"),
          downloadButton("download_processed_data", "Download")
        ),
        
        tabPanel(
          "Univariate",
          # Filter status indicator
          conditionalPanel(
            condition = "output.uv_filter_active",
            div(class = "alert alert-warning", style = "margin-bottom: 15px;",
                icon("filter"), 
                strong(" Filtered Mode: "),
                "Showing only features selected by ",
                textOutput("uv_filter_source", inline = TRUE),
                " (",
                textOutput("uv_filter_count", inline = TRUE),
                " features). P-values are still adjusted over all tested features. ",
                actionButton("clear_splsda_filter_top", "Show All Features", 
                             class = "btn-warning btn-xs")
            )
          ),
          div(class = "alert alert-info", style = "margin-bottom: 15px;",
              icon("info-circle"),
              textOutput("selected_groups_info_univ", inline = TRUE)
          ),
          div(class = "well well-sm",
              h5(icon("calculator"), " Statistics"),
              fluidRow(
                column(4, checkboxInput("use_nonparametric", "Non-parametric tests", value = FALSE)),
                column(4, selectInput("adj_method", "P-value correction:",
                                       choices = c("FDR" = "BH", "Bonferroni" = "bonferroni"))),
                column(4, numericInput("pval_threshold", "P-value threshold:", value = 0.05, min = 0.001, max = 0.5, step = 0.01))
              )
          ),
          tabsetPanel(
            tabPanel(
              "Tests",
              h4("Univariate test results"),
              fluidRow(
                column(4, numericInput("n_heatmap_features", "Nb features heatmap:", value = 40, min = 5, max = 200)),
                column(4, downloadButton("download_univ_table", "CSV")),
                column(4, downloadButton("download_univ_filtered_data", "Filtered data"))
              ),
              tags$hr(),
              DT::dataTableOutput("univariate_table"),
              tags$hr(),
              fluidRow(
                column(6, selectInput("univ_heatmap_dist", "Distance measure:",
                                       choices = c("Euclidean" = "euclidean",
                                                   "Correlation (Pearson)" = "correlation",
                                                   "Manhattan" = "manhattan",
                                                   "Maximum" = "maximum",
                                                   "Canberra" = "canberra",
                                                   "Binary" = "binary",
                                                   "Minkowski" = "minkowski"),
                                       selected = "euclidean")),
                column(6, selectInput("univ_heatmap_clust", "Clustering method:",
                                       choices = c("Complete" = "complete",
                                                   "Ward (D)" = "ward.D",
                                                   "Average (UPGMA)" = "average",
                                                   "Single" = "single",
                                                   "McQuitty (WPGMA)" = "mcquitty",
                                                   "Median (WPGMC)" = "median",
                                                   "Centroid (UPGMC)" = "centroid"),
                                       selected = "complete"))
              ),
              plotOutput("univariate_heatmap", height = "700px"),
              fluidRow(
                column(6, downloadButton("download_univ_heatmap", "Download Heatmap (PDF)")),
                column(6, downloadButton("download_univ_heatmap_data", "Download Heatmap Data (CSV)"))
              ),
              tags$hr(),
              verbatimTextOutput("filtering_stats"),
              tags$hr(),
              # Feature filtering controls
              fluidRow(
                column(12,
                       div(class = "well well-sm", style = "background-color: #f0f8ff;",
                           h5(icon("filter"), " Feature Filtering for Multivariate Analysis"),
                           p("Send these significant features to sPLS-DA for focused biomarker discovery."),
                           helpText("Caution: selecting features with the same samples that are then used for cross-validation makes the error rate over-optimistic (selection bias). Use it for exploration, and confirm biomarkers on independent samples."),
                           actionButton("send_to_splsda", 
                                        "Use Significant Features in sPLS-DA", 
                                        class = "btn-info",
                                        icon = icon("arrow-right")),
                           conditionalPanel(
                             condition = "output.univ_filter_active",
                             HTML("&nbsp;&nbsp;"),
                             actionButton("clear_univ_filter", 
                                          "Reset sPLS-DA to All Features", 
                                          class = "btn-warning btn-sm",
                                          icon = icon("times"))
                           )
                       )
                )
              )
            ),
            
            tabPanel(
              "Plots",
              fluidRow(
                column(6, radioButtons("plot_type", "Plot type:",
                                        choices = c("Boxplot" = "box", "Violin plot" = "violin"),
                                        selected = "box", inline = TRUE)),
                column(6, checkboxInput("add_points", "Add individual points", value = FALSE))
              ),
              plotOutput("univariate_boxplot", height = "500px"),
              downloadButton("download_univ_boxplot", "PDF")
            ),
            
            tabPanel(
              "Volcano",
              fluidRow(
                column(6, selectInput("volcano_group1", "Group 1:", choices = NULL)),
                column(6, selectInput("volcano_group2", "Group 2:", choices = NULL))
              ),
              fluidRow(
                column(6, sliderInput("volcano_pval", "Threshold -Log10(p):", min = 0, max = 10, value = 1.3, step = 0.1)),
                column(6, sliderInput("volcano_fc", "Threshold |Log2FC|:", min = 0, max = 5, value = 1, step = 0.1))
              ),
              checkboxInput("volcano_label", "Label features", value = TRUE),
              numericInput("volcano_n_labels", "Nb labels:", value = 15, min = 0, max = 50),
              conditionalPanel(
                condition = "output.metadata_loaded",
                fluidRow(
                  column(6, selectInput("volcano_size_col", "Volcano Point Size:", choices = NULL)),
                  column(6, sliderInput("volcano_size_range", "Size Range:",
                                         min = 0.5, max = 10, value = c(1, 5), step = 0.5))
                )
              ),
              tags$hr(),
              fluidRow(
                column(7, 
                       plotOutput("volcano_plot", height = "600px"),
                       downloadButton("download_volcano", "PDF")
                ),
                column(5,
                       plotOutput("volcano_boxplot", height = "400px"),
                       downloadButton("download_volcano_boxplot", "PDF")
                )
              ),
              tags$hr(),
              DT::dataTableOutput("volcano_table"),
              downloadButton("download_volcano_table", "CSV")
            )
          )
        ),
        
        tabPanel(
          "Multivariate",
          # Filter status indicator
          conditionalPanel(
            condition = "output.univ_filter_active",
            div(class = "alert alert-warning", style = "margin-bottom: 15px;",
                icon("filter"), 
                strong(" Filtered Mode: "),
                "Using only significant features from ANOVA (",
                textOutput("univ_filter_count", inline = TRUE),
                " features). ",
                actionButton("clear_univ_filter_top", "Use All Features", 
                             class = "btn-warning btn-xs")
            )
          ),
          div(class = "alert alert-info", style = "margin-bottom: 15px;",
              icon("info-circle"), 
              textOutput("selected_groups_info_multi", inline = TRUE)
          ),
          tabsetPanel(
            tabPanel(
              "PCA",
              h4("Principal Component Analysis (unsupervised)"),
              fluidRow(
                column(6,
                       h5("Scores"),
                       fluidRow(
                         column(6, numericInput("pca_x", "PC X:", value = 1, min = 1, max = 15)),
                         column(6, numericInput("pca_y", "PC Y:", value = 2, min = 1, max = 15))
                       ),
                       plotOutput("pca_scores", height = "500px"),
                       downloadButton("download_pca_scores", "PDF")
                ),
                column(6,
                       plotOutput("pca_scree", height = "250px"),
                       downloadButton("download_pca_scree", "PDF"),
                       tags$hr(),
                       plotOutput("pca_biplot", height = "400px"),
                       downloadButton("download_pca_biplot", "PDF")
                )
              ),
              tags$hr(),
              fluidRow(
                column(6,
                       h5("PERMANOVA"),
                       verbatimTextOutput("pca_permanova"),
                       tags$hr(),
                       h5("F-test per Component"),
                       plotOutput("pca_f_test", height = "350px"),
                       downloadButton("download_pca_ftest", "PDF")
                ),
                column(6,
                       plotOutput("pca_silhouette", height = "400px"),
                       downloadButton("download_pca_silhouette", "PDF")
                )
              )
            ),
            
            tabPanel(
              "PLS-DA",
              fluidRow(
                column(6,
                       fluidRow(
                         column(6, numericInput("plsda_x", "Comp X:", value = 1, min = 1, max = 15)),
                         column(6, numericInput("plsda_y", "Comp Y:", value = 2, min = 1, max = 15))
                       ),
                       plotOutput("plsda_scores", height = "500px"),
                       downloadButton("download_plsda_scores", "PDF")
                ),
                column(6,
                       selectInput("plsda_cv_method", "CV:",
                                   choices = c("M-fold" = "Mfold", "LOO" = "loo")),
                       conditionalPanel(
                         condition = "input.plsda_cv_method == 'Mfold'",
                         numericInput("plsda_folds", "Folds:", value = 5, min = 2, max = 10),
                         numericInput("plsda_nrepeat", "Repeats:", value = if (IS_WEBR) 10 else 50, min = 3, max = 200)
                       ),
                       actionButton("run_plsda_perf", "Performance", class = "btn-info"),
                       tags$hr(),
                       plotOutput("plsda_perf_plot", height = "400px"),
                       downloadButton("download_plsda_perf", "PDF"),
                       tags$hr(),
                       h5("Quality metrics (cross-validated)"),
                       verbatimTextOutput("plsda_perf_summary")
                )
              ),
              tags$hr(),
              fluidRow(
                column(6,
                       helpText("PERMANOVA on the PLS-DA scores. Caution: PLS-DA scores are built to separate the groups, so this p-value is circular (over-optimistic); rely on the cross-validated error rate above."),
                       verbatimTextOutput("plsda_permanova")),
                column(6,
                       helpText("AUROC is computed on the training data (optimistic); the cross-validated error rate/BER above is the honest estimate."),
                       numericInput("plsda_auc_comp", "Comp AUC:", value = 1, min = 1, max = 15),
                       plotOutput("plsda_auc", height = "400px"),
                       downloadButton("download_plsda_auc", "PDF")
                )
              )
            ),
            
            tabPanel(
              "sPLS-DA",
              conditionalPanel(
                condition = "input.splsda_mode == 'auto'",
                plotOutput("splsda_tune_plot", height = "400px"),
                downloadButton("download_splsda_tune", "PDF"),
                verbatimTextOutput("splsda_tune_summary"),
                tags$hr()
              ),
              fluidRow(
                column(6,
                       plotOutput("splsda_scores", height = "500px"),
                       downloadButton("download_splsda_scores", "PDF")
                ),
                column(6,
                       verbatimTextOutput("splsda_selected_features"),
                       downloadButton("download_splsda_features", "CSV")
                )
              ),
              tags$hr(),
              fluidRow(
                column(12,
                       h5("Final Model Performance (cross-validated)"),
                       actionButton("run_splsda_final_perf", "Evaluate Final Model", class = "btn-info"),
                       tags$hr(),
                       fluidRow(
                         column(6,
                                plotOutput("splsda_final_perf_plot", height = "350px"),
                                downloadButton("download_splsda_final_perf", "PDF")
                         ),
                         column(6, verbatimTextOutput("splsda_final_perf_summary"))
                       )
                )
              ),
              tags$hr(),
              # Feature export to univariate
              fluidRow(
                column(12,
                       div(class = "well well-sm", style = "background-color: #fff8f0;",
                           h5(icon("exchange-alt"), " Cross-Analysis with Univariate"),
                           p("Examine univariate statistics for sPLS-DA-selected features."),
                           actionButton("send_to_univariate", 
                                        "Show sPLS-DA Features in Univariate", 
                                        class = "btn-info",
                                        icon = icon("arrow-left")),
                           conditionalPanel(
                             condition = "output.splsda_filter_active",
                             HTML("&nbsp;&nbsp;"),
                             actionButton("clear_splsda_filter", 
                                          "Reset Univariate to All Features", 
                                          class = "btn-warning btn-sm",
                                          icon = icon("times"))
                           )
                       )
                )
              ),
              tags$hr(),
              fluidRow(
                column(6, 
                       plotOutput("splsda_loadings1", height = "400px"),
                       downloadButton("download_splsda_load1", "PDF")
                ),
                column(6,
                       plotOutput("splsda_loadings2", height = "400px"),
                       downloadButton("download_splsda_load2", "PDF")
                )
              ),
              tags$hr(),
              h5("Heatmap of Selected Biomarkers"),
              fluidRow(
                column(6, selectInput("splsda_heatmap_dist", "Distance measure:",
                                       choices = c("Euclidean" = "euclidean",
                                                   "Correlation (Pearson)" = "correlation",
                                                   "Manhattan" = "manhattan",
                                                   "Maximum" = "maximum",
                                                   "Canberra" = "canberra",
                                                   "Binary" = "binary",
                                                   "Minkowski" = "minkowski"),
                                       selected = "euclidean")),
                column(6, selectInput("splsda_heatmap_clust", "Clustering method:",
                                       choices = c("Complete" = "complete",
                                                   "Ward (D)" = "ward.D",
                                                   "Average (UPGMA)" = "average",
                                                   "Single" = "single",
                                                   "McQuitty (WPGMA)" = "mcquitty",
                                                   "Median (WPGMC)" = "median",
                                                   "Centroid (UPGMC)" = "centroid"),
                                       selected = "complete"))
              ),
              plotOutput("splsda_cim", height = "700px"),
              fluidRow(
                column(6, downloadButton("download_splsda_cim", "Download Heatmap (PDF)")),
                column(6, downloadButton("download_splsda_cim_data", "Download Heatmap Data (CSV)"))
              )
            ),

            tabPanel(
              "PLS-R",
              conditionalPanel(
                condition = "input.plsr_y_col == ''",
                div(class = "alert alert-warning",
                    icon("exclamation-triangle"),
                    " Select a numeric response variable (Y) in the sidebar (section 7) to run PLS Regression.")
              ),
              h4("PLS Regression (continuous response)"),
              fluidRow(
                column(6,
                       h5("Scores (colored by Y)"),
                       fluidRow(
                         column(6, numericInput("plsr_x", "Comp X:", value = 1, min = 1, max = 15)),
                         column(6, numericInput("plsr_y_comp", "Comp Y:", value = 2, min = 1, max = 15))
                       ),
                       plotOutput("plsr_scores", height = "500px"),
                       downloadButton("download_plsr_scores", "PDF")
                ),
                column(6,
                       h5("Arrow Plot (X vs Y latent space)"),
                       plotOutput("plsr_arrow", height = "500px"),
                       downloadButton("download_plsr_arrow", "PDF")
                )
              ),
              tags$hr(),
              fluidRow(
                column(6,
                       h5("Correlation Circle"),
                       plotOutput("plsr_varplot", height = "450px"),
                       downloadButton("download_plsr_varplot", "PDF")
                ),
                column(6,
                       h5("Fitted vs Observed (final model)"),
                       plotOutput("plsr_fitted", height = "450px"),
                       downloadButton("download_plsr_fitted", "PDF")
                )
              ),
              tags$hr(),
              fluidRow(
                column(6,
                       h5("Performance (cross-validated)"),
                       selectInput("plsr_cv_method", "CV:",
                                   choices = c("M-fold" = "Mfold", "LOO" = "loo")),
                       conditionalPanel(
                         condition = "input.plsr_cv_method == 'Mfold'",
                         numericInput("plsr_folds", "Folds:", value = 5, min = 2, max = 10),
                         numericInput("plsr_nrepeat", "Repeats:", value = if (IS_WEBR) 10 else 50, min = 3, max = 200)
                       ),
                       selectInput("plsr_perf_criterion", "Criterion:",
                                   choices = c("MSEP" = "MSEP", "R2" = "R2", "Q2.total" = "Q2.total")),
                       actionButton("run_plsr_perf", "Performance", class = "btn-info"),
                       tags$hr(),
                       plotOutput("plsr_perf_plot", height = "400px"),
                       downloadButton("download_plsr_perf", "PDF"),
                       tags$hr(),
                       h5("Quality metrics (cross-validated)"),
                       verbatimTextOutput("plsr_perf_summary")
                ),
                column(6,
                       h5("Loadings & VIP"),
                       numericInput("plsr_loadings_comp", "Component:", value = 1, min = 1, max = 15),
                       plotOutput("plsr_loadings", height = "350px"),
                       downloadButton("download_plsr_loadings", "PDF"),
                       tags$hr(),
                       plotOutput("plsr_vip_plot", height = "400px"),
                       fluidRow(
                         column(6, downloadButton("download_plsr_vip", "PDF")),
                         column(6, downloadButton("download_plsr_vip_csv", "VIP table (CSV)"))
                       )
                )
              )
            ),

            tabPanel(
              "sPLS-R",
              conditionalPanel(
                condition = "input.splsr_mode == 'auto'",
                plotOutput("splsr_tune_plot", height = "400px"),
                downloadButton("download_splsr_tune", "PDF"),
                verbatimTextOutput("splsr_tune_summary"),
                tags$hr()
              ),
              fluidRow(
                column(6,
                       plotOutput("splsr_scores", height = "500px"),
                       downloadButton("download_splsr_scores", "PDF")
                ),
                column(6,
                       verbatimTextOutput("splsr_selected_features"),
                       downloadButton("download_splsr_features", "CSV")
                )
              ),
              tags$hr(),
              fluidRow(
                column(12,
                       h5("Final Model Performance (cross-validated)"),
                       actionButton("run_splsr_final_perf", "Evaluate Final Model", class = "btn-info"),
                       tags$hr(),
                       fluidRow(
                         column(6,
                                plotOutput("splsr_final_perf_plot", height = "350px"),
                                downloadButton("download_splsr_final_perf", "PDF")
                         ),
                         column(6, verbatimTextOutput("splsr_final_perf_summary"))
                       )
                )
              ),
              tags$hr(),
              fluidRow(
                column(12,
                       div(class = "well well-sm", style = "background-color: #fff8f0;",
                           h5(icon("exchange-alt"), " Cross-Analysis with Univariate"),
                           p("Examine univariate statistics for sPLS-R-selected features."),
                           actionButton("send_to_univariate_splsr",
                                        "Show sPLS-R Features in Univariate",
                                        class = "btn-info",
                                        icon = icon("arrow-left")),
                           conditionalPanel(
                             condition = "output.splsr_filter_active",
                             HTML("&nbsp;&nbsp;"),
                             actionButton("clear_splsr_filter",
                                          "Reset Univariate to All Features",
                                          class = "btn-warning btn-sm",
                                          icon = icon("times"))
                           )
                       )
                )
              ),
              tags$hr(),
              fluidRow(
                column(6,
                       plotOutput("splsr_loadings1", height = "400px"),
                       downloadButton("download_splsr_load1", "PDF")
                ),
                column(6,
                       plotOutput("splsr_loadings2", height = "400px"),
                       downloadButton("download_splsr_load2", "PDF")
                )
              ),
              tags$hr(),
              h5("Heatmap of Selected Features (columns colored by Y)"),
              fluidRow(
                column(6, selectInput("splsr_heatmap_dist", "Distance measure:",
                                       choices = c("Euclidean" = "euclidean",
                                                   "Correlation (Pearson)" = "correlation",
                                                   "Manhattan" = "manhattan",
                                                   "Maximum" = "maximum",
                                                   "Canberra" = "canberra",
                                                   "Binary" = "binary",
                                                   "Minkowski" = "minkowski"),
                                       selected = "euclidean")),
                column(6, selectInput("splsr_heatmap_clust", "Clustering method:",
                                       choices = c("Complete" = "complete",
                                                   "Ward (D)" = "ward.D",
                                                   "Average (UPGMA)" = "average",
                                                   "Single" = "single",
                                                   "McQuitty (WPGMA)" = "mcquitty",
                                                   "Median (WPGMC)" = "median",
                                                   "Centroid (UPGMC)" = "centroid"),
                                       selected = "complete"))
              ),
              plotOutput("splsr_cim", height = "700px"),
              fluidRow(
                column(6, downloadButton("download_splsr_cim", "Download Heatmap (PDF)")),
                column(6, downloadButton("download_splsr_cim_data", "Download Heatmap Data (CSV)"))
              )
            )
          )
        )
      )
    )
  ),

  div(
    class = "text-muted text-center small",
    style = "margin: 25px 0 15px 0; padding-top: 10px; border-top: 1px solid #ddd;",
    HTML("&#9888; For educational purposes only. This software is provided \"as is\", without warranty of any kind: the authors do not guarantee the accuracy of the results and accept no liability for their use. Always verify results independently before using them in research, publications or regulated contexts."),
    tags$br(),
    "\u00a9 2025-2026 Guillaume Marti, University of Toulouse - ",
    tags$a("CeCILL-B licence", href = "https://cecill.info/licences/Licence_CeCILL-B_V1-en.html", target = "_blank"),
    if (IS_WEBR) " - Your data are processed locally in your browser and are never uploaded."
  )
)

#==============================================================================
# UTILS
#==============================================================================

# Separator detected from the header line: the candidate (tab, ;, ,) that occurs
# most often wins. The header is used on purpose: data lines may contain decimal
# commas ("1,5;2,5"), which would otherwise be mistaken for column separators.
detect_separator <- function(path) {
  header <- readLines(path, n = 1, encoding = "UTF-8", warn = FALSE)
  if (length(header) == 0) return(",")
  counts <- c("\t" = 0, ";" = 0, "," = 0)
  for (sep in names(counts)) {
    counts[sep] <- lengths(regmatches(header, gregexpr(sep, header, fixed = TRUE)))
  }
  if (max(counts) == 0) "," else names(counts)[which.max(counts)]
}

read_table_file <- function(path, sep) {
  fread(path, header = TRUE, sep = sep, data.table = FALSE, encoding = "UTF-8",
        check.names = FALSE, integer64 = "double", stringsAsFactors = FALSE)
}

# Grouping column to pre-select: a column named like a group / class / treatment,
# otherwise the second text column (the first one is usually the sample ID).
guess_group_col <- function(df) {
  text_cols <- names(df)[!vapply(df, is.numeric, logical(1))]
  if (length(text_cols) == 0) return(names(df)[1])
  by_name <- text_cols[grepl("^(class|group|treatment|condition|type|chemotype|genotype|phenotype)",
                             text_cols, ignore.case = TRUE)]
  if (length(by_name) > 0) return(by_name[1])
  if (length(text_cols) >= 2) text_cols[2] else text_cols[1]
}

# Numeric columns used as features: the grouping column, the PLS-R response and
# a leading sample-index column (unique integers) are NOT features - otherwise
# they would be normalized, scaled and analysed like metabolites.
get_feature_cols <- function(df, group_col = NULL, y_col = NULL) {
  numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  excluded <- c(group_col, y_col)
  if (ncol(df) > 0) {
    first <- df[[1]]
    if (is.numeric(first) && !anyNA(first) && !anyDuplicated(first) && all(first == round(first))) {
      excluded <- c(excluded, names(df)[1])
    }
  }
  setdiff(numeric_cols, excluded)
}

# Pick the metadata column whose values best match the data's feature (column)
# names, instead of always defaulting to the first metadata column. Datasets
# commonly carry several ID-like columns (e.g. "AnnoID" vs "Alignment ID") and
# only one of them lines up with the feature names actually used in the data.
guess_feature_id_col <- function(metadata_df, feature_names) {
  default_col <- names(metadata_df)[1]

  if (is.null(feature_names) || length(feature_names) == 0) return(default_col)

  match_rates <- sapply(metadata_df, function(col) {
    mean(feature_names %in% as.character(col))
  })

  best_col <- names(which.max(match_rates))
  if (length(best_col) > 0 && match_rates[best_col] > 0) best_col else default_col
}

# Replace missing values by min/10 of their column (detection-limit estimate)
impute_min10 <- function(X) {
  X <- as.data.frame(X)
  for (j in seq_along(X)) {
    v <- X[[j]]
    if (anyNA(v)) {
      col_min <- min(v, na.rm = TRUE)
      X[[j]][is.na(v)] <- if (col_min > 0) col_min / 10 else abs(col_min) / 10 + 0.01
    }
  }
  X
}

# Scaling of a samples x features table. UV / Pareto drop constant features.
scale_features <- function(X, method) {
  M <- as.matrix(X)
  if (method %in% c("uv", "pareto")) {
    sds <- apply(M, 2, sd, na.rm = TRUE)
    ok <- is.finite(sds) & sds > 1e-10
    if (sum(ok) < 2) stop("not enough variables with variance for ", toupper(method), " scaling")
    M <- M[, ok, drop = FALSE]
    sds <- sds[ok]
    M <- sweep(M, 2, colMeans(M, na.rm = TRUE), "-")
    # Pareto = mean-centering and division by sqrt(sd)
    M <- sweep(M, 2, if (method == "uv") sds else sqrt(sds), "/")
  } else if (method == "center") {
    M <- sweep(M, 2, colMeans(M, na.rm = TRUE), "-")
  }
  M[!is.finite(M)] <- 0
  as.data.frame(M)
}

# Number of cross-validation folds that every class (or sample) can support
cv_folds <- function(requested, groups = NULL, n = NULL) {
  limit <- if (!is.null(groups)) min(table(groups)) else floor(n / 2)
  max(2, min(requested, limit))
}

# Pulls the per-component mean/sd out of one measure of a mixOmics perf.pls
# result (perf_res$measures$R2, $Q2, $MSEP, ...) for a clean printed table.
format_pls_perf_measure <- function(measure_obj, digits = 4) {
  if (is.null(measure_obj) || is.null(measure_obj$summary)) return(NULL)
  df <- measure_obj$summary[, c("comp", "mean", "sd")]
  df$mean <- round(df$mean, digits)
  df$sd <- round(df$sd, digits)
  df
}

#==============================================================================
# SERVER
#==============================================================================

server <- function(input, output, session) {
  
  rv <- reactiveValues(
    raw_data = NULL,
    data_path = NULL,
    data_sep = NULL,
    metadata = NULL,
    meta_path = NULL,
    meta_sep = NULL,
    groups = NULL,
    group_col = NULL,
    data_validated = FALSE,
    metadata_validated = FALSE,
    plsda_perf = NULL,
    splsda_tune_result = NULL,
    splsda_final_perf = NULL,
    plsr_perf = NULL,
    splsr_tune_result = NULL,
    splsr_final_perf = NULL,
    selected_groups = NULL,
    groups_applied = FALSE,
    # Feature filtering. The two directions are kept separate so that the
    # "sPLS-DA -> Univariate" view never changes the multivariate models (and
    # vice versa).
    mv_features = NULL,   # restricts the multivariate analyses (set from the Univariate tab)
    uv_features = NULL,   # restricts the univariate view (set from sPLS-DA / sPLS-R)
    uv_source = NULL      # "sPLS-DA" or "sPLS-R"
  )

  # Clear every cached model / tuning / performance result (they are only valid
  # for the data, groups, preprocessing and feature set they were computed on).
  reset_results <- function() {
    rv$plsda_perf <- NULL
    rv$splsda_tune_result <- NULL
    rv$splsda_final_perf <- NULL
    rv$plsr_perf <- NULL
    rv$splsr_tune_result <- NULL
    rv$splsr_final_perf <- NULL
  }

  # Reproducible cross-validation / tuning / permutations
  seed_value <- function() {
    s <- suppressWarnings(as.integer(input$seed))
    if (length(s) == 1 && !is.na(s)) s else 123L
  }
  apply_seed <- function() set.seed(seed_value())

  # Output for conditional panel
  output$metadata_loaded <- reactive({
    !is.null(rv$metadata) && rv$metadata_validated
  })
  outputOptions(output, "metadata_loaded", suspendWhenHidden = FALSE)

  # Status indicators for filtering
  output$univ_filter_active <- reactive({ !is.null(rv$mv_features) })
  outputOptions(output, "univ_filter_active", suspendWhenHidden = FALSE)

  output$uv_filter_active <- reactive({ !is.null(rv$uv_features) })
  outputOptions(output, "uv_filter_active", suspendWhenHidden = FALSE)

  output$splsda_filter_active <- reactive({
    !is.null(rv$uv_features) && identical(rv$uv_source, "sPLS-DA")
  })
  outputOptions(output, "splsda_filter_active", suspendWhenHidden = FALSE)

  output$splsr_filter_active <- reactive({
    !is.null(rv$uv_features) && identical(rv$uv_source, "sPLS-R")
  })
  outputOptions(output, "splsr_filter_active", suspendWhenHidden = FALSE)

  output$univ_filter_count <- renderText(as.character(length(rv$mv_features)))
  output$uv_filter_count <- renderText(as.character(length(rv$uv_features)))
  output$uv_filter_source <- renderText(if (is.null(rv$uv_source)) "" else rv$uv_source)

  #----------------------------------------------------------------------------
  # File import (separator auto-detected from the header line)
  #----------------------------------------------------------------------------

  # Called whenever a new statistical data table is read (upload, example,
  # separator change). A new table always has to be validated again.
  set_main_data <- function(df, path, sep, group_hint = NULL) {
    names(df) <- make.unique(ifelse(is.na(names(df)) | names(df) == "", "V", names(df)))
    rv$raw_data <- df
    rv$data_path <- path
    rv$data_sep <- sep
    rv$data_validated <- FALSE
    rv$metadata_validated <- FALSE
    rv$groups <- NULL
    rv$group_col <- NULL
    rv$selected_groups <- NULL
    rv$mv_features <- NULL
    rv$uv_features <- NULL
    rv$uv_source <- NULL
    reset_results()

    num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
    group_guess <- if (!is.null(group_hint) && group_hint %in% names(df)) group_hint else guess_group_col(df)
    updateRadioButtons(session, "sep", selected = sep)
    updateSelectInput(session, "group_col", choices = names(df), selected = group_guess)
    updateSelectInput(session, "plsr_y_col", choices = c("None" = "", num_cols))

    # If metadata was already uploaded, re-guess the best-matching ID column
    # now that the data's feature names are known too.
    if (!is.null(rv$metadata)) {
      updateSelectInput(session, "feature_id_col", choices = names(rv$metadata),
                        selected = guess_feature_id_col(rv$metadata, num_cols))
    }
  }

  load_main_data <- function(path, group_hint = NULL) {
    tryCatch({
      sep <- detect_separator(path)
      df <- tryCatch(read_table_file(path, sep),
                     error = function(e) read_table_file(path, input$sep))

      if (ncol(df) < 2) {
        showNotification("File appears to have parsing issues. Please check the separator.", type = "warning")
      }
      if (nrow(df) < 5) {
        showNotification("File too small (at least 5 samples are required)", type = "error")
        return(invisible(FALSE))
      }

      set_main_data(df, path, sep, group_hint)
      showNotification(paste("Data file loaded:", nrow(df), "rows,", ncol(df), "columns"), type = "message")
      invisible(TRUE)
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      invisible(FALSE)
    })
  }

  set_metadata <- function(md, path, sep) {
    names(md) <- make.unique(ifelse(is.na(names(md)) | names(md) == "", "V", names(md)))
    rv$metadata <- md
    rv$meta_path <- path
    rv$meta_sep <- sep

    feature_names <- if (is.null(rv$raw_data)) NULL else names(rv$raw_data)[vapply(rv$raw_data, is.numeric, logical(1))]
    numeric_cols <- names(md)[vapply(md, is.numeric, logical(1))]
    character_cols <- names(md)[vapply(md, function(x) is.character(x) || is.factor(x), logical(1))]

    updateRadioButtons(session, "meta_sep", selected = sep)
    updateSelectInput(session, "feature_id_col", choices = names(md),
                      selected = guess_feature_id_col(md, feature_names))
    updateSelectInput(session, "heatmap_annotation_col", choices = c("None" = "", character_cols))
    updateSelectInput(session, "volcano_size_col", choices = c("None" = "", numeric_cols))
    # Any metadata column can serve as label (name, class, m/z, ...)
    updateSelectInput(session, "feature_label_col", choices = c("Feature ID" = "", names(md)),
                      selected = "")
  }

  load_metadata <- function(path) {
    tryCatch({
      sep <- detect_separator(path)
      md <- tryCatch(read_table_file(path, sep),
                     error = function(e) read_table_file(path, input$meta_sep))
      if (ncol(md) < 2) {
        showNotification("Metadata parsing issue detected. Please check the separator.", type = "warning")
      }
      set_metadata(md, path, sep)
      showNotification(paste("Metadata file loaded:", nrow(md), "features,", ncol(md), "columns"), type = "message")
      invisible(TRUE)
    }, error = function(e) {
      showNotification(paste("Metadata error:", e$message), type = "error")
      invisible(FALSE)
    })
  }

  observeEvent(input$file1, {
    req(input$file1)
    load_main_data(input$file1$datapath)
  })

  # Re-parse data when the user changes the separator by hand (the automatic
  # update of the radio button after detection is ignored: same separator).
  observeEvent(input$sep, {
    req(rv$data_path)
    if (identical(input$sep, rv$data_sep)) return(NULL)

    tryCatch({
      df <- read_table_file(rv$data_path, input$sep)
      if (ncol(df) < 2) {
        showNotification("Parsing issue detected. Try a different separator.", type = "warning")
      } else {
        set_main_data(df, rv$data_path, input$sep)
        showNotification(paste("Data re-loaded with new separator:", nrow(df), "rows,", ncol(df), "columns"), type = "message")
      }
    }, error = function(e) {
      showNotification(paste("Error with selected separator:", e$message), type = "error")
    })
  })

  observeEvent(input$metadata_file, {
    req(input$metadata_file)
    load_metadata(input$metadata_file$datapath)
  })

  observeEvent(input$meta_sep, {
    req(rv$meta_path)
    if (identical(input$meta_sep, rv$meta_sep)) return(NULL)

    tryCatch({
      md <- read_table_file(rv$meta_path, input$meta_sep)
      if (ncol(md) < 2) {
        showNotification("Metadata parsing issue. Try a different separator.", type = "warning")
      } else {
        set_metadata(md, rv$meta_path, input$meta_sep)
        showNotification(paste("Metadata re-loaded:", nrow(md), "features,", ncol(md), "columns"), type = "message")
      }
    }, error = function(e) {
      showNotification(paste("Error with selected separator:", e$message), type = "error")
    })
  })

  # Example datasets shipped next to the app
  output$example_ui <- renderUI({
    available <- vapply(EXAMPLES, function(e) file.exists(e$data), logical(1))
    if (!any(available)) return(NULL)
    tagList(
      selectInput("example_choice", "Example dataset:", choices = names(EXAMPLES)[available]),
      actionButton("load_example", "Load example data", class = "btn-default btn-sm",
                   style = "width: 100%; margin-bottom: 15px;")
    )
  })

  observeEvent(input$load_example, {
    ex <- EXAMPLES[[input$example_choice]]
    req(ex, file.exists(ex$data))

    if (isTRUE(load_main_data(ex$data, group_hint = ex$group))) {
      if (!is.null(ex$meta) && file.exists(ex$meta)) {
        load_metadata(ex$meta)
      } else {
        rv$metadata <- NULL
        rv$meta_path <- NULL
      }
      # Validate once the group selector has been updated in the browser
      shinyjs::delay(800, shinyjs::click("validate_data"))
    }
  })

  # Metadata matching check
  output$metadata_match_info <- renderPrint({
    req(rv$raw_data, rv$metadata, input$feature_id_col, rv$group_col)

    feature_names <- get_feature_cols(rv$raw_data, rv$group_col)
    metadata_ids <- trimws(as.character(rv$metadata[[input$feature_id_col]]))

    n_features <- length(feature_names)
    n_metadata <- length(metadata_ids)
    n_matched <- sum(trimws(feature_names) %in% metadata_ids)

    cat("Data features:", n_features, "\n")
    cat("Metadata entries:", n_metadata, "\n")
    cat("Matched features:", n_matched, "\n")
    cat("Match percentage:", round(n_matched/n_features * 100, 1), "%\n")
    if (anyDuplicated(metadata_ids)) {
      cat("Warning: duplicated IDs in the metadata (the first entry is used).\n")
    }

    if(n_matched < n_features) {
      unmatched <- feature_names[!trimws(feature_names) %in% metadata_ids]
      cat("\nFirst 10 unmatched features:\n")
      cat(head(unmatched, 10), sep = ", ")
    }
  })

  # Validation
  observeEvent(input$validate_data, {
    req(rv$raw_data, input$group_col)

    tryCatch({
      groups <- factor(rv$raw_data[[input$group_col]])

      if (length(levels(groups)) < 2) {
        showNotification("At least 2 groups required", type = "error")
        return(NULL)
      }

      feats <- get_feature_cols(rv$raw_data, input$group_col)
      if (length(feats) < 2) {
        showNotification("At least 2 numeric feature columns required", type = "error")
        return(NULL)
      }

      excluded_id <- setdiff(names(rv$raw_data)[vapply(rv$raw_data, is.numeric, logical(1))],
                             c(feats, input$group_col))
      if (length(excluded_id) > 0) {
        showNotification(paste0("Column '", excluded_id[1],
                                "' looks like a sample index (unique integers) and is not used as a feature."),
                         type = "message", duration = 8)
      }
      if (anyNA(groups)) {
        showNotification(paste(sum(is.na(groups)), "sample(s) without group are ignored."), type = "warning")
      }
      tab <- table(groups)
      small <- names(tab)[tab < 3]
      if (length(small) > 0) {
        showNotification(paste("Groups with fewer than 3 samples (cross-validation will be unreliable):",
                               paste(small, collapse = ", ")), type = "warning", duration = 10)
      }

      rv$groups <- groups
      rv$group_col <- input$group_col
      rv$selected_groups <- levels(groups)
      rv$data_validated <- TRUE
      rv$groups_applied <- FALSE
      rv$mv_features <- NULL
      rv$uv_features <- NULL
      rv$uv_source <- NULL
      reset_results()

      # Validate metadata if present
      if(!is.null(rv$metadata) && !is.null(input$feature_id_col)) {
        rv$metadata_validated <- TRUE
        showNotification("Data and metadata validated", type = "message")
      } else {
        rv$metadata_validated <- FALSE
        showNotification(paste("Validated:", nrow(rv$raw_data), "samples,",
                               length(feats), "variables"), type = "message")
      }

    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })

  # Get matched metadata (one row per data feature, NA when the feature has no
  # metadata; duplicated metadata IDs are resolved to the first entry)
  matched_metadata <- reactive({
    if(is.null(rv$metadata) || !rv$metadata_validated || is.null(rv$raw_data) ||
       is.null(rv$group_col) || is.null(input$feature_id_col) || input$feature_id_col == "") {
      return(NULL)
    }

    feature_names <- get_feature_cols(rv$raw_data, rv$group_col)
    metadata_df <- rv$metadata
    clean_metadata_ids <- trimws(as.character(metadata_df[[input$feature_id_col]]))

    idx <- match(trimws(feature_names), clean_metadata_ids)
    matched_df <- metadata_df[idx, , drop = FALSE]
    matched_df$feature_name <- feature_names
    rownames(matched_df) <- feature_names

    # Report matching statistics
    n_matched <- sum(!is.na(idx))
    if(n_matched == 0) {
      showNotification("Warning: No features matched between data and metadata. Check feature ID column selection.",
                       type = "warning", duration = 10)
    } else if(n_matched < length(feature_names) * 0.5) {
      showNotification(paste("Warning: Only", n_matched, "out of", length(feature_names),
                             "features matched with metadata."),
                       type = "warning", duration = 10)
    }

    return(matched_df)
  })

  # Display labels for feature IDs: the metadata column chosen in sidebar
  # section 5 (truncated), or the ID itself when no name is available.
  # Duplicated names get their ID appended so that labels stay unique.
  feature_labels <- function(ids) {
    ids <- as.character(ids)
    md <- matched_metadata()
    col <- input$feature_label_col
    if (is.null(md) || is.null(col) || col == "" || !(col %in% names(md))) return(ids)

    v <- md[ids, col]
    if (is.numeric(v)) v <- ifelse(is.na(v), NA, format(signif(v, 6), trim = TRUE))
    v <- trimws(as.character(v))
    ok <- !is.na(v) & v != "" & !(tolower(v) %in% c("unknown", "na", "nan", "null", "none"))
    max_chars <- suppressWarnings(as.integer(input$feature_label_max))
    if (length(max_chars) != 1 || is.na(max_chars) || max_chars < 5) max_chars <- 10
    long <- ok & nchar(v) > max_chars
    v[long] <- paste0(substr(v[long], 1, max_chars - 3), "...")

    lab <- ifelse(ok, v, ids)
    dup <- ok & (lab %in% lab[duplicated(lab)])
    lab[dup] <- paste0(lab[dup], " [", ids[dup], "]")
    lab
  }

  # Plot title: "Name (ID)" when a label is active
  feature_title <- function(id) {
    lab <- feature_labels(id)
    if (identical(lab, as.character(id))) id else paste0(lab, " (", id, ")")
  }

  # Adds a "Label" column after "Variable" when labels are active
  add_label_column <- function(df) {
    if (is.null(input$feature_label_col) || input$feature_label_col == "" || is.null(matched_metadata())) return(df)
    df <- as.data.frame(df)
    df$Label <- feature_labels(df$Variable)
    df[, c("Variable", "Label", setdiff(names(df), c("Variable", "Label")))]
  }

  # Dynamic UI for group selection
  output$group_selector_ui <- renderUI({
    req(rv$groups)
    
    all_groups <- levels(rv$groups)
    
    checkboxGroupInput(
      "selected_groups",
      "Select groups to include in analysis:",
      choices = all_groups,
      selected = if(is.null(rv$selected_groups)) all_groups else rv$selected_groups,
      inline = FALSE
    )
  })
  
  # Observer for group selection changes
  observeEvent(input$selected_groups, {
    req(input$selected_groups)
    if(length(input$selected_groups) >= 2) {
      rv$selected_groups <- input$selected_groups
    } else if(length(input$selected_groups) == 1) {
      showNotification("At least 2 groups must be selected", type = "warning")
    }
  }, ignoreNULL = FALSE)
  
  # Keep the volcano group selectors in sync with rv$selected_groups, whichever
  # code path changed it (validation, reset, live checkbox toggling, or apply) -
  # this is what populates them the first time, right after data validation,
  # instead of leaving them empty until "Apply" is clicked.
  observe({
    req(rv$selected_groups)

    selected1 <- if (isTRUE(input$volcano_group1 %in% rv$selected_groups)) {
      input$volcano_group1
    } else {
      rv$selected_groups[1]
    }
    selected2 <- if (isTRUE(input$volcano_group2 %in% rv$selected_groups) && input$volcano_group2 != selected1) {
      input$volcano_group2
    } else {
      rv$selected_groups[rv$selected_groups != selected1][1]
    }

    updateSelectInput(session, "volcano_group1", choices = rv$selected_groups, selected = selected1)
    updateSelectInput(session, "volcano_group2", choices = rv$selected_groups, selected = selected2)
  })

  # Reset button observer
  observeEvent(input$reset_groups, {
    req(rv$groups)
    rv$selected_groups <- levels(rv$groups)

    showNotification("Groups selection reset to all groups", type = "message")
  })

  # Apply button observer
  observeEvent(input$apply_groups, {
    req(rv$selected_groups)

    if(length(rv$selected_groups) < 2) {
      showNotification("At least 2 groups must be selected", type = "error")
      return(NULL)
    }

    # Force reactivity by invalidating all models
    rv$plsda_perf <- NULL
    rv$splsda_tune_result <- NULL
    rv$groups_applied <- !rv$groups_applied  # Toggle to force reactivity

    showNotification(paste("Analysis updated with", length(rv$selected_groups),
                           "groups (", paste(rv$selected_groups, collapse = ", "), ")"),
                     type = "message", duration = 5)
    
    # Auto-trigger sPLS-DA tuning if in auto mode and on Multivariate tab
    if(input$main_tabs == "Multivariate" && input$splsda_mode == "auto") {
      showNotification("Launching sPLS-DA tuning automatically...", type = "message", duration = 3)
      # Trigger the tuning button with a small delay to ensure data is ready
      shinyjs::delay(500, shinyjs::click("run_splsda_tuning"))
    }
  })
  
  # Filtered groups
  groups_filtered <- reactive({
    req(rv$groups, rv$selected_groups)
    
    selected_idx <- rv$groups %in% rv$selected_groups
    groups_filt <- droplevels(factor(rv$groups[selected_idx]))
    
    return(groups_filt)
  })

  # Sample names (first text column that is not the grouping column), used as
  # row names so that plots and exported tables show the real sample IDs.
  sample_ids <- reactive({
    req(rv$raw_data, rv$group_col)
    df <- rv$raw_data
    candidates <- setdiff(names(df)[!vapply(df, is.numeric, logical(1))], rv$group_col)
    ids <- if (length(candidates) > 0) as.character(df[[candidates[1]]]) else as.character(seq_len(nrow(df)))
    missing_id <- is.na(ids) | ids == ""
    ids[missing_id] <- paste0("Sample_", which(missing_id))
    make.unique(ids)
  })

  # Preprocessing, step 1 (sample-wise, independent of the group selection):
  # missing values -> normalization -> transformation. The grouping column, the
  # PLS-R response and any sample-index column are never treated as features.
  # Returns the transformed matrix plus the normalized (untransformed) abundances
  # used for the true fold change in the volcano plot.
  X_pre <- reactive({
    req(rv$data_validated, rv$raw_data, rv$group_col)

    tryCatch({
      feats <- get_feature_cols(rv$raw_data, rv$group_col, input$plsr_y_col)
      if (length(feats) < 2) {
        showNotification("At least 2 numeric feature columns are required", type = "error")
        return(NULL)
      }

      X <- rv$raw_data[, feats, drop = FALSE]
      rownames(X) <- sample_ids()

      na_prop <- colMeans(is.na(X))
      n_vars_before <- ncol(X)
      max_na_allowed <- if (!is.null(input$max_na_prop)) input$max_na_prop / 100 else 0.5
      impute_method <- if (!is.null(input$impute_method)) input$impute_method else "none"

      X <- X[, na_prop <= max_na_allowed, drop = FALSE]
      if (ncol(X) < 2) {
        showNotification("Too many missing values. Need at least 2 valid variables.", type = "error")
        return(NULL)
      }
      if (ncol(X) < n_vars_before) {
        showNotification(paste("Removed", n_vars_before - ncol(X), "variables with too many missing values"),
                         type = "message")
      }

      if (anyNA(X)) {
        if (impute_method == "none") {
          X[is.na(X)] <- 0
        } else {
          X_min10 <- impute_min10(X)
          if (impute_method == "min10") {
            X <- X_min10
          } else {
            # NIPALS: only the missing cells are replaced (observed values are
            # kept); non-positive imputed values fall back to min/10.
            X_matrix <- as.matrix(X)
            ncomp_nipals <- max(1, min(5, nrow(X_matrix) - 2, ncol(X_matrix) - 1))
            imputed <- withProgress(message = "NIPALS imputation in progress...", value = 0.5, {
              tryCatch(
                suppressWarnings(suppressMessages(mixOmics::impute.nipals(X_matrix, ncomp = ncomp_nipals))),
                error = function(e) NULL
              )
            })
            if (is.null(imputed)) {
              showNotification("NIPALS failed, using min/10 fallback", type = "warning")
              X <- X_min10
            } else {
              imputed <- as.matrix(imputed)
              bad <- is.na(X_matrix) & !(is.finite(imputed) & imputed > 0)
              imputed[bad] <- as.matrix(X_min10)[bad]
              dimnames(imputed) <- dimnames(X_matrix)
              X <- as.data.frame(imputed)
            }
          }
        }
      }

      M <- as.matrix(X)

      # Normalization / transformation need non-negative intensities
      if (any(M < 0) && (input$norm_method == "tss" || input$trans_method != "none")) {
        showNotification("Negative values found: they are set to 0 before normalization / transformation.",
                         type = "warning")
        M[M < 0] <- 0
      }

      # Normalization: TSS, rescaled to the median total signal so that the data
      # stay in intensity units (otherwise log10(x+1) of fractions ~ 1e-4 is
      # almost the identity and CLR pseudo-counts dominate the signal)
      if (input$norm_method == "tss") {
        row_sums <- rowSums(M)
        if (any(row_sums <= 0)) {
          showNotification("TSS skipped: at least one sample has a total signal of 0.", type = "warning")
        } else {
          M <- sweep(M, 1, row_sums, FUN = "/") * stats::median(row_sums)
        }
      }

      # Unscaled normalized abundances: true fold change in the volcano plot
      M_abund <- M

      # Transformation
      if (input$trans_method == "log10") {
        M <- log10(M + 1)
      } else if (input$trans_method == "clr") {
        positive <- M[M > 0]
        offset <- if (any(M <= 0) && length(positive) > 0) min(positive) / 2 else 0
        res <- mixOmics::logratio.transfo(M, logratio = "CLR", offset = offset)
        M <- matrix(as.numeric(res), nrow = nrow(M), dimnames = dimnames(M))
      }

      # Remove any column with NaN or Inf
      valid_cols <- apply(M, 2, function(col) all(is.finite(col)))
      M <- M[, valid_cols, drop = FALSE]
      M_abund <- M_abund[, valid_cols, drop = FALSE]

      if (ncol(M) < 2) {
        showNotification("Preprocessing resulted in too few valid variables", type = "error")
        return(NULL)
      }

      list(trans = as.data.frame(M), abund = as.data.frame(M_abund))

    }, error = function(e) {
      showNotification(paste("Preprocessing error:", e$message), type = "error")
      return(NULL)
    })
  })

  # Preprocessed data, all samples (Data tab preview / download only)
  X_processed <- reactive({
    pre <- X_pre()
    req(pre)
    tryCatch(
      scale_features(pre$trans, input$scale_method),
      error = function(e) {
        showNotification(paste("Scaling error:", e$message), type = "error")
        NULL
      }
    )
  })

  # Data used by every analysis: selected groups only, and scaled on THOSE
  # samples (UV / Pareto / centering must not be influenced by excluded groups).
  X_filtered <- reactive({
    pre <- X_pre()
    req(pre, rv$groups, rv$selected_groups)

    selected_idx <- rv$groups %in% rv$selected_groups

    if(sum(selected_idx) < 4) {
      showNotification("Not enough samples in selected groups (minimum 4)", type = "error")
      return(NULL)
    }

    tryCatch(
      scale_features(pre$trans[selected_idx, , drop = FALSE], input$scale_method),
      error = function(e) {
        showNotification(paste("Scaling error:", e$message), type = "warning")
        NULL
      }
    )
  })

  # Multivariate feature set (all features, or those sent from the Univariate tab)
  active_features <- reactive({
    req(X_filtered())

    if(!is.null(rv$mv_features) && length(rv$mv_features) > 0) {
      available_features <- intersect(rv$mv_features, colnames(X_filtered()))

      if(length(available_features) >= 2) {
        return(available_features)
      }
      showNotification("Not enough filtered features available. Using all features.",
                       type = "warning")
    }
    colnames(X_filtered())
  })

  # Feature-filtered data for the multivariate analyses
  X_feature_filtered <- reactive({
    req(X_filtered(), active_features())

    X_filtered()[, active_features(), drop = FALSE]
  })

  # Number of components that can really be computed on the current data
  ncomp_eff <- reactive({
    X <- X_feature_filtered()
    req(X, input$ncomp_mv)
    max(1, min(input$ncomp_mv, nrow(X) - 1, ncol(X)))
  })

  # Unscaled (normalized) abundances of the selected samples, for the true Fold-Change
  X_unscaled_filtered <- reactive({
    pre <- X_pre()
    req(pre, rv$groups, rv$selected_groups)
    selected_idx <- rv$groups %in% rv$selected_groups
    pre$abund[selected_idx, , drop = FALSE]
  })

  # Cached results are only valid for the data they were computed on
  observeEvent(X_feature_filtered(), { reset_results() }, ignoreInit = TRUE)
  observeEvent(input$ncomp_mv, { rv$plsda_perf <- NULL; rv$plsr_perf <- NULL }, ignoreInit = TRUE)
  observeEvent(list(rv$splsda_tune_result, input$splsda_mode, input$splsda_ncomp_manual,
                    input$splsda_keepx_manual), { rv$splsda_final_perf <- NULL }, ignoreInit = TRUE)
  observeEvent(list(rv$splsr_tune_result, input$splsr_mode, input$splsr_ncomp_manual,
                    input$splsr_keepx_manual, input$plsr_y_col), { rv$splsr_final_perf <- NULL }, ignoreInit = TRUE)

  #----------------------------------------------------------------------------
  # Feature filtering (cross-analysis)
  #----------------------------------------------------------------------------

  # Send ANOVA results to sPLS-DA
  observeEvent(input$send_to_splsda, {
    sig_features <- significant_features_all()

    if(length(sig_features) < 2) {
      showNotification("At least 2 significant features required", type = "error")
      return(NULL)
    }

    rv$mv_features <- sig_features

    showNotification(
      paste("The multivariate analyses will now use the", length(sig_features),
            "significant features from the univariate tests (p <", input$pval_threshold, ")"),
      type = "message",
      duration = 5
    )

    updateTabsetPanel(session, "main_tabs", selected = "Multivariate")
  })

  # Send sPLS-DA features to Univariate
  observeEvent(input$send_to_univariate, {
    req(final_splsda_model())
    
    tryCatch({
      model <- final_splsda_model()
      
      # Get all selected features from all components
      all_selected_features <- character()
      for (i in 1:model$ncomp) {
        selected_i <- mixOmics::selectVar(model, comp = i)
        all_selected_features <- c(all_selected_features, selected_i$name)
      }
      all_selected_features <- unique(all_selected_features)
      
      if(length(all_selected_features) < 2) {
        showNotification("At least 2 features required", type = "error")
        return(NULL)
      }
      
      rv$uv_features <- all_selected_features
      rv$uv_source <- "sPLS-DA"
      
      showNotification(
        paste("Univariate analysis will now focus on", length(all_selected_features), 
              "features selected by sPLS-DA"),
        type = "message",
        duration = 5
      )
      
      updateTabsetPanel(session, "main_tabs", selected = "Univariate")
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Send sPLS-R features to Univariate
  observeEvent(input$send_to_univariate_splsr, {
    req(final_splsr_model())

    tryCatch({
      model <- final_splsr_model()

      all_selected_features <- character()
      for (i in 1:model$ncomp) {
        selected_i <- mixOmics::selectVar(model, comp = i)$X
        all_selected_features <- c(all_selected_features, selected_i$name)
      }
      all_selected_features <- unique(all_selected_features)

      if(length(all_selected_features) < 2) {
        showNotification("At least 2 features required", type = "error")
        return(NULL)
      }

      rv$uv_features <- all_selected_features
      rv$uv_source <- "sPLS-R"

      showNotification(
        paste("Univariate analysis will now focus on", length(all_selected_features),
              "features selected by sPLS-R"),
        type = "message",
        duration = 5
      )

      updateTabsetPanel(session, "main_tabs", selected = "Univariate")

    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })

  # Clear filters (several buttons for the same action)
  clear_mv_filter <- function() {
    rv$mv_features <- NULL
    showNotification("Filter cleared. The multivariate analyses use all features.", type = "message")
  }
  clear_uv_filter <- function() {
    rv$uv_features <- NULL
    rv$uv_source <- NULL
    showNotification("Filter cleared. The univariate analysis shows all features.", type = "message")
  }

  observeEvent(input$clear_univ_filter, { clear_mv_filter() })
  observeEvent(input$clear_univ_filter_top, { clear_mv_filter() })
  observeEvent(input$clear_splsda_filter, { clear_uv_filter() })
  observeEvent(input$clear_splsda_filter_top, { clear_uv_filter() })
  observeEvent(input$clear_splsr_filter, { clear_uv_filter() })
  
  # Data preview
  output$data_summary <- renderPrint({
    req(rv$data_validated, rv$groups)
    
    X <- rv$raw_data[, get_feature_cols(rv$raw_data, rv$group_col, input$plsr_y_col), drop = FALSE]
    
    cat("Total samples:", nrow(X), "\n")
    cat("Total variables:", ncol(X), "\n")
    cat("All groups:", paste(levels(rv$groups), collapse = ", "), "\n")
    
    if(!is.null(rv$selected_groups) && length(rv$selected_groups) < length(levels(rv$groups))) {
      cat("\nSelected groups:", paste(rv$selected_groups, collapse = ", "), "\n")
      selected_idx <- rv$groups %in% rv$selected_groups
      cat("Selected samples:", sum(selected_idx), "\n")
    }
  })
  
  output$metadata_summary <- renderPrint({
    if(is.null(rv$metadata)) {
      cat("No metadata loaded\n")
    } else {
      cat("Metadata loaded:\n")
      cat("  Total features in metadata:", nrow(rv$metadata), "\n")
      cat("  Total columns:", ncol(rv$metadata), "\n")
      
      if(!is.null(matched_metadata())) {
        matched <- matched_metadata()
        n_matched <- sum(!is.na(matched[[input$feature_id_col]]))
        cat("  Matched features:", n_matched, "/", nrow(matched), "\n")
      }
    }
  })
  
  output$data_preview_raw <- DT::renderDataTable({
    req(rv$raw_data)
    head(rv$raw_data, 20)
  }, options = list(scrollX = TRUE))
  
  output$metadata_preview <- DT::renderDataTable({
    req(rv$metadata)
    head(rv$metadata, 20)
  }, options = list(scrollX = TRUE, pageLength = 10))
  
  output$data_preview_processed <- DT::renderDataTable({
    req(X_processed())
    head(X_processed(), 20)
  }, options = list(scrollX = TRUE))
  output$processed_note <- renderText(
    "Preview scaled over all samples. Each analysis re-scales the data on the selected groups only.")
  
  output$download_processed_data <- downloadHandler(
    filename = function() { "processed_data.csv" },
    content = function(file) {
      req(X_processed(), rv$groups)
      df_export <- cbind(Sample = rownames(X_processed()), Group = rv$groups, X_processed())
      write.csv(df_export, file, row.names = FALSE)
    },
    contentType = "text/csv"
  )
  
  # Group info displays
  output$selected_groups_info_univ <- renderText({
    req(groups_filtered())
    n_groups <- length(levels(groups_filtered()))
    n_samples <- length(groups_filtered())
    paste(" Analysis with", n_groups, "groups (", paste(levels(groups_filtered()), collapse = ", "), 
          ") and", n_samples, "samples")
  })
  
  output$selected_groups_info_multi <- renderText({
    req(groups_filtered())
    n_groups <- length(levels(groups_filtered()))
    n_samples <- length(groups_filtered())
    paste(" Analysis with", n_groups, "groups (", paste(levels(groups_filtered()), collapse = ", "), 
          ") and", n_samples, "samples")
  })
  
  # Univariate tests: always computed on ALL features of the selected samples, so
  # the multiple-testing correction is not weakened when the view is restricted
  # to features selected by sPLS-DA / sPLS-R.
  univariate_results <- reactive({
    req(X_filtered(), groups_filtered())
    
    withProgress(message = "Univariate tests...", value = 0, {
      X <- X_filtered()
      groups <- groups_filtered()
      
      results <- map_df(names(X), function(var_name) {
        incProgress(1/ncol(X))
        
        var_data <- X[[var_name]]
        test_data <- data.frame(value = var_data, group = groups)
        test_data <- test_data[complete.cases(test_data), ]
        
        if (nrow(test_data) < length(levels(groups)) * 2) {
          return(tibble(Variable = var_name, P_Value = NA))
        }
        
        pval <- tryCatch({
          if (input$use_nonparametric) {
            kruskal.test(value ~ group, data = test_data)$p.value
          } else {
            summary(aov(value ~ group, data = test_data))[[1]][["Pr(>F)"]][1]
          }
        }, error = function(e) NA)
        
        tibble(Variable = var_name, P_Value = pval)
      })
      
      results %>%
        filter(!is.na(P_Value)) %>%
        mutate(Adjusted_P_Value = p.adjust(P_Value, method = input$adj_method)) %>%
        arrange(Adjusted_P_Value)
    })
  })

  # Features below the p-value threshold (whatever the view filter)
  significant_features_all <- reactive({
    req(univariate_results())
    univariate_results() %>% filter(Adjusted_P_Value <= input$pval_threshold) %>% pull(Variable)
  })

  # Rows shown in the Univariate tab: significant features, or - when features
  # were sent from sPLS-DA / sPLS-R - those features whatever their p-value
  filtered_univariate_results <- reactive({
    req(univariate_results())
    if (!is.null(rv$uv_features)) {
      univariate_results() %>% filter(Variable %in% rv$uv_features)
    } else {
      univariate_results() %>% filter(Adjusted_P_Value <= input$pval_threshold)
    }
  })
  
  output$univariate_table <- DT::renderDataTable({
    req(filtered_univariate_results())
    
    results_df <- add_label_column(filtered_univariate_results())
    
    # Add metadata if available
    if(!is.null(matched_metadata()) && !is.null(input$heatmap_annotation_col) && 
       input$heatmap_annotation_col != "") {
      metadata_df <- matched_metadata()
      metadata_info <- metadata_df[results_df$Variable, input$heatmap_annotation_col, drop = FALSE]
      results_df <- cbind(results_df, metadata_info)
    }
    
    DT::formatSignif(
      DT::datatable(as.data.frame(results_df), selection = 'single', rownames = FALSE,
                    options = list(pageLength = 15, scrollX = TRUE)),
      columns = intersect(c("P_Value", "Adjusted_P_Value"), names(results_df)), digits = 3)
  })
  
  output$filtering_stats <- renderPrint({
    req(univariate_results(), filtered_univariate_results())
    cat("Variables tested:", nrow(univariate_results()), "\n")
    cat("Significant variables (adjusted p <=", input$pval_threshold, "):", length(significant_features_all()), "\n")
    if (!is.null(rv$uv_features)) {
      cat("Variables shown (selected by ", rv$uv_source, "): ", nrow(filtered_univariate_results()), "\n", sep = "")
    }
  })
  
  # Enhanced Heatmap with metadata annotation - single source of truth for the
  # cim() call args and the resulting clustered order, so the plot (drawn by
  # cim() itself, with real dendrograms) and the CSV/PDF downloads always agree.
  # cim() is run once headlessly here (to read back the rowInd/colInd it computed)
  # and re-run with the identical, deterministic arguments wherever it is drawn.
  univ_heatmap_ordered <- reactive({
    req(filtered_univariate_results(), X_filtered(), groups_filtered())

    if (nrow(filtered_univariate_results()) < 2) return(NULL)

    n_features <- min(input$n_heatmap_features, nrow(filtered_univariate_results()))
    top_vars <- filtered_univariate_results()$Variable[1:n_features]

    X_df <- X_filtered()
    mat <- as.matrix(X_df[, top_vars, drop = FALSE])
    mat <- scale(mat)
    mat[!is.finite(mat)] <- 0
    mat_t <- t(mat)                      # features (rows) x samples (cols)
    rownames(mat_t) <- feature_labels(top_vars)

    # Row (feature) side colors from metadata, in the SAME (pre-cluster) order
    # as mat_t's rows - cim() reorders these itself alongside the matrix rows.
    row_sideColors <- NULL
    feature_colors <- NULL

    if (!is.null(matched_metadata()) && !is.null(input$heatmap_annotation_col) &&
        input$heatmap_annotation_col != "") {
      metadata_df <- matched_metadata()
      annotation_data <- metadata_df[top_vars, input$heatmap_annotation_col]

      if (!all(is.na(annotation_data))) {
        unique_vals <- unique(na.omit(annotation_data))
        n_colors <- length(unique_vals)
        colors <- if (n_colors <= 12) {
          brewer.pal(min(n_colors, 12), "Set3")
        } else {
          viridis(n_colors)
        }
        color_map <- setNames(colors[1:n_colors], unique_vals)
        color_map["NA"] <- "grey80"

        row_sideColors <- sapply(annotation_data, function(x) {
          if (is.na(x)) "grey80" else color_map[as.character(x)]
        })

        feature_colors <- color_map
      }
    }

    # Column side colors are the sample Group - always shown in the legend,
    # same combining pattern as the sPLS-DA heatmap below.
    group_colors <- mixOmics::color.mixo(1:length(levels(groups_filtered())))
    names(group_colors) <- levels(groups_filtered())

    legend <- if (!is.null(row_sideColors) && !is.null(feature_colors)) {
      list(
        legend = c(names(group_colors), "", names(feature_colors)),
        col = c(group_colors, "white", feature_colors),
        title = paste("Groups /", input$heatmap_annotation_col),
        cex = 0.7
      )
    } else {
      list(
        legend = names(group_colors),
        col = group_colors,
        title = "Groups",
        cex = 0.8
      )
    }

    dist_method <- input$univ_heatmap_dist
    clust_method <- input$univ_heatmap_clust

    cim_args <- list(
      mat = mat_t,
      cluster = "both",
      dist.method = c(dist_method, dist_method),
      clust.method = c(clust_method, clust_method),
      col.sideColors = mixOmics::color.mixo(groups_filtered()),
      row.sideColors = row_sideColors,
      row.names = TRUE,
      col.names = FALSE,
      margins = if (!is.null(row_sideColors)) c(5, 12) else c(5, 10),
      legend = legend
    )

    # Run cim() into a null device just to read back its clustering order -
    # nothing is drawn to the user's plot from this call.
    grDevices::pdf(file = NULL)
    res <- tryCatch(do.call(mixOmics::cim, cim_args), error = function(e) NULL,
                     finally = grDevices::dev.off())
    validate(need(!is.null(res), "Could not compute heatmap clustering with the selected distance/clustering method"))

    groups_ordered <- groups_filtered()[res$colInd]
    sample_names_ordered <- if (!is.null(rownames(X_df))) {
      rownames(X_df)[res$colInd]
    } else {
      paste0("Sample_", seq_len(nrow(X_df)))[res$colInd]
    }

    list(
      cim_args = cim_args,
      mat_ordered = res$mat,
      top_vars_ordered = top_vars[res$rowInd],
      groups_ordered = groups_ordered,
      sample_names_ordered = sample_names_ordered
    )
  })

  output$univariate_heatmap <- renderPlot({
    hd <- univ_heatmap_ordered()
    validate(need(!is.null(hd), "At least 2 significant variables required for heatmap"))
    do.call(mixOmics::cim, hd$cim_args)
  })

  # Enhanced volcano plot with size based on metadata - UPDATED to use X_feature_filtered
  volcano_data <- reactive({
    req(X_filtered(), groups_filtered(), input$volcano_group1, input$volcano_group2)
    
    if (input$volcano_group1 == input$volcano_group2) {
      return(NULL)
    }
    
    # Check if selected groups are in filtered groups
    if(!(input$volcano_group1 %in% levels(groups_filtered())) || 
       !(input$volcano_group2 %in% levels(groups_filtered()))) {
      showNotification("Selected volcano groups not in current group selection", type = "warning")
      return(NULL)
    }
    
    withProgress(message = "Volcano plot...", value = 0, {
      X <- X_filtered()
      X_unscaled <- tryCatch(X_unscaled_filtered(), error = function(e) NULL)
      groups <- groups_filtered()
      
      idx1 <- which(groups == input$volcano_group1)
      idx2 <- which(groups == input$volcano_group2)
      
      if (length(idx1) < 2 || length(idx2) < 2) {
        return(NULL)
      }
      
      results <- map_df(names(X), function(var_name) {
        incProgress(1/ncol(X))
        
        g1_vals <- X[[var_name]][idx1]
        g2_vals <- X[[var_name]][idx2]
        
        if (sum(!is.na(g1_vals)) < 2 || sum(!is.na(g2_vals)) < 2) {
          return(tibble(Variable = var_name, log2FC = NA_real_, pvalue = NA_real_))
        }
        
        test_res <- tryCatch({
          if (input$use_nonparametric) {
            wilcox.test(g2_vals, g1_vals)
          } else {
            t.test(g2_vals, g1_vals)
          }
        }, error = function(e) NULL)
        
        if (is.null(test_res)) {
          return(tibble(Variable = var_name, log2FC = NA_real_, pvalue = NA_real_))
        }
        
        # Calculate true Log2 Fold-Change on unscaled positive abundance
        if (!is.null(X_unscaled) && var_name %in% colnames(X_unscaled)) {
          u1 <- X_unscaled[[var_name]][idx1]
          u2 <- X_unscaled[[var_name]][idx2]
          m1 <- mean(u1, na.rm = TRUE)
          m2 <- mean(u2, na.rm = TRUE)
          non_zeros <- c(u1[u1 > 0], u2[u2 > 0])
          min_val <- if (length(non_zeros) > 0) min(non_zeros, na.rm = TRUE) else NA
          eps <- if (is.finite(min_val) && min_val > 0) min_val / 10 else 1e-6
          log2fc <- log2((max(0, m2) + eps) / (max(0, m1) + eps))
        } else {
          log2fc <- mean(g2_vals, na.rm = TRUE) - mean(g1_vals, na.rm = TRUE)
        }
        
        tibble(Variable = var_name, log2FC = log2fc, pvalue = test_res$p.value)
      })
      
      results <- results %>%
        filter(!is.na(pvalue)) %>%
        mutate(
          adj_pvalue = p.adjust(pvalue, method = input$adj_method),
          negLog10P = -log10(adj_pvalue),
          negLog10P = ifelse(is.infinite(negLog10P), 
                             max(negLog10P[is.finite(negLog10P)], 10, na.rm = TRUE) * 1.2,
                             negLog10P),
          Significance = case_when(
            abs(log2FC) >= input$volcano_fc & negLog10P >= input$volcano_pval ~ 
              ifelse(log2FC > 0, "Up", "Down"),
            TRUE ~ "NS"
          )
        )

      # Features sent from sPLS-DA / sPLS-R: restrict the view (the p-values were
      # adjusted over all tested features above)
      if (!is.null(rv$uv_features)) {
        results <- results[results$Variable %in% rv$uv_features, , drop = FALSE]
      }
      
      # Add metadata for sizing if available
      if(!is.null(matched_metadata()) && !is.null(input$volcano_size_col) && 
         input$volcano_size_col != "") {
        metadata_df <- matched_metadata()
        size_data <- metadata_df[results$Variable, input$volcano_size_col]
        results$size_value <- as.numeric(size_data)
        
        # Normalize sizes to the specified range
        if(!all(is.na(results$size_value))) {
          min_size <- input$volcano_size_range[1]
          max_size <- input$volcano_size_range[2]
          results$point_size <- scales::rescale(results$size_value, 
                                                to = c(min_size, max_size),
                                                from = range(results$size_value, na.rm = TRUE))
          results$point_size[is.na(results$point_size)] <- min_size
        } else {
          results$point_size <- 2
        }
      } else {
        results$point_size <- 2
      }
      
      results
    })
  })
  
  volcano_gg <- function() {
    data <- volcano_data()
    
    p <- ggplot(data, aes(x = log2FC, y = negLog10P, color = Significance)) +
      geom_vline(xintercept = c(-input$volcano_fc, input$volcano_fc), linetype = "dashed") +
      geom_hline(yintercept = input$volcano_pval, linetype = "dashed") +
      scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
      labs(title = paste(input$volcano_group2, "vs", input$volcano_group1),
           x = "Log2 FC", y = "-Log10(p-adj)") +
      theme_minimal()
    
    # Add points with variable sizes if metadata is available
    if("point_size" %in% names(data) && length(unique(data$point_size)) > 1) {
      p <- p + geom_point(aes(size = point_size), alpha = 0.6) +
        scale_size_continuous(range = input$volcano_size_range, 
                              name = input$volcano_size_col,
                              guide = "legend")
    } else {
      p <- p + geom_point(alpha = 0.6, size = 2)
    }
    
    if (input$volcano_label && input$volcano_n_labels > 0) {
      sig_data <- data[data$Significance != "NS", ]
      sig_data <- sig_data[order(sig_data$adj_pvalue), ]
      sig_data <- sig_data[seq_len(min(input$volcano_n_labels, nrow(sig_data))), ]
      if (nrow(sig_data) > 0) {
        sig_data$Label <- feature_labels(sig_data$Variable)
        p <- p + geom_text_repel(data = sig_data, aes(label = Label), size = 3, max.overlaps = 30)
      }
    }
    
    p
  }

  output$volcano_plot <- renderPlot({
    req(volcano_data())
    volcano_gg()
  })
  
  significant_volcano_results <- reactive({
    req(volcano_data())
    volcano_data()[volcano_data()$Significance != "NS", c("Variable", "log2FC", "pvalue", "adj_pvalue", "Significance")]
  })
  
  output$volcano_table <- DT::renderDataTable({
    req(significant_volcano_results())
    
    results_df <- add_label_column(significant_volcano_results())
    
    # Add metadata if available
    if(!is.null(matched_metadata()) && !is.null(input$volcano_size_col) && 
       input$volcano_size_col != "") {
      metadata_df <- matched_metadata()
      size_info <- metadata_df[results_df$Variable, input$volcano_size_col, drop = FALSE]
      results_df <- cbind(results_df, size_info)
    }
    
    DT::formatSignif(
      DT::datatable(as.data.frame(results_df), selection = 'single', rownames = FALSE,
                    options = list(pageLength = 10, scrollX = TRUE)),
      columns = intersect(c("log2FC", "pvalue", "adj_pvalue"), names(results_df)), digits = 3)
  })
  
  selected_volcano_var <- reactive({
    req(input$volcano_table_rows_selected)
    as.character(significant_volcano_results()[input$volcano_table_rows_selected, "Variable"])
  })
  
  volcano_boxplot_gg <- function() {
    var_name <- selected_volcano_var()
    X_unscaled <- tryCatch(X_unscaled_filtered(), error = function(e) NULL)
    has_unscaled <- !is.null(X_unscaled) && var_name %in% colnames(X_unscaled)
    X_df <- if (has_unscaled) X_unscaled else X_filtered()
    plot_data <- data.frame(Group = groups_filtered(), Value = X_df[[var_name]])
    plot_data <- plot_data[plot_data$Group %in% c(input$volcano_group1, input$volcano_group2), ]
    plot_data <- plot_data[complete.cases(plot_data), ]
    plot_data$Group <- factor(plot_data$Group)
    
    y_lab <- if (has_unscaled) "Abundance (unscaled)" else "Value"
    
    ggboxplot(plot_data, x = "Group", y = "Value", fill = "Group",
              palette = c("blue", "red"), title = feature_title(var_name),
              ylab = y_lab) +
      stat_compare_means(method = if(input$use_nonparametric) "wilcox.test" else "t.test") +
      theme_minimal() + theme(legend.position = "none")
  }

  output$volcano_boxplot <- renderPlot({
    req(selected_volcano_var(), X_filtered(), groups_filtered())
    volcano_boxplot_gg()
  })
  
  # Boxplot / violin of the selected feature
  selected_univ_var <- reactive({
    req(input$univariate_table_rows_selected)
    as.character(filtered_univariate_results()[input$univariate_table_rows_selected, "Variable"])
  })

  univ_boxplot_gg <- function() {
    var_name <- selected_univ_var()
    X_df <- X_filtered()
    plot_data <- data.frame(
      Group = groups_filtered(),
      Value = X_df[[var_name]]
    )
    plot_data <- plot_data[complete.cases(plot_data), ]
    validate(need(nrow(plot_data) > 0, "No data available for selected variable"))
    
    n_groups <- length(levels(plot_data$Group))
    
    if (input$plot_type == "violin") {
      p <- ggplot(plot_data, aes(x = Group, y = Value, fill = Group)) +
        geom_violin(trim = FALSE, alpha = 0.7) +
        scale_fill_manual(values = mixOmics::color.mixo(1:n_groups)) +
        labs(title = paste("Violin plot:", feature_title(var_name)), x = "", y = "Value") +
        theme_minimal(base_size = 14) +
        theme(legend.position = "none")
    } else {
      p <- ggboxplot(plot_data, x = "Group", y = "Value", fill = "Group",
                     palette = mixOmics::color.mixo(1:n_groups),
                     title = paste("Boxplot:", feature_title(var_name)),
                     xlab = "", ylab = "Value")
    }
    if (input$add_points) {
      p <- p + geom_jitter(width = 0.2, alpha = 0.6, size = 2)
    }
    
    letters_added <- FALSE
    if (n_groups > 2 && !input$use_nonparametric) {
      # Tukey HSD letters. Groups are renamed G1..Gk for the test: multcompLetters
      # splits the pair names "A-B" on "-", which breaks group names containing "-".
      letters_df <- tryCatch({
        safe <- plot_data
        safe$Group <- factor(safe$Group, levels = levels(plot_data$Group),
                             labels = paste0("G", seq_len(n_groups)))
        tukey_res <- TukeyHSD(aov(Value ~ Group, data = safe))
        let <- multcompLetters(tukey_res$Group[, "p adj"])$Letters
        max_vals <- aggregate(Value ~ Group, data = plot_data, FUN = max)
        data.frame(
          Group = max_vals$Group,
          Letters = unname(let[paste0("G", match(max_vals$Group, levels(plot_data$Group)))]),
          MaxValue = max_vals$Value,
          stringsAsFactors = FALSE
        )
      }, error = function(e) NULL)
      
      if (!is.null(letters_df)) {
        y_span <- diff(range(plot_data$Value, na.rm = TRUE))
        if (!is.finite(y_span) || y_span <= 0) y_span <- 1
        letters_df$YPos <- letters_df$MaxValue + 0.08 * y_span
        letters_df$Group <- factor(letters_df$Group, levels = levels(plot_data$Group))
        p <- p + geom_text(
          data = letters_df,
          aes(x = Group, y = YPos, label = Letters),
          size = 6,
          inherit.aes = FALSE
        )
        letters_added <- TRUE
      }
    }
    if (!letters_added) {
      test_method <- if(input$use_nonparametric) "kruskal.test" else "anova"
      p <- p + stat_compare_means(method = test_method)
    }
    
    p + theme_minimal(base_size = 14) + theme(legend.position = "none")
  }
  
  output$univariate_boxplot <- renderPlot({
    req(selected_univ_var(), X_filtered(), groups_filtered())
    univ_boxplot_gg()
  })
  
  # CSV export: exact same order as the plot
  univariate_heatmap_data <- reactive({
    hd <- univ_heatmap_ordered()
    req(hd)

    # hd$mat_ordered is features (rows) x samples (cols), in clustered/plotted order.
    # Transpose back to samples x features for a human-friendly CSV, keeping order.
    df <- as.data.frame(t(hd$mat_ordered))
    colnames(df) <- hd$top_vars_ordered
    df <- cbind(Sample = hd$sample_names_ordered, Group = hd$groups_ordered, df)
    rownames(df) <- NULL
    df
  })
  
  # PCA and other multivariate analyses use X_feature_filtered (selected groups,
  # scaled on those samples, optional feature filter)
  pca_model <- reactive({
    req(X_feature_filtered(), ncomp_eff())
    tryCatch({
      # center = TRUE: without centering, PC1 mostly captures the mean of the
      # data whenever the scaling method is "None" (or groups were excluded)
      mixOmics::pca(X_feature_filtered(), ncomp = ncomp_eff(), center = TRUE, scale = FALSE)
    }, error = function(e) {
      showNotification(paste("PCA error:", e$message), type = "error")
      NULL
    })
  })
  
  # Component selectors must refer to computed components
  check_comps <- function(...) {
    comps <- c(...)
    validate(need(all(!is.na(comps)) && all(comps >= 1) && all(comps <= ncomp_eff()),
                  paste0("Component numbers must be between 1 and ", ncomp_eff(),
                         " (increase 'Number of components' in the sidebar for more).")))
  }

  output$pca_scores <- renderPlot({
    req(pca_model(), groups_filtered())
    check_comps(input$pca_x, input$pca_y)
    mixOmics::plotIndiv(pca_model(), group = groups_filtered(), 
                        comp = c(input$pca_x, input$pca_y),
                        ind.names = FALSE, legend = TRUE, ellipse = TRUE, 
                        star = TRUE, abline = TRUE,
                        title = "PCA Scores Plot")
  })
  
  output$pca_scree <- renderPlot({
    req(pca_model())
    plot(pca_model())
  })
  
  output$pca_biplot <- renderPlot({
    req(pca_model(), groups_filtered())
    check_comps(input$pca_x, input$pca_y)
    biplot(pca_model(), group = groups_filtered(), comp = c(input$pca_x, input$pca_y))
  })
  
  selected_pca_scores <- reactive({
    req(pca_model())
    check_comps(input$pca_x, input$pca_y)
    comps <- unique(c(input$pca_x, input$pca_y))
    pca_model()$x[, comps, drop = FALSE]
  })
  
  output$pca_permanova <- renderPrint({
    req(selected_pca_scores(), groups_filtered())
    df <- data.frame(group = groups_filtered())
    apply_seed()
    result <- vegan::adonis2(selected_pca_scores() ~ group, data = df, permutations = 999, method = "euclidean")
    cat("Selected groups:", paste(levels(groups_filtered()), collapse = ", "), "\n\n")
    print(result)
  })
  
  pca_f_test_data <- reactive({
    req(pca_model(), groups_filtered())
    scores <- pca_model()$x
    f_vals <- sapply(1:ncol(scores), function(i) {
      summary(aov(scores[, i] ~ groups_filtered()))[[1]][["F value"]][1]
    })
    data.frame(
      Component = paste0("PC", 1:ncol(scores)),
      F_Value = f_vals
    )
  })
  
  output$pca_f_test <- renderPlot({
    req(pca_f_test_data())
    ggplot(pca_f_test_data(), aes(x = Component, y = F_Value)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_minimal(base_size = 14) +
      labs(title = "ANOVA F-statistic per Component", 
           x = "Component", y = "F-value") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$pca_silhouette <- renderPlot({
    req(selected_pca_scores(), groups_filtered())
    sil <- cluster::silhouette(as.numeric(groups_filtered()), dist(selected_pca_scores()))
    plot(sil, col = mixOmics::color.mixo(1:length(levels(groups_filtered()))), 
         border = NA,
         main = paste("Silhouette -", length(levels(groups_filtered())), "groups"))
  })
  
  # PLS-DA - UPDATED to use X_feature_filtered
  plsda_model <- reactive({
    req(X_feature_filtered(), groups_filtered(), ncomp_eff())
    tryCatch({
      mixOmics::plsda(X_feature_filtered(), groups_filtered(), ncomp = ncomp_eff(), scale = FALSE)
    }, error = function(e) {
      showNotification(paste("PLS-DA error:", e$message), type = "error")
      NULL
    })
  })
  
  output$plsda_scores <- renderPlot({
    req(plsda_model(), groups_filtered())
    check_comps(input$plsda_x, input$plsda_y)
    mixOmics::plotIndiv(plsda_model(), comp = c(input$plsda_x, input$plsda_y),
                        group = groups_filtered(), ind.names = FALSE, 
                        ellipse = TRUE, legend = TRUE, star = TRUE, 
                        abline = TRUE,
                        title = "PLS-DA Scores Plot")
  })
  
  observeEvent(input$run_plsda_perf, {
    req(plsda_model())
    
    withProgress(message = "Computing performance...", value = 0.3, {
      tryCatch({
        args <- list(plsda_model(), validation = input$plsda_cv_method,
                     progressBar = FALSE, auc = TRUE, seed = seed_value())
        if (input$plsda_cv_method == "Mfold") {
          # Never more folds than the smallest group can fill
          folds <- cv_folds(input$plsda_folds, groups_filtered())
          if (folds < input$plsda_folds) {
            showNotification(paste("Folds reduced to", folds, "(smallest group size)"), type = "message")
          }
          args$folds <- folds
          args$nrepeat <- input$plsda_nrepeat
        }
        apply_seed()
        perf <- do.call(mixOmics::perf, args)
        
        rv$plsda_perf <- perf
        showNotification("Performance calculated", type = "message")
        
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
      })
    })
  })
  
  output$plsda_perf_plot <- renderPlot({
    req(rv$plsda_perf)
    plot(rv$plsda_perf, col = mixOmics::color.mixo(5:7), sd = TRUE)
  })

  output$plsda_perf_summary <- renderPrint({
    req(rv$plsda_perf)
    perf_res <- rv$plsda_perf

    cat("Classification error rate - overall (mean over repeats):\n")
    print(round(perf_res$error.rate$overall, 4))

    cat("\nBalanced Error Rate (BER):\n")
    print(round(perf_res$error.rate$BER, 4))

    if (!is.null(perf_res$choice.ncomp)) {
      cat("\nOptimal number of components (lowest BER per distance):\n")
      print(perf_res$choice.ncomp)
    }

    if (!is.null(perf_res$auc)) {
      cat("\nAUC available per component - see the AUC plot below for per-class detail.\n")
    }
  })

  output$plsda_permanova <- renderPrint({
    req(plsda_model(), groups_filtered())
    check_comps(input$plsda_x, input$plsda_y)
    scores <- plsda_model()$variates$X[, c(input$plsda_x, input$plsda_y), drop = FALSE]
    df <- data.frame(group = groups_filtered())
    apply_seed()
    result <- vegan::adonis2(scores ~ group, data = df, permutations = 999, method = "euclidean")
    cat("Selected groups:", paste(levels(groups_filtered()), collapse = ", "), "\n\n")
    print(result)
  })
  
  output$plsda_auc <- renderPlot({
    req(plsda_model())
    
    validate(
      need(!is.na(input$plsda_auc_comp) && input$plsda_auc_comp <= ncomp_eff() && input$plsda_auc_comp >= 1, 
           "Invalid component number")
    )
    
    tryCatch({
      mixOmics::auroc(plsda_model(), roc.comp = input$plsda_auc_comp)
    }, error = function(e) {
      validate(need(FALSE, paste("AUC plot error:", e$message)))
    })
  })
  
  # sPLS-DA Tuning - UPDATED to use X_feature_filtered
  observeEvent(input$run_splsda_tuning, {
    req(X_feature_filtered(), groups_filtered())
    
    if(length(levels(groups_filtered())) < 2) {
      showNotification("At least 2 groups required for sPLS-DA", type = "error")
      return(NULL)
    }
    
    # Show info about feature filtering
    if(!is.null(rv$mv_features)) {
      showNotification(
        paste("Tuning sPLS-DA with", length(active_features()), 
              "pre-filtered significant features from ANOVA"),
        type = "message",
        duration = 3
      )
    }
    
    withProgress(message = "sPLS-DA tuning (may take several minutes)...", value = 0, {
      tryCatch({
        if (input$splsda_strategy == "sparse") {
          test_keepX <- c(5, 10, 15, 20, 25, 30, 40, 50)
        } else {
          test_keepX <- c(seq(20, 100, by = 10), seq(120, 200, by = 20))
        }
        
        # Limit keepX to available features
        max_features <- ncol(X_feature_filtered())
        test_keepX <- test_keepX[test_keepX < max_features]
        
        if(length(test_keepX) == 0) {
          test_keepX <- c(min(5, max_features - 1))
        }
        
        n_folds <- cv_folds(5, groups_filtered())
        
        # scale = FALSE: the data are already scaled with the method chosen in the
        # sidebar (mixOmics would otherwise apply UV scaling on top of it)
        apply_seed()
        tune_result <- mixOmics::tune.splsda(
          X_feature_filtered(),
          groups_filtered(),
          ncomp = min(5, ncomp_eff()),
          validation = 'Mfold',
          folds = n_folds,
          nrepeat = if (IS_WEBR) 5 else 10,
          dist = 'max.dist',
          measure = "BER",
          test.keepX = test_keepX,
          scale = FALSE,
          progressBar = FALSE,
          seed = seed_value()
        )
        
        rv$splsda_tune_result <- tune_result
        showNotification("Tuning completed", type = "message")
        
      }, error = function(e) {
        showNotification(paste("Tuning error:", e$message), type = "error")
      })
    })
  })
  
  output$splsda_tune_plot <- renderPlot({
    req(rv$splsda_tune_result)
    plot(rv$splsda_tune_result, col = mixOmics::color.jet(5))
  })
  
  output$splsda_tune_summary <- renderPrint({
    req(rv$splsda_tune_result)
    cat("Optimal number of components:", rv$splsda_tune_result$choice.ncomp$ncomp, "\n")
    cat("Optimal keepX:\n")
    print(rv$splsda_tune_result$choice.keepX)
  })
  
  final_splsda_model <- reactive({
    req(X_feature_filtered(), groups_filtered())
    
    if (input$splsda_mode == "auto") {
      req(rv$splsda_tune_result)
      ncomp <- rv$splsda_tune_result$choice.ncomp$ncomp
      keepX <- rv$splsda_tune_result$choice.keepX[1:ncomp]
    } else {
      ncomp <- input$splsda_ncomp_manual
      keepX <- as.numeric(trimws(unlist(strsplit(input$splsda_keepx_manual, ","))))
      validate(
        need(!any(is.na(keepX)), "Invalid keepX"),
        need(length(keepX) == ncomp, "keepX length mismatch")
      )
    }
    
    tryCatch({
      mixOmics::splsda(X_feature_filtered(), groups_filtered(), ncomp = ncomp, keepX = keepX, scale = FALSE)
    }, error = function(e) {
      showNotification(paste("sPLS-DA error:", e$message), type = "error")
      NULL
    })
  })
  
  output$splsda_scores <- renderPlot({
    req(final_splsda_model(), groups_filtered())
    mixOmics::plotIndiv(final_splsda_model(), group = groups_filtered(),
                        ind.names = FALSE, ellipse = TRUE, 
                        legend = TRUE, star = TRUE, abline = TRUE,
                        title = "sPLS-DA Scores Plot")
  })
  
  output$splsda_selected_features <- renderPrint({
    req(final_splsda_model())
    
    tryCatch({
      model <- final_splsda_model()
      
      cat("Features per component:\n")
      for (i in 1:model$ncomp) {
        selected_i <- mixOmics::selectVar(model, comp = i)
        n_features <- length(selected_i$name)
        cat("Component", i, ":", n_features, "features\n")
      }
      
      # Get all unique features across all components
      all_features <- character()
      for (i in 1:model$ncomp) {
        selected_i <- mixOmics::selectVar(model, comp = i)
        all_features <- c(all_features, selected_i$name)
      }
      all_features <- unique(all_features)
      
      cat("\nTotal unique features:", length(all_features), "\n")
      
      # Add metadata info if available
      if(!is.null(matched_metadata()) && !is.null(input$heatmap_annotation_col) && 
         input$heatmap_annotation_col != "") {
        metadata_df <- matched_metadata()
        feature_info <- metadata_df[all_features, input$heatmap_annotation_col]
        if(!all(is.na(feature_info))) {
          cat("\nFeature annotations:\n")
          annotation_summary <- table(feature_info, useNA = "ifany")
          print(annotation_summary)
        }
      }
      
    }, error = function(e) {
      cat("Error extracting features:", e$message, "\n")
    })
  })

  # Performance of the FINAL (tuned or manual) sPLS-DA model - distinct from
  # the tuning curve above, which only reports the BER used to pick keepX/ncomp.
  observeEvent(input$run_splsda_final_perf, {
    req(final_splsda_model())

    withProgress(message = "Evaluating final sPLS-DA model...", value = 0, {
      tryCatch({
        n_folds <- cv_folds(5, groups_filtered())

        apply_seed()
        perf_res <- mixOmics::perf(
          final_splsda_model(),
          validation = "Mfold",
          folds = n_folds,
          nrepeat = if (IS_WEBR) 10 else 50,
          progressBar = FALSE,
          auc = TRUE,
          seed = seed_value()
        )

        rv$splsda_final_perf <- perf_res
        showNotification("Final model performance calculated", type = "message")

      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
      })
    })
  })

  output$splsda_final_perf_plot <- renderPlot({
    req(rv$splsda_final_perf)
    plot(rv$splsda_final_perf, col = mixOmics::color.mixo(5:7), sd = TRUE)
  })

  output$splsda_final_perf_summary <- renderPrint({
    req(rv$splsda_final_perf)
    perf_res <- rv$splsda_final_perf

    cat("Classification error rate - overall (mean over repeats):\n")
    print(round(perf_res$error.rate$overall, 4))

    cat("\nBalanced Error Rate (BER):\n")
    print(round(perf_res$error.rate$BER, 4))
  })

  splsda_features_table <- reactive({
    req(final_splsda_model())
    
    tryCatch({
      model <- final_splsda_model()
      
      results_list <- list()
      for (i in 1:model$ncomp) {
        loadings <- mixOmics::plotLoadings(model, comp = i, method = "median", 
                                           contrib = "max", plot = FALSE)
        df <- as.data.frame(loadings$X)
        df$Component <- paste0("Comp", i)
        df$Variable <- rownames(df)
        rownames(df) <- NULL
        df <- df[, c("Variable", "Component", "GroupContrib", "importance")]
        results_list[[i]] <- df
      }
      
      results_df <- do.call(rbind, results_list)
      
      # Add metadata if available
      if(!is.null(matched_metadata())) {
        metadata_df <- matched_metadata()
        # Add selected metadata columns
        meta_cols <- names(metadata_df)[!names(metadata_df) %in% c("feature_name", input$feature_id_col)]
        if(length(meta_cols) > 0) {
          meta_subset <- metadata_df[results_df$Variable, meta_cols[1:min(3, length(meta_cols))], drop = FALSE]
          results_df <- cbind(results_df, meta_subset)
        }
      }
      
      return(results_df)
      
    }, error = function(e) {
      NULL
    })
  })
  
  output$download_splsda_features <- downloadHandler(
    filename = function() { "splsda_features_with_metadata.csv" },
    content = function(file) {
      req(splsda_features_table())
      write.csv(splsda_features_table(), file, row.names = FALSE)
    },
    contentType = "text/csv"
  )
  
  output$splsda_loadings1 <- renderPlot({
    req(final_splsda_model())
    model <- final_splsda_model()
    if (model$ncomp >= 1) {
      tryCatch({
        mixOmics::plotLoadings(model, comp = 1, method = "median", contrib = "max", 
                               size.name = 0.5, name.var = feature_labels(colnames(model$X)))
      }, error = function(e) {
        validate(need(FALSE, "Loadings plot error: try enlarging plot window or reduce number of features"))
      })
    }
  })
  
  output$splsda_loadings2 <- renderPlot({
    req(final_splsda_model())
    model <- final_splsda_model()
    if (model$ncomp >= 2) {
      tryCatch({
        mixOmics::plotLoadings(model, comp = 2, method = "median", contrib = "max",
                               size.name = 0.5, name.var = feature_labels(colnames(model$X)))
      }, error = function(e) {
        validate(need(FALSE, "Loadings plot error: try enlarging plot window or reduce number of features"))
      })
    }
  })
  
  # Enhanced sPLS-DA heatmap - single source of truth for the cim() call args and
  # the resulting clustered order, so the plot (drawn by cim() itself, with real
  # dendrograms) and the CSV/PDF downloads always agree. cim() is run once
  # headlessly here (to read back the rowInd/colInd it computed) and re-run with
  # the identical, deterministic arguments wherever it is drawn.
  splsda_heatmap_ordered <- reactive({
    req(final_splsda_model(), groups_filtered(), X_feature_filtered())

    model <- final_splsda_model()
    X_data <- X_feature_filtered()

    all_selected_features <- character()
    for (i in 1:model$ncomp) {
      selected_i <- mixOmics::selectVar(model, comp = i)
      all_selected_features <- c(all_selected_features, selected_i$name)
    }
    all_selected_features <- unique(all_selected_features)
    validate(need(length(all_selected_features) > 0, "No features selected in sPLS-DA model"))

    mat <- as.matrix(X_data[, all_selected_features, drop = FALSE])
    mat_scaled <- scale(mat)
    mat_scaled[!is.finite(mat_scaled)] <- 0
    mat_t <- t(mat_scaled)   # features x samples
    rownames(mat_t) <- feature_labels(all_selected_features)

    # Row (feature) side colors from metadata, in the SAME (pre-cluster) order
    # as mat_t's rows - cim() reorders these itself alongside the matrix rows.
    row_sideColors <- NULL
    feature_colors <- NULL

    if (!is.null(matched_metadata()) && !is.null(input$heatmap_annotation_col) &&
        input$heatmap_annotation_col != "") {

      metadata_df <- matched_metadata()
      annotation_data <- trimws(as.character(
        metadata_df[all_selected_features, input$heatmap_annotation_col]
      ))
      valid_annotations <- annotation_data[!is.na(annotation_data) & annotation_data != ""]

      if (length(valid_annotations) > 0) {
        unique_vals <- unique(valid_annotations)
        n_colors <- length(unique_vals)
        colors <- if (n_colors <= 12) {
          if (n_colors < 3) c("#8DD3C7", "#FFFFB3", "#BEBADA")[1:n_colors]
          else brewer.pal(max(3, n_colors), "Set3")[1:n_colors]
        } else {
          viridis(n_colors)
        }
        color_map <- setNames(colors, unique_vals)

        row_sideColors <- vapply(annotation_data, function(val) {
          if (is.na(val) || val == "" || !(val %in% names(color_map))) "grey80"
          else color_map[[val]]
        }, character(1))

        if (length(row_sideColors) == length(all_selected_features)) {
          feature_colors <- color_map
        } else {
          row_sideColors <- NULL
        }
      }
    }

    group_colors <- mixOmics::color.mixo(1:length(levels(groups_filtered())))
    names(group_colors) <- levels(groups_filtered())

    legend <- if (!is.null(row_sideColors) && !is.null(feature_colors)) {
      list(
        legend = c(names(group_colors), "", names(feature_colors)),
        col = c(group_colors, "white", feature_colors),
        title = paste("Groups /", input$heatmap_annotation_col),
        cex = 0.7
      )
    } else {
      list(
        legend = names(group_colors),
        col = group_colors,
        title = "Groups",
        cex = 0.8
      )
    }

    dist_method <- input$splsda_heatmap_dist
    clust_method <- input$splsda_heatmap_clust

    cim_args <- list(
      mat = mat_t,
      cluster = "both",
      dist.method = c(dist_method, dist_method),
      clust.method = c(clust_method, clust_method),
      col.sideColors = mixOmics::color.mixo(as.numeric(groups_filtered())),
      row.sideColors = row_sideColors,
      row.names = TRUE,
      col.names = FALSE,
      margins = if (!is.null(row_sideColors)) c(5, 16) else c(5, 12),
      legend = legend
    )

    # Run cim() into a null device just to read back its clustering order -
    # nothing is drawn to the user's plot from this call.
    grDevices::pdf(file = NULL)
    res <- tryCatch(do.call(mixOmics::cim, cim_args), error = function(e) NULL,
                     finally = grDevices::dev.off())
    validate(need(!is.null(res), "Could not compute heatmap clustering with the selected distance/clustering method"))

    groups_ordered <- groups_filtered()[res$colInd]
    sample_names_ordered <- if (!is.null(rownames(X_data))) {
      rownames(X_data)[res$colInd]
    } else {
      paste0("Sample_", seq_len(nrow(X_data)))[res$colInd]
    }

    list(
      cim_args = cim_args,
      mat_ordered = res$mat,
      features_ordered = all_selected_features[res$rowInd],
      groups_ordered = groups_ordered,
      sample_names_ordered = sample_names_ordered
    )
  })

  output$splsda_cim <- renderPlot({
    hd <- splsda_heatmap_ordered()
    req(hd)
    tryCatch({
      do.call(mixOmics::cim, hd$cim_args)
    }, error = function(e) {
      showNotification(paste("Heatmap error:", e$message), type = "error", duration = 5)
    })
  })

    # CSV export: exact same order as the plot
    splsda_cim_data <- reactive({
      hd <- splsda_heatmap_ordered()
      req(hd)

      df <- as.data.frame(t(hd$mat_ordered))     # samples x features, in plotted order
      colnames(df) <- hd$features_ordered
      df <- cbind(Sample = hd$sample_names_ordered, Group = hd$groups_ordered, df)
      rownames(df) <- NULL
      df
    })

    #==============================================================================
    # PLS REGRESSION (continuous response variable)
    #==============================================================================

    # Response variable, restricted to the currently selected samples and kept
    # as a named 1-column matrix so mixOmics plots (plotVar, plotArrow, plotIndiv
    # legends) show the real variable name instead of a generic "Y".
    Y_plsr <- reactive({
      req(rv$raw_data, input$plsr_y_col, rv$groups, rv$selected_groups)
      validate(need(input$plsr_y_col != "", "Select a response variable (Y) in the sidebar (section 7)"))
      validate(need(is.numeric(rv$raw_data[[input$plsr_y_col]]), "Response variable must be numeric"))

      selected_idx <- rv$groups %in% rv$selected_groups
      y_vals <- rv$raw_data[[input$plsr_y_col]][selected_idx]
      validate(need(!anyNA(y_vals), "Response variable contains missing values in the selected samples"))

      matrix(y_vals, ncol = 1, dimnames = list(NULL, make.names(input$plsr_y_col)))
    })

    # Predictor matrix: the shared feature-filtered matrix with the response
    # column excluded (it is otherwise just another numeric column).
    X_plsr <- reactive({
      req(X_feature_filtered(), input$plsr_y_col)
      validate(need(input$plsr_y_col != "", "Select a response variable (Y) in the sidebar (section 7)"))

      X <- X_feature_filtered()
      X <- X[, setdiff(colnames(X), input$plsr_y_col), drop = FALSE]
      validate(need(ncol(X) >= 2, "At least 2 predictor variables required after excluding the response variable"))
      X
    })

    plsr_model <- reactive({
      req(X_plsr(), Y_plsr())
      tryCatch({
        mixOmics::pls(X_plsr(), Y_plsr(), ncomp = min(ncomp_eff(), ncol(X_plsr())), mode = "regression", scale = FALSE)
      }, error = function(e) {
        showNotification(paste("PLS Regression error:", e$message), type = "error")
        NULL
      })
    })

    plsr_scores_gg <- reactive({
      req(plsr_model())
      model <- plsr_model()
      validate(need(input$plsr_x >= 1 && input$plsr_y_comp >= 1 &&
                    input$plsr_x <= model$ncomp && input$plsr_y_comp <= model$ncomp,
                    "Invalid component number"))

      scores <- as.data.frame(model$variates$X[, c(input$plsr_x, input$plsr_y_comp), drop = FALSE])
      colnames(scores) <- c("CompX", "CompY")
      scores$Y <- as.numeric(Y_plsr())

      grp <- tryCatch(groups_filtered(), error = function(e) NULL)
      p <- if (!is.null(grp) && length(grp) == nrow(scores)) {
        scores$Group <- grp
        ggplot(scores, aes(x = CompX, y = CompY, color = Y, shape = Group))
      } else {
        ggplot(scores, aes(x = CompX, y = CompY, color = Y))
      }

      p +
        geom_point(size = 3, alpha = 0.85) +
        scale_color_viridis_c(name = colnames(model$Y)[1]) +
        theme_minimal(base_size = 14) +
        labs(title = "PLS Regression Scores",
             x = paste0("Component ", input$plsr_x),
             y = paste0("Component ", input$plsr_y_comp))
    })

    output$plsr_scores <- renderPlot({
      plsr_scores_gg()
    })

    output$plsr_arrow <- renderPlot({
      req(plsr_model())
      grp <- tryCatch(groups_filtered(), error = function(e) NULL)
      tryCatch({
        mixOmics::plotArrow(plsr_model(), group = grp, legend = !is.null(grp),
                            title = "PLS-R Arrow Plot (X vs Y variates)")
      }, error = function(e) {
        validate(need(FALSE, paste("Arrow plot error:", e$message)))
      })
    })

    output$plsr_varplot <- renderPlot({
      req(plsr_model())
      model <- plsr_model()
      validate(need(input$plsr_x <= model$ncomp && input$plsr_y_comp <= model$ncomp,
                    "Invalid component number"))
      tryCatch({
        mixOmics::plotVar(model, comp = c(input$plsr_x, input$plsr_y_comp),
                          title = "Correlation Circle Plot")
      }, error = function(e) {
        validate(need(FALSE, paste("Correlation circle error:", e$message)))
      })
    })

    # In-sample fitted values from the final model - a quick-look diagnostic,
    # deliberately labelled as such since it is NOT cross-validated (see the
    # Performance panel below for the honest, CV-based Q2/MSEP numbers).
    plsr_fitted_gg <- reactive({
      req(plsr_model(), Y_plsr())
      model <- plsr_model()
      pred <- tryCatch(predict(model, newdata = X_plsr()), error = function(e) NULL)
      validate(need(!is.null(pred), "Prediction failed"))

      fitted_vals <- as.numeric(pred$predict[, 1, model$ncomp])
      obs <- as.numeric(Y_plsr())
      r2 <- suppressWarnings(cor(fitted_vals, obs)^2)

      df <- data.frame(Observed = obs, Fitted = fitted_vals)
      grp <- tryCatch(groups_filtered(), error = function(e) NULL)
      if (!is.null(grp) && length(grp) == nrow(df)) df$Group <- grp

      p <- if ("Group" %in% names(df)) {
        ggplot(df, aes(x = Observed, y = Fitted, color = Group))
      } else {
        ggplot(df, aes(x = Observed, y = Fitted))
      }

      p +
        geom_point(size = 3) +
        geom_smooth(method = "lm", se = FALSE, color = "grey40", linetype = "dashed") +
        theme_minimal(base_size = 14) +
        labs(title = "Fitted vs Observed (final model, in-sample)",
             subtitle = paste0("R2 = ", round(r2, 3),
                               " (in-sample fit; not cross-validated - see Performance panel for Q2/MSEP)"),
             x = paste("Observed", input$plsr_y_col), y = "Fitted")
    })

    output$plsr_fitted <- renderPlot({
      plsr_fitted_gg()
    })

    observeEvent(input$run_plsr_perf, {
      req(plsr_model())

      withProgress(message = "Computing PLS-R performance...", value = 0.3, {
        tryCatch({
          args <- list(plsr_model(), validation = input$plsr_cv_method, progressBar = FALSE,
                       seed = seed_value())
          if (input$plsr_cv_method == "Mfold") {
            args$folds <- cv_folds(input$plsr_folds, n = nrow(X_plsr()))
            args$nrepeat <- input$plsr_nrepeat
          }
          apply_seed()
          perf_res <- do.call(mixOmics::perf, args)

          rv$plsr_perf <- perf_res
          showNotification("Performance calculated", type = "message")

        }, error = function(e) {
          showNotification(paste("Error:", e$message), type = "error")
        })
      })
    })

    output$plsr_perf_plot <- renderPlot({
      req(rv$plsr_perf)
      tryCatch({
        plot(rv$plsr_perf, criterion = input$plsr_perf_criterion, sd = TRUE, col = mixOmics::color.mixo(1))
      }, error = function(e) {
        validate(need(FALSE, paste("Performance plot error:", e$message)))
      })
    })

    output$plsr_perf_summary <- renderPrint({
      req(rv$plsr_perf)
      m <- rv$plsr_perf$measures
      cat("Cross-validated performance (mean +/- sd over repeats):\n\n")
      for (measure_name in c("MSEP", "R2", "Q2", "Q2.total")) {
        tbl <- format_pls_perf_measure(m[[measure_name]])
        if (!is.null(tbl)) {
          cat(measure_name, ":\n", sep = "")
          print(tbl, row.names = FALSE)
          cat("\n")
        }
      }
    })

    output$plsr_loadings <- renderPlot({
      req(plsr_model())
      model <- plsr_model()
      validate(need(input$plsr_loadings_comp >= 1 && input$plsr_loadings_comp <= model$ncomp,
                    "Invalid component number"))
      tryCatch({
        mixOmics::plotLoadings(model, block = "X", comp = input$plsr_loadings_comp,
                                 name.var = feature_labels(colnames(model$X)),
                               method = "median", size.name = 0.5,
                               title = paste("Loadings - Component", input$plsr_loadings_comp))
      }, error = function(e) {
        validate(need(FALSE, "Loadings plot error: try enlarging plot window or reduce number of features"))
      })
    })

    # VIP table built directly from mixOmics::vip() (proper feature names in
    # rownames); plotLoadings()'s own plot=FALSE return value loses feature
    # names for regression (non-DA) objects, so it is not used here for data export.
    plsr_vip_table <- reactive({
      req(plsr_model())
      model <- plsr_model()
      vip_mat <- tryCatch(mixOmics::vip(model), error = function(e) NULL)
      validate(need(!is.null(vip_mat), "VIP computation failed"))

      comp <- min(max(1, input$plsr_loadings_comp), ncol(vip_mat))
      df <- data.frame(Variable = rownames(vip_mat), VIP = vip_mat[, comp], stringsAsFactors = FALSE)
      df <- df[order(-df$VIP), ]
      rownames(df) <- NULL
      df
    })

    plsr_vip_gg <- reactive({
      df <- plsr_vip_table()
      req(df)

      plot_df <- head(df, 30)
      plot_df$Variable <- feature_labels(plot_df$Variable)
      plot_df$Variable <- factor(plot_df$Variable, levels = rev(plot_df$Variable))

      ggplot(plot_df, aes(x = Variable, y = VIP)) +
        geom_col(fill = "steelblue") +
        geom_hline(yintercept = 1, linetype = "dashed", color = "firebrick") +
        coord_flip() +
        theme_minimal(base_size = 12) +
        labs(title = paste("VIP scores - component", input$plsr_loadings_comp),
             x = NULL, y = "VIP")
    })

    output$plsr_vip_plot <- renderPlot({
      plsr_vip_gg()
    })

    #==============================================================================
    # sPLS REGRESSION (sparse, feature-selecting)
    #==============================================================================

    observeEvent(input$run_splsr_tuning, {
      req(X_plsr(), Y_plsr())

      if (!is.null(rv$mv_features)) {
        showNotification(
          paste("Tuning sPLS-R with", ncol(X_plsr()),
                "pre-filtered significant features from ANOVA"),
          type = "message",
          duration = 3
        )
      }

      withProgress(message = "sPLS-R tuning (may take several minutes)...", value = 0, {
        tryCatch({
          if (input$splsr_strategy == "sparse") {
            test_keepX <- c(5, 10, 15, 20, 25, 30, 40, 50)
          } else {
            test_keepX <- c(seq(20, 100, by = 10), seq(120, 200, by = 20))
          }

          max_features <- ncol(X_plsr())
          test_keepX <- test_keepX[test_keepX < max_features]
          if (length(test_keepX) == 0) test_keepX <- c(min(5, max_features - 1))

          n_folds <- cv_folds(5, n = nrow(X_plsr()))

          apply_seed()
          tune_result <- mixOmics::tune.spls(
            X_plsr(), Y_plsr(),
            ncomp = min(5, ncomp_eff()),
            validation = 'Mfold',
            folds = n_folds,
            nrepeat = if (IS_WEBR) 5 else 10,
            measure = input$splsr_tune_measure,
            test.keepX = test_keepX,
            scale = FALSE,
            progressBar = FALSE,
            seed = seed_value()
          )

          rv$splsr_tune_result <- tune_result
          showNotification("Tuning completed", type = "message")

        }, error = function(e) {
          showNotification(paste("Tuning error:", e$message), type = "error")
        })
      })
    })

    output$splsr_tune_plot <- renderPlot({
      req(rv$splsr_tune_result)
      tryCatch({
        plot(rv$splsr_tune_result)
      }, error = function(e) {
        validate(need(FALSE, paste("Tuning plot error:", e$message)))
      })
    })

    output$splsr_tune_summary <- renderPrint({
      req(rv$splsr_tune_result)
      cat("Optimal number of components:", rv$splsr_tune_result$choice.ncomp$ncomp, "\n")
      cat("Optimal keepX:\n")
      print(rv$splsr_tune_result$choice.keepX)
    })

    final_splsr_model <- reactive({
      req(X_plsr(), Y_plsr())

      if (input$splsr_mode == "auto") {
        req(rv$splsr_tune_result)
        ncomp <- rv$splsr_tune_result$choice.ncomp$ncomp
        keepX <- rv$splsr_tune_result$choice.keepX[1:ncomp]
      } else {
        ncomp <- input$splsr_ncomp_manual
        keepX <- as.numeric(trimws(unlist(strsplit(input$splsr_keepx_manual, ","))))
        validate(
          need(!any(is.na(keepX)), "Invalid keepX"),
          need(length(keepX) == ncomp, "keepX length mismatch")
        )
      }

      tryCatch({
        mixOmics::spls(X_plsr(), Y_plsr(), ncomp = ncomp, keepX = keepX, mode = "regression", scale = FALSE)
      }, error = function(e) {
        showNotification(paste("sPLS-R error:", e$message), type = "error")
        NULL
      })
    })

    splsr_scores_gg <- reactive({
      req(final_splsr_model(), Y_plsr())
      model <- final_splsr_model()

      comp_x <- 1
      comp_y <- if (model$ncomp >= 2) 2 else 1

      scores <- as.data.frame(model$variates$X[, c(comp_x, comp_y), drop = FALSE])
      colnames(scores) <- c("CompX", "CompY")
      scores$Y <- as.numeric(Y_plsr())

      grp <- tryCatch(groups_filtered(), error = function(e) NULL)
      p <- if (!is.null(grp) && length(grp) == nrow(scores)) {
        scores$Group <- grp
        ggplot(scores, aes(x = CompX, y = CompY, color = Y, shape = Group))
      } else {
        ggplot(scores, aes(x = CompX, y = CompY, color = Y))
      }

      p +
        geom_point(size = 3, alpha = 0.85) +
        scale_color_viridis_c(name = colnames(model$Y)[1]) +
        theme_minimal(base_size = 14) +
        labs(title = "sPLS Regression Scores",
             x = paste0("Component ", comp_x), y = paste0("Component ", comp_y))
    })

    output$splsr_scores <- renderPlot({
      splsr_scores_gg()
    })

    output$splsr_selected_features <- renderPrint({
      req(final_splsr_model())

      tryCatch({
        model <- final_splsr_model()

        cat("Features per component:\n")
        for (i in 1:model$ncomp) {
          sel_i <- mixOmics::selectVar(model, comp = i)$X
          cat("Component", i, ":", length(sel_i$name), "features\n")
        }

        all_features <- character()
        for (i in 1:model$ncomp) {
          all_features <- c(all_features, mixOmics::selectVar(model, comp = i)$X$name)
        }
        all_features <- unique(all_features)

        cat("\nTotal unique features:", length(all_features), "\n")

        if (!is.null(matched_metadata()) && !is.null(input$heatmap_annotation_col) &&
            input$heatmap_annotation_col != "") {
          metadata_df <- matched_metadata()
          feature_info <- metadata_df[all_features, input$heatmap_annotation_col]
          if (!all(is.na(feature_info))) {
            cat("\nFeature annotations:\n")
            print(table(feature_info, useNA = "ifany"))
          }
        }

      }, error = function(e) {
        cat("Error extracting features:", e$message, "\n")
      })
    })

    # Performance of the FINAL (tuned or manual) sPLS-R model - distinct from
    # the tuning curve above, which only reports the criterion used to pick keepX/ncomp.
    observeEvent(input$run_splsr_final_perf, {
      req(final_splsr_model())

      withProgress(message = "Evaluating final sPLS-R model...", value = 0, {
        tryCatch({
          n_folds <- cv_folds(5, n = nrow(X_plsr()))

          apply_seed()
          perf_res <- mixOmics::perf(
            final_splsr_model(),
            validation = "Mfold",
            folds = n_folds,
            nrepeat = if (IS_WEBR) 10 else 50,
            progressBar = FALSE,
            seed = seed_value()
          )

          rv$splsr_final_perf <- perf_res
          showNotification("Final model performance calculated", type = "message")

        }, error = function(e) {
          showNotification(paste("Error:", e$message), type = "error")
        })
      })
    })

    output$splsr_final_perf_plot <- renderPlot({
      req(rv$splsr_final_perf)
      tryCatch({
        plot(rv$splsr_final_perf, criterion = "MSEP", sd = TRUE, col = mixOmics::color.mixo(1))
      }, error = function(e) {
        validate(need(FALSE, paste("Performance plot error:", e$message)))
      })
    })

    output$splsr_final_perf_summary <- renderPrint({
      req(rv$splsr_final_perf)
      m <- rv$splsr_final_perf$measures
      cat("Cross-validated performance of the final sPLS-R model (mean +/- sd):\n\n")
      for (measure_name in c("MSEP", "R2", "Q2", "Q2.total")) {
        tbl <- format_pls_perf_measure(m[[measure_name]])
        if (!is.null(tbl)) {
          cat(measure_name, ":\n", sep = "")
          print(tbl, row.names = FALSE)
          cat("\n")
        }
      }
    })

    splsr_features_table <- reactive({
      req(final_splsr_model())

      tryCatch({
        model <- final_splsr_model()
        vip_mat <- mixOmics::vip(model)

        results_list <- list()
        for (i in 1:model$ncomp) {
          sel_i <- mixOmics::selectVar(model, comp = i)$X
          df <- data.frame(
            Variable = sel_i$name,
            Component = paste0("Comp", i),
            Loading = sel_i$value$value.var,
            VIP = vip_mat[sel_i$name, i],
            stringsAsFactors = FALSE
          )
          results_list[[i]] <- df
        }
        results_df <- do.call(rbind, results_list)

        if (!is.null(matched_metadata())) {
          metadata_df <- matched_metadata()
          meta_cols <- names(metadata_df)[!names(metadata_df) %in% c("feature_name", input$feature_id_col)]
          if (length(meta_cols) > 0) {
            meta_subset <- metadata_df[results_df$Variable, meta_cols[1:min(3, length(meta_cols))], drop = FALSE]
            results_df <- cbind(results_df, meta_subset)
          }
        }

        results_df

      }, error = function(e) {
        NULL
      })
    })

    output$splsr_loadings1 <- renderPlot({
      req(final_splsr_model())
      model <- final_splsr_model()
      if (model$ncomp >= 1) {
        tryCatch({
          mixOmics::plotLoadings(model, block = "X", comp = 1, method = "median",
                                 size.name = 0.5, name.var = feature_labels(colnames(model$X)))
        }, error = function(e) {
          validate(need(FALSE, "Loadings plot error: try enlarging plot window or reduce number of features"))
        })
      }
    })

    output$splsr_loadings2 <- renderPlot({
      req(final_splsr_model())
      model <- final_splsr_model()
      if (model$ncomp >= 2) {
        tryCatch({
          mixOmics::plotLoadings(model, block = "X", comp = 2, method = "median",
                                 size.name = 0.5, name.var = feature_labels(colnames(model$X)))
        }, error = function(e) {
          validate(need(FALSE, "Loadings plot error: try enlarging plot window or reduce number of features"))
        })
      }
    })

    # Heatmap of sPLS-R-selected features, columns colored by a continuous
    # gradient of the response variable (no discrete groups here, unlike the
    # sPLS-DA heatmap) - single source of truth shared by the plot and the
    # CSV/PDF downloads, same pattern as splsda_heatmap_ordered above.
    splsr_heatmap_ordered <- reactive({
      req(final_splsr_model(), Y_plsr(), X_plsr())

      model <- final_splsr_model()
      X_data <- X_plsr()
      y_vals <- as.numeric(Y_plsr())

      all_selected_features <- character()
      for (i in 1:model$ncomp) {
        all_selected_features <- c(all_selected_features, mixOmics::selectVar(model, comp = i)$X$name)
      }
      all_selected_features <- unique(all_selected_features)
      validate(need(length(all_selected_features) > 0, "No features selected in sPLS-R model"))

      mat <- as.matrix(X_data[, all_selected_features, drop = FALSE])
      mat_scaled <- scale(mat)
      mat_scaled[!is.finite(mat_scaled)] <- 0
      mat_t <- t(mat_scaled)
      rownames(mat_t) <- feature_labels(all_selected_features)

      row_sideColors <- NULL
      feature_colors <- NULL

      if (!is.null(matched_metadata()) && !is.null(input$heatmap_annotation_col) &&
          input$heatmap_annotation_col != "") {

        metadata_df <- matched_metadata()
        annotation_data <- trimws(as.character(
          metadata_df[all_selected_features, input$heatmap_annotation_col]
        ))
        valid_annotations <- annotation_data[!is.na(annotation_data) & annotation_data != ""]

        if (length(valid_annotations) > 0) {
          unique_vals <- unique(valid_annotations)
          n_colors <- length(unique_vals)
          colors <- if (n_colors <= 12) {
            if (n_colors < 3) c("#8DD3C7", "#FFFFB3", "#BEBADA")[1:n_colors]
            else brewer.pal(max(3, n_colors), "Set3")[1:n_colors]
          } else {
            viridis(n_colors)
          }
          color_map <- setNames(colors, unique_vals)

          row_sideColors <- vapply(annotation_data, function(val) {
            if (is.na(val) || val == "" || !(val %in% names(color_map))) "grey80"
            else color_map[[val]]
          }, character(1))

          if (length(row_sideColors) == length(all_selected_features)) {
            feature_colors <- color_map
          } else {
            row_sideColors <- NULL
          }
        }
      }

      y_ramp <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
      y_range <- range(y_vals, na.rm = TRUE)
      y_has_variance <- diff(y_range) > 0
      y_bins <- if (y_has_variance) {
        cut(y_vals, breaks = 100, labels = FALSE, include.lowest = TRUE)
      } else {
        rep(50L, length(y_vals))
      }
      col_sideColors <- y_ramp[y_bins]

      # Column side colors are a continuous gradient, so build a small
      # representative color key (5 values spanning the observed range) since
      # cim()'s legend only supports discrete swatches - always shown, same
      # "primary axis always in the legend" pattern as the other heatmaps.
      y_col_name <- if (!is.null(input$plsr_y_col) && input$plsr_y_col != "") input$plsr_y_col else "Y"
      y_legend_vals <- if (y_has_variance) seq(y_range[1], y_range[2], length.out = 5) else rep(y_range[1], 5)
      y_legend_idx <- if (y_has_variance) {
        pmin(100, pmax(1, round((y_legend_vals - y_range[1]) / diff(y_range) * 99) + 1))
      } else {
        rep(50L, 5)
      }
      y_legend_colors <- y_ramp[y_legend_idx]
      y_legend_labels <- formatC(y_legend_vals, format = "g", digits = 3)

      legend <- if (!is.null(row_sideColors) && !is.null(feature_colors)) {
        list(
          legend = c(y_legend_labels, "", names(feature_colors)),
          col = c(y_legend_colors, "white", feature_colors),
          title = paste(y_col_name, "/", input$heatmap_annotation_col),
          cex = 0.7
        )
      } else {
        list(
          legend = y_legend_labels,
          col = y_legend_colors,
          title = y_col_name,
          cex = 0.8
        )
      }

      dist_method <- input$splsr_heatmap_dist
      clust_method <- input$splsr_heatmap_clust

      cim_args <- list(
        mat = mat_t,
        cluster = "both",
        dist.method = c(dist_method, dist_method),
        clust.method = c(clust_method, clust_method),
        col.sideColors = col_sideColors,
        row.sideColors = row_sideColors,
        row.names = TRUE,
        col.names = FALSE,
        margins = if (!is.null(row_sideColors)) c(5, 16) else c(5, 12),
        legend = legend
      )

      grDevices::pdf(file = NULL)
      res <- tryCatch(do.call(mixOmics::cim, cim_args), error = function(e) NULL,
                       finally = grDevices::dev.off())
      validate(need(!is.null(res), "Could not compute heatmap clustering with the selected distance/clustering method"))

      y_ordered <- y_vals[res$colInd]
      sample_names_ordered <- if (!is.null(rownames(X_data))) {
        rownames(X_data)[res$colInd]
      } else {
        paste0("Sample_", seq_len(nrow(X_data)))[res$colInd]
      }

      list(
        cim_args = cim_args,
        mat_ordered = res$mat,
        features_ordered = all_selected_features[res$rowInd],
        y_ordered = y_ordered,
        sample_names_ordered = sample_names_ordered
      )
    })

    output$splsr_cim <- renderPlot({
      hd <- splsr_heatmap_ordered()
      req(hd)
      tryCatch({
        do.call(mixOmics::cim, hd$cim_args)
      }, error = function(e) {
        showNotification(paste("Heatmap error:", e$message), type = "error", duration = 5)
      })
    })

    # CSV export: exact same order as the plot
    splsr_cim_data <- reactive({
      hd <- splsr_heatmap_ordered()
      req(hd)

      col_name <- if (!is.null(input$plsr_y_col) && input$plsr_y_col != "") input$plsr_y_col else "Y"
      df <- as.data.frame(t(hd$mat_ordered))     # samples x features, in plotted order
      colnames(df) <- hd$features_ordered
      df <- cbind(Sample = hd$sample_names_ordered, df)
      df[[col_name]] <- hd$y_ordered
      df <- df[, c("Sample", col_name, setdiff(names(df), c("Sample", col_name)))]
      rownames(df) <- NULL
      df
    })

    #==============================================================================
    # DOWNLOAD HANDLERS
    #==============================================================================
    
    # Helper function for PDF downloads
    create_pdf_download <- function(plot_expr, default_name) {
      downloadHandler(
        filename = function() { default_name },
        content = function(file) {
          pdf(file, width = input$pdf_width, height = input$pdf_height)
          tryCatch({
            plot_expr
          }, finally = {
            dev.off()
          })
        },
        contentType = "application/pdf"
      )
    }
    
    # --- UNIVARIATE DOWNLOADS ---
    
    # Univariate results table
    output$download_univ_table <- downloadHandler(
      filename = function() { "univariate_results.csv" },
      content = function(file) {
        req(filtered_univariate_results())
        
        results_df <- filtered_univariate_results()
        
        # Add metadata if available
        if(!is.null(matched_metadata()) && !is.null(input$heatmap_annotation_col) && 
           input$heatmap_annotation_col != "") {
          metadata_df <- matched_metadata()
          metadata_info <- metadata_df[results_df$Variable, input$heatmap_annotation_col, drop = FALSE]
          results_df <- cbind(results_df, metadata_info)
        }
        
        write.csv(results_df, file, row.names = FALSE)
      },
      contentType = "text/csv"
    )
    
    # Filtered data (significant features only)
    output$download_univ_filtered_data <- downloadHandler(
      filename = function() { "filtered_data_significant.csv" },
      content = function(file) {
        req(filtered_univariate_results(), X_filtered(), groups_filtered())
        
        sig_vars <- filtered_univariate_results()$Variable
        X_df <- X_filtered()
        filtered_data <- X_df[, sig_vars, drop = FALSE]
        filtered_data <- cbind(Group = groups_filtered(), filtered_data)
        
        write.csv(filtered_data, file, row.names = TRUE)
      },
      contentType = "text/csv"
    )
    
    # Univariate heatmap PDF - uses the same frozen order as the plot/CSV
    output$download_univ_heatmap <- downloadHandler(
      filename = function() { "univariate_heatmap.pdf" },
      content = function(file) {
        hd <- univ_heatmap_ordered()

        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          if (is.null(hd)) {
            plot.new()
            text(0.5, 0.5, "At least 2 significant variables required for heatmap")
          } else {
            do.call(mixOmics::cim, hd$cim_args)
          }
        }, finally = {
          dev.off()
        })
      },
      contentType = "application/pdf"
    )

    # Univariate heatmap data CSV - exact same order as the plot
    output$download_univ_heatmap_data <- downloadHandler(
      filename = function() { "univariate_heatmap_data.csv" },
      content = function(file) {
        req(univariate_heatmap_data())
        write.csv(univariate_heatmap_data(), file, row.names = FALSE)
      },
      contentType = "text/csv"
    )

    # Univariate boxplot PDF (same figure as on screen)
    output$download_univ_boxplot <- downloadHandler(
      filename = function() { "univariate_boxplot.pdf" },
      content = function(file) {
        req(selected_univ_var(), X_filtered(), groups_filtered())
        p <- univ_boxplot_gg()
        
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch(print(p), finally = dev.off())
      },
      contentType = "application/pdf"
    )
    
    # --- VOLCANO PLOT DOWNLOADS ---
    
    output$download_volcano <- downloadHandler(
      filename = function() { "volcano_plot.pdf" },
      content = function(file) {
        req(volcano_data())
        p <- volcano_gg()
        
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch(print(p), finally = dev.off())
      },
      contentType = "application/pdf"
    )
    
    output$download_volcano_table <- downloadHandler(
      filename = function() { "volcano_significant_features.csv" },
      content = function(file) {
        req(significant_volcano_results())
        
        results_df <- significant_volcano_results()
        
        # Add metadata if available
        if(!is.null(matched_metadata()) && !is.null(input$volcano_size_col) && 
           input$volcano_size_col != "") {
          metadata_df <- matched_metadata()
          size_info <- metadata_df[results_df$Variable, input$volcano_size_col, drop = FALSE]
          results_df <- cbind(results_df, size_info)
        }
        
        write.csv(results_df, file, row.names = FALSE)
      },
      contentType = "text/csv"
    )
    
    output$download_volcano_boxplot <- downloadHandler(
      filename = function() { "volcano_boxplot.pdf" },
      content = function(file) {
        req(selected_volcano_var(), X_filtered(), groups_filtered())
        p <- volcano_boxplot_gg()
        
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch(print(p), finally = dev.off())
      },
      contentType = "application/pdf"
    )
    
    # --- PCA DOWNLOADS ---
    
    output$download_pca_scores <- downloadHandler(
      filename = function() { "pca_scores.pdf" },
      content = function(file) {
        req(pca_model(), groups_filtered())
        
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          mixOmics::plotIndiv(pca_model(), group = groups_filtered(), 
                              comp = c(input$pca_x, input$pca_y),
                              ind.names = FALSE, legend = TRUE, ellipse = TRUE, 
                              star = TRUE, abline = TRUE,
                              title = "PCA Scores Plot")
        }, finally = {
          dev.off()
        })
      },
      contentType = "application/pdf"
    )
    
    output$download_pca_scree <- downloadHandler(
      filename = function() { "pca_scree.pdf" },
      content = function(file) {
        req(pca_model())
        
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          plot(pca_model())
        }, finally = {
          dev.off()
        })
      },
      contentType = "application/pdf"
    )
    
    output$download_pca_biplot <- downloadHandler(
      filename = function() { "pca_biplot.pdf" },
      content = function(file) {
        req(pca_model(), groups_filtered())
        
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          biplot(pca_model(), group = groups_filtered(), comp = c(input$pca_x, input$pca_y))
        }, finally = {
          dev.off()
        })
      },
      contentType = "application/pdf"
    )
    
    output$download_pca_ftest <- downloadHandler(
      filename = function() { "pca_ftest.pdf" },
      content = function(file) {
        req(pca_f_test_data())
        
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          p <- ggplot(pca_f_test_data(), aes(x = Component, y = F_Value)) +
            geom_bar(stat = "identity", fill = "steelblue") +
            theme_minimal(base_size = 14) +
            labs(title = "ANOVA F-statistic per Component", 
                 x = "Component", y = "F-value") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
          print(p)
        }, finally = {
          dev.off()
        })
      },
      contentType = "application/pdf"
    )
    
    output$download_pca_silhouette <- downloadHandler(
      filename = function() { "pca_silhouette.pdf" },
      content = function(file) {
        req(selected_pca_scores(), groups_filtered())
        
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          sil <- cluster::silhouette(as.numeric(groups_filtered()), dist(selected_pca_scores()))
          plot(sil, col = mixOmics::color.mixo(1:length(levels(groups_filtered()))), 
               border = NA,
               main = paste("Silhouette -", length(levels(groups_filtered())), "groups"))
        }, finally = {
          dev.off()
        })
      },
      contentType = "application/pdf"
    )
    
    # --- PLS-DA DOWNLOADS ---
    
    output$download_plsda_scores <- downloadHandler(
      filename = function() { "plsda_scores.pdf" },
      content = function(file) {
        req(plsda_model(), groups_filtered())
        
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          mixOmics::plotIndiv(plsda_model(), comp = c(input$plsda_x, input$plsda_y),
                              group = groups_filtered(), ind.names = FALSE, 
                              ellipse = TRUE, legend = TRUE, star = TRUE, 
                              abline = TRUE,
                              title = "PLS-DA Scores Plot")
        }, finally = {
          dev.off()
        })
      },
      contentType = "application/pdf"
    )
    
    output$download_plsda_perf <- downloadHandler(
      filename = function() { "plsda_performance.pdf" },
      content = function(file) {
        req(rv$plsda_perf)
        
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          plot(rv$plsda_perf, col = mixOmics::color.mixo(5:7), sd = TRUE)
        }, finally = {
          dev.off()
        })
      },
      contentType = "application/pdf"
    )
    
    output$download_plsda_auc <- downloadHandler(
      filename = function() { "plsda_auc.pdf" },
      content = function(file) {
        req(plsda_model())
        
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          mixOmics::auroc(plsda_model(), roc.comp = input$plsda_auc_comp)
        }, finally = {
          dev.off()
        })
      },
      contentType = "application/pdf"
    )
    
    # --- sPLS-DA DOWNLOADS ---
    
    output$download_splsda_tune <- downloadHandler(
      filename = function() { "splsda_tuning.pdf" },
      content = function(file) {
        req(rv$splsda_tune_result)
        
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          plot(rv$splsda_tune_result, col = mixOmics::color.jet(5))
        }, finally = {
          dev.off()
        })
      },
      contentType = "application/pdf"
    )
    
    output$download_splsda_scores <- downloadHandler(
      filename = function() { "splsda_scores.pdf" },
      content = function(file) {
        req(final_splsda_model(), groups_filtered())
        
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          mixOmics::plotIndiv(final_splsda_model(), group = groups_filtered(),
                              ind.names = FALSE, ellipse = TRUE, 
                              legend = TRUE, star = TRUE, abline = TRUE,
                              title = "sPLS-DA Scores Plot")
        }, finally = {
          dev.off()
        })
      },
      contentType = "application/pdf"
    )

    output$download_splsda_final_perf <- downloadHandler(
      filename = function() { "splsda_final_model_performance.pdf" },
      content = function(file) {
        req(rv$splsda_final_perf)
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          plot(rv$splsda_final_perf, col = mixOmics::color.mixo(5:7), sd = TRUE)
        }, finally = { dev.off() })
      },
      contentType = "application/pdf"
    )

    output$download_splsda_load1 <- downloadHandler(
      filename = function() { "splsda_loadings_comp1.pdf" },
      content = function(file) {
        req(final_splsda_model())
        
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          model <- final_splsda_model()
          if (model$ncomp >= 1) {
            mixOmics::plotLoadings(model, comp = 1, method = "median", contrib = "max", 
                                   size.name = 0.5, name.var = feature_labels(colnames(model$X)))
          }
        }, finally = {
          dev.off()
        })
      },
      contentType = "application/pdf"
    )
    
    output$download_splsda_load2 <- downloadHandler(
      filename = function() { "splsda_loadings_comp2.pdf" },
      content = function(file) {
        req(final_splsda_model())
        
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          model <- final_splsda_model()
          if (model$ncomp >= 2) {
            mixOmics::plotLoadings(model, comp = 2, method = "median", contrib = "max",
                                   size.name = 0.5, name.var = feature_labels(colnames(model$X)))
          }
        }, finally = {
          dev.off()
        })
      },
      contentType = "application/pdf"
    )
    
    output$download_splsda_cim <- downloadHandler(
      filename = function() { "splsda_heatmap.pdf" },
      content = function(file) {
        hd <- splsda_heatmap_ordered()
        req(hd)

        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          do.call(mixOmics::cim, hd$cim_args)
        }, finally = {
          dev.off()
        })
      },
      contentType = "application/pdf"
    )

    output$download_splsda_cim_data <- downloadHandler(
      filename = function() { "splsda_heatmap_data.csv" },
      content = function(file) {
        hd <- splsda_heatmap_ordered()
        req(hd)

        data_export <- splsda_cim_data()

        if (!is.null(matched_metadata()) && !is.null(input$heatmap_annotation_col) &&
            input$heatmap_annotation_col != "") {

          metadata_df <- matched_metadata()
          annotation_info <- metadata_df[hd$features_ordered, input$heatmap_annotation_col, drop = FALSE]

          metadata_row <- data.frame(
            Sample = "Feature_Metadata",
            Group = input$heatmap_annotation_col,
            stringsAsFactors = FALSE
          )
          for (feat in hd$features_ordered) {
            metadata_row[[feat]] <- as.character(annotation_info[feat, 1])
          }

          data_export <- rbind(metadata_row, data_export)
        }

        write.csv(data_export, file, row.names = FALSE)
      },
      contentType = "text/csv"
    )

    # --- PLS REGRESSION DOWNLOADS ---

    output$download_plsr_scores <- downloadHandler(
      filename = function() { "plsr_scores.pdf" },
      content = function(file) {
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          print(plsr_scores_gg())
        }, finally = { dev.off() })
      },
      contentType = "application/pdf"
    )

    output$download_plsr_arrow <- downloadHandler(
      filename = function() { "plsr_arrow.pdf" },
      content = function(file) {
        req(plsr_model())
        grp <- tryCatch(groups_filtered(), error = function(e) NULL)
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          mixOmics::plotArrow(plsr_model(), group = grp, legend = !is.null(grp),
                              title = "PLS-R Arrow Plot (X vs Y variates)")
        }, finally = { dev.off() })
      },
      contentType = "application/pdf"
    )

    output$download_plsr_varplot <- downloadHandler(
      filename = function() { "plsr_correlation_circle.pdf" },
      content = function(file) {
        req(plsr_model())
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          mixOmics::plotVar(plsr_model(), comp = c(input$plsr_x, input$plsr_y_comp),
                            title = "Correlation Circle Plot")
        }, finally = { dev.off() })
      },
      contentType = "application/pdf"
    )

    output$download_plsr_fitted <- downloadHandler(
      filename = function() { "plsr_fitted_vs_observed.pdf" },
      content = function(file) {
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          print(plsr_fitted_gg())
        }, finally = { dev.off() })
      },
      contentType = "application/pdf"
    )

    output$download_plsr_perf <- downloadHandler(
      filename = function() { "plsr_performance.pdf" },
      content = function(file) {
        req(rv$plsr_perf)
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          plot(rv$plsr_perf, criterion = input$plsr_perf_criterion, sd = TRUE, col = mixOmics::color.mixo(1))
        }, finally = { dev.off() })
      },
      contentType = "application/pdf"
    )

    output$download_plsr_loadings <- downloadHandler(
      filename = function() { "plsr_loadings.pdf" },
      content = function(file) {
        req(plsr_model())
        model <- plsr_model()
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          mixOmics::plotLoadings(model, block = "X", comp = input$plsr_loadings_comp,
                                 name.var = feature_labels(colnames(model$X)),
                                 method = "median", size.name = 0.5,
                                 title = paste("Loadings - Component", input$plsr_loadings_comp))
        }, finally = { dev.off() })
      },
      contentType = "application/pdf"
    )

    output$download_plsr_vip <- downloadHandler(
      filename = function() { "plsr_vip.pdf" },
      content = function(file) {
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          print(plsr_vip_gg())
        }, finally = { dev.off() })
      },
      contentType = "application/pdf"
    )

    output$download_plsr_vip_csv <- downloadHandler(
      filename = function() { "plsr_vip_table.csv" },
      content = function(file) {
        req(plsr_vip_table())
        write.csv(plsr_vip_table(), file, row.names = FALSE)
      },
      contentType = "text/csv"
    )

    # --- sPLS REGRESSION DOWNLOADS ---

    output$download_splsr_tune <- downloadHandler(
      filename = function() { "splsr_tuning.pdf" },
      content = function(file) {
        req(rv$splsr_tune_result)
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          plot(rv$splsr_tune_result)
        }, finally = { dev.off() })
      },
      contentType = "application/pdf"
    )

    output$download_splsr_scores <- downloadHandler(
      filename = function() { "splsr_scores.pdf" },
      content = function(file) {
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          print(splsr_scores_gg())
        }, finally = { dev.off() })
      },
      contentType = "application/pdf"
    )

    output$download_splsr_features <- downloadHandler(
      filename = function() { "splsr_features_with_metadata.csv" },
      content = function(file) {
        req(splsr_features_table())
        write.csv(splsr_features_table(), file, row.names = FALSE)
      },
      contentType = "text/csv"
    )

    output$download_splsr_final_perf <- downloadHandler(
      filename = function() { "splsr_final_model_performance.pdf" },
      content = function(file) {
        req(rv$splsr_final_perf)
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          plot(rv$splsr_final_perf, criterion = "MSEP", sd = TRUE, col = mixOmics::color.mixo(1))
        }, finally = { dev.off() })
      },
      contentType = "application/pdf"
    )

    output$download_splsr_load1 <- downloadHandler(
      filename = function() { "splsr_loadings_comp1.pdf" },
      content = function(file) {
        req(final_splsr_model())
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          model <- final_splsr_model()
          if (model$ncomp >= 1) {
            mixOmics::plotLoadings(model, block = "X", comp = 1, method = "median",
                                   size.name = 0.5, name.var = feature_labels(colnames(model$X)))
          }
        }, finally = { dev.off() })
      },
      contentType = "application/pdf"
    )

    output$download_splsr_load2 <- downloadHandler(
      filename = function() { "splsr_loadings_comp2.pdf" },
      content = function(file) {
        req(final_splsr_model())
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          model <- final_splsr_model()
          if (model$ncomp >= 2) {
            mixOmics::plotLoadings(model, block = "X", comp = 2, method = "median",
                                   size.name = 0.5, name.var = feature_labels(colnames(model$X)))
          }
        }, finally = { dev.off() })
      },
      contentType = "application/pdf"
    )

    output$download_splsr_cim <- downloadHandler(
      filename = function() { "splsr_heatmap.pdf" },
      content = function(file) {
        hd <- splsr_heatmap_ordered()
        req(hd)
        pdf(file, width = input$pdf_width, height = input$pdf_height)
        tryCatch({
          do.call(mixOmics::cim, hd$cim_args)
        }, finally = { dev.off() })
      },
      contentType = "application/pdf"
    )

    output$download_splsr_cim_data <- downloadHandler(
      filename = function() { "splsr_heatmap_data.csv" },
      content = function(file) {
        hd <- splsr_heatmap_ordered()
        req(hd)

        data_export <- splsr_cim_data()

        if (!is.null(matched_metadata()) && !is.null(input$heatmap_annotation_col) &&
            input$heatmap_annotation_col != "") {

          metadata_df <- matched_metadata()
          annotation_info <- metadata_df[hd$features_ordered, input$heatmap_annotation_col, drop = FALSE]

          col_name <- if (!is.null(input$plsr_y_col) && input$plsr_y_col != "") input$plsr_y_col else "Y"
          metadata_row <- data.frame(
            Sample = "Feature_Metadata",
            stringsAsFactors = FALSE
          )
          metadata_row[[col_name]] <- input$heatmap_annotation_col
          for (feat in hd$features_ordered) {
            metadata_row[[feat]] <- as.character(annotation_info[feat, 1])
          }

          data_export <- rbind(metadata_row, data_export)
        }

        write.csv(data_export, file, row.names = FALSE)
      },
      contentType = "text/csv"
    )
    }

#==============================================================================
# LAUNCH APPLICATION
#==============================================================================

shinyApp(ui = ui, server = server)