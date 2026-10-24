## =============================================================================
## DiagProgMeta Shiny app
##
## Diagnostic test accuracy meta-analysis (SROC curve):
##   - reitsma.sroc(): bivariate random-effects model (Reitsma et al., 2005)
##   - srocsa.t():     likelihood-based sensitivity analysis for publication
##                      bias on the SROC curve / SAUC (Zhou et al., 2023)
##
## Launch with DiagProgMeta::run_app().
## =============================================================================

library(shiny)
library(DT)

## ---- helpers ---------------------------------------------------------------

## load a built-in dataset of the package without attaching the package
.get_builtin_data <- function(name) {
  e <- new.env()
  utils::data(list = name, package = "DiagProgMeta", envir = e)
  e[[name]]
}

## round only the numeric columns of a data frame (for display)
.round_df <- function(df, digits = 3) {
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], round, digits = digits)
  df
}

required_cols <- c("TP", "FN", "FP", "TN")

## =============================================================================
## UI
## =============================================================================

### fluidPage(): Bootstrapの12列グリッドに基づくレスポンシブなページのレイアウトを作成する関数
ui <- fluidPage(
  theme = bslib::bs_theme(version = 5, bootswatch = "flatly"), #Bootstrap 5 に Bootswatch の Flatly テーマを当てる指定
  titlePanel(
    "DiagProgMeta:Meta-analysis of diagnostic test accuracy (SROC curve)"
  ),
  p(
    "An interface for meta-analysis of diagnostic test accuracy",
    "(SROC curve and sensitivity analysis for publication bias) using ",
    tags$code("reitsma.sroc()"), " and ", tags$code("srocsa.t()"),
    " from the ", tags$code("DiagProgMeta"), " package."
  ),

  sidebarLayout(
    sidebarPanel(
      width = 4,

      h4("1. Data"),
      radioButtons(
        "data_source", "Data source",
        choices = c(
          "Built-in data" = "builtin",
          "Upload a CSV file" = "upload"
        ),
        selected = "builtin"
      ),

      conditionalPanel(
        condition = "input.data_source == 'builtin'",
        selectInput(
          "builtin_data", "Dataset",
          choices = c(
            "IVD (diagnosis of intravascular device-related infection, 33 studies)" = "IVD",
            "CD64 (neutrophil CD64 for diagnosing bacterial infection, 27 studies)" = "CD64"
          )
        )
      ),

      conditionalPanel(
        condition = "input.data_source == 'upload'",
        fileInput(
          "csv_file", "CSV file (with TP/FN/FP/TN columns)",
          accept = c(".csv", "text/csv")
        ),
        downloadLink("dl_sample_csv", "Download a sample CSV (IVD)"),
        uiOutput("col_mapping_ui")
      ),

      hr(),
      h4("2. Continuity correction"),
      numericInput(
        "cc_value", "cc.value(range:(0, 1])",
        value = 0.5, min = 0.001, max = 1, step = 0.1
      ),
      selectInput(
        "cc_type", "cc.type",
        choices = c(
          "single (only studies with zero cells)" = "single",
          "all (all studies)" = "all"
        ),
        selected = "single"
      ),

      hr(),
      h4("3. Estimation settings"),
      numericInput(
        "ci_level", "Confidence level (ci.level)",
        value = 0.95, min = 0.5, max = 0.999, step = 0.01
      ),
      checkboxInput(
        "run_sa",
        "Run sensitivity analysis for publication bias (srocsa.t)",
        value = FALSE
      ),
      conditionalPanel(
        condition = "input.run_sa == true",
        textInput(
          "sa_p", "sa.p (marginal selection probabilities, comma-separated)",
          value = "1, 0.8, 0.6"
        ),
        checkboxInput("estimate_c11", "Estimate c1^2 (recommended)", value = TRUE),
        conditionalPanel(
          condition = "input.estimate_c11 == false",
          numericInput(
            "c11_value", "Fixed value of c1^2 (0 to 1)",
            value = 0.5, min = 0, max = 1, step = 0.05
          )
        )
      ),

      hr(),
      actionButton("run", "Run analysis", class = "btn-primary", width = "100%"),
      br(), br(),
      helpText(
        "Reference: Reitsma et al. (2005) doi:10.1016/j.jclinepi.2005.02.022; ",
        "Zhou et al. (2023) doi:10.1002/sim.9643"
      )
    ),

    mainPanel(
      width = 8,
      tabsetPanel(
        id = "tabs",

        tabPanel(
          "Data preview",
          br(),
          textOutput("data_msg"),
          br(),
          DTOutput("data_table")
        ),

        tabPanel(
          "Results (no publication bias)",
          br(),
          h4("print()"),
          verbatimTextOutput("print_main"),
          h4("summary()"),
          verbatimTextOutput("summary_main"),
          hr(),
          fluidRow(
            column(3, checkboxInput("plot_ci_main", "Confidence intervals (plot.ci)", TRUE)),
            column(3, checkboxInput("spoint_region_main", "Summary point confidence regions", TRUE)),
            column(3, checkboxInput("add_study_main", "Show individual studies", TRUE)),
            column(
              3,
              selectInput(
                "study_size_main", "Size of study points",
                choices = c("fixed", "precision"), selected = "fixed"
              )
            )
          ),
          plotOutput("plot_main", height = "480px"),
          downloadButton("dl_plot_main", "Save plot as PNG")
        ),

        tabPanel(
          "Sensitivity analysis (publication bias)",
          br(),
          p(
            style = "color:#888;",
            "To see results here, check ‘Run sensitivity analysis for publication bias’ in the sidebar and click ‘Run analysis’."
          ),
          h4("print()"),
          verbatimTextOutput("print_sa"),
          h4("SAUC estimates"),
          DTOutput("sauc_table"),
          h4("Parameter estimates"),
          DTOutput("coef_table"),
          hr(),
          fluidRow(
            column(4, checkboxInput("add_study_sa", "Show individual studies", TRUE)),
            column(4, checkboxInput("plot_ci_sa", "Confidence intervals (plot.ci)", FALSE))
          ),
          plotOutput("plot_sa", height = "480px"),
          downloadButton("dl_plot_sa", "Save plot as PNG"),
          downloadButton("dl_sauc_csv", "Save SAUC table as CSV")
        )
      )
    )
  )
)

## =============================================================================
## Server
## =============================================================================

server <- function(input, output, session) {

  ## ---- 1. data -------------------------------------------------------------

  uploaded_df <- reactive({
    req(input$csv_file)
    df <- tryCatch(
      utils::read.csv(input$csv_file$datapath, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    validate(need(!is.null(df) && ncol(df) >= 1, "Could not read the CSV file."))
    df
  })

  output$col_mapping_ui <- renderUI({
    req(uploaded_df())
    cols <- names(uploaded_df())
    guess <- function(pat) {
      hit <- grep(pat, cols, ignore.case = TRUE, value = TRUE)
      if (length(hit)) hit[1] else cols[1]
    }
    tagList(
      hr(),
      strong("Column mapping"),
      selectInput("col_tp", "TP (True Positive) column", choices = cols, selected = guess("^tp$")),
      selectInput("col_fn", "FN (False Negative) column", choices = cols, selected = guess("^fn$")),
      selectInput("col_fp", "FP (False Positive) column", choices = cols, selected = guess("^fp$")),
      selectInput("col_tn", "TN (True Negative) column", choices = cols, selected = guess("^tn$"))
    )
  })

  output$dl_sample_csv <- downloadHandler(
    filename = function() "sample_dta_data.csv",
    content = function(file) {
      d <- .get_builtin_data("IVD")
      utils::write.csv(d[, c("study", required_cols)], file, row.names = FALSE)
    }
  )

  ## standardized data frame with exactly the columns TP, FN, FP, TN
  study_data <- reactive({
    if (input$data_source == "builtin") {
      d <- .get_builtin_data(input$builtin_data)
      validate(need(all(required_cols %in% names(d)), "Could not find the required columns (TP/FN/FP/TN) in the built-in data."))
      d[, required_cols, drop = FALSE]
    } else {
      req(uploaded_df(), input$col_tp, input$col_fn, input$col_fp, input$col_tn)
      df <- uploaded_df()
      out <- data.frame(
        TP = suppressWarnings(as.numeric(df[[input$col_tp]])),
        FN = suppressWarnings(as.numeric(df[[input$col_fn]])),
        FP = suppressWarnings(as.numeric(df[[input$col_fp]])),
        TN = suppressWarnings(as.numeric(df[[input$col_tn]]))
      )
      out
    }
  })

  output$data_table <- renderDT({
    DT::datatable(study_data(), options = list(pageLength = 10), rownames = FALSE)
  })

  output$data_msg <- renderText({
    df <- study_data()
    if (anyNA(df)) {
      return("Warning: There are values in TP/FN/FP/TN that cannot be read as numbers (missing values). Please check the column mapping.")
    }
    if (any(df < 0)) {
      return("Warning: There are negative values in TP/FN/FP/TN.")
    }
    sprintf("Number of studies: %d (The table below is the data passed directly to reitsma.sroc() / srocsa.t())", nrow(df))
  })

  ## ---- 2. fit: no publication bias (reitsma.sroc) ---------------------------

  fit_main <- eventReactive(input$run, {
    df <- study_data()
    validate(
      need(nrow(df) >= 2, "Too few studies: at least 2 are required."),
      need(!anyNA(df), "TP/FN/FP/TN contain missing values (NA)."),
      need(all(df >= 0), "TP/FN/FP/TN contain negative values."),
      need(input$cc_value > 0 && input$cc_value <= 1, "cc.value must be in the range (0, 1]."),
      need(input$ci_level > 0 && input$ci_level < 1, "ci.level must be in the range (0, 1).")
    )
    DiagProgMeta::reitsma.sroc(
      TP = TP, FN = FN, FP = FP, TN = TN, data = df,
      cc.value = input$cc_value, cc.type = input$cc_type, ci.level = input$ci_level
    )
  })

  output$print_main  <- renderPrint(print(fit_main()))
  output$summary_main <- renderPrint(summary(fit_main()))

  plot_main_expr <- function() {
    plot(
      fit_main(),
      plot.ci = input$plot_ci_main,
      spoint.region = input$spoint_region_main,
      add.study = input$add_study_main,
      study.size = input$study_size_main
    )
  }
  output$plot_main <- renderPlot(plot_main_expr())

  output$dl_plot_main <- downloadHandler(
    filename = function() "sroc_no_bias.png",
    content = function(file) {
      grDevices::png(file, width = 1800, height = 1350, res = 220)
      on.exit(grDevices::dev.off())
      plot_main_expr()
    }
  )

  ## ---- 3. fit: sensitivity analysis for publication bias (srocsa.t) --------

  fit_sa <- eventReactive(input$run, {
    req(input$run_sa)
    df <- study_data()
    sa.p <- suppressWarnings(as.numeric(trimws(strsplit(input$sa_p, ",")[[1]])))
    validate(
      need(nrow(df) >= 2, "Too few studies: at least 2 are required."),
      need(!anyNA(df), "TP/FN/FP/TN contain missing values (NA)."),
      need(length(sa.p) >= 1 && !anyNA(sa.p), "Please enter the values of sa.p correctly (e.g., 1, 0.8, 0.6)."),
      need(all(sa.p > 0 & sa.p <= 1), "sa.p must be in the range (0, 1].")
    )
    c11 <- if (isTRUE(input$estimate_c11)) NULL else input$c11_value
    DiagProgMeta::srocsa.t(
      TP = TP, FN = FN, TN = TN, FP = FP, data = df,
      cc.value = input$cc_value, cc.type = input$cc_type,
      sa.p = sa.p, c11 = c11, ci.level = input$ci_level
    )
  })

  output$print_sa <- renderPrint({
    req(fit_sa())
    print(fit_sa())
  })

  sa_summary <- reactive({
    req(fit_sa())
    utils::capture.output(res <- summary(fit_sa()))
    res
  })

  output$sauc_table <- renderDT({
    DT::datatable(.round_df(sa_summary()$sauc), options = list(dom = "t"), rownames = FALSE)
  })

  output$coef_table <- renderDT({
    DT::datatable(.round_df(sa_summary()$coef), options = list(scrollX = TRUE, dom = "t"), rownames = FALSE)
  })

  plot_sa_expr <- function() {
    plot(fit_sa(), add.study = input$add_study_sa, plot.ci = input$plot_ci_sa)
  }
  output$plot_sa <- renderPlot({
    req(fit_sa())
    plot_sa_expr()
  })

  output$dl_plot_sa <- downloadHandler(
    filename = function() "sroc_sensitivity_analysis.png",
    content = function(file) {
      req(fit_sa())
      grDevices::png(file, width = 1800, height = 1350, res = 220)
      on.exit(grDevices::dev.off())
      plot_sa_expr()
    }
  )

  output$dl_sauc_csv <- downloadHandler(
    filename = function() "sauc_sensitivity_analysis.csv",
    content = function(file) {
      utils::write.csv(sa_summary()$sauc, file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)
