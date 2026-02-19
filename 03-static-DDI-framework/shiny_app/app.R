# =============================================================================
# DDI Static Screening Calculator — Interactive Shiny App
# =============================================================================
# Single-file Shiny app using bslib for modern layout.
# Reuses all computation functions from analysis/00_setup.R.
# =============================================================================

library(shiny)
library(bslib)
library(DT)
library(plotly)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)

# --- Source computation functions from the analysis pipeline -----------------
APP_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
  error = function(e) getwd()
)
if (is.na(APP_DIR) || APP_DIR == "") APP_DIR <- getwd()
PROJECT_ROOT <- dirname(APP_DIR)

setup_path <- file.path(PROJECT_ROOT, "analysis", "00_setup.R")

# Minimal sourcing: load only computation functions, suppress side effects
local({
  # Temporarily override log_msg so sourcing doesn't create log files
  env <- new.env(parent = globalenv())
  env$log_msg <- function(...) invisible(NULL)

  # Load packages needed by 00_setup.R
  suppressPackageStartupMessages({
    library(tidyverse)
    library(readr)
    library(dplyr)
    library(scales)
  })

  # Source into our environment
  source(setup_path, local = env)

  # Export computation functions to global environment
  fns <- c("cmax_to_uM", "compute_Imax_u", "compute_Igut",
           "compute_R1", "compute_R1_gut", "compute_TDI_factor",
           "compute_R3", "compute_AUCR_hepatic", "compute_AUCR_gut",
           "compute_AUCR_total", "classify_ddi", "screen_transporter",
           "parse_inputs", "COLORS_RISK", "COLOR_FLAG", "COLOR_PASS")
  for (fn in fns) {
    if (exists(fn, envir = env)) assign(fn, get(fn, envir = env), envir = globalenv())
  }
})

# AUCR wrapper (from 02_sensitivity_uncertainty.R)
compute_aucr_from_params <- function(Ki, Imax_u, Igut, fm, fg, kinact, KI, kdeg) {
  R1      <- compute_R1(Imax_u, Ki)
  R1_gut  <- compute_R1_gut(Igut, Ki)
  tdi_fac <- compute_TDI_factor(kinact, Imax_u, KI, kdeg)
  aucr_h  <- compute_AUCR_hepatic(fm, R1, tdi_fac)
  aucr_g  <- compute_AUCR_gut(fg, R1_gut)
  list(total = aucr_h * aucr_g, hepatic = aucr_h, gut = aucr_g)
}

# --- Preset data paths ------------------------------------------------------
PRESET_FILES <- list(
  "Ketoconazole (Benchmark)" = file.path(PROJECT_ROOT, "data", "inputs_midazolam_case.csv"),
  "Sotorasib (Case B)"       = file.path(PROJECT_ROOT, "data", "inputs_sotorasib_case.csv")
)

# --- Badge styling helpers ---------------------------------------------------
classification_color <- function(cls) {
  switch(cls,
    "Strong"         = "danger",
    "Moderate"       = "warning",
    "Weak"           = "warning",
    "No interaction" = "success",
    "secondary"
  )
}

classification_action <- function(cls) {
  switch(cls,
    "Strong"         = "Clinical DDI study required per FDA guidance.",
    "Moderate"       = "Clinical DDI study or PBPK modeling recommended.",
    "Weak"           = "Evaluate with PBPK; clinical study may not be needed.",
    "No interaction" = "No further DDI evaluation required.",
    "Insufficient data to classify."
  )
}

# --- Tooltip helper ----------------------------------------------------------
tip <- function(text, tooltip_text) {
  tagList(
    text,
    tags$span(
      class = "ms-1",
      style = "cursor: help; color: #6c757d;",
      title = tooltip_text,
      icon("circle-question")
    )
  )
}

# =============================================================================
# UI
# =============================================================================

ui <- page_sidebar(
  title = "DDI Static Screening Calculator",
  theme = bs_theme(bootswatch = "flatly", version = 5),

  # --- Sidebar: Inputs -------------------------------------------------------

sidebar = sidebar(
    width = 350,

    # Preset selector
    selectInput("preset", "Load Preset",
      choices = c("Custom", names(PRESET_FILES)),
      selected = "Custom"
    ),
    conditionalPanel(
      condition = "input.preset != 'Custom'",
      actionButton("unlock_btn", "Unlock Fields for Editing",
                    class = "btn-sm btn-outline-secondary mb-2")
    ),

    hr(),
    h6("Perpetrator"),
    textInput("drug_name", "Drug Name", "My Compound"),
    numericInput("dose_mg", tip("Dose (mg)", "Oral dose in milligrams"), 400, min = 0),
    numericInput("mw", tip("MW (g/mol)", "Molecular weight"), 531.43, min = 1),
    numericInput("cmax", tip("Cmax (ug/mL)", "Total maximum plasma concentration"), 8.07, min = 0),
    numericInput("fu", tip("fu (fraction unbound)", "Fraction unbound in plasma"), 0.01, min = 0, max = 1, step = 0.01),

    hr(),
    h6("CYP Inhibition (Ki, uM)"),
    numericInput("ki_3a4", tip("CYP3A4 Ki", "R1 = 1 + [I]max,u / Ki; threshold >= 1.02"), 0.015, min = 0, step = 0.01),
    numericInput("ki_2c8", "CYP2C8 Ki", NA_real_, min = 0),
    numericInput("ki_2c9", "CYP2C9 Ki", 10.0, min = 0),
    numericInput("ki_2c19", "CYP2C19 Ki", 3.4, min = 0),
    numericInput("ki_2d6", "CYP2D6 Ki", 100, min = 0),
    numericInput("ki_1a2", "CYP1A2 Ki", 100, min = 0),

    hr(),
    h6("Time-Dependent Inhibition (CYP3A4)"),
    numericInput("kinact", tip("kinact (1/min)", "Maximum inactivation rate constant"), 0.048, min = 0, step = 0.001),
    numericInput("KI_tdi", tip("KI (uM)", "Inhibitor concentration producing 50% of kinact"), 0.86, min = 0, step = 0.01),
    numericInput("kdeg", tip("kdeg (1/min)", "Hepatic enzyme degradation rate (CYP3A4 t1/2 ~36h)"), 0.00032, min = 0, step = 0.0001),

    hr(),
    h6("Induction (CYP3A4)"),
    numericInput("emax", tip("Emax (fold)", "Maximum fold induction"), NA_real_, min = 0),
    numericInput("ec50_ind", tip("EC50 (uM)", "Inducer concentration producing 50% of Emax"), NA_real_, min = 0),
    sliderInput("d_scaling", tip("d-scaling factor", "R3 = 1/(1 + d*Emax*[I]/(EC50+[I])); threshold <= 0.8"),
                min = 0.1, max = 1.0, value = 1.0, step = 0.1),

    hr(),
    h6("Victim Parameters (for AUCR)"),
    numericInput("fm", tip("fm (CYP3A4)", "Fraction metabolized by CYP3A4 in liver"), 0.94, min = 0, max = 1, step = 0.01),
    numericInput("fg", tip("fg (gut extraction)", "Fraction metabolized by CYP3A4 in gut wall (= 1 - Fg)"), 0.43, min = 0, max = 1, step = 0.01),

    hr(),
    h6("Transporter IC50 (uM)"),
    numericInput("ic50_pgp", "P-gp", NA_real_, min = 0),
    numericInput("ic50_bcrp", "BCRP", NA_real_, min = 0),
    numericInput("ic50_oatp1b1", "OATP1B1", NA_real_, min = 0),
    numericInput("ic50_oatp1b3", "OATP1B3", NA_real_, min = 0),
    numericInput("ic50_mate1", "MATE1", NA_real_, min = 0),
    numericInput("ic50_mate2k", "MATE2-K", NA_real_, min = 0),
    numericInput("ic50_oat1", "OAT1", NA_real_, min = 0),
    numericInput("ic50_oat3", "OAT3", NA_real_, min = 0),
    numericInput("ic50_oct2", "OCT2", NA_real_, min = 0)
  ),

  # --- Main panel -------------------------------------------------------------
  layout_columns(
    col_widths = 12,

    # Top card: Regulatory Conclusion
    card(
      card_header("Regulatory Conclusion"),
      card_body(
        uiOutput("conclusion_badge"),
        textOutput("conclusion_action")
      )
    ),

    # Tabbed panels
    navset_card_tab(
      id = "main_tabs",

      # Tab 1: Screening Results (CORE)
      nav_panel(
        "Screening Results",
        h5("CYP Enzyme Screening"),
        DTOutput("cyp_table"),
        br(),
        h5("Induction Screening"),
        uiOutput("induction_card"),
        br(),
        h5("AUCR Breakdown"),
        uiOutput("aucr_card"),
        br(),
        downloadButton("download_report", "Download Report (CSV)", class = "btn-primary mt-2")
      ),

      # Tab 2: Sensitivity
      nav_panel(
        "Sensitivity",
        h5("Tornado Plot (One-at-a-Time Parameter Variation)"),
        plotlyOutput("tornado_plot", height = "400px"),
        br(),
        h5("Ki vs [I]max,u Heatmap"),
        plotlyOutput("heatmap_plot", height = "450px")
      ),

      # Tab 3: Monte Carlo
      nav_panel(
        "Monte Carlo",
        layout_columns(
          col_widths = c(4, 4, 4),
          numericInput("mc_n", "Number of Samples", 1000, min = 100, max = 10000, step = 100),
          sliderInput("mc_cv", "CV for Parameter Uncertainty (%)", min = 10, max = 60, value = 30, step = 5),
          actionButton("mc_run", "Run Simulation", class = "btn-primary mt-4")
        ),
        plotlyOutput("mc_histogram", height = "400px"),
        br(),
        uiOutput("mc_summary_card")
      ),

      # Tab 4: Transporter Dashboard
      nav_panel(
        "Transporter Dashboard",
        plotlyOutput("transporter_bar", height = "450px")
      )
    )
  )
)

# =============================================================================
# SERVER
# =============================================================================

server <- function(input, output, session) {

  # --- Preset loading --------------------------------------------------------
  preset_locked <- reactiveVal(FALSE)

  observeEvent(input$preset, {
    req(input$preset != "Custom")
    preset_locked(TRUE)

    csv_path <- PRESET_FILES[[input$preset]]
    if (!file.exists(csv_path)) return()

    params <- parse_inputs(csv_path)

    updateTextInput(session, "drug_name", value = params$perpetrator)
    updateNumericInput(session, "dose_mg", value = params$dose_mg)
    updateNumericInput(session, "mw", value = params$mw)
    updateNumericInput(session, "cmax", value = params$cmax_total)
    updateNumericInput(session, "fu", value = params$fu)
    updateNumericInput(session, "fm", value = ifelse(is.na(params$fm_cyp3a4), 0.94, params$fm_cyp3a4))
    updateNumericInput(session, "fg", value = ifelse(is.na(params$fg_cyp3a4), 0.43, params$fg_cyp3a4))

    # CYP Ki values
    cyp_map <- list(ki_3a4 = "CYP3A4", ki_2c8 = "CYP2C8", ki_2c9 = "CYP2C9",
                    ki_2c19 = "CYP2C19", ki_2d6 = "CYP2D6", ki_1a2 = "CYP1A2")
    for (field in names(cyp_map)) {
      enz <- cyp_map[[field]]
      val <- if (!is.null(params$cyp[[enz]])) params$cyp[[enz]]$ki else NA_real_
      updateNumericInput(session, field, value = val)
    }

    # TDI
    cyp3a4 <- params$cyp[["CYP3A4"]]
    updateNumericInput(session, "kinact", value = if (!is.null(cyp3a4)) cyp3a4$kinact else NA_real_)
    updateNumericInput(session, "KI_tdi", value = if (!is.null(cyp3a4)) cyp3a4$KI else NA_real_)
    updateNumericInput(session, "kdeg", value = ifelse(is.na(params$kdeg_hepatic), 0.00032, params$kdeg_hepatic))

    # Induction
    if (params$induction$available) {
      updateNumericInput(session, "emax", value = params$induction$emax)
      updateNumericInput(session, "ec50_ind", value = params$induction$ec50)
      updateSliderInput(session, "d_scaling", value = params$induction$d)
    } else {
      updateNumericInput(session, "emax", value = NA_real_)
      updateNumericInput(session, "ec50_ind", value = NA_real_)
      updateSliderInput(session, "d_scaling", value = 1.0)
    }

    # Transporters
    trans_map <- list(ic50_pgp = "P-gp", ic50_bcrp = "BCRP",
                      ic50_oatp1b1 = "OATP1B1", ic50_oatp1b3 = "OATP1B3",
                      ic50_mate1 = "MATE1", ic50_mate2k = "MATE2-K",
                      ic50_oat1 = "OAT1", ic50_oat3 = "OAT3", ic50_oct2 = "OCT2")
    for (field in names(trans_map)) {
      tn <- trans_map[[field]]
      val <- if (!is.null(params$transporters[[tn]])) params$transporters[[tn]] else NA_real_
      updateNumericInput(session, field, value = val)
    }
  })

  observeEvent(input$unlock_btn, {
    preset_locked(FALSE)
    updateSelectInput(session, "preset", selected = "Custom")
  })

  # --- Derived concentrations ------------------------------------------------
  Imax_u <- reactive({
    req(input$cmax, input$fu, input$mw)
    compute_Imax_u(input$cmax, input$fu, input$mw)
  })

  Igut <- reactive({
    req(input$dose_mg, input$mw)
    compute_Igut(input$dose_mg, input$mw)
  })

  # --- CYP Screening Table (Tab 1) ------------------------------------------
  cyp_data <- reactive({
    req(Imax_u(), Igut())
    imax <- Imax_u()
    igut <- Igut()

    enzymes <- list(
      list(name = "CYP3A4", ki = input$ki_3a4, has_tdi = TRUE),
      list(name = "CYP2C8", ki = input$ki_2c8, has_tdi = FALSE),
      list(name = "CYP2C9", ki = input$ki_2c9, has_tdi = FALSE),
      list(name = "CYP2C19", ki = input$ki_2c19, has_tdi = FALSE),
      list(name = "CYP2D6", ki = input$ki_2d6, has_tdi = FALSE),
      list(name = "CYP1A2", ki = input$ki_1a2, has_tdi = FALSE)
    )

    map_dfr(enzymes, function(e) {
      ki <- e$ki
      r1 <- if (!is.na(ki) && ki > 0) compute_R1(imax, ki) else NA_real_
      r1g <- if (!is.na(ki) && ki > 0) compute_R1_gut(igut, ki) else NA_real_
      tdi <- if (e$has_tdi) compute_TDI_factor(input$kinact, imax, input$KI_tdi, input$kdeg) else NA_real_

      tibble(
        Enzyme = e$name,
        `Ki (uM)` = ki,
        R1 = round(r1, 3),
        `R1 Flag` = ifelse(!is.na(r1) && r1 >= 1.02, "FLAGGED", "Pass"),
        `R1,gut` = round(r1g, 1),
        `R1,gut Flag` = ifelse(!is.na(r1g) && r1g >= 11, "FLAGGED", "Pass"),
        `TDI Factor` = if (!is.na(tdi)) round(tdi, 4) else NA_real_,
        `TDI Flag` = if (!is.na(tdi)) ifelse(tdi < 0.5, "FLAGGED", "Pass") else NA_character_
      )
    })
  })

  output$cyp_table <- renderDT({
    df <- cyp_data()
    datatable(df, rownames = FALSE, options = list(dom = 't', pageLength = 10)) %>%
      formatStyle("R1 Flag",
        backgroundColor = styleEqual(c("FLAGGED", "Pass"), c("#f8d7da", "#d1e7dd")),
        fontWeight = styleEqual("FLAGGED", "bold")) %>%
      formatStyle("R1,gut Flag",
        backgroundColor = styleEqual(c("FLAGGED", "Pass"), c("#f8d7da", "#d1e7dd")),
        fontWeight = styleEqual("FLAGGED", "bold")) %>%
      formatStyle("TDI Flag",
        backgroundColor = styleEqual(c("FLAGGED", "Pass"), c("#f8d7da", "#d1e7dd")),
        fontWeight = styleEqual("FLAGGED", "bold"))
  })

  # --- Induction card --------------------------------------------------------
  r3_val <- reactive({
    if (is.na(input$emax) || is.na(input$ec50_ind)) return(NA_real_)
    compute_R3(input$emax, Imax_u(), input$ec50_ind, input$d_scaling)
  })

  output$induction_card <- renderUI({
    r3 <- r3_val()
    if (is.na(r3)) {
      card(card_body("No induction data provided."), class = "bg-light")
    } else {
      flagged <- r3 <= 0.8
      card(
        card_body(
          tags$p(tags$strong("R3 = "), sprintf("%.4f", r3)),
          tags$p(
            if (flagged) tags$span(class = "badge bg-danger", "FLAGGED (R3 <= 0.8)")
            else tags$span(class = "badge bg-success", "Pass (R3 > 0.8)")
          ),
          tags$p(tags$small(sprintf("Emax=%.1f, EC50=%.1f uM, d=%.1f, [I]max,u=%.3f uM",
                                     input$emax, input$ec50_ind, input$d_scaling, Imax_u())))
        )
      )
    }
  })

  # --- AUCR Breakdown --------------------------------------------------------
  aucr_vals <- reactive({
    req(input$fm, Imax_u(), Igut(), input$ki_3a4)
    ki <- input$ki_3a4
    if (is.na(ki) || ki <= 0) return(list(hepatic = NA, gut = NA, total = NA, class = "Not assessable"))

    r1 <- compute_R1(Imax_u(), ki)
    r1g <- compute_R1_gut(Igut(), ki)
    tdi <- compute_TDI_factor(input$kinact, Imax_u(), input$KI_tdi, input$kdeg)
    aucr_h <- compute_AUCR_hepatic(input$fm, r1, tdi)
    aucr_g <- compute_AUCR_gut(input$fg, r1g)
    total <- aucr_h * aucr_g
    cls <- classify_ddi(total)
    list(hepatic = aucr_h, gut = aucr_g, total = total, class = cls)
  })

  output$aucr_card <- renderUI({
    vals <- aucr_vals()
    if (is.na(vals$total)) {
      return(card(card_body("Insufficient data to compute AUCR."), class = "bg-light"))
    }
    cls_color <- classification_color(vals$class)
    card(
      card_body(
        layout_columns(
          col_widths = c(4, 4, 4),
          value_box(
            title = "AUCR Hepatic",
            value = sprintf("%.2f", vals$hepatic),
            theme = "light"
          ),
          value_box(
            title = "AUCR Gut",
            value = sprintf("%.2f", vals$gut),
            theme = "light"
          ),
          value_box(
            title = "AUCR Total",
            value = sprintf("%.2f", vals$total),
            theme = cls_color
          )
        ),
        tags$p(
          tags$strong("FDA Classification: "),
          tags$span(class = paste0("badge bg-", cls_color, " fs-6"), vals$class)
        )
      )
    )
  })

  # --- Top Regulatory Conclusion Badge ---------------------------------------
  output$conclusion_badge <- renderUI({
    vals <- aucr_vals()
    cls <- vals$class
    cls_color <- classification_color(cls)

    # Also check induction
    r3 <- r3_val()
    r3_text <- if (!is.na(r3) && r3 <= 0.8) {
      tags$p(tags$span(class = "badge bg-warning", "CYP3A4 Induction Flagged"),
             sprintf(" R3 = %.4f", r3))
    }

    tagList(
      tags$div(
        style = "text-align: center; margin-bottom: 10px;",
        if (!is.na(vals$total)) {
          tagList(
            tags$span(class = paste0("badge bg-", cls_color, " fs-3 px-4 py-2"), cls),
            tags$p(class = "mt-2 fs-5", sprintf("Predicted AUCR = %.2f", vals$total))
          )
        } else {
          tags$span(class = "badge bg-secondary fs-4 px-4 py-2", "Enter Parameters")
        }
      ),
      r3_text
    )
  })

  output$conclusion_action <- renderText({
    vals <- aucr_vals()
    classification_action(vals$class)
  })

  # --- Tab 2: Sensitivity (Tornado) ------------------------------------------
  tornado_data <- reactive({
    req(Imax_u(), Igut(), input$ki_3a4, input$fm, input$fg)
    imax <- Imax_u()
    igut <- Igut()

    base_args <- list(Ki = input$ki_3a4, Imax_u = imax, Igut = igut,
                      fm = input$fm, fg = input$fg,
                      kinact = ifelse(is.na(input$kinact), 0, input$kinact),
                      KI = ifelse(is.na(input$KI_tdi), 1, input$KI_tdi),
                      kdeg = ifelse(is.na(input$kdeg), 0.00032, input$kdeg))

    base_aucr <- do.call(compute_aucr_from_params, base_args)$total
    if (is.na(base_aucr)) return(NULL)

    params <- list(
      list(name = "Ki", label = "Ki CYP3A4 (uM)", lo = base_args$Ki * 0.1, hi = base_args$Ki * 10),
      list(name = "Imax_u", label = "[I]max,u (uM)", lo = imax * 0.5, hi = imax * 1.5),
      list(name = "fm", label = "fm CYP3A4", lo = 0.5, hi = min(0.99, input$fm + 0.01)),
      list(name = "fg", label = "fg (gut)", lo = 0.0, hi = min(0.8, max(input$fg + 0.1, 0.7))),
      list(name = "kinact", label = "kinact (1/min)", lo = base_args$kinact * 0.5, hi = base_args$kinact * 1.5),
      list(name = "KI", label = "KI TDI (uM)", lo = base_args$KI * 0.5, hi = base_args$KI * 1.5),
      list(name = "kdeg", label = "kdeg (1/min)", lo = base_args$kdeg * 0.5, hi = base_args$kdeg * 1.5)
    )

    map_dfr(params, function(p) {
      args_lo <- base_args; args_lo[[p$name]] <- p$lo
      args_hi <- base_args; args_hi[[p$name]] <- p$hi
      aucr_lo <- do.call(compute_aucr_from_params, args_lo)$total
      aucr_hi <- do.call(compute_aucr_from_params, args_hi)$total
      tibble(
        param = p$name, label = p$label,
        aucr_lo = aucr_lo, aucr_hi = aucr_hi,
        base_aucr = base_aucr,
        delta = max(abs(aucr_lo - base_aucr), abs(aucr_hi - base_aucr))
      )
    }) %>% arrange(delta)
  })

  output$tornado_plot <- renderPlotly({
    td <- tornado_data()
    if (is.null(td) || nrow(td) == 0) return(plotly_empty())

    td$label <- factor(td$label, levels = td$label)

    plot_ly(td, y = ~label) %>%
      add_segments(x = ~aucr_lo, xend = ~aucr_hi, y = ~label, yend = ~label,
                   line = list(width = 14, color = "#3498db"),
                   hovertemplate = "%{y}<br>Range: %{x:.1f} - %{customdata:.1f}<extra></extra>",
                   customdata = ~aucr_hi) %>%
      add_markers(x = ~base_aucr, y = ~label, marker = list(size = 10, color = "black"),
                  name = "Base Case", hoverinfo = "text",
                  text = ~paste0("Base AUCR: ", round(base_aucr, 2))) %>%
      layout(
        title = "Tornado Sensitivity Analysis",
        xaxis = list(title = "Predicted AUCR"),
        yaxis = list(title = ""),
        showlegend = FALSE,
        margin = list(l = 150)
      )
  })

  # Heatmap
  output$heatmap_plot <- renderPlotly({
    req(Imax_u(), Igut(), input$fm, input$fg)
    imax <- Imax_u()

    ki_seq <- 10^seq(log10(0.001), log10(10), length.out = 50)
    imax_seq <- 10^seq(log10(0.01), log10(10), length.out = 50)

    grid <- expand.grid(Ki = ki_seq, Im = imax_seq)

    grid$AUCR <- mapply(function(ki, im) {
      compute_aucr_from_params(Ki = ki, Imax_u = im, Igut = Igut(),
                                fm = input$fm, fg = input$fg,
                                kinact = ifelse(is.na(input$kinact), 0, input$kinact),
                                KI = ifelse(is.na(input$KI_tdi), 1, input$KI_tdi),
                                kdeg = ifelse(is.na(input$kdeg), 0.00032, input$kdeg))$total
    }, grid$Ki, grid$Im)

    # Cap for display
    grid$AUCR_cap <- pmin(grid$AUCR, 100)

    z_mat <- matrix(log10(grid$AUCR_cap), nrow = 50, ncol = 50)

    plot_ly(x = log10(ki_seq), y = log10(imax_seq), z = t(z_mat),
            type = "heatmap",
            colorscale = list(
              list(0, "#2CA02C"), list(0.3, "#FFD700"),
              list(0.5, "#FF8C00"), list(1, "#D62728")
            ),
            colorbar = list(title = "log10(AUCR)"),
            hovertemplate = "Ki: %{x:.2f} (log10)<br>[I]max,u: %{y:.2f} (log10)<br>log10(AUCR): %{z:.2f}<extra></extra>"
    ) %>%
      add_markers(x = log10(input$ki_3a4), y = log10(imax),
                  marker = list(size = 14, color = "white", symbol = "x",
                                line = list(color = "black", width = 2)),
                  name = "Your Compound",
                  hoverinfo = "text",
                  text = paste0("Ki=", input$ki_3a4, " uM, [I]max,u=", round(imax, 3), " uM")) %>%
      layout(
        title = "Ki vs [I]max,u: AUCR Heatmap",
        xaxis = list(title = "log10(Ki, uM)"),
        yaxis = list(title = "log10([I]max,u, uM)"),
        showlegend = FALSE
      )
  })

  # --- Tab 3: Monte Carlo ----------------------------------------------------
  mc_results <- reactiveVal(NULL)

  observeEvent(input$mc_run, {
    req(Imax_u(), Igut(), input$ki_3a4, input$fm, input$fg)
    n <- input$mc_n
    cv <- input$mc_cv / 100
    imax <- Imax_u()
    igut <- Igut()

    set.seed(42)
    ki_samp <- rlnorm(n, log(input$ki_3a4), cv)
    imax_samp <- rlnorm(n, log(imax), cv)
    fm_samp <- pmin(0.99, pmax(0.3, rnorm(n, input$fm, 0.05)))
    fg_samp <- pmin(0.80, pmax(0.0, rnorm(n, input$fg, 0.10)))
    kinact_samp <- if (!is.na(input$kinact)) rlnorm(n, log(input$kinact), cv) else rep(0, n)
    KI_samp <- if (!is.na(input$KI_tdi)) rlnorm(n, log(input$KI_tdi), cv) else rep(1, n)
    kdeg_samp <- if (!is.na(input$kdeg)) rlnorm(n, log(input$kdeg), cv) else rep(0.00032, n)

    aucr_vec <- numeric(n)
    for (i in seq_len(n)) {
      res <- compute_aucr_from_params(ki_samp[i], imax_samp[i], igut,
                                       fm_samp[i], fg_samp[i],
                                       kinact_samp[i], KI_samp[i], kdeg_samp[i])
      aucr_vec[i] <- res$total
    }

    cls_vec <- sapply(aucr_vec, classify_ddi)
    mc_results(list(aucr = aucr_vec, class = cls_vec))
  })

  output$mc_histogram <- renderPlotly({
    res <- mc_results()
    if (is.null(res)) return(plotly_empty())

    df <- tibble(AUCR = res$aucr, Class = res$class)

    # Color map
    color_map <- c("No interaction" = "#2CA02C", "Weak" = "#FFD700",
                   "Moderate" = "#FF8C00", "Strong" = "#D62728")

    plot_ly(df, x = ~AUCR, color = ~Class, colors = color_map,
            type = "histogram", nbinsx = 50,
            hovertemplate = "AUCR: %{x:.1f}<br>Count: %{y}<extra></extra>") %>%
      layout(
        title = paste0("Monte Carlo AUCR Distribution (N=", length(res$aucr), ")"),
        xaxis = list(title = "Predicted AUCR"),
        yaxis = list(title = "Count"),
        barmode = "stack",
        shapes = list(
          list(type = "line", x0 = 1.25, x1 = 1.25, y0 = 0, y1 = 1, yref = "paper",
               line = list(color = "#FFD700", dash = "dot", width = 1)),
          list(type = "line", x0 = 2, x1 = 2, y0 = 0, y1 = 1, yref = "paper",
               line = list(color = "#FF8C00", dash = "dot", width = 1)),
          list(type = "line", x0 = 5, x1 = 5, y0 = 0, y1 = 1, yref = "paper",
               line = list(color = "#D62728", dash = "dot", width = 1))
        )
      )
  })

  output$mc_summary_card <- renderUI({
    res <- mc_results()
    if (is.null(res)) return(NULL)

    aucr <- res$aucr
    cls <- res$class
    med <- median(aucr)
    p5 <- quantile(aucr, 0.05)
    p95 <- quantile(aucr, 0.95)
    pct_strong <- mean(cls == "Strong") * 100
    pct_moderate <- mean(cls == "Moderate") * 100
    pct_weak <- mean(cls == "Weak") * 100
    pct_none <- mean(cls == "No interaction") * 100

    card(
      card_header("Monte Carlo Summary"),
      card_body(
        layout_columns(
          col_widths = c(4, 4, 4),
          value_box(title = "Median AUCR", value = sprintf("%.1f", med), theme = "primary"),
          value_box(title = "5th Percentile", value = sprintf("%.1f", p5), theme = "light"),
          value_box(title = "95th Percentile", value = sprintf("%.1f", p95), theme = "light")
        ),
        tags$table(
          class = "table table-sm mt-2",
          tags$thead(tags$tr(tags$th("Classification"), tags$th("Percentage"))),
          tags$tbody(
            tags$tr(tags$td(tags$span(class = "badge bg-danger", "Strong")),
                    tags$td(sprintf("%.1f%%", pct_strong))),
            tags$tr(tags$td(tags$span(class = "badge bg-warning", "Moderate")),
                    tags$td(sprintf("%.1f%%", pct_moderate))),
            tags$tr(tags$td(tags$span(class = "badge bg-warning text-dark", "Weak")),
                    tags$td(sprintf("%.1f%%", pct_weak))),
            tags$tr(tags$td(tags$span(class = "badge bg-success", "No interaction")),
                    tags$td(sprintf("%.1f%%", pct_none)))
          )
        )
      )
    )
  })

  # --- Tab 4: Transporter Dashboard ------------------------------------------
  transporter_data <- reactive({
    req(Imax_u(), Igut())
    imax <- Imax_u()
    igut <- Igut()

    transporters <- list(
      list(name = "P-gp", ic50 = input$ic50_pgp, I = igut, threshold = 10, type = "Intestinal"),
      list(name = "BCRP", ic50 = input$ic50_bcrp, I = igut, threshold = 10, type = "Intestinal"),
      list(name = "OATP1B1", ic50 = input$ic50_oatp1b1, I = imax, threshold = 0.1, type = "Hepatic"),
      list(name = "OATP1B3", ic50 = input$ic50_oatp1b3, I = imax, threshold = 0.1, type = "Hepatic"),
      list(name = "MATE1", ic50 = input$ic50_mate1, I = imax, threshold = 0.1, type = "Renal"),
      list(name = "MATE2-K", ic50 = input$ic50_mate2k, I = imax, threshold = 0.1, type = "Renal"),
      list(name = "OAT1", ic50 = input$ic50_oat1, I = imax, threshold = 0.1, type = "Renal"),
      list(name = "OAT3", ic50 = input$ic50_oat3, I = imax, threshold = 0.1, type = "Renal"),
      list(name = "OCT2", ic50 = input$ic50_oct2, I = imax, threshold = 0.1, type = "Renal")
    )

    map_dfr(transporters, function(t) {
      if (is.na(t$ic50) || t$ic50 <= 0) {
        return(tibble(Transporter = t$name, Type = t$type,
                      Ratio = NA_real_, Threshold = t$threshold,
                      Flag = "No data"))
      }
      res <- screen_transporter(t$I, t$ic50, t$threshold)
      tibble(Transporter = t$name, Type = t$type,
             Ratio = res$ratio, Threshold = t$threshold,
             Flag = ifelse(res$flag, "FLAGGED", "Pass"))
    })
  })

  output$transporter_bar <- renderPlotly({
    td <- transporter_data()
    td <- td %>% filter(!is.na(Ratio))
    if (nrow(td) == 0) return(plotly_empty())

    td$color <- ifelse(td$Flag == "FLAGGED", "#D62728", "#2CA02C")
    td$Transporter <- factor(td$Transporter, levels = rev(td$Transporter))

    # Use log scale for ratio
    td$log_ratio <- log10(td$Ratio)
    td$log_threshold <- log10(td$Threshold)

    plot_ly(td, y = ~Transporter, x = ~log_ratio, type = "bar",
            orientation = "h",
            marker = list(color = ~color),
            hovertemplate = "%{y}<br>Ratio: %{customdata:.3f}<extra></extra>",
            customdata = ~Ratio) %>%
      layout(
        title = "Transporter Screening: [I]/IC50 Ratio vs Threshold",
        xaxis = list(title = "log10([I]/IC50)"),
        yaxis = list(title = ""),
        shapes = lapply(unique(td$log_threshold), function(lt) {
          list(type = "line", x0 = lt, x1 = lt, y0 = -0.5, y1 = nrow(td) - 0.5,
               line = list(color = "black", dash = "dash", width = 1.5))
        }),
        margin = list(l = 100)
      )
  })

  # --- Download Report -------------------------------------------------------
  output$download_report <- downloadHandler(
    filename = function() {
      paste0("DDI_screening_", gsub(" ", "_", input$drug_name), "_",
             format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      cyp <- cyp_data()
      trans <- transporter_data()
      vals <- aucr_vals()
      r3 <- r3_val()

      # Build summary
      summary_rows <- tibble(
        Section = c("Drug", "Dose (mg)", "MW (g/mol)", "Cmax (ug/mL)", "fu",
                     "[I]max,u (uM)", "[I]gut (uM)",
                     "AUCR Hepatic", "AUCR Gut", "AUCR Total", "FDA Classification",
                     "R3 (Induction)"),
        Value = c(input$drug_name, input$dose_mg, input$mw, input$cmax, input$fu,
                  round(Imax_u(), 4), round(Igut(), 1),
                  round(vals$hepatic, 2), round(vals$gut, 2),
                  round(vals$total, 2), vals$class,
                  ifelse(is.na(r3), "Not assessed", round(r3, 4)))
      )

      # Combine
      out <- bind_rows(
        summary_rows %>% mutate(Table = "Summary"),
        cyp %>% mutate(Table = "CYP Screening") %>% rename_with(~"Value", .cols = 1),
        trans %>% mutate(Table = "Transporter Screening") %>% rename_with(~"Value", .cols = 1)
      )

      # Simpler: write multiple sections
      con <- file(file, "w")
      writeLines("=== DDI Screening Report ===", con)
      writeLines(paste0("Drug: ", input$drug_name), con)
      writeLines(paste0("Date: ", Sys.Date()), con)
      writeLines("", con)
      close(con)

      # Append tables
      write_csv(summary_rows, file, append = TRUE)
      write_csv(tibble(Section = ""), file, append = TRUE)
      write_csv(tibble(Section = "--- CYP Screening ---"), file, append = TRUE)
      write_csv(cyp, file, append = TRUE)
      write_csv(tibble(Section = ""), file, append = TRUE)
      write_csv(tibble(Section = "--- Transporter Screening ---"), file, append = TRUE)
      write_csv(trans, file, append = TRUE)
    }
  )
}

# =============================================================================
# Run App
# =============================================================================
shinyApp(ui = ui, server = server)
