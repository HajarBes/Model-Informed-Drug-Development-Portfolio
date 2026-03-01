# ==============================================================================
# Project 06: QC Dashboard - Interactive Dataset Inspection
# Shiny Application
# ==============================================================================

library(shiny)
library(tidyverse)
library(DT)
library(plotly)

# ---- Load data ----
proj_root <- dirname(getwd())
if (!file.exists(file.path(proj_root, "outputs", "nm_dataset.csv"))) {
  proj_root <- getwd()  # fallback
}

nm_data <- read_csv(file.path(proj_root, "outputs", "nm_dataset.csv"),
                    show_col_types = FALSE, na = ".")
qc_results <- read_csv(file.path(proj_root, "outputs", "qc_checklist.csv"),
                        show_col_types = FALSE)

obs_data <- nm_data %>% filter(EVID == 0)
dose_data <- nm_data %>% filter(EVID == 1)
subj_data <- nm_data %>% distinct(ID, .keep_all = TRUE)

# ---- UI ----
ui <- fluidPage(
  titlePanel("QC Dashboard - NONMEM PK Dataset (CDISCPILOT01)"),

  tabsetPanel(
    # Tab 1: QC Checklist
    tabPanel("QC Checklist",
      br(),
      h4("Automated QC Check Results"),
      DTOutput("qc_table"),
      br(),
      verbatimTextOutput("qc_summary")
    ),

    # Tab 2: Missingness
    tabPanel("Missingness",
      br(),
      h4("Missing Values by Variable"),
      plotlyOutput("missing_plot", height = "500px"),
      br(),
      h4("Missing Values by Subject"),
      plotlyOutput("missing_by_subj", height = "400px")
    ),

    # Tab 3: Subject Timelines
    tabPanel("Subject Timelines",
      br(),
      fluidRow(
        column(4, selectInput("subj_id", "Select Subject ID:",
                              choices = sort(unique(nm_data$ID)),
                              selected = sort(unique(nm_data$ID))[1])),
        column(4, checkboxInput("show_blq", "Highlight BLQ", TRUE))
      ),
      plotlyOutput("timeline_plot", height = "500px"),
      br(),
      DTOutput("subj_table")
    ),

    # Tab 4: BLQ Analysis
    tabPanel("BLQ Analysis",
      br(),
      fluidRow(
        column(6, plotlyOutput("blq_by_time", height = "400px")),
        column(6, plotlyOutput("blq_by_arm", height = "400px"))
      ),
      br(),
      h4("BLQ Summary"),
      verbatimTextOutput("blq_summary")
    ),

    # Tab 5: Covariate Distributions
    tabPanel("Covariates",
      br(),
      fluidRow(
        column(4, selectInput("cov_var", "Select Covariate:",
                              choices = c("AGE", "WT", "CREAT", "ALT", "AST", "BILI"),
                              selected = "AGE"))
      ),
      fluidRow(
        column(6, plotlyOutput("cov_hist", height = "400px")),
        column(6, plotlyOutput("cov_by_arm", height = "400px"))
      ),
      br(),
      verbatimTextOutput("cov_summary")
    ),

    # Tab 6: Concentration-Time
    tabPanel("PK Profiles",
      br(),
      fluidRow(
        column(4, selectInput("pk_time", "Time Variable:",
                              choices = c("TIME", "TAD"), selected = "TAD")),
        column(4, checkboxInput("pk_log", "Log Scale (DV)", TRUE)),
        column(4, sliderInput("pk_xlim", "Time Limit (hours):",
                              min = 0, max = 200, value = 50))
      ),
      plotlyOutput("pk_spaghetti", height = "500px")
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {

  # Tab 1: QC Checklist
  output$qc_table <- renderDT({
    qc_results %>%
      mutate(
        status = case_when(
          status == "PASS" ~ "\U2705 PASS",
          status == "FAIL" ~ "\U274C FAIL",
          status == "WARNING" ~ "\U26A0\UFE0F WARNING",
          status == "INFO" ~ "\U2139\UFE0F INFO",
          TRUE ~ status
        )
      ) %>%
      datatable(options = list(pageLength = 25, dom = 't'),
                rownames = FALSE, escape = FALSE)
  })

  output$qc_summary <- renderText({
    n_pass <- sum(qc_results$status == "PASS")
    n_fail <- sum(qc_results$status == "FAIL")
    n_warn <- sum(qc_results$status == "WARNING")
    n_info <- sum(qc_results$status == "INFO")
    overall <- if_else(n_fail == 0, "QC PASS", "QC FAIL")
    sprintf("Overall: %s | Checks: %d PASS, %d FAIL, %d WARNING, %d INFO",
            overall, n_pass, n_fail, n_warn, n_info)
  })

  # Tab 2: Missingness
  output$missing_plot <- renderPlotly({
    miss_summary <- nm_data %>%
      summarise(across(everything(), ~sum(is.na(.)))) %>%
      pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
      mutate(pct_missing = round(100 * n_missing / nrow(nm_data), 1)) %>%
      filter(n_missing > 0) %>%
      arrange(desc(n_missing))

    if (nrow(miss_summary) == 0) {
      plot_ly() %>% layout(title = "No missing values in dataset")
    } else {
      plot_ly(miss_summary, x = ~reorder(variable, n_missing), y = ~pct_missing,
              type = "bar", text = ~paste(n_missing, "missing"),
              marker = list(color = "#2E75B6")) %>%
        layout(title = "Missing Values by Variable",
               xaxis = list(title = ""), yaxis = list(title = "% Missing"),
               margin = list(b = 100))
    }
  })

  output$missing_by_subj <- renderPlotly({
    miss_per_subj <- nm_data %>%
      group_by(ID) %>%
      summarise(n_missing = sum(is.na(DV) & EVID == 0), .groups = "drop")

    plot_ly(miss_per_subj, x = ~ID, y = ~n_missing, type = "bar",
            marker = list(color = "#E67E22")) %>%
      layout(title = "Missing DV per Subject (observations only)",
             xaxis = list(title = "Subject ID"), yaxis = list(title = "N Missing"))
  })

  # Tab 3: Subject Timelines
  output$timeline_plot <- renderPlotly({
    sid <- as.integer(input$subj_id)
    subj <- nm_data %>% filter(ID == sid)

    p <- plot_ly()

    # Dose events
    doses_subj <- subj %>% filter(EVID == 1)
    if (nrow(doses_subj) > 0) {
      p <- p %>% add_markers(data = doses_subj, x = ~TIME, y = ~AMT,
                              name = "Dose", marker = list(color = "red", size = 12, symbol = "triangle-up"))
    }

    # Observations
    obs_subj <- subj %>% filter(EVID == 0)
    if (input$show_blq) {
      obs_q <- obs_subj %>% filter(BLQ == 0)
      obs_b <- obs_subj %>% filter(BLQ == 1)
      p <- p %>%
        add_markers(data = obs_q, x = ~TIME, y = ~DV, name = "Quantifiable",
                    marker = list(color = "#2E75B6", size = 8)) %>%
        add_markers(data = obs_b, x = ~TIME, y = ~DV, name = "BLQ",
                    marker = list(color = "grey", size = 8, symbol = "x"))
    } else {
      p <- p %>% add_markers(data = obs_subj, x = ~TIME, y = ~DV, name = "Observation",
                              marker = list(color = "#2E75B6", size = 8))
    }

    p %>% layout(title = sprintf("Subject %d - Dosing & PK Timeline", sid),
                 xaxis = list(title = "Time (hours)"),
                 yaxis = list(title = "Concentration (ug/mL) / Dose (mg)"))
  })

  output$subj_table <- renderDT({
    sid <- as.integer(input$subj_id)
    nm_data %>%
      filter(ID == sid) %>%
      select(TIME, TAD, AMT, DV, EVID, CMT, BLQ, MDV) %>%
      datatable(options = list(pageLength = 20, dom = 'tp'), rownames = FALSE)
  })

  # Tab 4: BLQ Analysis
  output$blq_by_time <- renderPlotly({
    blq_time <- obs_data %>%
      mutate(time_bin = cut(TAD, breaks = c(-Inf, 0, 2, 6, 12, 24, 48, Inf),
                            labels = c("Pre-dose", "0-2h", "2-6h", "6-12h", "12-24h", "24-48h", ">48h"))) %>%
      group_by(time_bin) %>%
      summarise(pct_blq = round(100 * mean(BLQ == 1, na.rm = TRUE), 1),
                n = n(), .groups = "drop")

    plot_ly(blq_time, x = ~time_bin, y = ~pct_blq, type = "bar",
            text = ~paste("n =", n), marker = list(color = "#95A5A6")) %>%
      layout(title = "% BLQ by Time After Dose",
             xaxis = list(title = ""), yaxis = list(title = "% BLQ"))
  })

  output$blq_by_arm <- renderPlotly({
    blq_arm <- obs_data %>%
      group_by(ARM) %>%
      summarise(pct_blq = round(100 * mean(BLQ == 1, na.rm = TRUE), 1),
                n = n(), .groups = "drop")

    plot_ly(blq_arm, x = ~ARM, y = ~pct_blq, type = "bar",
            text = ~paste("n =", n), marker = list(color = "#3498DB")) %>%
      layout(title = "% BLQ by Treatment Arm",
             xaxis = list(title = ""), yaxis = list(title = "% BLQ"))
  })

  output$blq_summary <- renderText({
    n_total <- nrow(obs_data)
    n_blq <- sum(obs_data$BLQ == 1, na.rm = TRUE)
    sprintf("Total observations: %d\nBLQ observations: %d (%.1f%%)\nQuantifiable: %d (%.1f%%)",
            n_total, n_blq, 100*n_blq/n_total, n_total-n_blq, 100*(n_total-n_blq)/n_total)
  })

  # Tab 5: Covariates
  output$cov_hist <- renderPlotly({
    var <- input$cov_var
    vals <- subj_data[[var]]
    if (all(is.na(vals))) {
      plot_ly() %>% layout(title = paste(var, "- all values missing"))
    } else {
      plot_ly(x = ~vals[!is.na(vals)], type = "histogram",
              marker = list(color = "#2E75B6")) %>%
        layout(title = paste(var, "Distribution (baseline)"),
               xaxis = list(title = var), yaxis = list(title = "Count"))
    }
  })

  output$cov_by_arm <- renderPlotly({
    var <- input$cov_var
    if (all(is.na(subj_data[[var]]))) {
      plot_ly() %>% layout(title = paste(var, "- all values missing"))
    } else {
      plot_ly(subj_data, x = ~ARM, y = ~get(var), type = "box",
              color = ~ARM) %>%
        layout(title = paste(var, "by Treatment Arm"),
               yaxis = list(title = var), showlegend = FALSE)
    }
  })

  output$cov_summary <- renderText({
    var <- input$cov_var
    vals <- subj_data[[var]]
    if (all(is.na(vals))) {
      sprintf("%s: All values missing (not available in source data)", var)
    } else {
      sprintf("%s: N=%d, Mean=%.1f, SD=%.1f, Median=%.1f, Range=[%.1f-%.1f], Missing=%d",
              var, sum(!is.na(vals)), mean(vals, na.rm=TRUE), sd(vals, na.rm=TRUE),
              median(vals, na.rm=TRUE), min(vals, na.rm=TRUE), max(vals, na.rm=TRUE),
              sum(is.na(vals)))
    }
  })

  # Tab 6: PK Profiles
  output$pk_spaghetti <- renderPlotly({
    time_var <- input$pk_time
    obs_plot <- obs_data %>%
      filter(BLQ == 0, .data[[time_var]] <= input$pk_xlim, .data[[time_var]] >= 0)

    p <- plot_ly(obs_plot, x = ~get(time_var), y = ~DV, color = ~factor(ID),
                 type = "scatter", mode = "lines+markers",
                 line = list(width = 0.5), marker = list(size = 3),
                 showlegend = FALSE, opacity = 0.4) %>%
      layout(title = "Individual PK Profiles (quantifiable only)",
             xaxis = list(title = paste(time_var, "(hours)")),
             yaxis = list(title = "Concentration (ug/mL)"))

    if (input$pk_log) {
      p <- p %>% layout(yaxis = list(type = "log", title = "Concentration (ug/mL, log)"))
    }
    p
  })
}

shinyApp(ui, server)
