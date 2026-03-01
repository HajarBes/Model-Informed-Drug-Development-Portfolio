# ==============================================================================
# Project 06: Clinical Dataset Engineering for Pharmacometrics
# Script: 07_generate_figures.R - Generate publication-quality PNG figures
# ==============================================================================

source("R/00_setup.R")

cat("=== Generating Portfolio Figures ===\n\n")

fig_dir <- file.path(proj_root, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Load datasets
nm_data <- read_csv(file.path(paths$outputs, "nm_dataset.csv"),
                    show_col_types = FALSE, na = ".")
qc_results <- read_csv(file.path(paths$outputs, "qc_checklist.csv"),
                       show_col_types = FALSE)

obs_data <- nm_data %>% filter(EVID == 0)
dose_data <- nm_data %>% filter(EVID == 1)
subj_data <- nm_data %>% distinct(ID, .keep_all = TRUE)

# Theme
theme_portfolio <- theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40", size = 11),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# ---- Figure 1: PK Spaghetti Plot (Hero) ----
cat("  Fig 1: PK concentration profiles...\n")

pk_plot <- obs_data %>%
  filter(BLQ == 0, TAD >= 0, TAD <= 50)

p1 <- ggplot(pk_plot, aes(x = TAD, y = DV, group = factor(ID), color = ARM)) +
  geom_line(alpha = 0.25, linewidth = 0.3) +
  geom_point(alpha = 0.25, size = 0.6) +
  scale_y_log10(labels = scales::label_number()) +
  labs(
    title = "Individual PK Concentration-Time Profiles",
    subtitle = "Quantifiable observations, TAD 0-50h | CDISCPILOT01 Xanomeline",
    x = "Time After Dose (hours)",
    y = "Concentration (ug/mL, log scale)",
    color = "Treatment Arm"
  ) +
  theme_portfolio +
  scale_color_brewer(palette = "Set1")

ggsave(file.path(fig_dir, "pk_concentration_profiles.png"), p1,
       width = 10, height = 6, dpi = 150, bg = "white")

# ---- Figure 2: BLQ by Time After Dose ----
cat("  Fig 2: BLQ analysis by time...\n")

blq_time <- obs_data %>%
  mutate(time_bin = cut(TAD,
    breaks = c(-Inf, 0, 2, 6, 12, 24, 48, Inf),
    labels = c("Pre-dose", "0-2h", "2-6h", "6-12h", "12-24h", "24-48h", ">48h"))) %>%
  group_by(time_bin) %>%
  summarise(
    pct_blq = round(100 * mean(BLQ == 1, na.rm = TRUE), 1),
    n = n(),
    .groups = "drop"
  )

p2 <- ggplot(blq_time, aes(x = time_bin, y = pct_blq)) +
  geom_col(fill = "#2E75B6", alpha = 0.85) +
  geom_text(aes(label = paste0("n=", n)), vjust = -0.5, size = 3.5) +
  labs(
    title = "BLQ Percentage by Time After Dose",
    subtitle = paste0("LLOQ = ", config$lloq, " ", config$conc_unit,
                      " | BLQ rule: ", config$blq_rule),
    x = "", y = "% Below LLOQ"
  ) +
  theme_portfolio +
  ylim(0, max(blq_time$pct_blq) * 1.15)

ggsave(file.path(fig_dir, "blq_by_time.png"), p2,
       width = 8, height = 5, dpi = 150, bg = "white")

# ---- Figure 3: QC Checklist Summary ----
cat("  Fig 3: QC checklist summary...\n")

qc_colors <- c("PASS" = "#27AE60", "FAIL" = "#E74C3C",
                "WARNING" = "#F39C12", "INFO" = "#3498DB")

p3 <- ggplot(qc_results, aes(x = reorder(check_name, rev(seq_len(nrow(qc_results)))),
                               y = 1, fill = status)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = status), size = 3.5, fontface = "bold") +
  scale_fill_manual(values = qc_colors) +
  coord_flip() +
  labs(
    title = "QC Checklist Results",
    subtitle = sprintf("%d checks | %d PASS | %d FAIL | %d INFO",
                       nrow(qc_results),
                       sum(qc_results$status == "PASS"),
                       sum(qc_results$status == "FAIL"),
                       sum(qc_results$status == "INFO")),
    x = "", y = "", fill = "Status"
  ) +
  theme_portfolio +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

ggsave(file.path(fig_dir, "qc_checklist_summary.png"), p3,
       width = 9, height = 6, dpi = 150, bg = "white")

# ---- Figure 4: Covariate Distributions ----
cat("  Fig 4: Covariate distributions...\n")

cov_vars <- c("AGE", "WT")
available_covs <- intersect(cov_vars, names(subj_data))

cov_long <- subj_data %>%
  select(ID, ARM, all_of(available_covs)) %>%
  pivot_longer(all_of(available_covs), names_to = "covariate", values_to = "value") %>%
  filter(!is.na(value))

p4 <- ggplot(cov_long, aes(x = ARM, y = value, fill = ARM)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  facet_wrap(~covariate, scales = "free_y") +
  labs(
    title = "Baseline Covariate Distributions by Treatment Arm",
    subtitle = sprintf("N = %d subjects", n_distinct(subj_data$ID)),
    x = "", y = ""
  ) +
  theme_portfolio +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 15, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")

ggsave(file.path(fig_dir, "covariate_distributions.png"), p4,
       width = 9, height = 5, dpi = 150, bg = "white")

# ---- Figure 5: Dataset Composition Waterfall ----
cat("  Fig 5: Dataset composition...\n")

composition <- tibble(
  category = c("Total Subjects\n(DM)", "Active Arms\nOnly", "With Dose\nRecords",
               "Dosing\nEvents", "PK\nObservations", "BLQ\nRecords"),
  count = c(254, 168, 168,
            nrow(dose_data), nrow(obs_data), sum(obs_data$BLQ == 1)),
  type = c("total", "filter", "filter", "component", "component", "component")
)

p5 <- ggplot(composition, aes(x = fct_inorder(category), y = count, fill = type)) +
  geom_col(alpha = 0.85, width = 0.7) +
  geom_text(aes(label = scales::comma(count)), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("total" = "#7F8C8D", "filter" = "#2E75B6",
                                "component" = "#27AE60")) +
  labs(
    title = "Dataset Composition: SDTM to NONMEM",
    subtitle = "Subject filtering and record breakdown",
    x = "", y = "Count"
  ) +
  theme_portfolio +
  theme(legend.position = "none") +
  ylim(0, max(composition$count) * 1.15)

ggsave(file.path(fig_dir, "dataset_composition.png"), p5,
       width = 10, height = 5, dpi = 150, bg = "white")

# ---- Figure 6: Missingness Heatmap ----
cat("  Fig 6: Missingness heatmap...\n")

miss_summary <- nm_data %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(pct_missing = round(100 * n_missing / nrow(nm_data), 1))

p6 <- ggplot(miss_summary, aes(x = reorder(variable, pct_missing), y = pct_missing)) +
  geom_col(fill = ifelse(miss_summary$pct_missing > 0, "#E74C3C", "#27AE60"),
           alpha = 0.8) +
  coord_flip() +
  labs(
    title = "Missingness by Variable",
    subtitle = sprintf("Total records: %s | Green = complete, Red = missing",
                       scales::comma(nrow(nm_data))),
    x = "", y = "% Missing"
  ) +
  theme_portfolio

ggsave(file.path(fig_dir, "missingness_heatmap.png"), p6,
       width = 8, height = 7, dpi = 150, bg = "white")

# ---- Figure 7: Pipeline Overview (text-based) ----
cat("  Fig 7: Pipeline overview...\n")

pipeline_data <- tibble(
  step = factor(1:6),
  label = c("Extract\nSDTM", "Build\nNONMEM", "Build\nADNCA",
            "QC\nChecks", "QC\nReport", "Submission\nPackage"),
  script = c("01_extract", "02_build_nm", "03_build_adnca",
             "04_qc_nm", "05_render_qc", "06_assemble"),
  y = rep(1, 6)
)

p7 <- ggplot(pipeline_data, aes(x = step, y = y)) +
  geom_tile(fill = "#2E75B6", alpha = 0.15, width = 0.85, height = 0.6) +
  geom_text(aes(label = label), size = 4, fontface = "bold", color = "#2E75B6") +
  geom_text(aes(label = script, y = 0.6), size = 2.8, color = "grey50") +
  # Arrows between steps
  annotate("segment", x = 1.45, xend = 1.55, y = 1, yend = 1,
           arrow = arrow(length = unit(0.15, "cm")), color = "grey40") +
  annotate("segment", x = 2.45, xend = 2.55, y = 1, yend = 1,
           arrow = arrow(length = unit(0.15, "cm")), color = "grey40") +
  annotate("segment", x = 3.45, xend = 3.55, y = 1, yend = 1,
           arrow = arrow(length = unit(0.15, "cm")), color = "grey40") +
  annotate("segment", x = 4.45, xend = 4.55, y = 1, yend = 1,
           arrow = arrow(length = unit(0.15, "cm")), color = "grey40") +
  annotate("segment", x = 5.45, xend = 5.55, y = 1, yend = 1,
           arrow = arrow(length = unit(0.15, "cm")), color = "grey40") +
  labs(title = "Pipeline: SDTM to Submission Package",
       subtitle = "6 modular R scripts | Config-driven | Fully reproducible") +
  theme_void(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(color = "grey40", size = 11, hjust = 0.5)
  ) +
  ylim(0.4, 1.4)

ggsave(file.path(fig_dir, "pipeline_overview.png"), p7,
       width = 12, height = 3.5, dpi = 150, bg = "white")

# ---- Done ----
cat(sprintf("\n=== %d figures saved to %s ===\n",
            length(list.files(fig_dir, pattern = "\\.png$")), fig_dir))
