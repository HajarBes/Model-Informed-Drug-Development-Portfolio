# =============================================================================
# 03_make_figures.R — Generate all publication-quality figures
# =============================================================================

script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(script_dir) || !file.exists(file.path(script_dir, "00_setup.R"))) {
  script_dir <- file.path(getwd(), "analysis")
}
source(file.path(script_dir, "00_setup.R"))

log_msg("=== Starting Figure Generation ===")

# Load results
ddi_res <- readRDS(file.path(DIR_PROCESSED, "ddi_results.rds"))
sens    <- readRDS(file.path(DIR_PROCESSED, "sensitivity_results.rds"))

SAVE_W <- 8
SAVE_H <- 6

# =============================================================================
# FIGURE 1: MIDAZOLAM TORNADO PLOT
# =============================================================================

tornado <- sens$tornado %>%
  mutate(label = factor(label, levels = rev(label)))

fig_tornado <- ggplot(tornado) +
  geom_segment(aes(x = aucr_low, xend = aucr_high, y = label, yend = label),
               linewidth = 6, color = "#4472C4", alpha = 0.7) +
  geom_point(aes(x = aucr_base, y = label), size = 3, color = "black") +
  geom_vline(xintercept = tornado$aucr_base[1], linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 11.2, linetype = "solid", color = COLOR_FLAG, linewidth = 0.8) +
  annotate("text", x = 11.2, y = 0.5, label = "Observed\nAUCR = 11.2",
           color = COLOR_FLAG, size = 3, hjust = -0.1, fontface = "bold") +
  scale_x_continuous(trans = "log10",
                     breaks = c(1, 2, 5, 10, 20, 50, 100, 500),
                     labels = c("1", "2", "5", "10", "20", "50", "100", "500")) +
  labs(
    title = "Tornado Sensitivity: Predicted AUCR (Midazolam + Ketoconazole)",
    subtitle = "One-at-a-time variation | Black dot = base case | Red line = observed",
    x = "Predicted AUCR (log scale)",
    y = NULL
  )

ggsave(file.path(DIR_FIGURES, "midazolam_tornado.png"), fig_tornado,
       width = SAVE_W, height = SAVE_H, dpi = 300, bg = "white")
log_msg("Saved: midazolam_tornado.png")

# =============================================================================
# FIGURE 2: MIDAZOLAM HEATMAP (Ki vs [I]max,u)
# =============================================================================

hm <- sens$heatmap %>%
  mutate(AUCR_capped = pmin(AUCR, 100))

fig_heatmap <- ggplot(hm, aes(x = Ki_uM, y = Imax_u_uM, fill = AUCR_capped)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("#2CA02C", "#FFD700", "#FF8C00", "#D62728", "#800000"),
    values = scales::rescale(c(1, 1.25, 2, 5, 100)),
    name = "Predicted\nAUCR",
    limits = c(1, 100),
    trans = "log10",
    breaks = c(1, 1.25, 2, 5, 10, 50, 100),
    labels = c("1", "1.25", "2", "5", "10", "50", "100+")
  ) +
  scale_x_log10(name = expression(K[i] ~ "(uM)"),
                breaks = c(0.001, 0.01, 0.1, 1, 10)) +
  scale_y_log10(name = expression("[I]"[max*","*u] ~ "(uM)"),
                breaks = c(0.01, 0.1, 1, 10)) +
  annotate("point", x = 0.015, y = 0.152, shape = 4, size = 4,
           color = "white", stroke = 2) +
  annotate("text", x = 0.015, y = 0.25, label = "Ketoconazole",
           color = "white", size = 3.5, fontface = "bold") +
  labs(
    title = expression("Decision Boundary: " * K[i] ~ "vs" ~ "[I]"[max*","*u]),
    subtitle = "Midazolam as victim (fm = 0.94) | X = ketoconazole base case"
  )

ggsave(file.path(DIR_FIGURES, "midazolam_heatmap_Ki_vs_Imaxu.png"), fig_heatmap,
       width = SAVE_W, height = SAVE_H, dpi = 300, bg = "white")
log_msg("Saved: midazolam_heatmap_Ki_vs_Imaxu.png")

# =============================================================================
# FIGURE 3: MIDAZOLAM PREDICTED VS OBSERVED BAND
# =============================================================================

mc <- sens$mc_samples
obs_aucr <- 11.2

fig_predobs <- ggplot(mc, aes(x = AUCR)) +
  geom_histogram(aes(fill = class), bins = 60, alpha = 0.8, color = "white",
                 linewidth = 0.2) +
  scale_fill_manual(values = COLORS_RISK, name = "Classification") +
  geom_vline(xintercept = obs_aucr, color = COLOR_FLAG, linewidth = 1.2,
             linetype = "solid") +
  annotate("text", x = obs_aucr, y = Inf, label = sprintf("Observed = %.1f", obs_aucr),
           vjust = 2, hjust = -0.1, color = COLOR_FLAG, fontface = "bold", size = 3.5) +
  geom_vline(xintercept = c(1.25, 2, 5), linetype = "dotted", color = "grey40") +
  annotate("text", x = c(1.25, 2, 5), y = rep(Inf, 3),
           label = c("Weak", "Moderate", "Strong"),
           vjust = 4, color = "grey40", size = 2.5) +
  scale_x_continuous(trans = "log10",
                     breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500),
                     labels = c("1", "2", "5", "10", "20", "50", "100", "200", "500"),
                     limits = c(1, 500)) +
  labs(
    title = "Monte Carlo: Predicted AUCR Distribution (N = 10,000)",
    subtitle = sprintf("Median = %.1f | 90%% CI: [%.1f, %.1f] | Strong in %.0f%% of samples",
                       sens$mc_summary$median, sens$mc_summary$p5,
                       sens$mc_summary$p95, sens$mc_summary$pct_strong),
    x = "Predicted AUCR (log scale)",
    y = "Count"
  )

ggsave(file.path(DIR_FIGURES, "midazolam_pred_obs_band.png"), fig_predobs,
       width = SAVE_W, height = 5, dpi = 300, bg = "white")
log_msg("Saved: midazolam_pred_obs_band.png")

# =============================================================================
# FIGURE 4: SOTORASIB MECHANISM DASHBOARD
# =============================================================================

sot <- ddi_res[["sotorasib_perpetrator"]]

# Panel A: R1 by CYP enzyme
cyp_plot_data <- sot$cyp %>%
  mutate(
    enzyme = factor(enzyme, levels = rev(enzyme)),
    color  = ifelse(R1_flag, COLOR_FLAG, COLOR_PASS)
  )

panel_a <- ggplot(cyp_plot_data, aes(x = R1, y = enzyme, fill = R1_flag)) +
  geom_col(width = 0.6) +
  geom_vline(xintercept = 1.02, linetype = "dashed", color = "grey30") +
  annotate("text", x = 1.02, y = 0.4, label = "Threshold\n(R1 = 1.02)",
           size = 2.5, hjust = -0.1, color = "grey30") +
  scale_fill_manual(values = c("FALSE" = COLOR_PASS, "TRUE" = COLOR_FLAG),
                    guide = "none") +
  scale_x_continuous(limits = c(0, max(cyp_plot_data$R1) * 1.15)) +
  labs(title = "A. CYP Reversible Inhibition (R1)", x = "R1", y = NULL)

# Panel B: R3 Induction
induct_data <- tibble(
  metric = c("R3 (CYP3A4)"),
  value  = sot$induction$R3,
  flag   = sot$induction$R3_flag
)

panel_b <- ggplot(induct_data, aes(x = value, y = metric, fill = flag)) +
  geom_col(width = 0.4) +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "grey30") +
  annotate("text", x = 0.8, y = 0.6, label = "Threshold\n(R3 = 0.8)",
           size = 2.5, hjust = 1.1, color = "grey30") +
  scale_fill_manual(values = c("FALSE" = COLOR_PASS, "TRUE" = COLOR_FLAG),
                    guide = "none") +
  scale_x_continuous(limits = c(0, 1)) +
  labs(title = "B. CYP3A4 Induction (R3)", x = "R3", y = NULL)

# Panel C: Transporter screening
trans_plot <- sot$transporter %>%
  filter(!is.na(ratio)) %>%
  mutate(
    transporter = factor(transporter, levels = rev(transporter)),
    log_ratio   = log10(pmax(ratio, 0.001)),
    log_thresh  = log10(threshold)
  )

panel_c <- ggplot(trans_plot, aes(x = ratio, y = transporter, fill = flag)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("FALSE" = COLOR_PASS, "TRUE" = COLOR_FLAG),
                    labels = c("Below threshold", "Flagged"),
                    name = "Status") +
  scale_x_continuous(trans = "log10",
                     breaks = c(0.01, 0.1, 1, 10, 100, 1000, 10000),
                     labels = c("0.01", "0.1", "1", "10", "100", "1K", "10K")) +
  labs(title = "C. Transporter Screening Ratios", x = "Ratio (log scale)", y = NULL)

fig_dashboard <- (panel_a / panel_b / panel_c) +
  plot_annotation(
    title = "Sotorasib DDI Mechanism Dashboard",
    subtitle = "Static screening: CYP inhibition, induction, and transporter risk",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

ggsave(file.path(DIR_FIGURES, "sotorasib_mechanism_dashboard.png"), fig_dashboard,
       width = SAVE_W, height = 10, dpi = 300, bg = "white")
log_msg("Saved: sotorasib_mechanism_dashboard.png")

# =============================================================================
# FIGURE 5: SOTORASIB NET EFFECT SCENARIOS
# =============================================================================

scen <- sens$scenarios %>%
  mutate(
    scenario = factor(scenario, levels = rev(scenario)),
    bar_color = case_when(
      observed ~ "#4472C4",
      AUCR > 1.25 ~ COLOR_FLAG,
      AUCR < 0.80 ~ "#4472C4",
      TRUE ~ "#FFD700"
    )
  )

fig_scenarios <- ggplot(scen, aes(x = AUCR, y = scenario, fill = direction)) +
  geom_col(width = 0.6) +
  geom_vline(xintercept = 1, linetype = "solid", color = "black", linewidth = 0.6) +
  geom_vline(xintercept = c(0.5, 0.8, 1.25, 2, 5), linetype = "dotted",
             color = "grey50") +
  annotate("text", x = c(0.5, 0.8, 1.25, 2, 5), y = rep(6.8, 5),
           label = c("Strong\nind.", "Mod.\nind.", "Weak\ninhib.", "Mod.\ninhib.", "Strong\ninhib."),
           size = 2, color = "grey40") +
  scale_fill_manual(
    values = c("Inhibition dominant" = COLOR_FLAG,
               "Induction dominant"  = "#4472C4",
               "Mixed / balanced"    = "#FFD700"),
    name = "Net Direction"
  ) +
  scale_x_continuous(
    trans = "log10",
    breaks = c(0.1, 0.25, 0.5, 1, 2, 5, 10),
    labels = c("0.1", "0.25", "0.5", "1", "2", "5", "10"),
    limits = c(0.1, 15)
  ) +
  labs(
    title = "Sotorasib: Net CYP3A4 Effect Depends on Induction Scaling",
    subtitle = "Hypothetical CYP3A4 substrate (fm=0.94) | Static model cannot resolve net direction",
    x = "Predicted AUCR (log scale) | <1 = induction | >1 = inhibition",
    y = NULL
  ) +
  theme(legend.position = "right")

ggsave(file.path(DIR_FIGURES, "sotorasib_net_effect_scenarios.png"), fig_scenarios,
       width = SAVE_W + 1, height = 5, dpi = 300, bg = "white")
log_msg("Saved: sotorasib_net_effect_scenarios.png")

log_msg("=== Figure Generation Complete ===")
log_msg("All figures saved to: ", DIR_FIGURES)
