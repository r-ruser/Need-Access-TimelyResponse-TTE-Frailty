suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(patchwork); library(scales); library(grid)
})

ink <- "#222222"
blue <- "#0072B2"
orange <- "#D55E00"
green <- "#009E73"
purple <- "#7B4AB5"
grey <- "#8A8A8A"
light_grey <- "#E6E6E6"
wrap_text <- function(x, width = 120) paste(strwrap(x, width = width), collapse = "\n")

theme_nature <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      plot.title = element_text(size = base_size + 1, face = "bold", colour = ink, hjust = 0),
      plot.subtitle = element_text(size = base_size, colour = "#4D4D4D", hjust = 0, margin = margin(b = 5)),
      plot.caption = element_text(size = base_size - 1, colour = "#555555", hjust = 0, lineheight = 1.05),
      axis.title = element_text(size = base_size, colour = ink),
      axis.text = element_text(size = base_size - 0.5, colour = ink),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold", colour = ink),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.height = unit(3.5, "mm"),
      panel.spacing = unit(5, "mm"),
      plot.margin = margin(5, 7, 5, 5)
    )
}

save_nature <- function(plot, stem, width = 7.2, height = 4.8) {
  svg <- file.path(paths$figures, paste0(stem, ".svg"))
  pdf <- file.path(paths$figures, paste0(stem, ".pdf"))
  tif <- file.path(paths$figures, paste0(stem, ".tiff"))
  png <- file.path(paths$figures, paste0(stem, ".png"))
  svglite::svglite(svg, width = width, height = height, bg = "white"); print(plot); dev.off()
  grDevices::cairo_pdf(pdf, width = width, height = height, family = "Arial", bg = "white"); print(plot); dev.off()
  ragg::agg_tiff(tif, width = width, height = height, units = "in", res = 600,
                 compression = "lzw", background = "white"); print(plot); dev.off()
  ragg::agg_png(png, width = width, height = height, units = "in", res = 300,
                background = "white"); print(plot); dev.off()
  invisible(c(svg, pdf, tif, png))
}

table_path <- function(x) file.path(paths$tables, x)

# Figure 1: HRS 2 x 2 construct decomposition.
decomp <- read_csv(table_path("Table_construct_decomposition.csv"), show_col_types = FALSE) %>%
  mutate(
    construct = factor(
      paste0(analysis, ". ", need_definition, " / ", outcome_definition),
      levels = rev(paste0(analysis, ". ", need_definition, " / ", outcome_definition))
    ),
    rd = 100 * rd_early_persistent, low = 100 * rd_ep_low, high = 100 * rd_ep_high,
    outcome_group = ifelse(outcome_definition == "H-mPFI-28", "H-mPFI-28", "common-FI20")
  )
shifts <- read_csv(table_path("Table_construct_shifts.csv"), show_col_types = FALSE) %>%
  mutate(
    shift_short = case_when(
      shift == "B - A" ~ "B - A\nNeed shift / H-mPFI-28",
      shift == "D - C" ~ "D - C\nNeed shift / common-FI20",
      shift == "C - A" ~ "C - A\nOutcome shift / PA+SMK+ALC",
      TRUE ~ "D - B\nOutcome shift / PA+SMK"
    ),
    shift_label = factor(shift_short, levels = rev(shift_short)),
    rd = 100 * estimate, low = 100 * conf_low, high = 100 * conf_high,
    shift_type = ifelse(grepl("Need", interpretation), "Need definition", "Outcome definition")
  )
write_csv(decomp, file.path(paths$source, "Figure_construct_decomposition_estimates.csv"))
write_csv(shifts, file.path(paths$source, "Figure_construct_decomposition_shifts.csv"))

p1a <- ggplot(decomp, aes(rd, construct, colour = outcome_group)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, colour = grey) +
  geom_errorbar(aes(xmin = low, xmax = high), orientation = "y", width = 0, linewidth = 0.55) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c("H-mPFI-28" = blue, "common-FI20" = orange), name = "Outcome") +
  scale_x_continuous(labels = label_number(accuracy = 1)) +
  labs(title = "a  Construct-specific estimates", x = "Early vs Persistent risk difference (percentage points)", y = NULL) +
  theme_nature() + theme(legend.position = "none")

p1b <- ggplot(shifts, aes(rd, shift_label, colour = shift_type)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, colour = grey) +
  geom_errorbar(aes(xmin = low, xmax = high), orientation = "y", width = 0, linewidth = 0.55) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c("Need definition" = green, "Outcome definition" = purple), name = "Shift") +
  labs(title = "b  Estimate shifts", x = "Difference in risk difference (percentage points)", y = NULL) +
  theme_nature() + theme(legend.position = "none")

fig_construct <- p1a + p1b + plot_layout(widths = c(1, 1.15)) +
  plot_annotation(
    title = "HRS construct decomposition",
    subtitle = "Outcome definition, rather than included Need domains, accounted for most of the directional shift",
    caption = wrap_text("Points show estimates; bars show 95% confidence intervals. Shift intervals use 300 person-level bootstrap samples. Shifts are sensitivity contrasts, not causal interactions.", 125),
    theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = 10),
                  plot.subtitle = element_text(family = "Arial", size = 8),
                  plot.caption = element_text(family = "Arial", size = 7, hjust = 0))
  )
save_nature(fig_construct, "Figure_construct_decomposition", 7.2, 4.7)

# Figure 2: Need-component specification forest plot.
need <- read_csv(table_path("Table_need_components.csv"), show_col_types = FALSE) %>%
  filter(population == "Primary") %>%
  mutate(
    definition = factor(need_definition, levels = rev(c("PA", "SMK", "ALC", "PA+SMK", "PA+ALC", "SMK+ALC", "PA+SMK+ALC"))),
    rd = 100 * rd_early_persistent, low = 100 * rd_low, high = 100 * rd_high,
    all_gates = support_pass & balance_pass & positivity_pass,
    gate_label = ifelse(all_gates, "All gates passed", "At least one gate failed")
  )
write_csv(need, file.path(paths$source, "Figure_need_components_source.csv"))
p_need <- ggplot(need, aes(rd, definition)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, colour = grey) +
  geom_errorbar(aes(xmin = low, xmax = high), orientation = "y", width = 0, linewidth = 0.6, colour = ink) +
  geom_point(aes(fill = gate_label), shape = 21, size = 2.6, stroke = 0.65, colour = ink) +
  scale_fill_manual(values = c("All gates passed" = blue, "At least one gate failed" = "white"), name = NULL) +
  labs(
    title = "Sensitivity to actionable-Need specification",
    subtitle = "HRS primary population; Early versus Persistent strategy",
    x = "Risk difference (percentage points)", y = "Need definition",
    caption = wrap_text("Bars show 95% clustered robust confidence intervals. Open points indicate insufficient strategy support or balance; estimates are shown for transparency and not treated as primary causal results.", 90)
  ) + theme_nature() + theme(legend.position = "top")
save_nature(p_need, "Figure_need_components", 5.2, 4.2)

# Figure 3: Need x Access matrix.
matrix_data <- read_csv(table_path("Table_need_access_matrix.csv"), show_col_types = FALSE) %>%
  mutate(
    need_burden = factor(need_burden, levels = c(">=3 needs", "2 needs", "1 need")),
    access_level = factor(access_level, levels = c("Limited access", "Favorable access")),
    cell_label = paste0("n=", n_people, "\nE ", number(early_pct, accuracy = 0.1),
                        "%\nO ", number(frailty_or_death_pct, accuracy = 0.1), "%")
  )
write_csv(matrix_data, file.path(paths$source, "Figure_need_access_heatmap_source.csv"))
p_heat <- ggplot(matrix_data, aes(access_level, need_burden, fill = early_pct)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = cell_label), family = "Arial", size = 1.7, lineheight = 0.9, colour = ink) +
  facet_wrap(~access_indicator, nrow = 1) +
  scale_fill_gradient(low = "#F2F2F2", high = blue, na.value = "white", name = "Early\nresponse (%)") +
  labs(
    title = "Need burden and Access context",
    subtitle = "HRS weighted policy-management matrix",
    x = NULL, y = "Baseline Need burden",
    caption = wrap_text("Cell labels show unique people, weighted Early-response percentage (E) and weighted frailty-or-death percentage (O). Blank cells indicate absent or unsupported Access categories.", 125)
  ) + theme_nature() +
  theme(axis.text.x = element_text(angle = 28, hjust = 1), legend.position = "right")
save_nature(p_heat, "Figure_need_access_heatmap", 7.2, 3.7)

# Figure 4: standardized Access differences in Early response.
access <- read_csv(table_path("Table_access_equity.csv"), show_col_types = FALSE) %>%
  filter(tier == "Primary") %>%
  mutate(
    cohort = factor(cohort, levels = c("HRS", "CLHLS")),
    indicator = recode(exposure, insured = "Insurance", wealth_high = "Wealth", urban = "Residence", living_alone = "Living arrangement"),
    pd_pp = 100 * pd, low_pp = 100 * pd_low, high_pp = 100 * pd_high,
    comparison = paste0(level1, " minus ", level0)
  )
write_csv(access, file.path(paths$source, "Figure_access_response_source.csv"))
access_supported <- access %>% filter(support)
access_missing <- access %>% filter(!support) %>%
  mutate(pd_pp = -4.8, label = "Not estimable")
p_access <- ggplot(access_supported, aes(pd_pp, indicator, colour = cohort)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, colour = grey) +
  geom_errorbar(aes(xmin = low_pp, xmax = high_pp), orientation = "y", width = 0, linewidth = 0.55,
                 position = position_dodge(width = 0.42)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.42)) +
  geom_text(data = access_missing, aes(x = pd_pp, y = indicator, label = label),
            inherit.aes = FALSE, family = "Arial", size = 2.3, colour = grey, hjust = 0) +
  facet_wrap(~cohort, nrow = 1) +
  scale_colour_manual(values = c(HRS = blue, CLHLS = orange), guide = "none") +
  coord_cartesian(xlim = c(-5, 9)) +
  labs(
    title = "Need-standardized probability of an Early response",
    subtitle = "Absolute difference for level 1 versus reference Access context",
    x = "Probability difference (percentage points)", y = NULL,
    caption = wrap_text("Bars show 95% person-bootstrap intervals. Wealth was not harmonized in CLHLS; insurance lacked variation in the HRS analytic sample. Cohort contrasts are descriptive.", 115)
  ) + theme_nature()
save_nature(p_access, "Figure_access_response", 6.6, 3.6)

# Figure 5: HRS policy-coverage scenarios.
policy <- read_csv(table_path("Table_policy_scenarios.csv"), show_col_types = FALSE) %>%
  mutate(scenario = factor(scenario, levels = c("Observed", "50%", "70%", "90%")))
write_csv(policy, file.path(paths$source, "Figure_policy_coverage_source.csv"))
p5a <- ggplot(policy, aes(scenario, 100 * expected_risk, group = 1)) +
  geom_errorbar(aes(ymin = 100 * conf_low, ymax = 100 * conf_high),
                width = 0.12, linewidth = 0.5, colour = blue) +
  geom_line(colour = blue, linewidth = 0.65) + geom_point(colour = blue, size = 2.1) +
  labs(title = "a  Expected frailty-or-death risk", x = "Early-response coverage", y = "Expected risk (%)") +
  theme_nature()
p5b_dat <- bind_rows(
  policy %>% transmute(scenario, outcome = "Frailty or death", estimate = event_difference_per_1000,
                       low = event_diff_low, high = event_diff_high),
  policy %>% transmute(scenario, outcome = "Hospitalization", estimate = hospitalization_difference_per_1000,
                       low = hosp_diff_low, high = hosp_diff_high)
)
p5b <- ggplot(p5b_dat, aes(scenario, estimate, colour = outcome, group = outcome)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35, colour = grey) +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.12, linewidth = 0.5,
                position = position_dodge(width = 0.22)) +
  geom_line(position = position_dodge(width = 0.22), linewidth = 0.55) +
  geom_point(position = position_dodge(width = 0.22), size = 2) +
  scale_colour_manual(values = c("Frailty or death" = blue, "Hospitalization" = orange), name = NULL) +
  labs(title = "b  Difference from observed coverage", x = "Early-response coverage", y = "Events per 1,000") +
  theme_nature() + theme(legend.position = "top")
fig_policy <- p5a + p5b + plot_annotation(
  title = "Observational policy-impact scenarios",
  subtitle = "HRS Need-eligible population",
    caption = wrap_text("Intervals use 500 person-level bootstrap samples. Scenarios reallocate the observed non-Early mixture and do not represent a policy implementation trial.", 125),
  theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = 10),
                plot.subtitle = element_text(family = "Arial", size = 8),
                plot.caption = element_text(family = "Arial", size = 7, hjust = 0))
)
save_nature(fig_policy, "Figure_policy_coverage", 7.2, 4.2)

# Figure 6: construct pathway and common-construct external transport.
transport <- read_csv(table_path("Table_cross_cohort_transportability.csv"), show_col_types = FALSE) %>%
  mutate(rd = 100 * rd_early_persistent, low = 100 * rd_low, high = 100 * rd_high)
pathway <- bind_rows(
  decomp %>% filter(analysis %in% c("A", "B", "D")) %>%
    transmute(stage = recode(analysis, A = "HRS full", B = "Need harmonization", D = "Outcome harmonization"),
              detail = paste0(need_definition, "\n", outcome_definition), rd, low, high, order = c(1, 2, 3)),
  transport %>% filter(cohort == "CLHLS") %>%
    transmute(stage = "External transport", detail = "PA+SMK\ncommon-FI20", rd, low, high, order = 4)
) %>% arrange(order) %>% mutate(stage = factor(stage, levels = stage))
write_csv(pathway, file.path(paths$source, "Figure_cross_cohort_transportability_pathway.csv"))
write_csv(transport, file.path(paths$source, "Figure_cross_cohort_transportability_common.csv"))

p6a <- ggplot(pathway, aes(order, rd)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35, colour = grey) +
  geom_line(colour = grey, linewidth = 0.6) +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.08, linewidth = 0.5, colour = ink) +
  geom_point(aes(fill = stage), shape = 21, size = 3, stroke = 0.6, colour = ink, show.legend = FALSE) +
  geom_text(aes(label = detail), y = min(pathway$low) - 4.8, family = "Arial", size = 2.25, lineheight = 0.95) +
  scale_fill_manual(values = c(blue, "#56B4E9", green, orange)) +
  scale_x_continuous(breaks = pathway$order, labels = pathway$stage, expand = expansion(mult = c(.08, .08))) +
  coord_cartesian(ylim = c(min(pathway$low) - 8, max(pathway$high) + 1.5), clip = "off") +
  labs(title = "a  Construct transport pathway", x = NULL, y = "Risk difference (percentage points)") +
  theme_nature() + theme(axis.text.x = element_text(face = "bold"))

p6b <- ggplot(transport, aes(rd, factor(cohort, levels = rev(c("HRS", "CLHLS"))), colour = cohort)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, colour = grey) +
  geom_errorbar(aes(xmin = low, xmax = high), orientation = "y", width = 0, linewidth = 0.6) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c(HRS = blue, CLHLS = orange), guide = "none") +
  labs(title = "b  Common-construct estimates", subtitle = "PA+SMK / common-FI20",
       x = "Early vs Persistent risk difference (percentage points)", y = NULL) +
  theme_nature()

fig_transport <- p6a / p6b + plot_layout(heights = c(1.35, 1)) +
  plot_annotation(
    title = "Transportability of policy conclusions",
    subtitle = "Specification changes within HRS precede external replication with the common construct",
    caption = wrap_text("Points show risk differences; bars show 95% confidence intervals. The HRS-CLHLS comparison is descriptive heterogeneity, not a causal health-system effect.", 125),
    theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = 10),
                  plot.subtitle = element_text(family = "Arial", size = 8),
                  plot.caption = element_text(family = "Arial", size = 7, hjust = 0))
  )
save_nature(fig_transport, "Figure_cross_cohort_transportability", 7.2, 6.0)

cat("Nature-style figures exported in SVG, PDF, TIFF (600 dpi), and PNG.\n")
