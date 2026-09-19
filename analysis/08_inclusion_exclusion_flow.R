suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork); library(grid)
})

hrs_raw <- readRDS(input_paths$hrs_candidates) %>% filter(cohort == "HRS")
clhls_raw <- readRDS(input_paths$clhls_candidates)

fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

build_flow_counts <- function(d, cohort) {
  if (cohort == "HRS") {
    responded <- d$responded == 1
    age_ok <- responded & d$age >= 60
    dep_ok <- age_ok & d$dep_trigger == 1
    need_observed <- dep_ok & d$need_primary_observed_domains == 3
    any_need <- need_observed & d$any_need_primary == 1
    eligible <- any_need & d$h_mpfi28_n >= 23
  } else {
    responded <- d$responded == 1
    age_ok <- responded & d$age >= 60
    dep_ok <- age_ok & d$dep_trigger == 1
    need_observed <- dep_ok & d$need_observed_domains == 2
    any_need <- need_observed & d$any_need == 1
    eligible <- any_need & d$common_fi20_n >= 16
  }
  t1_ok <- eligible & !is.na(d$adequate_t1)
  t2_ok <- t1_ok & !is.na(d$adequate_t2)
  strategy_ok <- t2_ok & !is.na(d$observed_strategy)
  outcome_ok <- strategy_ok & d$outcome_observed == 1

  conditions <- list(
    candidate = rep(TRUE, nrow(d)),
    responded = responded,
    age = age_ok,
    depression = dep_ok,
    need_observed = need_observed,
    any_need = any_need,
    eligible = eligible,
    t1 = t1_ok,
    t2 = t2_ok,
    strategy = strategy_ok,
    outcome = outcome_ok
  )
  count_condition <- function(ok) {
    z <- d[ok %in% TRUE, , drop = FALSE]
    c(trials = nrow(z), people = n_distinct(z$person_id))
  }
  counts <- lapply(conditions, count_condition)
  get_n <- function(name, field = "trials") unname(counts[[name]][[field]])
  need_label <- if (identical(cohort, "HRS")) {
    "Actionable Need and valid H-mPFI-28"
  } else {
    "Actionable Need and valid common-FI20"
  }

  main <- tibble(
    cohort = cohort,
    node_type = "included",
    node = c("candidate", "age", "depression", "eligible", "strategy", "outcome"),
    n_trials = c(get_n("candidate"), get_n("age"), get_n("depression"),
                 get_n("eligible"), get_n("strategy"), get_n("outcome")),
    n_people = c(get_n("candidate", "people"), get_n("age", "people"),
                 get_n("depression", "people"), get_n("eligible", "people"),
                 get_n("strategy", "people"), get_n("outcome", "people")),
    y = c(95, 79, 63, 47, 31, 15),
    label = c(
      "Candidate person-trials",
      "Responded at t0 and age >=60 years",
      "Elevated depressive symptoms",
      need_label,
      "Classifiable response strategy",
      "Outcome observed"
    )
  ) %>%
    mutate(display = paste0(label, "\n", fmt(n_trials), " trials | ", fmt(n_people), " people"))

  excluded <- tibble(
    cohort = cohort,
    node_type = "excluded",
    node = c("baseline", "depression", "need", "strategy", "outcome"),
    y = c(87, 71, 55, 39, 23),
    n_trials = c(
      get_n("candidate") - get_n("age"),
      get_n("age") - get_n("depression"),
      get_n("depression") - get_n("eligible"),
      get_n("eligible") - get_n("strategy"),
      get_n("strategy") - get_n("outcome")
    ),
    n_people = NA_integer_,
    display = c(
      paste0("Excluded: ", fmt(get_n("candidate") - get_n("age")), " trials\n",
             "No t0 response: ", fmt(get_n("candidate") - get_n("responded")), "\n",
             "Age <60: ", fmt(get_n("responded") - get_n("age"))),
      paste0("Below/missing depression threshold\n", fmt(get_n("age") - get_n("depression")), " trials"),
      paste0("Need/FI criteria: ", fmt(get_n("depression") - get_n("eligible")), " trials\n",
             "Need domains missing: ", fmt(get_n("depression") - get_n("need_observed")), "\n",
             "No actionable Need: ", fmt(get_n("need_observed") - get_n("any_need")), "\n",
             "Invalid baseline FI: ", fmt(get_n("any_need") - get_n("eligible"))),
      paste0("Response path not classifiable: ", fmt(get_n("eligible") - get_n("strategy")), "\n",
             "t1 response missing: ", fmt(get_n("eligible") - get_n("t1")), "\n",
             "t2 response missing: ", fmt(get_n("t1") - get_n("t2")), "\n",
             "Incompatible path: ", fmt(get_n("t2") - get_n("strategy"))),
      paste0("Outcome unobserved\n", fmt(get_n("strategy") - get_n("outcome")), " trials")
    )
  )

  strategies <- d %>%
    filter(outcome_ok) %>%
    group_by(observed_strategy) %>%
    summarise(n_trials = n(), n_people = n_distinct(person_id), .groups = "drop") %>%
    complete(observed_strategy = c("Persistent", "Delayed", "Early"),
             fill = list(n_trials = 0, n_people = 0)) %>%
    mutate(
      cohort = cohort, node_type = "strategy", node = as.character(observed_strategy),
      y = 4.5,
      x = recode(as.character(observed_strategy), Persistent = 0.14, Delayed = 0.34, Early = 0.54),
      display = paste0(observed_strategy, "\n", fmt(n_trials), " | ", fmt(n_people))
    ) %>%
    select(cohort, node_type, node, n_trials, n_people, y, x, display)

  list(main = main, excluded = excluded, strategies = strategies,
       source = bind_rows(main, excluded, strategies))
}

hrs_flow <- build_flow_counts(hrs_raw, "HRS")
clhls_flow <- build_flow_counts(clhls_raw, "CLHLS")
flow_source <- bind_rows(hrs_flow$source, clhls_flow$source)
write_csv(flow_source, file.path(paths$source, "Figure_inclusion_exclusion_flow_source.csv"))

arrow_closed <- arrow(type = "closed", length = unit(1.5, "mm"))

draw_panel <- function(flow, panel_title) {
  main <- flow$main %>% mutate(xmin = 0.04, xmax = 0.60, ymin = y - 4.4, ymax = y + 4.4)
  excl <- flow$excluded %>% mutate(xmin = 0.66, xmax = 0.99, ymin = y - 5.4, ymax = y + 5.4)
  strat <- flow$strategies %>% mutate(xmin = x - 0.085, xmax = x + 0.085, ymin = 1, ymax = 8)
  vertical <- tibble(x = 0.32, xend = 0.32,
                     y = main$ymin[-nrow(main)], yend = main$ymax[-1])
  branches <- tibble(y = excl$y, x = 0.32, xend = 0.655)
  strategy_vertical <- tibble(x = strat$x, y = 9.2, yend = strat$ymax)

  ggplot() +
    geom_segment(data = vertical, aes(x = x, xend = xend, y = y, yend = yend),
                 linewidth = 0.42, colour = "black", arrow = arrow_closed) +
    geom_segment(data = branches, aes(x = x, xend = xend, y = y, yend = y),
                 linewidth = 0.36, colour = "black", arrow = arrow_closed) +
    geom_segment(aes(x = 0.32, xend = 0.32, y = 10.6, yend = 9.2), linewidth = 0.42) +
    geom_segment(aes(x = 0.14, xend = 0.54, y = 9.2, yend = 9.2), linewidth = 0.42) +
    geom_segment(data = strategy_vertical, aes(x = x, xend = x, y = y, yend = yend),
                 linewidth = 0.36, arrow = arrow_closed) +
    geom_rect(data = main, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "white", colour = "black", linewidth = 0.55) +
    geom_rect(data = excl, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "#F2F2F2", colour = "black", linewidth = 0.42, linetype = "22") +
    geom_rect(data = strat, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "white", colour = "black", linewidth = 0.5) +
    geom_text(data = main, aes(x = 0.32, y = y, label = display),
              family = "Arial", fontface = "plain", size = 2.55, lineheight = 0.95) +
    geom_text(data = excl, aes(x = 0.825, y = y, label = display),
              family = "Arial", size = 2.15, lineheight = 0.92) +
    geom_text(data = strat, aes(x = x, y = y, label = display),
              family = "Arial", size = 2.1, lineheight = 0.92) +
    annotate("text", x = 0.34, y = 0.15, label = "person-trials | unique people",
             family = "Arial", size = 1.9, colour = "#444444") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 101), clip = "off") +
    labs(title = panel_title) +
    theme_void(base_family = "Arial", base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", size = 9, hjust = 0.04, margin = margin(b = 2)),
      plot.margin = margin(2, 3, 2, 3)
    )
}

p_hrs <- draw_panel(hrs_flow, "a  HRS primary analysis")
p_clhls <- draw_panel(clhls_flow, "b  CLHLS external replication")
fig_flow <- p_hrs / p_clhls + plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Participant and person-trial inclusion flow",
    subtitle = "Need-eligible sequential target trial emulation",
    caption = paste(strwrap(paste(
      "Counts inside main and strategy boxes are person-trials | unique people; exclusion counts are sequential person-trial counts.",
      "HRS uses PA+smoking+risky alcohol with H-mPFI-28; CLHLS uses PA+smoking with common-FI20."
    ), width = 125), collapse = "\n"),
    theme = theme(
      plot.title = element_text(family = "Arial", face = "bold", size = 11),
      plot.subtitle = element_text(family = "Arial", size = 8.5),
      plot.caption = element_text(family = "Arial", size = 7, hjust = 0)
    )
  )

stem <- file.path(paths$figures, "Figure_inclusion_exclusion_flow")
svglite::svglite(paste0(stem, ".svg"), width = 7.2, height = 10.0, bg = "white"); print(fig_flow); dev.off()
grDevices::cairo_pdf(paste0(stem, ".pdf"), width = 7.2, height = 10.0, family = "Arial", bg = "white"); print(fig_flow); dev.off()
ragg::agg_tiff(paste0(stem, ".tiff"), width = 7.2, height = 10.0, units = "in", res = 600,
               compression = "lzw", background = "white"); print(fig_flow); dev.off()
ragg::agg_png(paste0(stem, ".png"), width = 7.2, height = 10.0, units = "in", res = 300,
              background = "white"); print(fig_flow); dev.off()

cat("Black-and-white inclusion/exclusion flow figure exported.\n")
