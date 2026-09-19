suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(purrr)})

hrs_raw <- readRDS(input_paths$hrs_candidates)

specs <- tibble(
  analysis_id = c("A", "B", "C", "D"),
  needs = list(c("pa", "smk", "alc"), c("pa", "smk"),
               c("pa", "smk", "alc"), c("pa", "smk")),
  outcome = c("fi28", "fi28", "common20", "common20")
)

fits <- list()
rows <- list()
for (i in seq_len(nrow(specs))) {
  id <- specs$analysis_id[[i]]
  dat <- derive_hrs_spec(hrs_raw, specs$needs[[i]], specs$outcome[[i]], paste0("HRS_", id))
  fit <- fit_spec(dat)
  fits[[id]] <- fit
  if (!fit$success) {
    rows[[id]] <- tibble(analysis = id, status = paste("FAILED:", fit$error))
    next
  }
  saveRDS(fit$weighted, file.path(paths$models, paste0("construct_", id, "_weighted.rds")), compress = "xz")
  s <- fit$summary
  risk_w <- s$risks %>% select(strategy, risk, conf_low, conf_high) %>%
    pivot_wider(names_from = strategy, values_from = c(risk, conf_low, conf_high), names_sep = "_")
  get_eff <- function(measure_value, contrast_value) {
    s$effects %>%
      filter(.data$measure == .env$measure_value,
             .data$contrast == .env$contrast_value)
  }
  ep <- get_eff("RD", "Early vs Persistent"); ed <- get_eff("RD", "Early vs Delayed")
  dp <- get_eff("RD", "Delayed vs Persistent"); rrep <- get_eff("RR", "Early vs Persistent")
  rows[[id]] <- bind_cols(
    tibble(
      row_type = "2x2 estimate", analysis = id, status = "OK",
      need_definition = paste(toupper(specs$needs[[i]]), collapse = "+"),
      outcome_definition = ifelse(specs$outcome[[i]] == "fi28", "H-mPFI-28", "common-FI20"),
      eligible_person_trials = nrow(dat), eligible_people = n_distinct(dat$person_id),
      rd_early_persistent = ep$estimate, rd_ep_low = ep$conf_low, rd_ep_high = ep$conf_high,
      rd_early_delayed = ed$estimate, rd_ed_low = ed$conf_low, rd_ed_high = ed$conf_high,
      rd_delayed_persistent = dp$estimate, rd_dp_low = dp$conf_low, rd_dp_high = dp$conf_high,
      rr_early_persistent = rrep$estimate, rr_ep_low = rrep$conf_low, rr_ep_high = rrep$conf_high,
      min_strategy_n = s$gate$min_strategy_people, min_strategy_ess = s$gate$min_strategy_ess,
      max_weighted_abs_smd = s$gate$max_weighted_abs_smd,
      p1_treatment_probability = s$gate$p1_treatment_probability,
      support_pass = s$gate$support_pass, balance_pass = s$gate$balance_pass,
      positivity_pass = s$gate$positivity_pass
    ), risk_w
  )
  write_csv(s$risks %>% mutate(analysis = id), file.path(paths$source, paste0("construct_", id, "_risks.csv")))
  write_csv(s$effects %>% mutate(analysis = id), file.path(paths$source, paste0("construct_", id, "_effects.csv")))
}

stopifnot(all(vapply(fits, function(x) x$success, logical(1))))

shift_defs <- tibble::tribble(
  ~shift, ~from, ~to, ~interpretation,
  "B - A", "A", "B", "Need-definition shift within H-mPFI-28",
  "D - C", "C", "D", "Need-definition shift within common-FI20",
  "C - A", "A", "C", "Outcome-definition shift within PA+SMK+ALC",
  "D - B", "B", "D", "Outcome-definition shift within PA+SMK"
)
shift_rows <- pmap_dfr(shift_defs, function(shift, from, to, interpretation) {
  z <- bootstrap_rd_shift(fits[[from]]$weighted, fits[[to]]$weighted, B = 300,
                          seed = analysis_seed + utf8ToInt(from) + utf8ToInt(to))
  bind_cols(tibble(shift, from, to, interpretation), z)
})
write_csv(shift_rows, file.path(paths$tables, "Table_construct_shifts.csv"))

decomp <- bind_rows(rows) %>%
  left_join(
    shift_rows %>% filter(shift %in% c("B - A", "D - C")) %>%
      transmute(analysis = to, need_shift = estimate, need_shift_low = conf_low, need_shift_high = conf_high),
    by = "analysis"
  ) %>%
  left_join(
    shift_rows %>% filter(shift %in% c("C - A", "D - B")) %>%
      transmute(analysis = to, outcome_shift = estimate, outcome_shift_low = conf_low, outcome_shift_high = conf_high),
    by = "analysis"
  )
write_csv(decomp, file.path(paths$tables, "Table_construct_decomposition.csv"))
saveRDS(fits, file.path(paths$models, "construct_decomposition_fits.rds"), compress = "xz")

cat("Priority 1 complete: HRS 2x2 construct decomposition\n")
print(decomp %>% select(analysis, need_definition, outcome_definition, rd_early_persistent,
                        rd_ep_low, rd_ep_high, need_shift, outcome_shift,
                        min_strategy_n, min_strategy_ess, max_weighted_abs_smd,
                        p1_treatment_probability))
