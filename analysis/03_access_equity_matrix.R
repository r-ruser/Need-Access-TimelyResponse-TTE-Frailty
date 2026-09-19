suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(purrr)})

hrs_raw <- readRDS(input_paths$hrs_candidates)
clhls_raw <- readRDS(input_paths$clhls_candidates)
hrs <- derive_hrs_spec(hrs_raw, c("pa", "smk", "alc"), "fi28", "HRS_A") %>% add_access_derivations()
clhls <- derive_clhls_common(clhls_raw) %>% add_access_derivations()

access_specs <- list(
  HRS = list(
    primary = list(
      insured = c("Not insured", "Insured"),
      wealth_high = c("Lower wealth", "Higher wealth"),
      urban = c("Rural", "Urban"),
      living_alone = c("Not living alone", "Living alone")
    ),
    secondary = list(
      income_high = c("Lower income", "Higher income"),
      education_high = c("Lower education", "Higher education"),
      social_active = c("No social activity", "Socially active"),
      prior_doctor = c("No prior physician contact", "Prior physician contact"),
      prior_hosp = c("No prior hospitalization", "Prior hospitalization")
    )
  ),
  CLHLS = list(
    primary = list(
      insured = c("Not insured", "Insured"),
      urban = c("Rural", "Urban"),
      living_alone = c("Not living alone", "Living alone")
    ),
    secondary = list(
      income_high = c("Lower income", "Higher income"),
      education_high = c("Lower education", "Higher education"),
      adequate_medical_access = c("Limited medical access", "Adequate medical access"),
      prior_hosp = c("No prior hospitalization", "Prior hospitalization")
    )
  )
)

run_access <- function(dat, cohort, tier, specs, B) {
  bind_rows(imap(specs, function(labels, exposure) {
    standardized_access(dat, exposure, labels[[1]], labels[[2]], B = B,
                        seed = analysis_seed + sum(utf8ToInt(paste(cohort, exposure)))) %>%
      mutate(cohort = cohort, tier = tier)
  }))
}

access_results <- bind_rows(
  run_access(hrs, "HRS", "Primary", access_specs$HRS$primary, 200),
  run_access(hrs, "HRS", "Secondary", access_specs$HRS$secondary, 120),
  run_access(clhls, "CLHLS", "Primary", access_specs$CLHLS$primary, 200),
  run_access(clhls, "CLHLS", "Secondary", access_specs$CLHLS$secondary, 120),
  tibble(
    exposure = "wealth_high", level0 = "Lower wealth", level1 = "Higher wealth", n = nrow(clhls),
    p0 = NA_real_, p1 = NA_real_, pd = NA_real_, pr = NA_real_,
    p0_low = NA_real_, p0_high = NA_real_, p1_low = NA_real_, p1_high = NA_real_,
    pd_low = NA_real_, pd_high = NA_real_, pr_low = NA_real_, pr_high = NA_real_,
    support = FALSE, cohort = "CLHLS", tier = "Primary"
  )
) %>%
  mutate(note = case_when(
    cohort == "CLHLS" & exposure == "wealth_high" ~ "Wealth not harmonized in CLHLS",
    !support ~ "Insufficient support",
    TRUE ~ "Need-standardized marginal probability"
  ))
write_csv(access_results, file.path(paths$tables, "Table_access_equity.csv"))

# Need x Access policy-management matrix, based on the HRS primary full-construct analysis.
hrs_w <- readRDS(file.path(paths$models, "construct_A_weighted.rds")) %>%
  add_access_derivations() %>%
  mutate(
    need_burden = factor(case_when(
      need_count == 1 ~ "1 need", need_count == 2 ~ "2 needs", need_count >= 3 ~ ">=3 needs"
    ), levels = c("1 need", "2 needs", ">=3 needs")),
    low_wealth = ifelse(is.na(wealth_high), NA_real_, 1 - wealth_high),
    limited_insurance = ifelse(is.na(insured), NA_real_, 1 - insured)
  )

matrix_vars <- c(low_wealth = "Wealth", limited_insurance = "Insurance",
                 rural = "Residence", living_alone = "Living arrangement")
matrix_data <- bind_rows(imap(matrix_vars, function(access_name, v) {
  hrs_w %>% filter(.data[[v]] %in% c(0, 1), !is.na(need_burden), !is.na(strategy)) %>%
    mutate(access_level = ifelse(.data[[v]] == 1, "Limited access", "Favorable access")) %>%
    group_by(need_burden, access_level) %>%
    summarise(
      access_indicator = access_name,
      n_person_trials = n(), n_people = n_distinct(person_id),
      early_pct = 100 * weighted_mean(as.numeric(strategy == "Early"), sw_treatment * survey_norm),
      delayed_pct = 100 * weighted_mean(as.numeric(strategy == "Delayed"), sw_treatment * survey_norm),
      persistent_pct = 100 * weighted_mean(as.numeric(strategy == "Persistent"), sw_treatment * survey_norm),
      frailty_or_death_pct = 100 * weighted_mean(outcome, final_weight),
      .groups = "drop"
    )
}))
write_csv(matrix_data, file.path(paths$tables, "Table_need_access_matrix.csv"))
write_csv(matrix_data, file.path(paths$source, "Figure_need_access_heatmap_source.csv"))

cat("Priority 3 complete: Need -> Access -> Early-response equity and policy matrix\n")
print(access_results %>% filter(tier == "Primary") %>%
        select(cohort, exposure, level0, level1, p0, p1, pd, pd_low, pd_high, pr, support, note))
print(matrix_data)

