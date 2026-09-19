suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(purrr); library(survey)})

hrs_raw <- readRDS(input_paths$hrs_candidates)
clhls_raw <- readRDS(input_paths$clhls_candidates)

hrs_common_base <- derive_hrs_spec(hrs_raw, c("pa", "smk"), "common20", "HRS_D")
hrs_den <- attr(hrs_common_base, "need_prevalence_den")
hrs_num <- attr(hrs_common_base, "need_prevalence_num")
hrs_common_dat <- hrs_common_base %>% add_access_derivations()
hrs_common_w <- readRDS(file.path(paths$models, "construct_D_weighted.rds")) %>% add_access_derivations()
clhls_dat <- derive_clhls_common(clhls_raw) %>% add_access_derivations()
clhls_fit <- fit_spec(clhls_dat)
if (!clhls_fit$success) stop("CLHLS common-construct fit failed: ", clhls_fit$error)
clhls_w <- clhls_fit$weighted %>% add_access_derivations()
saveRDS(clhls_w, file.path(paths$models, "CLHLS_common_weighted.rds"), compress = "xz")

transport_row <- function(dat, weighted, label, need_den, need_num) {
  s <- summarise_fit(weighted)
  ep_rd <- s$effects %>% filter(measure == "RD", contrast == "Early vs Persistent")
  ep_rr <- s$effects %>% filter(measure == "RR", contrast == "Early vs Persistent")
  risks <- s$risks %>% select(strategy, risk) %>% pivot_wider(names_from = strategy, values_from = risk, names_prefix = "risk_")
  shares <- weighted %>% filter(!is.na(strategy)) %>% group_by(strategy) %>%
    summarise(p = weighted_mean(rep(1, n()), sw_treatment * survey_norm) * sum(sw_treatment * survey_norm) / sum(weighted$sw_treatment * weighted$survey_norm, na.rm = TRUE),
              n = n(), n_people = n_distinct(person_id), .groups = "drop")
  # Recompute normalized strategy shares directly.
  shares <- weighted %>% filter(!is.na(strategy)) %>% group_by(strategy) %>%
    summarise(w = sum(sw_treatment * survey_norm, na.rm = TRUE), n = n(), n_people = n_distinct(person_id), .groups = "drop") %>%
    mutate(p = w / sum(w))
  early_cov <- shares$p[shares$strategy == "Early"]
  access <- dat %>% summarise(
    insured_pct = 100 * mean(insured == 1, na.rm = TRUE),
    higher_wealth_pct = 100 * mean(wealth_high == 1, na.rm = TRUE),
    urban_pct = 100 * mean(urban == 1, na.rm = TRUE),
    living_alone_pct = 100 * mean(living_alone == 1, na.rm = TRUE)
  )
  bind_cols(
    tibble(
      cohort = label, need_definition = "PA+SMK", outcome_definition = "common-FI20",
      eligible_person_trials = nrow(dat), eligible_people = n_distinct(dat$person_id),
      baseline_need_prevalence = need_num / need_den,
      early_response_coverage = early_cov,
      rd_early_persistent = ep_rd$estimate, rd_low = ep_rd$conf_low, rd_high = ep_rd$conf_high,
      rr_early_persistent = ep_rr$estimate, rr_low = ep_rr$conf_low, rr_high = ep_rr$conf_high,
      min_strategy_people = s$gate$min_strategy_people,
      max_weighted_abs_smd = s$gate$max_weighted_abs_smd,
      p1_treatment_probability = s$gate$p1_treatment_probability
    ), risks, access
  )
}

clhls_den <- sum(clhls_raw$eligible_trigger & clhls_raw$need_observed_domains == 2 & clhls_raw$common_fi20_n >= 16, na.rm = TRUE)
clhls_num <- sum(clhls_raw$eligible_need, na.rm = TRUE)
transport <- bind_rows(
  transport_row(hrs_common_dat, hrs_common_w, "HRS", hrs_den, hrs_num),
  transport_row(clhls_dat, clhls_w, "CLHLS", clhls_den, clhls_num)
)
rd_diff <- with(transport, rd_early_persistent[cohort == "HRS"] - rd_early_persistent[cohort == "CLHLS"])
se_hrs <- with(transport, (rd_high[cohort == "HRS"] - rd_low[cohort == "HRS"]) / (2 * 1.96))
se_clhls <- with(transport, (rd_high[cohort == "CLHLS"] - rd_low[cohort == "CLHLS"]) / (2 * 1.96))
transport <- transport %>% mutate(
  descriptive_rd_difference_hrs_minus_clhls = rd_diff,
  descriptive_rd_difference_low = rd_diff - 1.96 * sqrt(se_hrs^2 + se_clhls^2),
  descriptive_rd_difference_high = rd_diff + 1.96 * sqrt(se_hrs^2 + se_clhls^2),
  interpretation = "Descriptive heterogeneity; not a health-system causal effect"
)
write_csv(transport, file.path(paths$tables, "Table_cross_cohort_transportability.csv"))

# Outcome-specific observation weighting for post-response health-system outcomes.
fit_secondary_outcome <- function(weighted, outcome_var, outcome_label, type = c("binary", "continuous")) {
  type <- match.arg(type)
  d <- weighted %>% prepare_covariates()
  d$secondary_outcome <- d[[outcome_var]]
  d$secondary_observed <- as.numeric(!is.na(d$secondary_outcome))
  wave_term <- if (n_distinct(d$wave) > 1) "+ factor(wave)" else ""
  denom <- as.formula(paste(
    "secondary_observed ~ strategy + age_z + female_i + education_z + fi0_z + dep_z + need_count +",
    "insured_i + income_z + rural_i + alone_i + prior_hosp_i", wave_term
  ))
  num <- as.formula(paste("secondary_observed ~ strategy", wave_term))
  if (n_distinct(d$secondary_observed) > 1) {
    md <- suppressWarnings(glm(denom, d, family = binomial()))
    mn <- suppressWarnings(glm(num, d, family = binomial()))
    wd <- clamp(predict(md, type = "response")); wn <- clamp(predict(mn, type = "response"))
    d$sw_secondary <- ifelse(d$secondary_observed == 1, wn / wd, NA_real_)
  } else d$sw_secondary <- ifelse(d$secondary_observed == 1, 1, NA_real_)
  d$secondary_weight_raw <- d$sw_treatment * d$sw_secondary
  q <- quantile(d$secondary_weight_raw[d$secondary_observed == 1], c(.01, .99), na.rm = TRUE)
  d$secondary_weight <- ifelse(d$secondary_observed == 1,
                               pmin(q[[2]], pmax(q[[1]], d$secondary_weight_raw)) * d$survey_norm,
                               NA_real_)
  dd <- d %>% filter(secondary_observed == 1, is.finite(secondary_weight), secondary_weight > 0)
  des <- svydesign(ids = ~person_id, weights = ~secondary_weight, data = dd)
  fam <- gaussian()
  fit <- svyglm(secondary_outcome ~ strategy, des, family = fam)
  risks <- svyby(~secondary_outcome, ~strategy, des, svymean, vartype = c("se", "ci"), na.rm = TRUE) %>%
    as_tibble() %>% transmute(strategy = as.character(strategy), estimate = secondary_outcome,
                              se, conf_low = ci_l, conf_high = ci_u)
  eff <- bind_rows(
    contrast_from_fit(fit, "Early vs Persistent", c(strategyEarly = 1)),
    contrast_from_fit(fit, "Early vs Delayed", c(strategyEarly = 1, strategyDelayed = -1))
  ) %>% transmute(contrast, estimate, conf_low, conf_high, p_value)
  list(
    data = dd %>% mutate(outcome = secondary_outcome, final_weight = secondary_weight,
                         outcome_observed = secondary_observed),
    table = eff %>% mutate(outcome = outcome_label, outcome_type = type),
    means = risks %>% mutate(outcome = outcome_label, outcome_type = type)
  )
}

hrs_primary_w <- readRDS(file.path(paths$models, "construct_A_weighted.rds"))
health_fits <- list(
  HRS_hospital = fit_secondary_outcome(hrs_primary_w, "t3_hospital", "Subsequent hospitalization", "binary"),
  HRS_physician = fit_secondary_outcome(hrs_primary_w, "t3_physician_count", "Subsequent physician visits", "continuous"),
  HRS_oop = fit_secondary_outcome(hrs_primary_w, "t3_oop_cost", "Subsequent out-of-pocket expenditure", "continuous"),
  CLHLS_hospital = fit_secondary_outcome(clhls_w, "t3_hospital", "Subsequent hospitalization", "binary")
)
health_effects <- bind_rows(imap(health_fits, function(x, nm) x$table %>% mutate(cohort_analysis = nm)))
health_means <- bind_rows(imap(health_fits, function(x, nm) x$means %>% mutate(cohort_analysis = nm)))
write_csv(health_effects, file.path(paths$tables, "Table_health_system_outcomes.csv"))
write_csv(health_means, file.path(paths$tables, "Table_health_system_outcome_means.csv"))

# HRS primary policy coverage scenarios. Increases in Early coverage proportionally
# reallocate the observed Delayed/Persistent mixture; this is a policy-impact scenario.
scenarios <- policy_scenarios(hrs_primary_w, B = 500, seed = analysis_seed + 50)
hosp_scenarios <- policy_scenarios(health_fits$HRS_hospital$data, B = 500, seed = analysis_seed + 51) %>%
  select(scenario, expected_hospitalization_risk = expected_risk,
         hosp_conf_low = conf_low, hosp_conf_high = conf_high,
         hospitalization_difference_per_1000 = event_difference_per_1000,
         hosp_diff_low = event_diff_low, hosp_diff_high = event_diff_high)
scenarios <- scenarios %>% left_join(hosp_scenarios, by = "scenario") %>%
  mutate(interpretation = "Observational policy-impact scenario; not a policy implementation trial")
write_csv(scenarios, file.path(paths$tables, "Table_policy_scenarios.csv"))
write_csv(scenarios, file.path(paths$source, "Figure_policy_coverage_source.csv"))

cat("Priorities 5-7 complete: transportability, health-system outcomes, and policy scenarios\n")
print(transport)
print(health_effects)
print(scenarios)
