suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(purrr); library(survey)})

hrs_raw <- readRDS(input_paths$hrs_candidates)
need_sets <- list(
  PA = "pa", SMK = "smk", ALC = "alc",
  `PA+SMK` = c("pa", "smk"), `PA+ALC` = c("pa", "alc"),
  `SMK+ALC` = c("smk", "alc"), `PA+SMK+ALC` = c("pa", "smk", "alc")
)

component_rows <- list()
component_fits <- list()
for (pop in c("Primary", "Stable-health")) {
  for (nm in names(need_sets)) {
    dat <- derive_hrs_spec(
      hrs_raw, need_sets[[nm]], "fi28",
      analysis_id = paste0("COMP_", gsub("[+]", "_", nm), "_", gsub("-", "_", pop)),
      stable_only = pop == "Stable-health"
    )
    fit <- fit_spec(dat)
    key <- paste(nm, pop, sep = "__")
    component_fits[[key]] <- fit
    den <- attr(dat, "need_prevalence_den"); num <- attr(dat, "need_prevalence_num")
    if (!fit$success) {
      component_rows[[key]] <- tibble(
        need_definition = nm, population = pop, status = paste("FAILED:", fit$error),
        eligible_person_trials = nrow(dat), eligible_people = n_distinct(dat$person_id),
        baseline_need_prevalence = ifelse(den > 0, num / den, NA_real_)
      )
      next
    }
    s <- fit$summary
    ep_rd <- s$effects %>% filter(measure == "RD", contrast == "Early vs Persistent")
    ep_rr <- s$effects %>% filter(measure == "RR", contrast == "Early vs Persistent")
    counts <- fit$weighted %>% count(strategy, name = "n_person_trials") %>%
      left_join(fit$weighted %>% group_by(strategy) %>% summarise(n_people = n_distinct(person_id), .groups = "drop"), by = "strategy") %>%
      pivot_wider(names_from = strategy, values_from = c(n_person_trials, n_people), names_sep = "_")
    component_rows[[key]] <- bind_cols(
      tibble(
        need_definition = nm, population = pop, status = "OK",
        eligible_person_trials = nrow(dat), eligible_people = n_distinct(dat$person_id),
        baseline_need_prevalence = ifelse(den > 0, num / den, NA_real_),
        rd_early_persistent = ep_rd$estimate, rd_low = ep_rd$conf_low, rd_high = ep_rd$conf_high,
        rr_early_persistent = ep_rr$estimate, rr_low = ep_rr$conf_low, rr_high = ep_rr$conf_high,
        min_strategy_n = s$gate$min_strategy_people, min_strategy_ess = s$gate$min_strategy_ess,
        max_weighted_abs_smd = s$gate$max_weighted_abs_smd,
        p1_treatment_probability = s$gate$p1_treatment_probability,
        support_pass = s$gate$support_pass, balance_pass = s$gate$balance_pass,
        positivity_pass = s$gate$positivity_pass
      ), counts
    )
    saveRDS(fit$weighted, file.path(paths$models, paste0("component_", gsub("[+]", "_", nm), "_", gsub("-", "_", pop), ".rds")), compress = "xz")
  }
}

components <- bind_rows(component_rows) %>%
  group_by(need_definition) %>%
  mutate(
    stable_health_rd_shift = ifelse(
      population == "Stable-health" & all(c("Primary", "Stable-health") %in% population),
      rd_early_persistent - rd_early_persistent[population == "Primary"], NA_real_
    )
  ) %>% ungroup()
write_csv(components, file.path(paths$tables, "Table_need_components.csv"))
saveRDS(component_fits, file.path(paths$models, "need_component_fits.rds"), compress = "xz")

# Need-severity gradient in the HRS full-need primary estimand.
full_dat <- derive_hrs_spec(hrs_raw, c("pa", "smk", "alc"), "fi28", "HRS_A") %>%
  mutate(need_burden = factor(case_when(
    need_count == 1 ~ "1 need", need_count == 2 ~ "2 needs", need_count >= 3 ~ ">=3 needs"
  ), levels = c("1 need", "2 needs", ">=3 needs")))
full_w <- component_fits[["PA+SMK+ALC__Primary"]]$weighted %>%
  mutate(need_burden = factor(case_when(
    need_count == 1 ~ "1 need", need_count == 2 ~ "2 needs", need_count >= 3 ~ ">=3 needs"
  ), levels = c("1 need", "2 needs", ">=3 needs")))

dd <- full_dat %>% filter(!is.na(a1), !is.na(need_burden)) %>% prepare_covariates()
wave_term <- if (n_distinct(dd$wave) > 1) "+ factor(wave)" else ""
sev_formula <- as.formula(paste(
  "a1 ~ need_burden + age_z + female_i + education_z + fi0_z + dep_z +",
  "multimorbidity_z + cognition_z", wave_term
))
sev_mod <- glm(sev_formula, dd, family = binomial())
sev_levels <- levels(dd$need_burden)
predict_sev <- function(mod, data, mult = rep(1, nrow(data))) {
  sapply(sev_levels, function(lv) {
    nd <- data; nd$need_burden <- factor(lv, levels = sev_levels)
    weighted_mean(predict(mod, nd, type = "response"), mult)
  })
}
sev_point <- predict_sev(sev_mod, dd)
set.seed(analysis_seed + 20)
ids <- unique(dd$person_id)
sev_boot <- replicate(300, {
  draw <- sample(ids, length(ids), replace = TRUE); tab <- table(draw)
  mult <- as.numeric(tab[match(dd$person_id, names(tab))]); mult[is.na(mult)] <- 0
  tryCatch({
    db <- dd
    db$.boot_weight <- mult
    m <- suppressWarnings(glm(sev_formula, data = db, family = binomial(),
                              weights = .boot_weight))
    predict_sev(m, db, db$.boot_weight)
  }, error = function(e) rep(NA_real_, length(sev_levels)))
})
sev_q <- t(apply(sev_boot, 1, quantile, c(.025, .975), na.rm = TRUE))
severity_response <- tibble(
  need_burden = sev_levels, n = as.integer(table(dd$need_burden)[sev_levels]),
  standardized_early_probability = sev_point,
  conf_low = sev_q[, 1], conf_high = sev_q[, 2]
)
write_csv(severity_response, file.path(paths$tables, "Table_need_severity_response.csv"))

# Additive-scale effect modification by baseline need burden using the primary CCW weights.
dw <- full_w %>% filter(outcome_observed == 1, !is.na(outcome), !is.na(need_burden), is.finite(final_weight), final_weight > 0)
des <- svydesign(ids = ~person_id, weights = ~final_weight, data = dw)
# A survey-weighted linear probability model targets additive interaction
# directly and avoids boundary failures from binomial identity-link fitting.
int_fit <- svyglm(outcome ~ strategy * need_burden, des, family = gaussian())
b <- coef(int_fit); V <- vcov(int_fit)
linear_est <- function(L) {
  L2 <- setNames(rep(0, length(b)), names(b)); L2[names(L)] <- L
  est <- sum(L2 * b); se <- sqrt(as.numeric(t(L2) %*% V %*% L2))
  c(est = est, low = est - 1.96 * se, high = est + 1.96 * se, p = 2 * pnorm(abs(est / se), lower.tail = FALSE))
}
severity_effects <- bind_rows(lapply(sev_levels, function(lv) {
  L <- c(strategyEarly = 1)
  if (lv == "2 needs") L["strategyEarly:need_burden2 needs"] <- 1
  if (lv == ">=3 needs") L["strategyEarly:need_burden>=3 needs"] <- 1
  rd <- linear_est(L)
  interaction <- if (lv == "1 need") c(est = 0, low = NA, high = NA, p = NA) else {
    nm <- if (lv == "2 needs") "strategyEarly:need_burden2 needs" else "strategyEarly:need_burden>=3 needs"
    linear_est(setNames(1, nm))
  }
  sub <- dw %>% filter(need_burden == lv)
  cnt <- sub %>% group_by(strategy) %>% summarise(n_people = n_distinct(person_id), .groups = "drop")
  risk_e <- weighted_mean(sub$outcome[sub$strategy == "Early"], sub$final_weight[sub$strategy == "Early"])
  risk_p <- weighted_mean(sub$outcome[sub$strategy == "Persistent"], sub$final_weight[sub$strategy == "Persistent"])
  tibble(
    need_burden = lv, n = nrow(sub), min_strategy_people = ifelse(nrow(cnt) == 3, min(cnt$n_people), 0),
    risk_early = risk_e, risk_persistent = risk_p,
    rd_early_persistent = rd[["est"]], rd_low = rd[["low"]], rd_high = rd[["high"]],
    additive_interaction_vs_1_need = interaction[["est"]], interaction_low = interaction[["low"]],
    interaction_high = interaction[["high"]], interaction_p = interaction[["p"]],
    support_pass = nrow(cnt) == 3 && min(cnt$n_people) >= 100
  )
}))
write_csv(severity_effects, file.path(paths$tables, "Table_need_severity_TTE.csv"))

cat("Priority 2 complete: Need-component and need-severity decompositions\n")
print(components %>% filter(population == "Primary") %>%
        select(need_definition, eligible_people, baseline_need_prevalence,
               rd_early_persistent, rd_low, rd_high, rr_early_persistent,
               min_strategy_n, min_strategy_ess, support_pass))
print(severity_response)
print(severity_effects)
