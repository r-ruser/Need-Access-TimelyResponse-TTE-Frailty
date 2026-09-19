suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(survey)
})

common20_stems <- c(
  "adl_dress", "adl_bath", "adl_eat", "adl_bed", "adl_toilet",
  "iadl_shop", "iadl_meals", "mob_long_walk", "mob_chair", "mob_stoop", "mob_lift",
  "dis_hypertension", "dis_diabetes", "dis_cancer", "dis_lung", "dis_heart",
  "dis_stroke", "dis_arthritis", "oth_srh", "oth_underweight"
)

need_labels <- c(pa = "Physical activity", smk = "Smoking", alc = "Risky alcohol")

calc_common20 <- function(d, prefix = "") {
  nms <- paste0(prefix, common20_stems)
  miss <- setdiff(nms, names(d))
  if (length(miss)) stop("Missing common-FI20 columns: ", paste(miss, collapse = ", "))
  mat <- as.matrix(d[, nms, drop = FALSE])
  n <- rowSums(!is.na(mat))
  value <- rowSums(mat, na.rm = TRUE) / n
  value[n < 16] <- NA_real_
  list(value = value, n = n)
}

calc_adequate <- function(d, needs, prefix) {
  need_cols <- paste0("need_", needs)
  need_mat <- as.matrix(d[, need_cols, drop = FALSE])
  target_list <- list(
    pa = d[[paste0(prefix, "pa_active")]],
    smk = ifelse(is.na(d[[paste0(prefix, "current_smoking")]]), NA_real_,
                 as.numeric(d[[paste0(prefix, "current_smoking")]] == 0)),
    alc = ifelse(is.na(d[[paste0(prefix, "alcohol_risk")]]), NA_real_,
                 as.numeric(d[[paste0(prefix, "alcohol_risk")]] == 0))
  )
  target_mat <- do.call(cbind, target_list[needs])
  if (length(needs) == 1) target_mat <- matrix(target_mat, ncol = 1)
  miss_required <- rowSums(need_mat == 1 & is.na(target_mat), na.rm = TRUE) > 0
  need_count <- rowSums(need_mat, na.rm = TRUE)
  success <- rowSums(ifelse(need_mat == 1, target_mat, 0), na.rm = TRUE)
  prop <- success / need_count
  prop[miss_required | need_count < 1] <- NA_real_
  adequate <- ifelse(
    is.na(prop), NA_real_,
    ifelse(need_count == 1, as.numeric(prop == 1), as.numeric(prop >= 0.5))
  )
  list(value = adequate, prop = prop, count = need_count,
       observed = rowSums(!is.na(need_mat)) == length(needs))
}

make_need_pattern <- function(d, needs) {
  flags <- as.matrix(d[, paste0("need_", needs), drop = FALSE])
  apply(flags, 1, function(z) {
    if (any(is.na(z))) return(NA_character_)
    active <- needs[z == 1]
    if (!length(active)) "None" else paste(toupper(active), collapse = "+")
  })
}

derive_stable_health <- function(d, fi0, fi1) {
  no_recent_hosp <- ifelse(is.na(d$prior_hosp), NA, d$prior_hosp == 0)
  no_new_cancer <- ifelse(is.na(d$dis_cancer) | is.na(d$t1_dis_cancer), NA,
                          !(d$dis_cancer == 0 & d$t1_dis_cancer == 1))
  no_new_cvd <- ifelse(
    is.na(d$dis_heart) | is.na(d$t1_dis_heart) | is.na(d$dis_stroke) | is.na(d$t1_dis_stroke),
    NA,
    !((d$dis_heart == 0 & d$t1_dis_heart == 1) | (d$dis_stroke == 0 & d$t1_dis_stroke == 1))
  )
  no_weight_loss <- ifelse(
    is.na(d$bmi) | is.na(d$t1_bmi) | d$bmi <= 0, NA,
    d$t1_bmi >= 0.95 * d$bmi
  )
  no_rapid_fi <- ifelse(is.na(fi0) | is.na(fi1), NA, fi1 - fi0 <= 0.05)
  no_mobility_decline <- ifelse(
    is.na(d$mobility_count) | is.na(d$t1_mobility_count), NA,
    d$t1_mobility_count <= d$mobility_count
  )
  no_srh_decline <- ifelse(
    is.na(d$oth_srh) | is.na(d$t1_oth_srh), NA,
    d$t1_oth_srh - d$oth_srh <= 0.25
  )
  checks <- cbind(no_recent_hosp, no_new_cancer, no_new_cvd, no_weight_loss,
                  no_rapid_fi, no_mobility_decline, no_srh_decline)
  n_observed <- rowSums(!is.na(checks))
  stable <- n_observed >= 5 & rowSums(checks == FALSE, na.rm = TRUE) == 0
  list(stable = stable, n_observed = n_observed)
}

derive_hrs_spec <- function(raw, needs, outcome = c("fi28", "common20"),
                            analysis_id, stable_only = FALSE) {
  outcome <- match.arg(outcome)
  d <- raw %>% filter(cohort == "HRS")
  if (outcome == "fi28") {
    fi0 <- d$h_mpfi28; fi0_n <- d$h_mpfi28_n
    fi1 <- d$t1_h_mpfi28; fi2 <- d$t2_h_mpfi28; fi3 <- d$t3_h_mpfi28
    valid0 <- fi0_n >= 23
  } else {
    f0 <- calc_common20(d, ""); f1 <- calc_common20(d, "t1_")
    f2 <- calc_common20(d, "t2_"); f3 <- calc_common20(d, "t3_")
    fi0 <- f0$value; fi0_n <- f0$n; fi1 <- f1$value; fi2 <- f2$value; fi3 <- f3$value
    valid0 <- fi0_n >= 16
  }
  a1 <- calc_adequate(d, needs, "t1_")
  a2 <- calc_adequate(d, needs, "t2_")
  a1_value <- a1$value; a2_value <- a2$value
  a1_prop <- a1$prop; a2_prop <- a2$prop
  n_need <- a1$count; need_observed <- a1$observed
  pattern <- make_need_pattern(d, needs)
  stable <- derive_stable_health(d, fi0, fi1)
  mm <- rowSums(as.matrix(d[, paste0("dis_", c("hypertension", "diabetes", "cancer", "lung", "heart", "stroke", "arthritis"))]), na.rm = TRUE)
  mm[rowSums(!is.na(as.matrix(d[, paste0("dis_", c("hypertension", "diabetes", "cancer", "lung", "heart", "stroke", "arthritis"))]))) < 4] <- NA_real_
  need_prevalence_den <- d$responded == 1 & d$age >= 60 & d$dep_trigger == 1 & valid0 & need_observed
  eligible <- need_prevalence_den & n_need >= 1
  if (stable_only) eligible <- eligible & stable$stable

  out <- d %>% transmute(
    person_id, cohort = "HRS", trial_id, wave = t0_wave, year = t0_year,
    analysis_id = analysis_id, need_definition = paste(toupper(needs), collapse = "+"),
    outcome_definition = ifelse(outcome == "fi28", "H-mPFI-28", "common-FI20"),
    population = ifelse(stable_only, "Stable-health", "Primary"),
    age, female, education, married, fi0 = fi0, fi0_n = fi0_n,
    dep = cesd, need_count = n_need, need_pattern = pattern,
    need_pa, need_smk, need_alc,
    insured, wealth, income, rural, living_alone, prior_hosp, prior_doctor,
    social_active, cognition, multimorbidity = mm, survey_weight,
    a1 = a1_value, a2 = a2_value,
    response_prop_t1 = a1_prop, response_prop_t2 = a2_prop,
    t1_fi = fi1, t1_insured, t1_income, t2_fi = fi2,
    t3_fi = fi3, t3_responded,
    died = as.numeric(t3_death_by_wave == 1),
    t3_hospital = t3_prior_hosp, t3_physician_count = t3_doctor_count,
    t3_oop_cost = t3_oop_cost,
    stable_health = stable$stable, stable_checks_observed = stable$n_observed,
    eligible = eligible,
    strategy = case_when(
      a1_value == 1 & a2_value == 1 ~ "Early",
      a1_value == 0 & a2_value == 1 ~ "Delayed",
      a1_value == 0 & a2_value == 0 ~ "Persistent",
      TRUE ~ NA_character_
    ),
    outcome_observed = as.numeric(died == 1 | (t3_responded == 1 & !is.na(fi3))),
    outcome = case_when(
      died == 1 ~ 1,
      t3_responded == 1 & !is.na(fi3) ~ as.numeric(fi3 >= 0.25),
      TRUE ~ NA_real_
    )
  )
  result <- out %>% filter(eligible)
  attr(result, "need_prevalence_den") <- sum(need_prevalence_den, na.rm = TRUE)
  attr(result, "need_prevalence_num") <- sum(need_prevalence_den & n_need >= 1, na.rm = TRUE)
  result
}

derive_clhls_common <- function(raw, analysis_id = "CLHLS_COMMON") {
  mm_cols <- paste0("dis_", c("hypertension", "diabetes", "cancer", "lung", "heart", "stroke", "arthritis"))
  mm_mat <- as.matrix(raw[, mm_cols])
  mm <- rowSums(mm_mat, na.rm = TRUE)
  mm[rowSums(!is.na(mm_mat)) < 4] <- NA_real_
  raw %>% transmute(
    person_id, cohort = "CLHLS", trial_id, wave = t0_wave, year = t0_year,
    analysis_id, need_definition = "PA+SMK", outcome_definition = "common-FI20",
    population = "Primary", age, female, education, married = NA_real_,
    fi0 = common_fi20, fi0_n = common_fi20_n, dep = dep_proxy,
    need_count, need_pattern = case_when(
      need_pa == 1 & need_smk == 1 ~ "PA+SMK",
      need_pa == 1 ~ "PA",
      need_smk == 1 ~ "SMK",
      TRUE ~ "None"
    ),
    need_pa, need_smk, need_alc = NA_real_,
    insured, wealth = NA_real_, income, rural, living_alone, prior_hosp,
    prior_doctor = NA_real_, social_active = NA_real_, cognition = NA_real_,
    adequate_medical_access, multimorbidity = mm, survey_weight,
    a1 = adequate_t1, a2 = adequate_t2,
    response_prop_t1, response_prop_t2,
    t1_fi = t1_common_fi20, t1_insured, t1_income,
    t2_fi = t2_common_fi20, t3_fi = t3_common_fi20, t3_responded,
    died = as.numeric(died_during_followup),
    t3_hospital = t3_prior_hosp, t3_physician_count = NA_real_, t3_oop_cost = NA_real_,
    stable_health = NA, stable_checks_observed = NA_integer_,
    eligible = eligible_need, strategy = observed_strategy,
    outcome_observed, outcome = frailty_or_death
  ) %>% filter(eligible)
}

impute_numeric <- function(x) {
  med <- median(x[is.finite(x)], na.rm = TRUE)
  if (!is.finite(med)) med <- 0
  ifelse(is.finite(x), x, med)
}
zscore <- function(x) {
  z <- impute_numeric(x)
  s <- sd(z)
  if (!is.finite(s) || s == 0) return(rep(0, length(z)))
  as.numeric((z - mean(z)) / s)
}
impute01 <- function(x) {
  p <- mean(x %in% 1, na.rm = TRUE)
  fill <- ifelse(is.finite(p) && p >= 0.5, 1, 0)
  ifelse(x %in% c(0, 1), x, fill)
}
clamp <- function(p, lo = 0.01) pmin(1 - lo, pmax(lo, p))
actual_prob <- function(p, a) ifelse(a == 1, p, 1 - p)

prepare_covariates <- function(d) {
  d %>% mutate(
    age_m = as.numeric(is.na(age)), age_z = zscore(age),
    female_m = as.numeric(is.na(female)), female_i = impute01(female),
    education_m = as.numeric(is.na(education)), education_z = zscore(education),
    fi0_m = as.numeric(is.na(fi0)), fi0_z = zscore(fi0),
    dep_m = as.numeric(is.na(dep)), dep_z = zscore(dep),
    insured_m = as.numeric(is.na(insured)), insured_i = impute01(insured),
    income_m = as.numeric(is.na(income)), income_z = zscore(log1p(pmax(impute_numeric(income), 0))),
    rural_m = as.numeric(is.na(rural)), rural_i = impute01(rural),
    alone_m = as.numeric(is.na(living_alone)), alone_i = impute01(living_alone),
    prior_hosp_m = as.numeric(is.na(prior_hosp)), prior_hosp_i = impute01(prior_hosp),
    t1_fi_m = as.numeric(is.na(t1_fi)), t1_fi_z = zscore(t1_fi),
    t1_insured_m = as.numeric(is.na(t1_insured)), t1_insured_i = impute01(t1_insured),
    t1_income_m = as.numeric(is.na(t1_income)), t1_income_z = zscore(log1p(pmax(impute_numeric(t1_income), 0))),
    cognition_m = as.numeric(is.na(cognition)), cognition_z = zscore(cognition),
    multimorbidity_m = as.numeric(is.na(multimorbidity)), multimorbidity_z = zscore(multimorbidity)
  )
}

safe_glm <- function(formula, data, family) {
  suppressWarnings(glm(formula, data = data, family = family,
                       control = glm.control(maxit = 100)))
}

fit_ccw <- function(d) {
  d <- d %>% filter(!is.na(a1), !is.na(a2), !is.na(strategy)) %>% prepare_covariates()
  if (!nrow(d) || n_distinct(d$strategy) < 3) stop("Fewer than three supported strategies.")
  wave_term <- if (n_distinct(d$wave) > 1) "factor(wave)" else NULL
  denom_terms <- c(
    wave_term, "age_z", "female_i", "education_z", "fi0_z", "dep_z", "need_count",
    "insured_i", "income_z", "rural_i", "alone_i", "prior_hosp_i",
    "age_m", "female_m", "education_m", "fi0_m", "dep_m", "insured_m",
    "income_m", "rural_m", "alone_m", "prior_hosp_m"
  )
  f1d <- as.formula(paste("a1 ~", paste(denom_terms, collapse = " + ")))
  f2d <- as.formula(paste(
    "a2 ~ a1 + t1_fi_z + t1_insured_i + t1_income_z + t1_fi_m + t1_insured_m + t1_income_m +",
    paste(denom_terms, collapse = " + ")
  ))
  m1d <- safe_glm(f1d, d, binomial())
  m1n <- safe_glm(a1 ~ 1, d, binomial())
  m2d <- safe_glm(f2d, d, binomial())
  m2n <- safe_glm(a2 ~ a1, d, binomial())
  d$p1d <- actual_prob(clamp(predict(m1d, type = "response")), d$a1)
  d$p1n <- actual_prob(clamp(predict(m1n, type = "response")), d$a1)
  d$p2d <- actual_prob(clamp(predict(m2d, type = "response")), d$a2)
  d$p2n <- actual_prob(clamp(predict(m2n, type = "response")), d$a2)
  d$sw_treatment <- (d$p1n / d$p1d) * (d$p2n / d$p2d)

  if (n_distinct(d$outcome_observed) == 1) {
    d$sw_outcome <- ifelse(d$outcome_observed == 1, 1, NA_real_)
  } else {
    mod_formula <- as.formula(paste("outcome_observed ~ strategy +", paste(denom_terms, collapse = " + ")))
    num_formula <- as.formula(paste("outcome_observed ~ strategy", if (!is.null(wave_term)) "+ factor(wave)" else ""))
    mo_d <- safe_glm(mod_formula, d, binomial())
    mo_n <- safe_glm(num_formula, d, binomial())
    po_d <- clamp(predict(mo_d, type = "response"))
    po_n <- clamp(predict(mo_n, type = "response"))
    d$sw_outcome <- ifelse(d$outcome_observed == 1, po_n / po_d, NA_real_)
  }
  d$sw_raw <- d$sw_treatment * d$sw_outcome
  q <- quantile(d$sw_raw[d$outcome_observed == 1], c(.01, .99), na.rm = TRUE)
  d$sw_trim <- ifelse(d$outcome_observed == 1, pmin(q[[2]], pmax(q[[1]], d$sw_raw)), NA_real_)
  survey_med <- median(d$survey_weight[d$outcome_observed == 1], na.rm = TRUE)
  if (!is.finite(survey_med) || survey_med <= 0) survey_med <- 1
  # Preserve design-zero survey weights exactly as in the audited v5 analysis;
  # only genuinely missing weights receive the neutral value 1.
  d$survey_norm <- ifelse(is.na(d$survey_weight), 1, d$survey_weight / survey_med)
  d$final_weight <- d$sw_trim * d$survey_norm
  d$strategy <- factor(d$strategy, levels = c("Persistent", "Delayed", "Early"))
  d
}

weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

contrast_from_fit <- function(fit, label, L, transform = identity) {
  b <- coef(fit); V <- vcov(fit)
  L2 <- setNames(rep(0, length(b)), names(b))
  L2[names(L)] <- L
  est <- sum(L2 * b)
  se <- sqrt(as.numeric(t(L2) %*% V %*% L2))
  tibble(
    contrast = label,
    estimate_link = est, se_link = se,
    estimate = transform(est), conf_low = transform(est - 1.96 * se),
    conf_high = transform(est + 1.96 * se),
    p_value = 2 * pnorm(abs(est / se), lower.tail = FALSE)
  )
}

analyse_binary <- function(d) {
  dd <- d %>% filter(outcome_observed == 1, !is.na(outcome), is.finite(final_weight), final_weight > 0)
  des <- svydesign(ids = ~person_id, weights = ~final_weight, data = dd)
  risks <- svyby(~outcome, ~strategy, des, svymean, vartype = c("se", "ci"), na.rm = TRUE) %>%
    as_tibble() %>% transmute(strategy = as.character(strategy), risk = outcome, se,
                              conf_low = ci_l, conf_high = ci_u)
  # Linear-probability and modified-Poisson survey models are numerically stable
  # for common outcomes while retaining the target RD and RR contrasts.
  rd_fit <- svyglm(outcome ~ strategy, des, family = gaussian())
  rr_fit <- svyglm(outcome ~ strategy, des, family = quasipoisson("log"))
  L_ep <- c(strategyEarly = 1)
  L_dp <- c(strategyDelayed = 1)
  L_ed <- c(strategyEarly = 1, strategyDelayed = -1)
  rd <- bind_rows(
    contrast_from_fit(rd_fit, "Early vs Persistent", L_ep),
    contrast_from_fit(rd_fit, "Delayed vs Persistent", L_dp),
    contrast_from_fit(rd_fit, "Early vs Delayed", L_ed)
  ) %>% mutate(measure = "RD")
  rr <- bind_rows(
    contrast_from_fit(rr_fit, "Early vs Persistent", L_ep, exp),
    contrast_from_fit(rr_fit, "Delayed vs Persistent", L_dp, exp),
    contrast_from_fit(rr_fit, "Early vs Delayed", L_ed, exp)
  ) %>% mutate(measure = "RR")
  list(risks = risks, effects = bind_rows(rd, rr) %>% select(contrast, measure, estimate, conf_low, conf_high, p_value))
}

balance_variables <- c("age_z", "female_i", "education_z", "fi0_z", "dep_z",
                       "need_count", "insured_i", "income_z", "rural_i",
                       "alone_i", "prior_hosp_i")

pair_balance <- function(d, use_weights = FALSE) {
  pairs <- list(c("Early", "Persistent"), c("Delayed", "Persistent"), c("Early", "Delayed"))
  bind_rows(lapply(pairs, function(pp) {
    a <- d %>% filter(strategy == pp[1]); b <- d %>% filter(strategy == pp[2])
    bind_rows(lapply(balance_variables, function(v) {
      pool_sd <- sd(d[[v]], na.rm = TRUE)
      if (!is.finite(pool_sd) || pool_sd == 0) pool_sd <- 1
      ma <- if (use_weights) weighted_mean(a[[v]], a$sw_treatment) else mean(a[[v]], na.rm = TRUE)
      mb <- if (use_weights) weighted_mean(b[[v]], b$sw_treatment) else mean(b[[v]], na.rm = TRUE)
      tibble(contrast = paste(pp, collapse = " vs "), variable = v, smd = (ma - mb) / pool_sd)
    }))
  }))
}

summarise_fit <- function(weighted, analysis_meta = NULL) {
  res <- analyse_binary(weighted)
  wd <- weighted %>% filter(outcome_observed == 1, !is.na(outcome)) %>%
    group_by(strategy) %>% summarise(
      n = n(), n_people = n_distinct(person_id), events = sum(outcome == 1),
      ess = sum(final_weight)^2 / sum(final_weight^2),
      sw_p01 = quantile(sw_trim, .01), sw_p99 = quantile(sw_trim, .99), .groups = "drop"
    )
  bal <- pair_balance(weighted, TRUE)
  gate <- tibble(
    min_strategy_people = min(wd$n_people),
    min_strategy_ess = min(wd$ess),
    max_weighted_abs_smd = max(abs(bal$smd), na.rm = TRUE),
    p1_treatment_probability = min(quantile(weighted$p1d, .01, na.rm = TRUE),
                                   quantile(weighted$p2d, .01, na.rm = TRUE)),
    support_pass = min_strategy_people >= 100,
    balance_pass = max_weighted_abs_smd <= .10,
    positivity_pass = p1_treatment_probability >= .01
  )
  list(risks = res$risks, effects = res$effects, weight_diag = wd,
       balance = bal, gate = gate)
}

fit_spec <- function(d) {
  tryCatch({
    w <- fit_ccw(d)
    s <- summarise_fit(w)
    list(success = TRUE, weighted = w, summary = s, error = NA_character_)
  }, error = function(e) {
    list(success = FALSE, weighted = NULL, summary = NULL, error = conditionMessage(e))
  })
}

bootstrap_rd_shift <- function(fit_from, fit_to, B = 300, seed = 20260919) {
  set.seed(seed)
  point_rd <- function(d) {
    dd <- d %>% filter(outcome_observed == 1, !is.na(outcome), strategy %in% c("Early", "Persistent"))
    weighted_mean(dd$outcome[dd$strategy == "Early"], dd$final_weight[dd$strategy == "Early"]) -
      weighted_mean(dd$outcome[dd$strategy == "Persistent"], dd$final_weight[dd$strategy == "Persistent"])
  }
  ids <- union(unique(fit_from$person_id), unique(fit_to$person_id))
  boot <- replicate(B, {
    draw <- sample(ids, length(ids), replace = TRUE)
    mult <- table(draw)
    calc <- function(d) {
      m <- as.numeric(mult[match(d$person_id, names(mult))]); m[is.na(m)] <- 0
      d2 <- d; d2$final_weight <- d2$final_weight * m
      point_rd(d2)
    }
    calc(fit_to) - calc(fit_from)
  })
  tibble(
    estimate = point_rd(fit_to) - point_rd(fit_from),
    conf_low = quantile(boot, .025, na.rm = TRUE),
    conf_high = quantile(boot, .975, na.rm = TRUE),
    bootstrap_success = sum(is.finite(boot)), bootstrap_replicates = B
  )
}

add_access_derivations <- function(d) {
  d %>% group_by(cohort, wave) %>% mutate(
    wealth_high = ifelse(is.na(wealth), NA_real_, as.numeric(wealth > median(wealth, na.rm = TRUE))),
    income_high = ifelse(is.na(income), NA_real_, as.numeric(income > median(income, na.rm = TRUE))),
    education_high = ifelse(is.na(education), NA_real_, as.numeric(education > median(education, na.rm = TRUE))),
    urban = ifelse(is.na(rural), NA_real_, 1 - rural)
  ) %>% ungroup()
}

standardized_access <- function(d, exposure, label0, label1, B = 200, seed = 20260919) {
  dd <- d %>% filter(!is.na(a1), .data[[exposure]] %in% c(0, 1)) %>%
    prepare_covariates() %>% mutate(need_pattern = factor(need_pattern))
  exposure_counts <- table(dd[[exposure]])
  if (nrow(dd) < 100 || length(exposure_counts) < 2 || min(exposure_counts) < 20) {
    return(tibble(exposure, level0 = label0, level1 = label1, n = nrow(dd),
                  p0 = NA_real_, p1 = NA_real_, pd = NA_real_, pr = NA_real_,
                  p0_low = NA_real_, p0_high = NA_real_, p1_low = NA_real_, p1_high = NA_real_,
                  pd_low = NA_real_, pd_high = NA_real_, pr_low = NA_real_, pr_high = NA_real_,
                  support = FALSE))
  }
  wave_term <- if (n_distinct(dd$wave) > 1) "+ factor(wave)" else ""
  f <- as.formula(paste(
    "a1 ~", exposure,
    "+ need_count + need_pattern + age_z + female_i + fi0_z + dep_z + multimorbidity_z + cognition_z",
    wave_term
  ))
  fit_once <- function(data, mult = NULL) {
    if (is.null(mult)) mult <- rep(1, nrow(data))
    data$.analysis_weight <- mult
    mod <- suppressWarnings(glm(f, data = data, family = binomial(), weights = .analysis_weight,
                                control = glm.control(maxit = 100)))
    d0 <- data; d1 <- data; d0[[exposure]] <- 0; d1[[exposure]] <- 1
    p0i <- predict(mod, newdata = d0, type = "response")
    p1i <- predict(mod, newdata = d1, type = "response")
    p0 <- weighted_mean(p0i, mult); p1 <- weighted_mean(p1i, mult)
    c(p0 = p0, p1 = p1, pd = p1 - p0, pr = p1 / p0)
  }
  point <- fit_once(dd)
  set.seed(seed)
  ids <- unique(dd$person_id)
  boot <- replicate(B, {
    draw <- sample(ids, length(ids), replace = TRUE)
    tab <- table(draw)
    mult <- as.numeric(tab[match(dd$person_id, names(tab))]); mult[is.na(mult)] <- 0
    tryCatch(fit_once(dd, mult), error = function(e) rep(NA_real_, 4))
  })
  qs <- apply(boot, 1, quantile, c(.025, .975), na.rm = TRUE)
  tibble(
    exposure, level0 = label0, level1 = label1, n = nrow(dd),
    p0 = point[["p0"]], p1 = point[["p1"]], pd = point[["pd"]], pr = point[["pr"]],
    p0_low = qs[1, "p0"], p0_high = qs[2, "p0"],
    p1_low = qs[1, "p1"], p1_high = qs[2, "p1"],
    pd_low = qs[1, "pd"], pd_high = qs[2, "pd"],
    pr_low = qs[1, "pr"], pr_high = qs[2, "pr"], support = TRUE
  )
}

weighted_strategy_risk <- function(d, strategy_value, outcome_var = "outcome", weight_var = "final_weight") {
  z <- d %>% filter(strategy == strategy_value)
  weighted_mean(z[[outcome_var]], z[[weight_var]])
}

policy_scenarios <- function(weighted, coverage = c(NA, .50, .70, .90), B = 500, seed = 20260919) {
  dd <- weighted %>% filter(outcome_observed == 1, !is.na(outcome), is.finite(final_weight), final_weight > 0)
  shares <- dd %>% group_by(strategy) %>% summarise(w = sum(final_weight), .groups = "drop") %>%
    mutate(p = w / sum(w))
  pE <- shares$p[shares$strategy == "Early"]
  pD <- shares$p[shares$strategy == "Delayed"]
  pP <- shares$p[shares$strategy == "Persistent"]
  rE <- weighted_strategy_risk(dd, "Early")
  rD <- weighted_strategy_risk(dd, "Delayed")
  rP <- weighted_strategy_risk(dd, "Persistent")
  non_early_risk <- (pD * rD + pP * rP) / (pD + pP)
  cov_values <- coverage; cov_values[is.na(cov_values)] <- pE
  labels <- c("Observed", paste0(round(coverage[-1] * 100), "%"))
  calc <- function(dat, covs) {
    sh <- dat %>% group_by(strategy) %>% summarise(w = sum(final_weight), .groups = "drop") %>% mutate(p = w / sum(w))
    pe <- sh$p[sh$strategy == "Early"]; pd <- sh$p[sh$strategy == "Delayed"]; pp <- sh$p[sh$strategy == "Persistent"]
    re <- weighted_strategy_risk(dat, "Early"); rd <- weighted_strategy_risk(dat, "Delayed"); rp <- weighted_strategy_risk(dat, "Persistent")
    nr <- (pd * rd + pp * rp) / (pd + pp)
    covs2 <- covs; covs2[1] <- pe
    covs2 * re + (1 - covs2) * nr
  }
  point <- calc(dd, cov_values)
  set.seed(seed)
  ids <- unique(dd$person_id)
  boot <- replicate(B, {
    draw <- sample(ids, length(ids), replace = TRUE); tab <- table(draw)
    mult <- as.numeric(tab[match(dd$person_id, names(tab))]); mult[is.na(mult)] <- 0
    db <- dd; db$final_weight <- db$final_weight * mult
    tryCatch(calc(db, cov_values), error = function(e) rep(NA_real_, length(cov_values)))
  })
  qs <- t(apply(boot, 1, quantile, c(.025, .975), na.rm = TRUE))
  diff_boot <- 1000 * sweep(boot, 2, boot[1, ], FUN = "-")
  diff_qs <- t(apply(diff_boot, 1, quantile, c(.025, .975), na.rm = TRUE))
  observed_risk <- point[1]
  tibble(
    scenario = labels, early_coverage = c(pE, coverage[-1]),
    expected_risk = point, conf_low = qs[, 1], conf_high = qs[, 2],
    event_difference_per_1000 = 1000 * (point - observed_risk),
    event_diff_low = diff_qs[, 1],
    event_diff_high = diff_qs[, 2],
    bootstrap_replicates = B
  )
}
