suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(purrr); library(survey)})

hrs_w <- readRDS(file.path(paths$models, "construct_A_weighted.rds")) %>% add_access_derivations()

access_modifiers <- list(
  wealth_high = c("Lower wealth", "Higher wealth"),
  insured = c("Not insured", "Insured"),
  urban = c("Rural", "Urban"),
  living_alone = c("Not living alone", "Living alone")
)

estimate_modifier <- function(d, variable, labels) {
  all_dat <- d %>% filter(.data[[variable]] %in% c(0, 1), !is.na(strategy)) %>%
    mutate(stratum = as.numeric(.data[[variable]]), stratum_f = factor(stratum, levels = c(0, 1)))
  if (n_distinct(all_dat$stratum) < 2) {
    return(tibble(
      access_indicator = variable, access_stratum = labels, stratum_value = 0:1,
      risk_persistent = NA_real_, risk_persistent_low = NA_real_, risk_persistent_high = NA_real_,
      risk_early = NA_real_, risk_early_low = NA_real_, risk_early_high = NA_real_,
      rd_early_persistent = NA_real_, rd_low = NA_real_, rd_high = NA_real_,
      min_strategy_people = NA_integer_, p1_treatment_probability = NA_real_,
      support_pass = FALSE, positivity_pass = FALSE,
      status = "Insufficient exposure variation", additive_interaction = NA_real_,
      interaction_low = NA_real_, interaction_high = NA_real_, interaction_p = NA_real_
    ))
  }
  gates <- all_dat %>% group_by(stratum) %>% summarise(
    min_strategy_people = {
      z <- pick(strategy, person_id) %>%
        group_by(strategy) %>% summarise(n = n_distinct(person_id), .groups = "drop")
      ifelse(nrow(z) == 3, min(z$n), 0)
    },
    p1_treatment_probability = min(quantile(p1d, .01, na.rm = TRUE), quantile(p2d, .01, na.rm = TRUE)),
    support_pass = min_strategy_people >= 100,
    positivity_pass = p1_treatment_probability >= .01,
    .groups = "drop"
  )
  dd <- all_dat %>% filter(outcome_observed == 1, !is.na(outcome), is.finite(final_weight), final_weight > 0)
  des <- svydesign(ids = ~person_id, weights = ~final_weight, data = dd)
  fit <- svyglm(outcome ~ strategy * stratum_f, des, family = gaussian())
  b <- coef(fit); V <- vcov(fit)
  lin <- function(L) {
    L2 <- setNames(rep(0, length(b)), names(b)); L2[names(L)] <- L
    est <- sum(L2 * b); se <- sqrt(as.numeric(t(L2) %*% V %*% L2))
    c(est = est, low = est - 1.96 * se, high = est + 1.96 * se,
      p = 2 * pnorm(abs(est / se), lower.tail = FALSE))
  }
  int_name <- grep("strategyEarly:stratum_f1|stratum_f1:strategyEarly", names(b), value = TRUE)
  if (!length(int_name)) int_name <- "strategyEarly:stratum_f1"
  rows <- bind_rows(lapply(0:1, function(s) {
    Lp <- c(`(Intercept)` = 1)
    Le <- c(`(Intercept)` = 1, strategyEarly = 1)
    Lrd <- c(strategyEarly = 1)
    if (s == 1) {
      Lp["stratum_f1"] <- 1
      Le["stratum_f1"] <- 1; Le[int_name] <- 1
      Lrd[int_name] <- 1
    }
    rp <- lin(Lp); re <- lin(Le); rd <- lin(Lrd)
    g <- gates %>% filter(stratum == s)
    ok <- g$support_pass && g$positivity_pass
    tibble(
      access_indicator = variable, access_stratum = labels[[s + 1]], stratum_value = s,
      risk_persistent = ifelse(ok, rp[["est"]], NA_real_),
      risk_persistent_low = ifelse(ok, rp[["low"]], NA_real_),
      risk_persistent_high = ifelse(ok, rp[["high"]], NA_real_),
      risk_early = ifelse(ok, re[["est"]], NA_real_),
      risk_early_low = ifelse(ok, re[["low"]], NA_real_),
      risk_early_high = ifelse(ok, re[["high"]], NA_real_),
      rd_early_persistent = ifelse(ok, rd[["est"]], NA_real_),
      rd_low = ifelse(ok, rd[["low"]], NA_real_), rd_high = ifelse(ok, rd[["high"]], NA_real_),
      min_strategy_people = g$min_strategy_people,
      p1_treatment_probability = g$p1_treatment_probability,
      support_pass = g$support_pass, positivity_pass = g$positivity_pass,
      status = ifelse(ok, "Estimated", "Insufficient support")
    )
  }))
  interaction <- lin(setNames(1, int_name))
  rows %>% mutate(
    additive_interaction = ifelse(all(gates$support_pass & gates$positivity_pass), interaction[["est"]], NA_real_),
    interaction_low = ifelse(all(gates$support_pass & gates$positivity_pass), interaction[["low"]], NA_real_),
    interaction_high = ifelse(all(gates$support_pass & gates$positivity_pass), interaction[["high"]], NA_real_),
    interaction_p = ifelse(all(gates$support_pass & gates$positivity_pass), interaction[["p"]], NA_real_)
  )
}

stratified <- bind_rows(imap(access_modifiers, function(labels, v) estimate_modifier(hrs_w, v, labels)))
write_csv(stratified, file.path(paths$tables, "Table_access_stratified_TTE.csv"))

cat("Priority 4 complete: additive-scale Access-stratified TTE\n")
print(stratified)
