suppressPackageStartupMessages({library(dplyr); library(readr)})

required_tables <- c(
  "Table_construct_decomposition.csv", "Table_need_components.csv",
  "Table_access_equity.csv", "Table_access_stratified_TTE.csv",
  "Table_policy_scenarios.csv", "Table_cross_cohort_transportability.csv",
  "Table_health_system_outcomes.csv"
)
required_figure_stems <- c(
  "Figure_construct_decomposition", "Figure_need_components",
  "Figure_need_access_heatmap", "Figure_access_response",
  "Figure_policy_coverage", "Figure_cross_cohort_transportability"
)
required_figures <- as.vector(outer(required_figure_stems, c("svg", "pdf", "tiff", "png"), paste, sep = "."))

checks <- bind_rows(
  tibble(
    check = paste0("table_exists:", required_tables),
    pass = file.exists(file.path(paths$tables, required_tables)),
    detail = file.path(paths$tables, required_tables)
  ),
  tibble(
    check = paste0("figure_exists:", required_figures),
    pass = file.exists(file.path(paths$figures, required_figures)),
    detail = file.path(paths$figures, required_figures)
  )
)

decomp <- read_csv(file.path(paths$tables, "Table_construct_decomposition.csv"), show_col_types = FALSE)
transport <- read_csv(file.path(paths$tables, "Table_cross_cohort_transportability.csv"), show_col_types = FALSE)
checks <- bind_rows(
  checks,
  tibble(check = "construct_rows_A_to_D", pass = identical(decomp$analysis, c("A", "B", "C", "D")), detail = paste(decomp$analysis, collapse = ",")),
  tibble(check = "construct_all_primary_gates", pass = all(decomp$support_pass & decomp$balance_pass & decomp$positivity_pass), detail = "support, balance, positivity"),
  tibble(check = "v5_HRS_full_anchor_reproduced", pass = abs(decomp$rd_early_persistent[decomp$analysis == "A"] - 0.00656235621) < 1e-8, detail = format(decomp$rd_early_persistent[decomp$analysis == "A"], digits = 12)),
  tibble(check = "CLHLS_external_anchor_reproduced", pass = abs(transport$rd_early_persistent[transport$cohort == "CLHLS"] + 0.02446135504) < 1e-8, detail = format(transport$rd_early_persistent[transport$cohort == "CLHLS"], digits = 12)),
  tibble(check = "no_country_causal_label", pass = all(grepl("not a health-system causal effect", transport$interpretation)), detail = unique(transport$interpretation))
)

checks <- checks %>% mutate(checked_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
write_csv(checks, file.path(paths$audit, "validation_manifest.csv"))
capture.output(sessionInfo(), file = file.path(paths$audit, "sessionInfo.txt"))

if (!all(checks$pass)) {
  print(checks %>% filter(!pass))
  stop("Validation failed.")
}
cat("All ", nrow(checks), " validation checks passed.\n", sep = "")
