options(stringsAsFactors = FALSE, scipen = 999)

repo_root <- normalizePath(
  Sys.getenv("NATR_REPO_ROOT", unset = getwd()),
  winslash = "/", mustWork = TRUE
)

source(file.path(repo_root, "R", "config.R"))
source(file.path(repo_root, "R", "ccw_core.R"))

scripts <- c(
  "01_construct_decomposition.R",
  "02_need_components_severity.R",
  "03_access_equity_matrix.R",
  "04_access_stratified_tte.R",
  "05_transport_policy_health.R",
  "06_figures.R",
  "07_validate.R"
)

for (script in scripts) {
  message("Running ", script)
  source(file.path(repo_root, "analysis", script), local = new.env(parent = globalenv()))
}

message("All v6 analyses completed.")
