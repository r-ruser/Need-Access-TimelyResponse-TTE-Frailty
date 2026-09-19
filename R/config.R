options(stringsAsFactors = FALSE, scipen = 999)

repo_root <- normalizePath(
  Sys.getenv("NATR_REPO_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
output_root <- normalizePath(
  Sys.getenv("NATR_OUTPUT_ROOT", unset = file.path(repo_root, "outputs")),
  winslash = "/", mustWork = FALSE
)

input_paths <- list(
  hrs_candidates = Sys.getenv(
    "NATR_HRS_CANDIDATES",
    unset = file.path(repo_root, "data", "02_all_person_trial_candidates.rds")
  ),
  clhls_candidates = Sys.getenv(
    "NATR_CLHLS_CANDIDATES",
    unset = file.path(repo_root, "data", "02_clhls_all_trial_candidates.rds")
  )
)

paths <- list(
  root = output_root,
  tables = file.path(output_root, "output", "tables"),
  audit = file.path(output_root, "output", "audit"),
  models = file.path(output_root, "output", "models"),
  figures = file.path(output_root, "figures"),
  source = file.path(output_root, "source_data"),
  logs = file.path(output_root, "logs")
)
invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))

required_packages <- c(
  "dplyr", "tidyr", "purrr", "readr", "survey", "sandwich",
  "ggplot2", "patchwork", "scales", "svglite", "ragg"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

if (!file.exists(input_paths$hrs_candidates)) {
  stop("HRS derived input not found. Set NATR_HRS_CANDIDATES.")
}
if (!file.exists(input_paths$clhls_candidates)) {
  stop("CLHLS derived input not found. Set NATR_CLHLS_CANDIDATES.")
}

analysis_seed <- 20260919L
set.seed(analysis_seed)

