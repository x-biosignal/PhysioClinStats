.inference_repo <- function() {
  candidates <- c(
    Sys.getenv("PHYSIO_REPO_ROOT"),
    file.path(testthat::test_path(), "..", "..", "..", ".."),
    getwd()
  )
  for (candidate in candidates[nzchar(candidates)]) {
    candidate <- normalizePath(candidate, mustWork = FALSE)
    if (file.exists(file.path(
      candidate, "physio-ecosystem", "validation", "inference", "surface.csv"
    ))) return(candidate)
  }
  NULL
}

.require_inference_repo <- function() {
  repo <- .inference_repo()
  if (is.null(repo)) {
    skip("central inference fixtures are not shipped with the package")
  }
  repo
}

test_that("WS8 PhysioClinStats exports pass the offline parity gate", {
  repo <- .require_inference_repo()
  script <- file.path(
    repo, "physio-ecosystem", "validation", "inference", "run_parity.R"
  )
  output <- tempfile("physioclinstats-parity-")
  log <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      shQuote(script), "--packages", "PhysioClinStats",
      "--output", shQuote(output), "--fail-fast"
    ),
    stdout = TRUE, stderr = TRUE
  )
  expect_null(attr(log, "status"), info = paste(log, collapse = "\n"))
  results <- utils::read.csv(
    file.path(output, "inference_parity.csv"),
    stringsAsFactors = FALSE
  )
  expect_true(nrow(results) >= 80L)
  expect_true(all(results$status %in% c("PASS", "STRUCTURAL")))
})

test_that("survival fixture factor direction is explicit and detectable", {
  repo <- .require_inference_repo()
  root <- file.path(
    repo, "physio-ecosystem", "validation", "inference"
  )
  input <- readRDS(file.path(
    root, "fixtures", "survival", "survival-v1", "input.rds"
  ))
  forward <- PhysioCore::resultValue(coxModel(
    input$survival, time = "time", event = "event", covariates = "group"
  ))$coefficients$hr[1]
  reversed <- input$survival
  reversed$group <- factor(reversed$group, levels = c("B", "A"))
  backward <- PhysioCore::resultValue(coxModel(
    reversed, time = "time", event = "event", covariates = "group"
  ))$coefficients$hr[1]
  expect_equal(forward * backward, 1, tolerance = 1e-10)
  expect_gt(abs(log(forward)), 0.01)
})
