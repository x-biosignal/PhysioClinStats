test_that("defineEstimand validates the strategy and round-trips through the slot", {
  e <- defineEstimand("A vs B", "post-op", "gait",
                      intercurrent = list(event = "discontinuation",
                                          strategy = "hypothetical"),
                      summary_measure = "difference in means")
  expect_s3_class(e, "estimand")
  expect_equal(e$strategy, "hypothetical")
  expect_true(!is.null(e$recipe))
  # invalid strategy is rejected
  expect_error(defineEstimand("A", "P", "E",
                              intercurrent = list(event = "d", strategy = "made-up")),
               "invalid.*strategy")
  # round-trips through the PhysioCore AnalysisResult estimand slot
  ar <- PhysioCore::AnalysisResult(type = "x", estimand = e)
  rt <- PhysioCore::estimandOf(ar)
  expect_equal(rt$strategy, e$strategy)
  expect_equal(rt$endpoint, e$endpoint)
})

test_that("poolEstimates reproduces mice::pool on the nhanes example", {
  skip_if_not_installed("mice")
  suppressWarnings(suppressMessages({
    imp <- mice::mice(mice::nhanes, m = 5, seed = 23109, printFlag = FALSE)
    fit <- with(imp, lm(bmi ~ age + hyp + chl))
  }))
  mp <- summary(mice::pool(fit))
  our <- poolEstimates(fit$analyses)@result$estimates
  expect_equal(as.character(our$term), as.character(mp$term))
  expect_lt(max(abs(our$estimate - mp$estimate)), 1e-6)
  expect_lt(max(abs(our$std.error - mp$std.error)), 1e-6)
  expect_lt(max(abs(our$df - mp$df)), 1e-6)
  expect_lt(max(abs(our$p.value - mp$p.value)), 1e-6)
})

test_that("poolEstimates matches mice::pool for glm fits and is NA-robust", {
  skip_if_not_installed("mice")
  suppressWarnings(suppressMessages({
    imp <- mice::mice(mice::nhanes2, m = 6, seed = 55, printFlag = FALSE)
    fit <- with(imp, glm(hyp ~ age + bmi, family = binomial))
  }))
  mp <- summary(mice::pool(fit))
  our <- suppressWarnings(poolEstimates(fit$analyses)@result$estimates)
  expect_lt(max(abs(our$df - mp$df)), 1e-6)               # no lambda-floor drift
  expect_lt(max(abs(our$estimate - mp$estimate)), 1e-6)
  # an aliased (rank-deficient) coefficient yields an NA row, not a crash
  dup <- mtcars; dup$hp2 <- dup$hp
  f <- lm(mpg ~ hp + hp2, data = dup)
  p <- poolEstimates(list(f, f))@result$estimates
  expect_true(is.na(p$estimate[p$term == "hp2"]))
  expect_true(is.finite(p$estimate[p$term == "hp"]))
})

test_that("poolEstimates rejects a tidy data frame missing std.error clearly", {
  td <- data.frame(term = c("a", "b"), estimate = c(1, 2))
  expect_error(poolEstimates(list(td, td)), "missing required column")
})

test_that("multipleImputation is reproducible and honours m/method; reference-based needs rbmi", {
  skip_if_not_installed("mice")
  a <- multipleImputation(mice::nhanes, m = 3, method = "pmm", seed = 42)
  b <- multipleImputation(mice::nhanes, m = 3, method = "pmm", seed = 42)
  expect_identical(mice::complete(a, 1), mice::complete(b, 1))
  expect_identical(mice::complete(a, 3), mice::complete(b, 3))
  expect_equal(a$m, 3L)
  expect_equal(unname(a$method[["bmi"]]), "pmm")
  expect_error(multipleImputation(mice::nhanes, reference_based = TRUE), "rbmi")
  expect_error(multipleImputation(mice::nhanes, m = 0), "positive integer")
})

test_that("analyseEstimand recovers a known treatment effect under MAR dropout", {
  skip_if_not_installed("mice")
  skip_on_cran()
  sim <- function(n, delta, seed) {
    set.seed(seed); arm <- rep(c("A", "B"), each = n / 2); bi <- rnorm(n, 0, 2)
    d <- do.call(rbind, lapply(seq_len(n), function(i) {
      t <- 1:4
      data.frame(subject = i, time = t, arm = arm[i],
                 y = 10 + bi[i] + delta * (arm[i] == "B") * (t - 1) + rnorm(4, 0, 1.5))
    }))
    for (i in unique(d$subject)) {
      rows <- which(d$subject == i)
      for (k in 2:4) {
        prev <- d$y[rows[k - 1]]
        if (!is.na(prev) && stats::plogis(-3.5 + 0.15 * prev) > runif(1)) {
          d$y[rows[k:4]] <- NA; break
        }
      }
    }
    d
  }
  e <- defineEstimand("B vs A", "sim", "y at t4",
                      intercurrent = list(event = "dropout", strategy = "hypothetical"))
  d <- sim(120, 1.0, seed = 42)                          # true t4 effect = 3.0
  res <- suppressWarnings(analyseEstimand(d, e, "y", "arm", "time", "subject",
                                          m = 15, seed = 1))
  pooled <- res@result$pooled
  row <- pooled[pooled$term == ".trtB:.time4", ]
  # the pooled estimate is close to the truth and its CI covers it
  expect_lt(abs(row$estimate - 3.0), 0.8)
  expect_true(row$conf.low <= 3.0 && row$conf.high >= 3.0)
  expect_equal(PhysioCore::estimandOf(res)$strategy, "hypothetical")
})
