.sim_exp_panel <- function(ns = 30, seed = 101) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(ns), function(s) {
    A <- 60 + rnorm(1, 0, 5); rate <- 0.30 * exp(rnorm(1, 0, 0.15))
    t <- 0:8
    data.frame(subject = factor(s), time = t,
               y = A * (1 - exp(-rate * t)) + rnorm(9, 0, 2))
  }))
}

test_that("recoveryTrajectoryLME matches a hand-built nlme fit and recovers params", {
  skip_if_not_installed("nlme")
  df <- .sim_exp_panel(30)
  our <- recoveryTrajectoryLME(df, "subject", "time", "y")
  gd <- nlme::groupedData(y ~ t | subj,
                          data = data.frame(subj = df$subject, t = df$time, y = df$y))
  ref <- nlme::nlme(y ~ SSasymp(t, Asym, R0, lrc), data = gd,
                    fixed = Asym + R0 + lrc ~ 1, random = Asym + R0 + lrc ~ 1 | subj,
                    control = nlme::nlmeControl(maxIter = 200, returnObject = TRUE))
  expect_lt(max(abs(our@estimate - nlme::fixef(ref))), 1e-4)
  # recovers the known asymptote (~60) and rate (~0.30) within 5%
  A_hat <- our@estimate[["Asym"]]; r_hat <- exp(our@estimate[["lrc"]])
  expect_lt(abs(A_hat - 60) / 60, 0.05)
  expect_lt(abs(r_hat - 0.30) / 0.30, 0.05)
  expect_equal(nrow(our@result$per_subject), 30L)
  expect_true(all(c("asymptote", "rate", "time_to_90") %in%
                    names(our@result$per_subject)))
})

test_that("recoveryTrajectoryLME validates its input", {
  df <- .sim_exp_panel(4)
  expect_error(recoveryTrajectoryLME(df, "nope", "time", "y"), "not in data")
  one <- df[df$subject == 1, ]
  expect_error(recoveryTrajectoryLME(one, "subject", "time", "y"),
               "at least 2 subjects")
})

test_that("latentClassGrowth recovers a 2-class mixture with high accuracy", {
  skip_if_not_installed("flexmix")
  set.seed(202)
  truth <- rep(c(1L, 2L), c(30, 30))
  df <- do.call(rbind, lapply(seq_along(truth), function(s) {
    sl <- if (truth[s] == 1L) 5 else 0.5
    t <- 0:6; data.frame(subject = s, time = t, y = sl * t + rnorm(7, 0, 1.5))
  }))
  r <- latentClassGrowth(df, "subject", "time", "y", n_classes = 1:3, seed = 1)
  expect_equal(r@result$n_classes, 2L)                   # BIC selects 2
  asg <- r@result$assignment$class
  acc <- max(mean(asg == truth), mean(asg == (3L - truth)))
  expect_gt(acc, 0.8)
  expect_gt(r@result$entropy, 0.8)
})

test_that("proportionalRecoveryRule reproduces ~0.7 and flags the RTM artefact", {
  set.seed(11); maxs <- 66
  init <- runif(80, 5, 55)
  fu <- init + 0.7 * (maxs - init) + rnorm(80, 0, 3)      # genuine 0.7 rule
  g <- proportionalRecoveryRule(init, fu, max_score = maxs, seed = 1)@result
  expect_lt(abs(g$slope - 0.7), 0.05)
  expect_lt(abs(g$followup_on_initial_slope - 0.30), 0.08)  # ~ 1 - 0.7
  expect_false(g$artefact_suspected)

  # a follow-up unrelated to the rule -> the slope is a coupling artefact
  set.seed(22); init2 <- runif(80, 5, 55)
  fu2 <- rnorm(80, 45, 8)
  s <- proportionalRecoveryRule(init2, fu2, max_score = maxs, seed = 2)@result
  expect_true(s$artefact_suspected)
})

test_that("proportionalRecoveryRule mixture diagnostic separates one line from two", {
  skip_if_not_installed("flexmix")
  set.seed(1); maxs <- 66
  init <- runif(80, 5, 55); fu <- init + 0.7 * (maxs - init) + rnorm(80, 0, 3)
  expect_false(proportionalRecoveryRule(init, fu, max_score = maxs, seed = 1)@result$mixture_preferred)
  # fitters (0.7) + non-fitters (0.1)
  set.seed(3); init2 <- runif(90, 5, 50)
  prop <- rep(c(0.7, 0.1), c(50, 40))
  fu2 <- init2 + prop * (maxs - init2) + rnorm(90, 0, 2.5)
  expect_true(proportionalRecoveryRule(init2, fu2, max_score = maxs, seed = 3)@result$mixture_preferred)
})

test_that("proportionalRecoveryRule validates its input", {
  expect_error(proportionalRecoveryRule(1:5, 1:6, 10), "same length")
  expect_error(proportionalRecoveryRule(1:3, 1:3, 10), "at least 5")
  # the RTM control cannot be silently disabled
  expect_error(proportionalRecoveryRule(runif(20, 5, 55), runif(20, 5, 55), 66,
                                        n_shuffle = 0), "positive integer")
})
