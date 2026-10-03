test_that("clinStatsBackends reports availability of each optional backend", {
  b <- clinStatsBackends()
  expect_type(b, "logical")
  expect_setequal(names(b), c("mmrm", "emmeans", "lme4", "SingleCaseES"))
})

test_that("requireBackend errors helpfully when a backend is absent", {
  expect_error(requireBackend("a-package-that-does-not-exist"), "optional")
})
