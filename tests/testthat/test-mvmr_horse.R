
# testthat library must be loaded


## Basic tests ####

testthat::test_that("Basic multivariable JAGS model runs successfully and produces expected results", {
  set.seed(100)
  # Run the model
  result = mvmr_horse(D = data_mv_ex, n.iter = 1000, n.burnin = 500)

  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  # Check estimates match expected
  expect_equal(result$MR_Estimate, data.frame(
    Parameter = c("theta[1]", "theta[2]"),
    Estimate = c(0.100, 0.104),
    SD = c(0.018, 0.018),
    `2.5% quantile` = c(0.065, 0.069),
    `97.5% quantile` = c(0.135, 0.140),
    Rhat = c(1.007, 1.001),
    check.names = FALSE
  ))
})


testthat::test_that("Basic multivariable Stan model runs successfully and produces expected output type", {
  set.seed(100)
  warnings = capture_warnings({
    result = mvmr_horse(D = data_mv_ex, n.iter = 1000, n.burnin = 500, stan = TRUE)
  })
  expect_length(warnings, 3)
  expect_match(warnings[1], "divergent transitions after warmup")
  expect_match(warnings[2], "transitions after warmup that exceeded the maximum treedepth")
  expect_match(warnings[3], "Examine the pairs")

  # Ensure the model returns a list with expected elements
  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  # Check estimates match expected
  expect_equal(result$MR_Estimate, data.frame(
    Parameter = c("theta[1]", "theta[2]"),
    Estimate = c(0.100, 0.106),
    SD = c(0.019, 0.017),
    `2.5% quantile` = c(0.065, 0.071),
    `97.5% quantile` = c(0.138, 0.140),
    Rhat = c(1.005, 1.000),
    check.names = FALSE
  ))
})
