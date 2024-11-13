
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
  print(result$MR_Estimate)
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


## Fixed tau tests ####

testthat::test_that("Multivariable JAGS model with fixed tau runs successfully and produces expected results", {
  set.seed(100)
  # Run the model
  result = mvmr_horse(D = data_mv_ex, n.iter = 1000, n.burnin = 500, variable.names = 'tau', fixed_tau = 0.01)

  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  expect_equal(summary(result$MR_Coda)$statistics['tau',1], 0.01)
  expect_equal(summary(result$MR_Coda)$statistics['tau',2], 0.00)

  # Check estimates match expected
  expect_equal(result$MR_Estimate, data.frame(
    Parameter = c("theta[1]", "theta[2]"),
    Estimate = c(0.099, 0.104),
    SD = c(0.018, 0.018),
    `2.5% quantile` = c(0.063, 0.07),
    `97.5% quantile` = c(0.134, 0.138),
    Rhat = c(1.003, 1.000),
    check.names = FALSE
  ))
  print(result$MR_Estimate)
})


testthat::test_that("Multivariable Stan model with fixed tau runs successfully and produces expected output type", {
  set.seed(100)
  warnings = capture_warnings({
    result = mvmr_horse(D = data_mv_ex, n.iter = 1000, n.burnin = 500, stan = TRUE, variable.names = 'tau', fixed_tau = 0.01)
  })
  expect_length(warnings, 3)
  expect_match(warnings[1], "divergent transitions after warmup")
  expect_match(warnings[2], "transitions after warmup that exceeded the maximum treedepth")
  expect_match(warnings[3], "Examine the pairs")

  # Ensure the model returns a list with expected elements
  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  expect_equal(summary(result$MR_Coda)$statistics['tau',1], 0.01)
  expect_equal(summary(result$MR_Coda)$statistics['tau',2], 0.00)

  # Check estimates match expected
  expect_equal(result$MR_Estimate, data.frame(
    Parameter = c("theta[1]", "theta[2]"),
    Estimate = c(0.100, 0.105),
    SD = c(0.018, 0.018),
    `2.5% quantile` = c(0.064, 0.070),
    `97.5% quantile` = c(0.136, 0.142),
    Rhat = c(1.002, 1.000),
    check.names = FALSE
  ))
})

## Non zero omega tests ####

testthat::test_that("Multivariable JAGS model runs successfully and produces expected results with non zero omega", {
  set.seed(100)
  # Run the model
  result = mvmr_horse(D = data_mv_ex, n.iter = 1000, n.burnin = 500, omega = 0.1)

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
  print(result$MR_Estimate)
})
