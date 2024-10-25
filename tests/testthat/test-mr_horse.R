
# testthat library must be loaded

## Basic tests ####

testthat::test_that("Basic JAGS model runs successfully and produces expected results", {
  set.seed(100)
  # Run the model
  result = mr_horse(D = data_ex, n.iter = 1000, n.burnin = 500)

  # Ensure the model returns a list with expected elements
  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  # Check estimates match expected
  expect_equal(result$MR_Estimate, data.frame("Estimate"=0.098, "SD"=0.018, "2.5% quantile"=0.063, "97.5% quantile"=0.135, "Rhat"=1, check.names = FALSE))
})


testthat::test_that("Basic Stan model runs successfully and produces expected output type", {
  set.seed(100)
  warnings = capture_warnings({
    result = mr_horse(D = data_ex, n.iter = 1000, n.burnin = 500, stan = TRUE)
  })
  expect_length(warnings, 4)
  expect_match(warnings[1], "divergent transitions after warmup")
  expect_match(warnings[2], "transitions after warmup that exceeded the maximum treedepth")
  expect_match(warnings[3], "Examine the pairs")
  expect_match(warnings[4], "Tail Effective Samples Size")

  # Ensure the model returns a list with expected elements
  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  # Check estimates match expected
  expect_equal(result$MR_Estimate, data.frame("Estimate"=0.098, "SD"=0.018, "2.5% quantile"=0.064, "97.5% quantile"=0.133, "Rhat"=1.009, check.names = FALSE))
})


## Fixed tau value tests ####

testthat::test_that("JAGS model runs successfully and produces expected results with fixed tau value", {
  set.seed(100)
  # Run the model
  result = mr_horse(D = data_ex, n.iter = 1000, n.burnin = 500, fixed_tau = 0.01, variable.names = 'tau')

  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  # Check estimates match expected
  expect_equal(summary(result$MR_Coda)$statistics['tau',1], 0.01)
  expect_equal(summary(result$MR_Coda)$statistics['theta',1], 0.0976, tolerance = 0.001)

  expect_equal(result$MR_Estimate, data.frame("Estimate"=0.098, "SD"=0.018, "2.5% quantile"=0.062, "97.5% quantile"=0.134, "Rhat"=1.005, check.names = FALSE))
})


testthat::test_that("Stan model runs successfully and produces expected results with fixed tau value", {
  set.seed(100)
  # Run the model
  warnings = capture_warnings({
    result = mr_horse(D = data_ex, n.iter = 1000, n.burnin = 500, stan = TRUE, fixed_tau = 0.01, variable.names = 'tau')
  })
  expect_length(warnings, 5)
  expect_match(warnings[1], "divergent transitions after warmup")
  expect_match(warnings[2], "transitions after warmup that exceeded the maximum treedepth")
  expect_match(warnings[3], "Examine the pairs")
  expect_match(warnings[4], "Bulk Effective Samples Size")
  expect_match(warnings[5], "Tail Effective Samples Size")

  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  # Check estimates match expected
  expect_equal(summary(result$MR_Coda)$statistics['tau',1], 0.01)
  expect_equal(summary(result$MR_Coda)$statistics['theta',1], 0.0976, tolerance = 0.001)
  expect_equal(result$MR_Estimate, data.frame("Estimate"=0.098, "SD"=0.018, "2.5% quantile"=0.064, "97.5% quantile"=0.133, "Rhat"=1.009, check.names = FALSE))
})
