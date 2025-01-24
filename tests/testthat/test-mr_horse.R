
# testthat library must be loaded
library(testthat)
# These tests are not necessarily proving the model produces correct results, but more so to tell
# if the results change if a new feature is added or something is modified


## Basic tests ####

testthat::test_that("Basic JAGS model runs successfully and produces expected results format", {
  # Note this test is checking against Andrew's original results with same seed
  set.seed(20230531)
  # Run the model
  result = mr_horse(D = data_ex)

  # Ensure the model returns a list with expected elements
  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  # Compare to last result
  cat('\nOriginal results:\n')
  print(data.frame("Estimate"=0.097, "SD"=0.018, "2.5% quantile"=0.063, "97.5% quantile"=0.132, "Rhat"=1, check.names = F))
  cat('\nCurrent results:\n')
  print(result$MR_Estimate)
})


testthat::test_that("Basic Stan model runs successfully and produces expected output type", {
  set.seed(20230531)
  warnings = capture_warnings({
    result = mr_horse(D = data_ex, n.iter = 1000, n.burnin = 500, stan = TRUE)
  })

  cat('\nExpected warning types:\n')
  print(c('Divergent transitions', 'maximum treedepth', 'Examine the pairs', 'Bulk effective samples size'))
  cat('\nCurrent warnings:\n')
  print(warnings)

  # Ensure the model returns a list with expected elements
  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  # Compare to original result
  cat('\nOriginal results (JAGS):\n')
  print(data.frame("Estimate"=0.097, "SD"=0.018, "2.5% quantile"=0.063, "97.5% quantile"=0.132, "Rhat"=1, check.names = F))
  cat('\nCurrent results:\n')
  print(result$MR_Estimate)
  # large rhat likely due to low iterations
})


## Test for MRInput ####
testthat::test_that("Basic JAGS model runs successfully and produces expected output type with MRInput", {
  data = MendelianRandomization::mr_input(bx = data_ex$betaX, by = data_ex$betaY,
                                        bxse = data_ex$betaXse, byse = data_ex$betaYse)
  set.seed(100)
  # Run the model
  result = mr_horse(D = data, n.iter = 1000, n.burnin = 500)

  # Ensure the model returns a list with expected elements
  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  # Compare to past result
  cat('\nOriginal results:\n')
  print(data.frame("Estimate"=0.098, "SD"=0.018, "2.5% quantile"=0.064, "97.5% quantile"=0.134, "Rhat"=1.002, check.names = F))
  cat('\nCurrent results:\n')
  print(result$MR_Estimate)
})



## Fixed tau value tests ####
testthat::test_that("JAGS model runs successfully and produces expected results type with fixed tau value", {
  set.seed(100)
  # Run the model
  result = mr_horse(D = data_ex, n.iter = 1000, n.burnin = 500, fixed_tau = 0.01, variable.names = 'tau')

  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  # Check estimates match expected
  expect_equal(summary(result$MR_Coda)$statistics['tau',1], 0.01)
  expect_equal(summary(result$MR_Coda)$statistics['tau',2], 0.00)

  # Compare to past result
  cat('\nOriginal results:\n')
  print(data.frame("Estimate"=0.098, "SD"=0.018, "2.5% quantile"=0.062, "97.5% quantile"=0.132, "Rhat"=1.002, check.names = F))
  cat('\nCurrent results:\n')
  print(result$MR_Estimate)
})


testthat::test_that("Stan model runs successfully and produces expected results with fixed tau value", {
  set.seed(100)
  # Run the model
  warnings = capture_warnings({
    result = mr_horse(D = data_ex, n.iter = 1000, n.burnin = 1000, stan = TRUE, fixed_tau = 0.01, variable.names = 'tau')
  })

  cat('\nExpected warning types:\n')
  print(c('Divergent transitions', 'maximum treedepth', 'Examine the pairs', 'Bulk effective samples size', 'Tail effective samples size'))
  cat('\nCurrent warnings:\n')
  print(warnings)

  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  # Check estimates match expected
  expect_equal(summary(result$MR_Coda)$statistics['tau',1], 0.01)
  expect_equal(summary(result$MR_Coda)$statistics['tau',2], 0.00)

  # Compare to past result
  cat('\nOriginal results:\n')
  print(data.frame("Estimate"=0.098, "SD"=0.017, "2.5% quantile"=0.065, "97.5% quantile"=0.131, "Rhat"=1.001, check.names = F))
  cat('\nCurrent results:\n')
  print(result$MR_Estimate)
})

## Sample overlap with non zero omega tests ######

testthat::test_that("JAGS model runs successfully and produces expected results nonzero omega", {
  set.seed(100)
  # Run the model
  result = mr_horse(D = data_ex, n.iter = 1000, n.burnin = 500, omega = 0.1)

  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  # Compare to past result
  cat('\nOriginal results:\n')
  print(data.frame("Estimate"=0.096, "SD"=0.017, "2.5% quantile"=0.062, "97.5% quantile"=0.129, "Rhat"=1.031, check.names = F))
  cat('\nCurrent results:\n')
  print(result$MR_Estimate)
})


testthat::test_that("Basic Stan model runs successfully and produces expected output type with nonzero omega", {
  set.seed(100)
  warnings = capture_warnings({
    result = mr_horse(D = data_ex, n.iter = 1000, n.burnin = 500, stan = TRUE, omega = 0.1)
  })

  cat('\nExpected warning types:\n')
  print(c('Divergent transitions', 'maximum treedepth', 'Examine the pairs', 'Bulk effective samples size'))
  cat('\nCurrent warnings:\n')
  print(warnings)

  # Ensure the model returns a list with expected elements
  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  # Compare to past result
  cat('\nOriginal results:\n')
  print(data.frame("Estimate"=0.097, "SD"=0.018, "2.5% quantile"=0.061, "97.5% quantile"=0.133, "Rhat"=1.000, check.names = F))

  cat('\nCurrent results:\n')
  print(result$MR_Estimate)
})

testthat::test_that("JAGS model runs successfully and produces expected results with fixed tau value and nonzero omega", {
  set.seed(100)
  # Run the model
  result = mr_horse(D = data_ex, n.iter = 1000, n.burnin = 500, fixed_tau = 0.01, variable.names = 'tau', omega = 0.1)

  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  # Check estimates match expected
  expect_equal(summary(result$MR_Coda)$statistics['tau',1], 0.01)
  expect_equal(summary(result$MR_Coda)$statistics['tau',2], 0.00)

  # Compare to past result
  cat('\nOriginal results:\n')
  print(data.frame("Estimate"=0.096, "SD"=0.018, "2.5% quantile"=0.060, "97.5% quantile"=0.131, "Rhat"=1.025, check.names = F))

  cat('\nCurrent results:\n')
  print(result$MR_Estimate)
})


