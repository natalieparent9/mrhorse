
# testthat library must be loaded
library(testthat)

## Basic tests ####

testthat::test_that("Basic multivariable JAGS model runs successfully and produces expected results", {
  # Note this test is checking against Andrew's original results - slightly different
  set.seed(20230531)
  # Run the model
  result = mvmr_horse(D = data_mv_ex)

  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  # Compare to last result
  cat('\nOriginal results:\n')
  print(data.frame(
      Parameter = c("theta[1]", "theta[2]"),
      Estimate = c(0.100, 0.105),
      SD = c(0.019, 0.018),
      `2.5% quantile` = c(0.065, 0.069),
      `97.5% quantile` = c(0.137, 0.139),
      Rhat = c(1.000, 1.001),
      check.names = FALSE
      ))
  cat('\nCurrent results:\n')
  print(result$MR_Estimate)
})


testthat::test_that("Basic multivariable Stan model runs successfully and produces expected output type", {
  set.seed(20230531)
  warnings = capture_warnings({
    result = mvmr_horse(D = data_mv_ex, n.iter = 1000, n.burnin = 500, stan = TRUE)
  })

  cat('\nExpected warning types:\n')
  print(c('Divergent transitions', 'maximum treedepth', 'Examine the pairs'))
  cat('\nCurrent warnings:\n')
  print(warnings)

  # Ensure the model returns a list with expected elements
  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  # Compare to last result
  cat('\nOriginal results (JAGS):\n')
  print(data.frame(
    Parameter = c("theta[1]", "theta[2]"),
    Estimate = c(0.100, 0.105),
    SD = c(0.019, 0.018),
    `2.5% quantile` = c(0.065, 0.069),
    `97.5% quantile` = c(0.137, 0.139),
    Rhat = c(1.000, 1.001),
    check.names = FALSE
    )
  )
  cat('\nCurrent results:\n')
  print(result$MVMR_Estimate)
})


## Test for MRMVInput ####

testthat::test_that("Basic multivariable JAGS model with MVMR_Input object runs successfully and produces expected results", {
  data = MendelianRandomization::mr_mvinput(bx = as.matrix(data_mv_ex[,c('betaX1', 'betaX2')], ncol=2), by = data_mv_ex$betaY,
                                          bxse = as.matrix(data_mv_ex[,c('betaX1se', 'betaX2se')], ncol=2), byse = data_mv_ex$betaYse)

  set.seed(100)
  # Run the model
  result = mvmr_horse(D = data, n.iter = 1000, n.burnin = 500)

  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  # Compare to last result
  cat('\nOriginal results:\n')
  print(data.frame(
    Parameter = c("theta[1]", "theta[2]"),
    Estimate = c(0.100, 0.105),
    SD = c(0.019, 0.018),
    `2.5% quantile` = c(0.064, 0.069),
    `97.5% quantile` = c(0.136, 0.140),
    Rhat = c(1.027, 1.031),
    check.names = FALSE
  ))
  cat('\nCurrent results:\n')
  print(result$MR_Estimate)
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

  # Compare to last result
  cat('\nOriginal results:\n')
  print(data.frame(
    Parameter = c("theta[1]", "theta[2]"),
    Estimate = c(0.100, 0.103),
    SD = c(0.019, 0.018),
    `2.5% quantile` = c(0.062, 0.068),
    `97.5% quantile` = c(0.135, 0.139),
    Rhat = c(1.008, 1.008),
    check.names = FALSE
  ))
  cat('\nCurrent results:\n')
  print(result$MR_Estimate)
})


testthat::test_that("Multivariable Stan model with fixed tau runs successfully and produces expected output type", {
  set.seed(100)
  warnings = capture_warnings({
    result = mvmr_horse(D = data_mv_ex, n.iter = 1000, n.burnin = 500, stan = TRUE, variable.names = 'tau', fixed_tau = 0.01)
  })

  cat('\nExpected warning types:\n')
  print(c('Divergent transitions', 'maximum treedepth', 'Examine the pairs'))
  cat('\nCurrent warnings:\n')
  print(warnings)

  # Ensure the model returns a list with expected elements
  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  expect_equal(summary(result$MR_Coda)$statistics['tau',1], 0.01)
  expect_equal(summary(result$MR_Coda)$statistics['tau',2], 0.00)

  # Compare to last result
  cat('\nOriginal results:\n')
  print(data.frame(
    Parameter = c("theta[1]", "theta[2]"),
    Estimate = c(0.099, 0.104),
    SD = c(0.019, 0.017),
    `2.5% quantile` = c(0.063, 0.069),
    `97.5% quantile` = c(0.136, 0.139),
    Rhat = c(1.005, 1.01),
    check.names = FALSE
  ))
  cat('\nCurrent results:\n')
  print(result$MR_Estimate)
})

## Non zero omega tests ####

testthat::test_that("Multivariable JAGS model runs successfully and produces expected results with non zero omega", {
  set.seed(100)
  # Run the model
  result = mvmr_horse(D = data_mv_ex, n.iter = 1000, n.burnin = 500, omega = 0.2)

  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  # Compare to last result
  cat('\nOriginal results:\n')
  print(data.frame(
    Parameter = c("theta[1]", "theta[2]"),
    Estimate = c(0.098, 0.103),
    SD = c(0.017, 0.017),
    `2.5% quantile` = c(0.065, 0.069),
    `97.5% quantile` = c(0.133, 0.137),
    Rhat = c(1.008, 1.005),
    check.names = FALSE
  ))
  cat('\nCurrent results:\n')
  print(result$MR_Estimate)
})


testthat::test_that("Multivariable Stan model runs successfully and produces expected output type with non zero omega", {
  set.seed(100)
  warnings = capture_warnings({
    result = mvmr_horse(D = data_mv_ex, n.iter = 1000, n.burnin = 500, stan = TRUE, omega = 0.1)
  })

  cat('\nExpected warning types:\n')
  print(c('Divergent transitions', 'maximum treedepth', 'Examine the pairs', 'Largest R-hat', 'Bulk effective samples size', 'Tail effective samples size'))
  cat('\nCurrent warnings:\n')
  print(warnings)

  # Ensure the model returns a list with expected elements
  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  # Compare to last result
  cat('\nOriginal results:\n')
  print(data.frame(
    Parameter = c("theta[1]", "theta[2]"),
    Estimate = c(0.096, 0.104),
    SD = c(0.022, 0.017),
    `2.5% quantile` = c(0.050, 0.068),
    `97.5% quantile` = c(0.135, 0.139),
    Rhat = c(1.006, 1.01),
    check.names = FALSE
  ))
  cat('\nCurrent results:\n')
  print(result$MR_Estimate)
})


testthat::test_that("Multivariable JAGS model runs successfully and produces expected results with omega matrix", {
  set.seed(100)
  # Run the model
  omega_mat = matrix(c(1, 0.2, 0.2, 0.2, 1 , 0.5, 0.2, 0.5, 1), 3,3)
  result = mvmr_horse(D = data_mv_ex, n.iter = 1000, n.burnin = 500, omega = omega_mat)

  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  # Compare to last result
  cat('\nOriginal results:\n')
  print(data.frame(
    Parameter = c("theta[1]", "theta[2]"),
    Estimate = c(0.099, 0.101),
    SD = c(0.018, 0.017),
    `2.5% quantile` = c(0.064, 0.067),
    `97.5% quantile` = c(0.134, 0.135),
    Rhat = c(1.026, 1.021),
    check.names = FALSE
  ))
  cat('\nCurrent results:\n')
  print(result$MR_Estimate)
})

testthat::test_that("Multivariable JAGS model runs successfully and produces expected results with omega matrix", {
  set.seed(100)
  # Run the model
  omega_mat = matrix(c(1, 0, 0.1, 0, 1 , 0.5, 0.1, 0.5, 1), 3,3)
  result = mvmr_horse(D = data_mv_ex, n.iter = 1000, n.burnin = 500, omega = omega_mat)

  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  # Compare to last result
  cat('\nOriginal results:\n')
  print(data.frame(
    Parameter = c("theta[1]", "theta[2]"),
    Estimate = c(0.099, 0.101),
    SD = c(0.018, 0.017),
    `2.5% quantile` = c(0.064, 0.067),
    `97.5% quantile` = c(0.134, 0.140),
    Rhat = c(1.03, 1.02),
    check.names = FALSE
  ))
  cat('\nCurrent results:\n')
  print(result$MR_Estimate)
})


testthat::test_that("Multivariable Stan model runs successfully and produces expected output type with omega matrix", {
  set.seed(100)
  omega_mat = matrix(c(1, 0.2, 0.2, 0.2, 1 , 0.5, 0.2, 0.5, 1), 3,3)
  warnings = capture_warnings({
    result = mvmr_horse(D = data_mv_ex, n.iter = 1000, n.burnin = 500, stan = TRUE, omega = omega_mat)
  })

  cat('\nExpected warning types:\n')
  print(c('Divergent transitions', 'maximum treedepth', 'Examine the pairs'))
  cat('\nCurrent warnings:\n')
  print(warnings)

  # Ensure the model returns a list with expected elements
  expect_type(result, "list")
  expect_named(result, c("MR_Estimate", "MR_Coda"))
  expect_s3_class(result$MR_Coda, "mcmc.list")

  # Compare to last result
  cat('\nOriginal results:\n')
  print(data.frame(
    Parameter = c("theta[1]", "theta[2]"),
    Estimate = c(0.099, 0.101),
    SD = c(0.017, 0.017),
    `2.5% quantile` = c(0.065, 0.069),
    `97.5% quantile` = c(0.131, 0.132),
    Rhat = c(1.006, 1.001),
    check.names = FALSE
  ))
  cat('\nCurrent results:\n')
  print(result$MR_Estimate)
})
