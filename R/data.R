
#' Univariable example dataset
#'
#' Example dataset containing association estimates between 100 genetic instruments and a single exposure variable and an outcome, as well as their corresponding standard errors.
#' This is taken from the first replication of a simulation study, where 20% of variants are pleiotropic, and pleiotropy is balanced.
#'
#' @name data_ex
#' @docType data
#' @format A data frame with 100 rows and 4 columns:
#' \describe{
#'   \item{betaY}{Estimate of genetic association with the outcome}
#'   \item{betaYse}{Standard error of estimate of genetic association with the outcome}
#'   \item{betaX}{Estimate of genetic association with exposure}
#'   \item{betaXse}{Standard error of estimate of genetic association with exposure}
#' }
#'
"data_ex"

#' Multivariable example dataset
#'
#' Example dataset containing associations between 100 genetic instruments and two exposure variables and an outcome, as well as their corresponding standard errors.
#' This is taken from a multivariable simulation study, where 20% of variants are pleiotropic, and pleiotropy is balanced.
#'
#' @name data_mv_ex
#' @docType data
#' @format A data frame with 100 rows and 6 columns:
#' \describe{
#'   \item{betaY}{Estimate of genetic association with the outcome}
#'   \item{betaYse}{Standard error of estimate of genetic association with the outcome}
#'   \item{betaX1}{Estimate of genetic association with exposure 1}
#'   \item{betaX1se}{Standard error of estimate of genetic association with exposure 1}
#'   \item{betaX2}{Estimate of genetic association with exposure 2}
#'   \item{betaX2se}{Standard error of estimate of genetic association with exposure 2}
#' }
#'
"data_mv_ex"
