##############################################################################
# UMIT Tirol -  Private University for Health Sciences and Health Technology
#   Institute of Psychology, Statistics and Psychometrics Working Group
#
# tcl_scorefun
#
# Part of R/tlc - Testing in Conditional Likelihood Context package
#
# This file contains a routine that computes the score function
# using function eRm_cml
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2025, Last Modified 04/01/2025
######################################################################
#' Computation of score function.
#'
#' Uses function \code{jacobian()} from numDeriv package to compute (approximate numerically) score function
#'   (first order partial derivatives of conditional log likelihood function)
#'   evaluated at arbitrary values of item easiness parameters.
#'
#' @param data data matrix.
#' @param eta numeric vector of item easiness parameters.
#' @param W design matrix.
#' @param model RM, PCM, RSM, LLTM. Default is set to "RM".
#' @return Score function evaluated at eta.
#'
#' @references{
#' Gilbert, P., Gilbert, M. P., & Varadhan, R. (2016). numDeriv: Accurate Numerical Derivatives. R package
#' version 2016.8-1.1. url: https://CRAN.R-project.org/package=numDeriv
#'  }
#' @keywords htest
#' @export
#' @examples
#' \dontrun{
#' # Rasch model with beta_1 restricted to 0
#' y <- eRm::raschdat1
#' res <- eRm::RM(X = y, sum0 = FALSE)
#' scorefun <- tcl_scorefun(data = y, eta = res$etapar, model = "RM")
#' }

tcl_scorefun <- function(data, eta, W, model = "RM") {
# X = observed data matrix
# eta - numeric vector of item easiness parameters
# W = design matrix
# model = RM, PCM, RSM, LLTM

  if (!is.numeric(eta)) stop("eta needs to be a numeric vector!")

  e <- eRm_cml(X = data, eta = eta, W = W, model = model)$scorefun

  return(e)
}
