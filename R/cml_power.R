##############################################################################
# UMIT Tirol -  Private University for Health Sciences and Health Technology
#   Institute of Psychology, Statistics and Psychometrics Working Group
#
# Power and Power Curve Functions
#
# Part of R/tlc - Testing in Conditional Likelihood Context package
#
# This file contains a routine that computes Functions to compute power
# for statistical tests and plot power curves as functions of
# effect size and sample size.
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2025, Last Modified 16/09/2025
######################################################################
#' @name cml_power
#' @aliases cml_power p_curve p_ncurve
#'
#' @title Power and Power Curve Functions
#'
#' @description Functions to compute the power of \eqn{\chi^2} tests, i.e., Wald (W), likelihood ratio (LR), Rao score (RS)
#' and gradient (GR) test, and to plot power curves as functions of effect size and sample size.
#'
#' @details
#'  The effect is interpreted as a pseudo \eqn{R^2}-like measure of explained variance (as in linear models).
#' It is 0 when persons with the same person parameters yield the same response probabilities.
#' If two persons with the same person parameter but different covariate values yield different response probabilities,
#' an additional variance component is introduced and the effect is greater than 0
#'
#' The power of the tests is computed from the cumulative distribution function of the non-central \eqn{\chi^2} distribution,
#' where the respective non-centrality parameter is obtained by multiplying the effect with the informative sample size.
#' This is only an approximation based on results of asymptotic theory. The approximation may be poor when the
#' informative sample size is small and/or the effect is large.
#' @export

#' @rdname cml_power
#' @title Power Function
#' @description \code{cml_power()} computes the power of the tests given a specified effect size, type I error prob. alpha,
#' informative sample size, and degrees of freedom.
#'
#' @param obj An object of class 'tcl_sa_size', typically containing information such as degrees of freedom (\code{df}) and
#' informative sample size (\code{n}). If missing, values for \code{df} and \code{n} need to be set manually.
#' @param effect Numeric value representing the effect size. A real number between 0 and 1, interpreted as a proportion of
#' pseudo-variance between persons with different covariate values (but the same person parameter). Default is 0.03.
#' @param alpha Type I error probability. Default is 0.05.
#' @param n Informative sample size (excluding persons with a score of 0 or highest possible score).
#' Default is \code{"auto"}, in which case the value is extracted from \code{obj}.
#' @param df Degrees of freedom. Default is \code{"auto"}, in which case the value is extracted from \code{obj}.
#'
# #' @return The power of the tests as a numeric value.
#' @return
#' \itemize{
#'   \item \code{cml_power()}: Numeric vector of power values.
#'   \item \code{p_curve()}, \code{p_ncurve()}: A power curve plotted to the active graphics device.
#' }
#'
#' @references{
#' Draxler, C., & Kurz, A. (2025). Testing measurement invariance in a conditional likelihood framework by considering
#' multiple covariates simultaneously. Behavior Research Methods, 57(1), 50.
#'   }
#' @export

cml_power <- function(obj, effect = 0.03, alpha = 0.05, n = "auto", df = "auto") {
  # Check if df is set to "auto" and obj is missing
  if (is.character(df) && df == "auto" || is.character(n) && n == "auto") {
    if (missing(obj)) {
      stop("The argument 'obj' is required when 'df' or 'n' is set to 'auto'.")
    }
    # Fetch degrees of freedom using a helper function
    if (is.character(df) && df == "auto") {
      df <- get_tcl_arg(obj = obj, arg = "df")
      message("Degrees of freedom (df) extracted from obj: ", df)
    }
    if (is.character(n) && n == "auto") {
      n <- get_tcl_arg(obj =  obj, arg = "sample_size_informative")
      message("Informative sample size extracted from obj: ",
              paste(names(n), n, sep = " = ", collapse = ", "))
    }

  }
  q <- qchisq(1 - alpha, df)
  p <- 1 - pchisq(q, df, effect * n)
  return(p)
}

#' @rdname cml_power
#' @title Power Curve Function
#' @description \code{p_curve()} generates a power curve as a function of effect size.
#'
#' @param from Lower bound of the effect or sample size range (default is 0).
#' @param to Upper bound of the effect or sample size range (default is 0.2 for effect size, and 600 for sample site).
#' @param ... Additional graphical arguments passed to \code{\link[graphics]{plot}} (e.g., \code{col}, \code{lwd}, \code{ylim}) via \code{\link{p_curve}} and \code{\link{p_ncurve}}.
#'
#'
# #' @return A plot of the power curve as a function of effect size.
#' @export

p_curve <- function(obj, alpha = 0.05, n = 300, df = "auto", from = 0, to = 0.2,...) {
  p <- function(effect) {
    q <- qchisq(1 - alpha, df)
    p <- 1 - pchisq(q, df, effect * n)
    return(p)
  }
  # Check if df is set to "auto" and obj is missing
  if (is.character(df) && df == "auto") {
    if (missing(obj)) {
      stop("The argument 'obj' is required when 'df' is set to 'auto'.")
    }
    # Fetch degrees of freedom using a helper function
    df <- get_tcl_arg(obj = obj, arg = "df")
    message("Degrees of freedom (df) extracted from obj:  ", df)
  }

 curve(
   expr = p,
   from = from,
   to = to,
   xlab = 'Effect',
   ylab = 'Power',
   main = "Power curve as a function of effect size",
   ...
   )
}

#' @rdname cml_power
#' @title Power Sample Size Curve Function
#' @description \code{p_ncurve()} generates a power curve as a function of sample size.
#'
# #' @return A plot of the power curve as a function of sample size.
#' @seealso \code{\link{sa_sizeRM}}, \code{\link{sa_sizePCM}}, and \code{\link{sa_sizeChange}}
#' @keywords sample_size_planning
#' @examples
#' \dontrun{
#'##### Sample size of Rasch Model #####
#'
#' res <-  sa_sizeRM(local_dev = list( c(0, -0.5, 0, 0.5, 1) , c(0, 0.5, 0, -0.5, 1)))
#'
#' cml_power(obj = res)
#'
#' p_curve(obj = res)
#' p_curve(obj = res, col = "red", lwd = 2, ylim = c(0, 1))
#'
#' p_ncurve(obj = res)
#' p_ncurve(obj = res, col = "red", lwd = 2, ylim = c(0, 1))
#'
#'}
#' @export
p_ncurve <- function(obj, effect = 0.03, alpha = 0.05, df = "auto",  from = 0, to = 600,...) {
  # Check if df is set to "auto" and obj is missing
  if (is.character(df) && df == "auto") {
    if (missing(obj)) {
      stop("The argument 'obj' is required when 'df' is set to 'auto'.")
    }
    # Fetch degrees of freedom using a helper function
    df <- get_tcl_arg(obj, arg = "df")
    message("Degrees of freedom (df) extracted from obj: ", df)
  }

  # Define the power function as a function of sample size n
  p <- function(n) {
    q <- qchisq(1 - alpha, df)  # Critical value of chi-squared distribution
    1 - pchisq(q, df, effect * n)  # Power calculation
  }

  # Plot the power curve
  curve(
    expr = p,
    from = from,
    to = to,
    xlab = "Informative sample size",
    ylab = "Power",
    main = "Power curve as a function of sample size",
    ...
  )
}
