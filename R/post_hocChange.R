######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# post_hocChange
#
# Part of R/tlc - Testing in Conditional Likelihood Context package
#
# This file contains a routine that returns post hoc power of W, LR,
# RS, and  GR Test of hypothesis about parameter quantifying change
# between two time points as modeled by an LLTM given data and
# probability of error of first kind alpha
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2021, Last Modified 19/04/2021
######################################################################
#' Power analysis of tests in context of measurement of change using LLTM
#'
#' Returns post hoc power of Wald (W), likelihood ratio (LR), Rao score (RS)
#' and gradient (GR) test given data and probability of error of first kind \eqn{\alpha}.
#' The hypothesis to be tested states that the shift parameter quantifying the constant change
#' for all items between time points 1 and 2 equals 0. The alternative states that the
#' shift parameter is not equal to 0. It is assumed that the same items are presented at both
#' time points. See function \code{\link{change_test}}.
#'
#'The power of the tests (Wald, LR, score, and gradient) is determined from the assumption that the
#'approximate distributions of the four test statistics are from the family of noncentral \eqn{\chi^2}
#'distributions with \eqn{df = 1} and noncentrality parameter \eqn{\lambda}. In case of evaluating the post hoc power,
#'\eqn{\lambda} is assumed to be given by the observed value of the test statistic. Given the probability of the
#'error of the first kind \eqn{\alpha} the post hoc power of the tests can be determined from \eqn{\lambda}.
#'More details about the distributions of the test statistics and the relationship between \eqn{\lambda}, power, and
#'sample size can be found in Draxler and Alexandrowicz (2015).
#'
#'In particular, let \eqn{q_{\alpha}} be the \eqn{1- \alpha} quantile of the central \eqn{\chi^2} distribution with df = 1. Then,
#'
#'\deqn{power = 1 - F_{df, \lambda} (q_{\alpha}),}
#'
#'where \eqn{F_{df, \lambda}} is the cumulative distribution function of the noncentral \eqn{\chi^2} distribution with
#'\eqn{df = 1} and \eqn{\lambda} equal to the observed value of the test statistic.
#'
#'
#' @param alpha Probability of error of first kind.
#' @param data Data matrix as required for function \code{\link{change_test}}.
#'
#'@return A list of results.
#' \item{test}{A numeric vector of Wald (W), likelihood ratio (LR), Rao score (RS), and gradient (GR)  test statistics.}
#'  \item{power}{Posthoc power value for each test.}
#'  \item{observed deviation}{CML estimate of shift parameter expressing observed deviation from hypothesis to be tested.}
#'  \item{person score distribution}{Relative frequencies of person scores. Uninformative scores, i.e., minimum and maximum score,
#'  are omitted. Note that the person score distribution does also have an influence on the power of the tests.}
#'  \item{degrees of freedom}{Degrees of freedom \eqn{df}.}
#'  \item{noncentrality parameter}{Noncentrality parameter \eqn{\lambda} of \eqn{\chi^2} distribution from which power is determined.
#'  It equals observed value of test statistic.}
#'  \item{call}{The matched call.}
#'
#' @references{
#'  Draxler, C., & Alexandrowicz, R. W. (2015). Sample size determination within the scope of conditional
#'  maximum likelihood estimation with special focus on testing the Rasch model. Psychometrika, 80(4), 897-919.
#'
#'  Fischer, G. H. (1995). The Linear Logistic Test Model. In G. H. Fischer & I. W. Molenaar (Eds.),
#'  Rasch models: Foundations, Recent Developments, and Applications (pp. 131-155). New York: Springer.
#'
#'  Fischer, G. H. (1983). Logistic Latent Trait Models with Linear Constraints. Psychometrika, 48(1), 3-26.
#'
#'  }
#' @keywords sample_size_planning
#' @export
#' @seealso \code{\link{sa_sizeChange}}, and \code{\link{powerChange}}.
#' @examples
#' \dontrun{
#'# Numerical example with 200 persons and 4 items
#'# presented twice, thus 8 virtual items
#'
#'# Data y generated under the assumption that shift parameter equals 0.5
#'# (change from time point 1 to 2)
#'
#'# design matrix W used only for exmaple data generation
#'#     (not used for estimating in change_test function)
#' W <- rbind(c(1,0,0,0,0), c(0,1,0,0,0), c(0,0,1,0,0), c(0,0,0,1,0),
#'            c(1,0,0,0,1), c(0,1,0,0,1), c(0,0,1,0,1), c(0,0,0,1,1))
#'
#'# eta parameter vector, first 4 are nuisance, i.e., item parameters at time point 1.
#'# (easiness parameters of the 4 items at time point 1),
#'# last one is the shift parameter
#' eta <- c(-2,-1,1,2,0.5)
#'
#' y <- eRm::sim.rasch(persons=rnorm(150), items=colSums(-eta*t(W)))
#'
#' res <- post_hocChange(alpha = 0.05, data = y)
#'
#'# > res
#'# $test
#'#     W     LR     RS     GR
#'# 9.822 10.021  9.955 10.088
#'#
#'# $power
#'#     W    LR    RS    GR
#'# 0.880 0.886 0.884 0.888
#'#
#'# $`observed deviation (estimate of shift parameter)`
#'# [1] 0.504
#'#
#'# $`person score distribution`
#'#
#'#     1     2     3     4     5     6     7
#'# 0.047 0.047 0.236 0.277 0.236 0.108 0.047
#'#
#'# $`degrees of freedom`
#'# [1] 1
#'#
#'# $`noncentrality parameter`
#'#     W     LR     RS     GR
#'# 9.822 10.021  9.955 10.088
#'#
#'# $call
#'# post_hocChange(alpha = 0.05, data = y)
#' }




post_hocChange <- function(alpha = 0.05, data){

  subfunc <- function(stats) {
    nc <- stats
    beta <- pchisq(qchisq(1-alpha, 1), 1, nc)
    power <- 1 - beta
    return(round(power, digits = 3))
  }

  call<-match.call()

  t <- table(factor(rowSums(data),levels=1:(ncol(data)-1)))

  est <- eRm::LLTM(X=data, mpoints = 2,se = TRUE, sum0 = FALSE) # unrestricted CML estimates of eta Parameters
  dev_obs <- unname(est$etapar[ncol(data)/2])

  vstats <- change_test(X = data)$test[c(4,2,3,1)] # W, LR, RS, GR test

  res <- lapply(vstats, subfunc)

  results <- list(  'test' = round(vstats, digits = 3),
                    'power' = do.call(c, res),
                    'observed deviation (estimate of shift parameter)' = round(dev_obs, digits = 3),
                    'person score distribution' = round(t/sum(t), digits = 3),
                    'degrees of freedom' = 1,
                    'noncentrality parameter' = round(vstats, digits = 3),
                    "call" = call)
  return(results)
}
