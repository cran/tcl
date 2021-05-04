######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# post_hocPCM
#
# Part of R/tlc - Testing in Conditional Likelihood Context package
#
# This file contains a routine that returns post hoc power of W, LR,
# RS and GR given data and probability of error of first kind alpha.
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2021, Last Modified 19/04/2021
######################################################################
#' Power analysis of tests of invariance of item parameters between two groups of
#' persons in partial credit model
#'
#' Returns post hoc power of Wald (W), likelihood ratio (LR), Rao score (RS)
#' and gradient (GR) test  given data and probability of error of first kind \eqn{\alpha}.
#' The hypothesis to be tested assumes equal item-category parameters of the partial
#' credit model between two predetermined groups of persons. The alternative states that
#' at least one of the parameters differs between the two groups.
#'
#' The power of the tests (Wald, LR, score, and gradient) is determined from the assumption
#' that the approximate distributions of the four test statistics are from the family of
#' noncentral \eqn{\chi^2}  distributions with \eqn{df} equal to the number of free item-category
#' parameters in the partial credit model and noncentrality parameter \eqn{\lambda}. In case of evaluating
#' the post hoc power, \eqn{\lambda} is assumed to be given by the observed value of the test statistic.
#' Given the probability of the error of the first kind \eqn{\alpha} the post hoc power of the tests
#' can be determined from \eqn{\lambda}. More details about the distributions of the test statistics and the
#' relationship between \eqn{\lambda}, power, and sample size can be found in Draxler and Alexandrowicz (2015).
#'
#' In particular, let \eqn{q_{\alpha}} be the \eqn{1- \alpha} quantile of the central \eqn{\chi^2} distribution
#' with \eqn{df} equal to the number of free item-category parameters. Then,
#'
#' \deqn{power = 1 - F_{df, \lambda} (q_{\alpha}),}
#'
#' where \eqn{F_{df, \lambda}} is the cumulative distribution function of the noncentral \eqn{\chi^2}
#' distribution with \eqn{df} equal to the number of free item-category parameters and \eqn{\lambda} equal to the
#' observed value of the test statistic.
#'
#' @param alpha Probability of error of first kind.
#' @param data Data matrix with item responses (in ordered categories starting from 0).
#' @param x A numeric vector of length equal to number of persons that contains zeros and ones indicating group membership of the persons.
#'
#'@return A list of results.
#'  \item{test}{A numeric vector of Wald (W), likelihood ratio (LR), Rao score (RS), and gradient (GR)  test statistics.}
#'  \item{power}{Post hoc power value for each test.}
#'  \item{observed global deviation}{Observed global deviation from hypothesis to be tested represented by a single number.
#'  It is obtained by dividing the test statistic by the informative sample size.
#'  The latter does not include persons with minimum or maximum person score.}
#'  \item{observed local deviation}{CML estimates of free item-category parameters in both groups of persons representing observed deviation
#'  from hypothesis to be tested locally per item and response category.}
#'  \item{person score distribution in group 1}{Relative frequencies of person scores in group 1. Uninformative scores, i.e., minimum and maximum score, are omitted.
#'  Note that the person score distribution does also have an influence on the power of the tests.}
#'  \item{person score distribution in group 2}{Relative frequencies of person scores in group 2. Uninformative scores, i.e., minimum and maximum score, are omitted.
#'  Note that the person score distribution does also have an influence on the power of the tests.}
#'  \item{degrees of freedom}{Degrees of freedom \eqn{df}.}
#'  \item{noncentrality parameter}{Noncentrality parameter \eqn{\lambda} of \eqn{\chi^2} distribution from which power is determined.
#'  It equals observed value of test statistic.}
#'  \item{call}{The matched call.}
#'
#'
#' @references{
#' Draxler, C. (2010). Sample Size Determination for Rasch Model Tests. Psychometrika, 75(4), 708–724.
#'
#' Draxler, C., & Alexandrowicz, R. W. (2015). Sample Size Determination Within the Scope of Conditional Maximum Likelihood Estimation
#' with Special Focus on Testing the Rasch Model. Psychometrika, 80(4), 897–919.
#'
#' Draxler, C., Kurz, A., & Lemonte, A. J. (2020). The Gradient Test and its Finite Sample Size Properties in a Conditional Maximum Likelihood
#' and Psychometric Modeling Context. Communications in Statistics-Simulation and Computation, 1-19.
#'
#' Glas, C. A. W., & Verhelst, N. D. (1995a). Testing the Rasch Model. In G. H. Fischer & I. W. Molenaar (Eds.),
#' Rasch Models: Foundations, Recent Developments, and Applications (pp. 69–95). New York: Springer.
#'
#' Glas, C. A. W., & Verhelst, N. D. (1995b). Tests of Fit for Polytomous Rasch Models. In G. H. Fischer & I. W. Molenaar (Eds.),
#' Rasch Models: Foundations, Recent Developments, and Applications (pp. 325-352). New York: Springer.
#'
#'  }
#' @keywords sample_size_planning
#' @export
#' @seealso \code{\link{sa_sizePCM}}, and \code{\link{powerPCM}}.
#' @examples
#' \dontrun{
#'# Numerical example for post hoc power analysis for PCM
#'
#' y <- eRm::pcmdat2
#' n <- nrow(y) # sample size
#' x <- c( rep(0,n/2), rep(1,n/2) ) # binary covariate
#'
#' res <- post_hocPCM(alpha = 0.05, data = y, x = x)
#'
#'# > res
#'# $test
#'#      W     LR     RS     GR
#'# 11.395 11.818 11.628 11.978
#'#
#'# $power
#'#     W    LR    RS    GR
#'# 0.683 0.702 0.694 0.709
#'#
#'# $`observed global deviation`
#'#     W    LR    RS    GR
#'# 0.045 0.046 0.045 0.047
#'#
#'# $`observed local deviation`
#'#        I1-C2 I2-C1 I2-C2  I3-C1  I3-C2  I4-C1  I4-C2
#'# group1 2.556 0.503 2.573 -2.573 -2.160 -1.272 -0.683
#'# group2 2.246 0.878 3.135 -1.852 -0.824 -0.494  0.941
#'#
#'# $`person score distribution in group 1`
#'#
#'#     1     2     3     4     5     6     7
#'# 0.016 0.097 0.137 0.347 0.121 0.169 0.113
#'#
#'# $`person score distribution in group 2`
#'#
#'#     1     2     3     4     5     6     7
#'# 0.015 0.083 0.136 0.280 0.152 0.227 0.106
#'#
#'# $`degrees of freedom`
#'# [1] 7
#'#
#'# $`noncentrality parameter`
#'#      W     LR     RS     GR
#'# 11.395 11.818 11.628 11.978
#'#
#'# $call
#'# post_hocPCM(alpha = 0.05, data = y, x = x)
#' }

post_hocPCM <- function(alpha = 0.05, data, x){

  subfunc <- function(stats) {
    e <- stats/(sum(t1)+sum(t2))
    nc <- stats
    beta <- pchisq(qchisq(1-alpha, df), df, nc)
    power <- 1 - beta
    return(list('power' = round(power, digits = 3),
                'observed global deviation' = round(e, digits = 3)))
    }

  call<-match.call()

  e <- tcl_splitcr(X = data, splitcr = x, model = "PCM")
  y1 <- e$X.list[[1]]
  y2 <- e$X.list[[2]]

  # y1 <- x * data
  # y2 <- abs(x - 1) * data

  r1 <- psychotools::pcmodel(y1)
  r2 <- psychotools::pcmodel(y2)
  r3 <- psychotools::pcmodel(rbind(y1,y2))

  df <- length(r1$coefficients)

  t1 <- table(factor(rowSums(y1),levels=1:df))
  t2 <- table(factor(rowSums(y2),levels=1:df))

  vstats <- do.call(c,invar_test_obj(r1=r1, r2=r2, r3=r3, model = "PCM"))

  res <- lapply(vstats, subfunc)

  results <- list(  'test' = round(vstats, digits = 3),
                    'power' = unlist(do.call(cbind, res)[1,]),
                    'observed global deviation' = unlist(do.call(cbind, res)[2,]),
                    'observed local deviation' = round( rbind("group1"=r1$coefficients, "group2"= r2$coefficients), digits = 3 ),
                    'person score distribution in group 1' = round(t1/sum(t1), digits=3),
                    'person score distribution in group 2' = round(t2/sum(t2), digits=3),
                    'degrees of freedom' = df,
                    'noncentrality parameter' = round(vstats, digits = 3),
                    "call" = call)
  return(results)
}
