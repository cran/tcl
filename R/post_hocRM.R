######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# post_hocRM
#
# Part of R/tlc - Testing in Conditional Likelihood Context package
#
# This file contains a routine that returns post hoc power of W, LRT, RS, GR
# given data and probability of error of first kind alpha
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2021, Last Modified 19/04/2021
######################################################################
#' Power analysis of tests of invariance of item parameters between two groups of persons in binary Rasch model
#'
#' Returns post hoc power of Wald (W), likelihood ratio (LR), Rao score (RS)
#' and gradient (GR) test given data and probability of error of first kind \eqn{\alpha}.
#' The hypothesis to be tested assumes equal item parameters  between two predetermined groups
#' of persons. The alternative states that at least one of the parameters differs between the two
#' groups.
#'
#' The power of the tests (Wald, LR, score, and gradient) is determined from the assumption that the
#' approximate distributions of the four test statistics are from the family of noncentral \eqn{\chi^2}
#' distributions with \eqn{df} equal to the number of items minus 1 and noncentrality parameter \eqn{\lambda}. In case
#' of evaluating the post hoc power, \eqn{\lambda} is assumed to be given by the observed value of the test statistic.
#' Given the probability of the error of the first kind \eqn{\alpha} the post hoc power of the tests can be
#' determined from \eqn{\lambda}. More details about the distributions of the test statistics and the relationship
#' between \eqn{\lambda}, power, and sample size can be found in Draxler and Alexandrowicz (2015).
#'
#' In particular, let \eqn{q_{\alpha}} be the \eqn{1- \alpha} quantile of the central \eqn{\chi^2} distribution
#' with df equal to the number of items minus 1. Then,
#'
#' \deqn{power = 1 - F_{df, \lambda} (q_{\alpha}),}
#'
#' where \eqn{F_{df, \lambda}} is the cumulative distribution function of the noncentral \eqn{\chi^2} distribution
#' with \eqn{df} equal to the number of items reduced by 1 and \eqn{\lambda} equal to the observed value of the test statistic.
#'
#' @param alpha Probability of error of first kind.
#' @param data Binary data matrix.
#' @param x A numeric vector of length equal to number of persons containing zeros and ones indicating group membership of the persons.
#'
#'@return A list of results.
#'  \item{test}{A numeric vector of Wald (W), likelihood ratio (LR), Rao score (RS), and gradient (GR) test statistics.}
#'  \item{power}{Post hoc power value for each test.}
#'  \item{global deviation}{Observed global deviation from hypothesis to be tested represented by a single number.
#'  It is obtained by dividing the test statistic by the informative sample size. The latter does not include persons
#'  with minimum or maximum person score. }
#'  \item{local deviation}{CML estimates of free item parameters in both groups of persons (first item parameter set
#'  to 0 in both groups) representing observed deviation from hypothesis to be tested locally per item.}
#'  \item{person score distribution in group 1}{Relative frequencies of person scores in group 1. Uninformative scores,
#'  i.e., minimum and maximum score, are omitted. Note that the person score distribution does also have an influence on
#'  the power of the tests.}
#'  \item{person score distribution in group 2}{Relative frequencies of person scores in group 2. Uninformative scores,
#'  i.e., minimum and maximum score, are omitted. Note that the person score distribution does also have an influence on
#'  the power of the tests.}
#'  \item{degrees of freedom}{Degrees of freedom \eqn{df}.}
#'  \item{noncentrality parameter}{Noncentrality parameter \eqn{\lambda} of \eqn{\chi^2} distribution from which power is determined.
#'  It equals observed value of test statistic.}
#'  \item{call}{The matched call.}
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
#'
#'  }
#' @keywords sample_size_planning
#' @export
#' @seealso \code{\link{sa_sizeRM}}, and \code{\link{powerRM}}.
#' @examples
#' \dontrun{
#'# Numerical example for post hoc power analysis for Rasch Model
#'
#' y <- eRm::raschdat1
#' n <- nrow(y) # sample size
#' x <- c( rep(0,n/2), rep(1,n/2) ) # binary covariate
#'
#' res <-  post_hocRM(alpha = 0.05, data = y, x = x)
#'
#'# > res
#'# $test
#'#      W     LR     RS     GR
#'# 29.241 29.981 29.937 30.238
#'#
#'# $power
#'#     W    LR    RS    GR
#'# 0.890 0.900 0.899 0.903
#'#
#'# $`observed global deviation`
#'#     W    LR    RS    GR
#'# 0.292 0.300 0.299 0.302
#'#
#'# $`observed local deviation`
#'#           I2    I3    I4    I5    I6    I7    I8    I9   I10   I11
#'# group1 1.039 0.693 2.790 2.404 1.129 1.039 0.864 1.039 2.790 2.244
#'# group2 2.006 0.945 2.006 3.157 1.834 0.690 0.822 1.061 2.689 2.260
#'#          I12   I13   I14   I15   I16   I17   I18   I19   I20   I21
#'# group1 1.412 3.777 3.038 1.315 2.244 1.039 1.221 2.404 0.608 0.608
#'# group2 0.945 2.962 4.009 1.171 2.175 1.472 2.091 2.344 1.275 0.690
#'#          I22   I23   I24   I25   I26   I27   I28   I29   I30
#'# group1 0.438 0.608 1.617 3.038 0.438 1.617 2.100 2.583 0.864
#'# group2 0.822 1.275 1.565 2.175 0.207 1.746 1.746 2.260 0.822
#'#
#'# $`person score distribution in group 1`
#'#
#'#    1    2    3    4    5    6    7    8    9   10   11   12   13
#'# 0.02 0.02 0.02 0.06 0.02 0.10 0.10 0.06 0.10 0.12 0.08 0.12 0.12
#'#   14   15   16   17   18   19   20   21   22   23   24   25   26
#'# 0.06 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
#'#   27   28   29
#'# 0.00 0.00 0.00
#'#
#'# $`person score distribution in group 2`
#'#
#'#    1    2    3    4    5    6    7    8    9   10   11   12   13
#'# 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
#'#   14   15   16   17   18   19   20   21   22   23   24   25   26
#'# 0.08 0.12 0.10 0.16 0.06 0.04 0.10 0.12 0.08 0.02 0.02 0.02 0.08
#'#   27   28   29
#'# 0.00 0.00 0.00
#'#
#'# $`degrees of freedom`
#'# [1] 29
#'#
#'# $`noncentrality parameter`
#'#      W     LR     RS     GR
#'# 29.241 29.981 29.937 30.238
#'#
#'# $call
#'# post_hocRM(alpha = 0.05, data = y, x = x)
#'
#'}



post_hocRM <- function(alpha = 0.05, data, x){

  subfunc <- function(stats) {
    e <- stats/(sum(t1)+sum(t2))
    nc <- stats
    beta <- pchisq(qchisq(1-alpha, df), df, nc)
    power <- 1 - beta
    return(list('power' = round(power, digits = 3),
                'observed global deviation' = round(e, digits = 3)))
  }
  call<-match.call()

  # Check for inappropriate response patterns within subgroups
  e <- tcl_splitcr(X = data, splitcr = x, model = "RM")
  y1 <- e$X.list[[1]]
  y2 <- e$X.list[[2]]

  r1 <- psychotools::raschmodel(y1)
  r2 <- psychotools::raschmodel(y2)
  r3 <- psychotools::raschmodel(rbind(y1,y2))

  df <- ncol(y1)-1

  t1 <- table(factor(rowSums(y1),levels=1:(ncol(y1)-1)))
  t2 <- table(factor(rowSums(y2),levels=1:(ncol(y2)-1)))

  vstats <- do.call(c, invar_test_obj(r1=r1, r2=r2, r3=r3, model = "RM"))

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
