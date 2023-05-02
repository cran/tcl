######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# sa_sizeChange
#
# Part of R/tlc - Testing in Conditional Likelihood Context package
#
# This file contains a routine that returns sample size for W, LR, RS,
# and GR Test for hypothesis about parameter  quantifying change
# between two time points as modeled by an LLTM given probabilities
# of errors of first and second kinds alpha and beta, and a given
# deviation from the hypothesis to be tested
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2021, Last Modified 25/10/2021
######################################################################
#' Sample size planning for tests in context of measurement of change using LLTM
#'
#' Returns sample size for Wald (W), likelihood ratio (LR), Rao score (RS)
#' and gradient (GR) test given probabilities of errors of first and second kinds \eqn{\alpha} and \eqn{\beta}
#' as well as a deviation from the hypothesis to be tested. The hypothesis to be tested states that
#' the shift parameter quantifying the constant change for all items between time points 1 and 2
#' equals 0. The alternative states that the shift parameter is not equal to 0. It is assumed that the same
#' items are presented at both time points. See function \code{\link{change_test}}.
#'
#' In general, the sample size is determined from the assumption that the approximate distributions
#' of the four test statistics are from the family of noncentral \eqn{\chi^2} distributions with \eqn{df = 1}
#' and noncentrality parameter \eqn{\lambda}. The latter is, inter alia, a function of the sample size. Hence,
#' the sample size can be determined from the condition \eqn{\lambda = \lambda_0}, where \eqn{\lambda_0} is
#' a predetermined constant which depends on the probabilities of the errors of the first and second kinds
#' \eqn{\alpha} and  \eqn{\beta} (or power). More details about the distributions of the test statistics
#' and the relationship between \eqn{\lambda}, power, and sample size can be found in Draxler and Alexandrowicz (2015).
#'
#'In particular, the determination of \eqn{\lambda} and the sample size, respectively, is based on a simple Monte
#'Carlo approach. As regards the concept of sample size a distinction between informative and total
#'sample size has to be made. In the conditional maximum likelihood context, the responses of
#'persons with minimum or maximum person score are completely uninformative. They do not contribute
#'to the value of the test statistic. Thus, the informative sample size does not include these persons.
#'The total sample size is composed of all persons. The Monte Carlo approach used in the present
#'problem to determine \eqn{\lambda} and informative (and total) sample size can briefly be described as follows.
#'Data (responses of a large number of persons to a number of items presented at two time points) are
#'generated given a user-specified scenario of a deviation from the hypothesis to be tested. The
#'hypothesis to be tested assumes no change between time points 1 and 2. A scenario of a deviation
#'is given by a choice of the item parameters at time point 1 and the shift parameter, i.e., the
#'LLTM eta parameters, as well as the person parameters (to be drawn randomly from a specified distribution).
#'The shift parameter represents a constant change of all item parameters from time point 1 to time point 2.
#'A test statistic \eqn{T} (Wald, LR, score, or gradient) is computed from the simulated data. The observed
#'value \eqn{t} of the test statistic is then divided by the informative sample size \eqn{n_{infsim}} observed
#'in the simulated data. This yields the so-called global deviation \eqn{e = t / n_{infsim}}, i.e.,
#'the chosen scenario of a deviation from the hypothesis to be tested being represented by a
#'single number. Let the informative sample size sought be denoted by \eqn{n_{inf}} (thus, this is not
#'the informative sample size observed in the sim. data). The noncentrality parameter \eqn{\lambda} can
#'be expressed by the product \eqn{n_{inf} * e}. Then, it follows from the condition \eqn{\lambda = \lambda_0} that
#'
#'\deqn{n_{inf} * e = \lambda_0}
#'
#'and
#'
#'\deqn{n_{inf} = \lambda_0 / e.}
#'
#'Note that the sample of size \eqn{n_{inf}} is assumed to be composed only of persons with informative person scores, where
#'the relative frequency distribution of these informative scores is considered to be equal to the
#'observed relative frequency distribution of the informative scores in the simulated data. The total sample size
#'\eqn{n_{total}} is then obtained from the relation \eqn{n_{inf} = n_{total} * pr}, where \eqn{pr} is the proportion
#'or relative frequency of persons observed in the simulated data with a minimum or maximum score. Basing
#'the tests given a level \eqn{\alpha} on an informative sample of size \eqn{n_{inf}} the probability of rejecting
#'the hypothesis to be tested will be at least \eqn{1 - \beta} if the true global deviation \eqn{\geq e}.
#'
#'Note that in this approach the data have to be generated only once. There are no replications
#'needed. Thus, the procedure is computationally not very time-consuming.
#'
#'Since e is determined from the value of the test statistic observed in the simulated data it has
#'to be treated as a realized value of a random variable \eqn{E}. Consequently, \eqn{n_{inf}} is also a
#'realization of a random variable \eqn{N_{inf}}. Thus, the (realized) value \eqn{n_{inf}} need not be
#'equal to the exact value of the informative sample size that follows from the user-specified
#'(predetermined) \eqn{\alpha}, \eqn{\beta}, and scenario of a deviation from the hypothesis to be
#'tested, i.e., the selected item parameters and shift parameter used for the simulation of the data.
#'If the CML estimates of these parameters computed from the simulated data are close to the
#'predetermined parameters \eqn{n_{inf}} will be close to the exact value. This will generally be the
#'case if the number of person parameters used for simulating the data is large, e.g.,
#'\eqn{10^5} or even \eqn{10^6} persons. In such cases, the possible random error of the computation procedure
#'of \eqn{n_{inf}} based on the sim. data may not be of practical relevance any more. That is why a
#'large number (of persons for the simulation process) is generally recommended.
#'
#'For theoretical reasons, the random error involved in computing \eqn{n_{inf}} can be pretty well approximated.
#'A suitable approach is the well-known delta method. Basically, it is a Taylor polynomial of first order,
#'i.e., a linear approximation of a function. According to it the variance of a function of a random
#'variable can be linearly approximated by multiplying the variance of this random variable with the square
#'of the first derivative of the respective function. In the present problem, the variance of the test
#'statistic \eqn{T} is (approximately) given by the variance of a noncentral \eqn{\chi^2} distribution. Thus, \eqn{Var(T) = 2 (df + 2 \lambda)},
#'with \eqn{df = 1} and \eqn{\lambda = t}. Since the global deviation
#'\eqn{e = (1 / n_{infsim}) * t} it follows for the variance of the corresponding random variable \eqn{E} that
#'\eqn{Var(E) = (1 / n_{infsim})^2 * Var(T)}. Since \eqn{n_{inf} = f(e) = \lambda_0 / e} one obtains by the
#'delta method (for the variance of the corresponding random variable \eqn{N_{inf}})
#'
#'\deqn{Var(N_{inf}) = Var(E) * (f'(e))^2,}
#'
#'where \eqn{f'(e) = - \lambda_0 / e^2} is the derivative of \eqn{f(e)}. The square root of \eqn{Var(N_{inf})}
#'is then used to quantify the random error of the suggested Monte Carlo computation procedure. It is called
#'Monte Carlo error of informative sample size.
#'
#' @param alpha Probability of error of first kind.
#' @param beta Probability of error of second kind.
#' @param eta A vector of eta parameters of the LLTM. The last element represents the constant change or shift for all items
#' between time points 1 and 2. The other elements of the vector are the item parameters at time point 1. A choice of the eta
#' parameters constitutes a scenario of deviation from the hypothesis of no change.
#' @param persons A vector of person parameters (drawn from a specified distribution). By default \eqn{10^6} parameters
#' are drawn at random from the standard normal distribution. The larger this number the more accurate are the computations.
#' See Details.
#'
#'@return A list results.
#'  \item{informative sample size}{Informative sample size for each test, omitting persons with min. and max score.}
#'  \item{MC error of sample size}{Monte Carlo error of sample size computation for each test.}
#'  \item{deviation}{Shift parameter estimated from the simulated data representing the constant shift of item
#'  parameters between time points 1 and 2.}
#'  \item{person score distribution}{Relative frequencies of person scores observed in simulated data. Uninformative scores,
#'  i.e., minimum and maximum score, are omitted.
#'  Note that the person score distribution does also have an influence on the sample size.}
#'  \item{degrees of freedom}{Degrees of freedom \eqn{df}.}
#'  \item{noncentrality parameter}{Noncentrality parameter \eqn{\lambda} of \eqn{\chi^2} distribution from which sample size is determined.}
#'  \item{total sample size}{Total sample size for each test. See Details.}
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
#'  }
#' @keywords sample_size_planning
#' @export
#' @seealso \code{\link{powerChange}}, and \code{\link{post_hocChange}}.
#' @examples
#' \dontrun{
#'# Numerical example 4 items presented twice, thus 8 virtual items
#'
#'# eta Parameter, first 4 are nuisance
#'# (easiness parameters of the 4 items at time point 1),
#'# last one is the shift parameter
#' eta <- c(-2,-1,1,2,0.5)
#'
#'res <- sa_sizeChange(alpha = 0.05, beta = 0.05, eta=eta, persons = rnorm(10^6))
#'
#'# > res
#'# $`informative sample size`
#'#   W  LR  RS  GR
#'# 177 174 175 173
#'#
#'# $`MC error of sample size`
#'#     W    LR    RS    GR
#'# 1.321 1.287 1.299 1.276
#'#
#'# $`deviation (estimate of shift parameter)`
#'# [1] 0.501
#'#
#'# $`person score distribution`
#'#
#'#     1     2     3     4     5     6     7
#'# 0.034 0.094 0.181 0.249 0.227 0.147 0.068
#'#
#'# $`degrees of freedom`
#'# [1] 1
#'#
#'# $`noncentrality parameter`
#'# [1] 12.995
#'#
#'# $`total sample size`
#'#   W  LR  RS  GR
#'# 182 179 180 178
#'#
#'# $call
#'# sa_sizeChange(alpha = 0.05, beta = 0.05, eta = eta, persons = rnorm(10^6))
#' }

sa_sizeChange <- function(alpha = 0.05, beta = 0.05, eta, persons = rnorm(10^6)){

  subfunc <- function(stats) {
    e <- stats/sum(t)
    n <- ceiling(lambda0 / e)
    se <- sqrt((2*df + 4*stats) * lambda0^2 * sum(t)^2 * {stats}^-4)
    n1 <- ceiling((n)/(sum(t)/nrow(y)))

    return(list('informative sample size' = n,
                'MC of sample size' = round(se, digits = 3),
                'total sample size' = n1) )
  }

  call<-match.call()

  # design matrix W used only for data generation
  #    (not used for estimating in change_test() function)
  k <- length(eta) - 1
  W <- cbind (rbind( diag(x=1, nrow = k, ncol = k), diag(x=1, nrow = k, ncol = k)),
              c(rep(0,k), rep(1,k) ) )

  y <- eRm::sim.rasch(persons=persons, items=colSums(-eta*t(W)))

  df <- 1

  func <- function (x) {beta - pchisq(qchisq(1-alpha,1), 1, ncp = x)}
  lambda0 <- uniroot(f = func, interval = c(0,10^3), tol = .Machine$double.eps^0.5)$root
  t <- table(factor(rowSums(y),levels=1:(ncol(y)-1)))

  est <- eRm::LLTM(X=y, mpoints = 2,se = FALSE, sum0 = FALSE) # unrestricted CML estimates of eta Parameters
  dev <- unname(est$etapar[ncol(y)/2])

  vstats <- change_test(X = y)$test[c(4,2,3,1)] # W, LR, RS, GR test

  res <- lapply(vstats, subfunc)

  results <- list(  'informative sample size' = unlist(do.call(cbind, res)[1,]),
                    'MC error of sample size' = unlist(do.call(cbind, res)[2,]),
                    'deviation (estimate of shift parameter)' =  round(dev, digits = 3),
                    'person score distribution' = round(t/sum(t), digits=3),
                    'degrees of freedom' = df,
                    'noncentrality parameter' = round(lambda0, digits = 3),
                    'total sample size' = unlist(do.call(cbind, res)[3,]),
                    "call" = call)
  return(results)
}
