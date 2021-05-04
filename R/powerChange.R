######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# powerChange
#
# Part of R/tlc - Testing in Conditional Likelihood Context package
#
# This file contains a routine that returns power of for W, LR, RS,
# and GR Test for hypothesis about parameter quantifying change
# between two time points as modeled by an LLTM given probability
# of error of first kind alpha, sample size, and a deviation from
# the hypothesis to be tested
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2021, Last Modified 19/04/2021
######################################################################
#' Power analysis of tests in context of measurement of change using LLTM
#'
#' Returns power of Wald (W), likelihood ratio (LR), Rao score (RS)
#' and gradient (GR) test given probability of error of first kind \eqn{\alpha}, sample size, and
#' a deviation from the hypothesis to be tested. The latter states that the shift parameter
#' quantifying the constant change for all items between time points 1 and 2 equals 0.
#' The alternative states that the shift parameter is not equal to 0.
#' It is assumed that the same items are presented at both time points. See function \code{\link{change_test}}.
#'
#'In general, the power of the tests is determined from the assumption that the approximate distributions of
#'the four test statistics are from the family of noncentral \eqn{\chi^2} distributions with \eqn{df = 1} and noncentrality
#'parameter \eqn{\lambda}. The latter depends on a scenario of deviation from the hypothesis to be tested and a specified sample size.
#'Given the probability of the error of the first kind \eqn{\alpha} the power of the tests can be determined from \eqn{\lambda}.
#'More details about the distributions of the test statistics and the relationship between \eqn{\lambda}, power, and sample size can be found
#'in Draxler and Alexandrowicz (2015).
#'
#'As regards the concept of sample size a distinction between informative and total sample size has to be made since the power
#'of the tests depends only on the informative sample size. In the conditional maximum likelihood context, the responses of
#'persons with minimum or maximum person score are completely uninformative. They do not contribute to the value of the test
#'statistic. Thus, the informative sample size does not include these persons. The total sample size is composed of all persons.
#'
#'In particular, the determination of \eqn{\lambda} and the power of the tests, respectively, is based on a simple Monte Carlo approach.
#'Data (responses of a large number of persons to a number of items presented at two time points) are generated given a
#'user-specified scenario of a deviation from the hypothesis to be tested. The hypothesis to be tested assumes no change
#'between time points 1 and 2. A scenario of a deviation is given by a choice of the item parameters at time point 1 and
#'the shift parameter, i.e., the LLTM eta parameters, as well as the person parameters (to be drawn randomly from a specified
#'distribution). The shift parameter represents a constant change of all item parameters from time point 1 to time point 2.
#'A test statistic \eqn{T} (Wald, LR, score, or gradient) is computed from the simulated data. The observed value \eqn{t} of the test
#'statistic is then divided by the informative sample size \eqn{n_{infsim}} observed in the simulated data. This yields the so-called
#'global deviation \eqn{e = t / n_{infsim}}, i.e., the chosen scenario of a deviation from the hypothesis to be tested being represented
#'by a single number. The power of the tests can be determined given a user-specified total sample size denoted by \eqn{n_{total}}.
#'The noncentrality parameter \eqn{\lambda} can then be expressed by \eqn{\lambda = n_{total}* (n_{infsim} / n_{totalsim}) * e},
#'where \eqn{n_{totalsim}} denotes the total number of persons in the simulated data and \eqn{n_{infsim} / n_{totalsim}} is the proportion of
#'informative persons in the sim. data. Let \eqn{q_{\alpha}} be the \eqn{1 - \alpha} quantile of the central \eqn{\chi^2} distribution with \eqn{df = 1}.
#'Then,
#'
#'\deqn{power = 1 - F_{df, \lambda} (q_{\alpha}),}
#'
#'where \eqn{F_{df, \lambda}} is the cumulative distribution function of the noncentral \eqn{\chi^2} distribution with \eqn{df = 1} and
#'\eqn{\lambda = n_{total} * (n_{infsim} / n_{totalsim}) * e}. Thereby, it is assumed that \eqn{n_{total}} is composed of a frequency distribution
#'of person scores that is proportional to the observed distribution of person scores in the simulated data.
#'
#'Note that in this approach the data have to be generated only once. There are no replications needed. Thus, the procedure is
#'computationally not very time-consuming.
#'
#'Since \eqn{e} is determined from the value of the test statistic observed in the simulated data it has to be treated as a realized
#'value of a random variable \eqn{E}. The same holds true for \eqn{\lambda} as well as the power of the tests. Thus, the power is a realized
#'value of a random variable that shall be denoted by \eqn{P}. Consequently, the (realized) value of the power of the tests need
#'not be equal to the exact power that follows from the user-specified \eqn{n_{total}}, \eqn{\alpha}, and the chosen item parameters and shift
#'parameter used for the simulation of the data. If the CML estimates of these parameters computed from the simulated data are
#'close to the predetermined parameters the power of the tests will be close to the exact value. This will generally be the
#'case if the number of person parameters used for simulating the data is large, e.g., \eqn{10^5} or even \eqn{10^6} persons. In such
#'cases, the possible random error of the computation procedure based on the sim. data may not be of practical relevance
#'any more. That is why a large number (of persons for the simulation process) is generally recommended.
#'
#'For theoretical reasons, the random error involved in computing the power of the tests can be pretty well approximated.
#'A suitable approach is the well-known delta method. Basically, it is a Taylor polynomial of first order, i.e., a linear
#'approximation of a function. According to it the variance of a function of a random variable can be linearly approximated
#'by multiplying the variance of this random variable with the square of the first derivative of the respective function.
#'In the present problem, the variance of the test statistic \eqn{T} is (approximately) given by the variance of a noncentral
#'\eqn{\chi^2} distribution. Thus, \eqn{Var(T) = 2 (df + 2 \lambda)},
#'with \eqn{df = 1} and \eqn{\lambda = t}.
#'Since the global deviation \eqn{e = (1 / n_{infsim})* t} it follows for the variance of the corresponding random variable \eqn{E}
#'that \eqn{Var(E) = (1 / n_{infsim})^2 * Var(T)}. The power of the tests is a function of \eqn{e} which is given by
#'\eqn{F_{df, \lambda} (q_{\alpha})}, where \eqn{\lambda = n_{total} * (n_{infsim} / n_{totalsim}) * e} and \eqn{df = 1}.
#'Then, by the delta method one obtains (for the variance of \eqn{P})
#'
#'\deqn{Var(P) = Var(E) * (F'_{df, \lambda} (q_{\alpha}))^2,}
#'
#'where \eqn{F'_{df, \lambda}} is the derivative of \eqn{F_{df, \lambda}} with respect to \eqn{e}. This derivative is determined
#'numerically and evaluated at \eqn{e} using the package numDeriv. The square root of \eqn{Var(P)} is then used to quantify the random
#'error of the suggested Monte Carlo computation procedure. It is called Monte Carlo error of power.
#'
#' @param alpha Probability of the error of first kind.
#' @param n_total Total sample size for which power shall be determined.
#' @param eta A vector of eta parameters of the LLTM. The last element represents the constant change or shift for all items
#' between time points 1 and 2. The other elements of the vector are the item parameters at time point 1. A choice of the eta
#' parameters constitutes a scenario of deviation from the hypothesis of no change.
#' @param persons A vector of person parameters (drawn from a specified distribution). By default \eqn{10^6} parameters are drawn at
#' random from the standard normal distribution. The larger this number the more accurate are the computations. See Details.
#'
#'@return A list of results.
#'  \item{power}{Power value for each test.}
#'  \item{MC error of power}{Monte Carlo error of power computation for each test.}
#'  \item{deviation}{Shift parameter estimated from the simulated data representing the constant shift of item parameters between time points 1 and 2.}
#'  \item{person score distribution}{Relative frequencies of person scores observed in simulated data. Uninformative scores,
#'  i.e., minimum and maximum score, are omitted. Note that the person score distribution does also have an influence on the power of the tests.}
#'  \item{degrees of freedom}{Degrees of freedom \eqn{df}.}
#'  \item{noncentrality parameter}{Noncentrality parameter \eqn{\lambda} of \eqn{\chi^2} distribution from which power is determined.}
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
#' @seealso \code{\link{sa_sizeChange}}, and \code{\link{post_hocChange}}.
#' @examples
#' \dontrun{
#'
#'# Numerical example: 4 items presented twice, thus 8 virtual items
#'
#'# eta Parameter, first 4 are nuisance
#'# (easiness parameters of the 4 items at time point 1),
#'# last one is the shift parameter
#' eta <- c(-2,-1,1,2,0.5)
#' res <- powerChange(alpha = 0.05, n_total=150, eta=eta, persons=rnorm(10^6))
#'
#'# > res
#'# $power
#'#     W    LR    RS    GR
#'# 0.905 0.910 0.908 0.911
#'#
#'# $`MC error of power`
#'#     W    LR    RS    GR
#'# 0.002 0.002 0.002 0.002
#'#
#'# $`deviation (estimate of shift parameter)`
#'# [1] 0.499
#'#
#'# $`person score distribution`
#'#
#'#     1     2     3     4     5     6     7
#'# 0.034 0.093 0.181 0.249 0.228 0.147 0.068
#'#
#'# $`degrees of freedom`
#'# [1] 1
#'#
#'# $`noncentrality parameter`
#'#      W     LR     RS     GR
#'# 10.692 10.877 10.815 10.939
#'#
#'# $call
#'# powerChange(alpha = 0.05, n_total = 150, eta = eta, persons = rnorm(10^6))
#'#
#' }



powerChange <- function(alpha = 0.05, n_total, eta, persons = rnorm(10^6)){

  subfunc <- function(stats) {
    e <- stats/sum(t)
    nc <- e * n_total * (sum(t) / nrow(y))
    f <- function(nc){pchisq(qchisq(1-alpha, 1), 1, nc)}
    beta <- f(nc)
    se <- sqrt((2 + 4*stats) * (1/sum(t))^2 * (numDeriv::grad(f, nc) * n_total * (sum(t) / nrow(y)))^2)
    power <- 1 - beta
    return(list('power' = round(power, digits = 3),
                'MC error of power' = round(se, digits = 3),
                'noncentrality parameter' = round(nc, digits = 3) ))
  }

  call<-match.call()

  # design matrix W used only for data generation
  #    (not used for estimating in change_test() function)
  k <- length(eta) - 1
  W <- cbind (rbind( diag(x=1, nrow = k, ncol = k), diag(x=1, nrow = k, ncol = k)),
              c(rep(0,k), rep(1,k) ) )

  y <- eRm::sim.rasch(persons=persons, items=colSums(-eta*t(W)))

  t <- table(factor(rowSums(y),levels=1:(ncol(y)-1)))

  est <- eRm::LLTM(X=y, mpoints = 2,se = FALSE, sum0 = FALSE) # unrestricted CML estimates of eta Parameters
  dev <- unname(est$etapar[ncol(y)/2])

  vstats <- change_test(X = y)$test[c(4,2,3,1)] # W, LR, RS, GR test

  res <- lapply(vstats, subfunc)

  results <- list(  'power' = unlist(do.call(cbind, res)[1,]),
                    'MC error of power' = unlist(do.call(cbind, res)[2,]),
                    'deviation (estimate of shift parameter)' =  round(dev, digits = 3),
                    'person score distribution' = round(t/sum(t), digits=3),
                    'degrees of freedom' = 1,
                    'noncentrality parameter' = unlist(do.call(cbind, res)[3,]),
                    "call" = call)
  return(results)
}
