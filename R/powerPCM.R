######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# powerPCM
#
# Part of R/tlc - Testing in Conditional Likelihood Context package
#
# This file contains a routine that returns power of W, LR, RS
# and GR given probability of error of first kind alpha, sample size,
# and a deviation from the hypothesis to be tested,
# i.e. that the difference of two group parameters equals zero,
# where one assumes all items presented to all persons.
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2021, Last Modified 25/10/2021
######################################################################
#' Power analysis of tests of invariance of item parameters between two groups
#' of persons in partial credit model
#'
#' Returns power of Wald (W), likelihood ratio (LR), Rao score (RS)
#' and gradient (GR) test given probability of error of first kind
#' \eqn{\alpha}, sample size, and a deviation from the hypothesis to be tested.
#' The hypothesis to be tested assumes equal item-category parameters of the
#' partial credit model between two predetermined groups of persons. The alternative
#' states that at least one of the parameters differs between the two groups.
#'
#' In general, the power of the tests is determined from the assumption that the
#' approximate distributions of the four test statistics are from the family of
#' noncentral \eqn{\chi^2} distributions with \eqn{df} equal to the number of
#' free item-category parameters and noncentrality parameter \eqn{\lambda}.
#' The latter depends on a scenario of deviation from the hypothesis to be tested
#' and a specified sample size. Given the probability of the error of the first
#' kind \eqn{\alpha} the power of the tests can be determined from \eqn{\lambda}.
#' More details about the distributions of the test statistics and the relationship
#' between \eqn{\lambda}, power, and sample size can be found in Draxler and
#' Alexandrowicz (2015).
#'
#'As regards the concept of sample size a distinction between informative and total
#'sample size has to be made since the power of the tests depends only on the informative
#'sample size. In the conditional maximum likelihood context, the responses of persons
#'with minimum or maximum person score are completely uninformative. They do not contribute
#'to the value of the test statistic. Thus, the informative sample size does not include
#'these persons. The total sample size is composed of all persons.
#'
#'In particular, the determination of \eqn{\lambda} and the power of the tests, respectively,
#'is based on a simple Monte Carlo approach. Data (responses of a large number of persons
#'to a number of items) are generated given a user-specified scenario of a deviation from
#'the hypothesis to be tested. A scenario of a deviation is given by a choice of the
#'item-cat. parameters and the person parameters (to be drawn randomly from a specified
#'distribution) for each of the two groups. Such a scenario may be called local deviation
#'since deviations can be specified locally for each item-category. The relative group
#'sizes are determined by the choice of the number of person parameters for each of the
#'two groups. For instance, by default \eqn{10^6} person parameters are selected randomly for
#'each group. In this case, it is implicitly assumed that the two groups of persons are
#'of equal size. The user can specify the relative group sizes by choosing the length of
#'the arguments persons1 and persons2 appropriately. Note that the relative group sizes
#'do have an impact on power and sample size of the tests. The next step is to compute a
#'test statistic \eqn{T} (Wald, LR, score, or gradient) from the simulated data. The observed
#'value \eqn{t} of the test statistic is then divided by the informative sample size
#'\eqn{n_{infsim}} observed in the simulated data. This yields the so-called global deviation
#'\eqn{e = t / n_{infsim}}, i.e., the chosen scenario of a deviation from the hypothesis to
#'be tested being represented by a single number. The power of the tests can be determined
#'given a user-specified total sample size denoted by \code{n_total}. The noncentrality
#'parameter \eqn{\lambda} can then be expressed by
#'\eqn{\lambda = n_{total}* (n_{infsim} / n_{totalsim}) * e}, where \eqn{n_{totalsim}} denotes
#'the total number of persons in the simulated data and \eqn{n_{infsim} / n_{totalsim}} is
#'the proportion of informative persons in the sim. data. Let \eqn{q_{\alpha}} be the
#'\eqn{1 - \alpha} quantile of the central \eqn{\chi^2} distribution with df equal to the
#'number of free item-category parameters. Then,
#'
#'\deqn{power = 1 - F_{df, \lambda} (q_{\alpha}),}
#'
#'where \eqn{F_{df, \lambda}} is the cumulative distribution function of the noncentral
#'\eqn{\chi^2} distribution with \eqn{df} equal to the number of free item-category parameters
#'and \eqn{\lambda = n_{total}  (n_{infsim} / n_{totalsim}) * e}. Thereby, it is assumed that
#'\eqn{n_{total}} is composed of a frequency distribution of person scores that is proportional
#'to the observed distribution of person scores in the simulated data. The same holds
#'true in respect of the relative group sizes, i.e., the relative frequencies of the two
#'person groups in a sample of size \eqn{n_{total}} are assumed to be equal to the relative frequencies of the two
#'groups in the simulated data.
#'
#'Note that in this approach the data have to be generated only once. There are no
#'replications needed. Thus, the procedure is computationally not very time-consuming.
#'
#'Since \eqn{e} is determined from the value of the test statistic observed in the simulated
#'data it has to be treated as a realized value of a random variable \eqn{E}. The same holds
#'true for \eqn{\lambda} as well as the power of the tests. Thus, the power is a realized
#'value of a random variable that shall be denoted by \eqn{P}. Consequently, the (realized)
#'value of the power of the tests need not be equal to the exact power that follows from the
#'user-specified \eqn{n_{total}}, \eqn{\alpha}, and the chosen item-category parameters used
#'for the simulation of the data. If the CML estimates of these parameters computed from the
#'simulated data are close to the predetermined parameters the power of the tests will be
#'close to the exact value. This will generally be the case if the number of person parameters
#'used for simulating the data is large, e.g., \eqn{10^5} or even \eqn{10^6} persons. In such cases,
#'the possible random error of the computation procedure based on the sim. data may not be of
#'practical relevance any more. That is why a large number (of persons for the simulation process)
#'is generally recommended.
#'
#'For theoretical reasons, the random error involved in computing the power of the tests can
#'be pretty well approximated. A suitable approach is the well-known delta method. Basically,
#'it is a Taylor polynomial of first order, i.e., a linear approximation of a function.
#'According to it the variance of a function of a random variable can be linearly approximated
#'by multiplying the variance of this random variable with the square of the first derivative
#'of the respective function. In the present problem, the variance of the test statistic \eqn{T}
#'is (approximately) given by the variance of a noncentral \eqn{\chi^2} distribution with \eqn{df}
#'equal to the number of free item-category parameters and noncentrality parameter \eqn{\lambda}.
#'Thus,  \eqn{Var(T) = 2 (df + 2 \lambda)}, with \eqn{\lambda = t}. Since the global
#'deviation \eqn{e = (1 / n_{infsim}) * t} it follows for the variance of the corresponding random variable \eqn{E}
#'that \eqn{Var(E) = (1 / n_{infsim})^2 * Var(T)}.
#'The power of the tests is a function of \eqn{e} which is given by \eqn{F_{df, \lambda} (q_{\alpha})},
#'where \eqn{\lambda = n_{total} * (n_{infsim} / n_{totalsim}) * e} and \eqn{df} equal to the
#'number of free item-category parameters. Then, by the delta method one obtains (for the variance of P).
#'
#'\deqn{Var(P) = Var(E) * (F'_{df, \lambda} (q_{\alpha}))^2,}
#'
#'where \eqn{F'_{df, \lambda}} is the derivative of \eqn{F_{df, \lambda}} with respect to \eqn{e}.
#'This derivative is determined numerically and evaluated at \eqn{e} using the package numDeriv. The square root of
#'\eqn{Var(P)} is then used to quantify the random error of the suggested Monte Carlo computation
#'procedure. It is called Monte Carlo error of power.
#'
#' @param alpha Probability of error of first kind.
#' @param n_total Total sample size for which power shall be determined.
#' @param persons1 A vector of person parameters in group 1 (drawn from a specified distribution).
#' By default \eqn{10^6} parameters are drawn at random from the standard normal distribution. The larger
#' this number the more accurate are the computations. See Details.
#' @param persons2 A vector of person parameters in group 2 (drawn from a specified distribution).
#' By default \eqn{10^6} parameters are drawn at random from the standard normal distribution. The larger
#' this number the more accurate are the computations. See Details.
#' @param local_dev A list consisting of two lists. One list refers to group 1, the other to group 2.
#' Each of the two lists contains a numeric vector per item, i.e., each list contains as many vectors as items.
#' Each vector contains the free item-cat. parameters of the respective item. The number of free item-cat.
#' parameters per item equals the number of categories of the item minus 1.
#'
#'@return A list of results.
#'  \item{power}{Power value for each test.}
#'  \item{MC error of power}{Monte Carlo error of power computation for each test.}
#'  \item{global deviation}{Global deviation computed from simulated data for each test. See Details.}
#'  \item{local deviation}{CML estimates of free item-category parameters in both groups of persons
#'  obtained from the simulated data expressing a deviation from the hypothesis to be tested locally
#'  per item and response category.}
#'  \item{person score distribution in group 1}{Relative frequencies of person scores in group 1 observed in
#'  simulated data. Uninformative scores, i.e., minimum and maximum score, are omitted.
#'  Note that the person score distribution does also have an influence on the power of the tests.}
#'  \item{person score distribution in group 2}{Relative frequencies of person scores in group 2 observed in
#'  simulated data. Uninformative scores, i.e., minimum and maximum score, are omitted.
#'  Note that the person score distribution does also have an influence on the power of the tests.}
#'  \item{degrees of freedom}{Degrees of freedom \eqn{df}.}
#'  \item{noncentrality parameter}{Noncentrality parameter \eqn{\lambda} of \eqn{\chi^2} distribution from which power is determined.}
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
#' @seealso \code{\link{sa_sizePCM}}, and \code{\link{post_hocPCM}}.
#' @examples
#' \dontrun{
#' # Numerical example
#'
#' # free item-category parameters for group 1 and 2  with 5 items, with 3 categories each
#' local_dev <-  list (  list(c( 0, 0), c( -1, 0), c( 0, 0),  c( 1, 0), c( 1, 0.5)) ,
#'                       list(c( 0, 0), c( -1, 0), c( 0, 0),  c( 1, 0), c( 0, -0.5))  )
#'
#' res <-  powerPCM(alpha = 0.05, n_total = 200, persons1 = rnorm(10^6),
#'                   persons2 = rnorm(10^6), local_dev = local_dev)
#'
#'# > res
#'# $power
#'#     W    LR    RS    GR
#'# 0.863 0.885 0.876 0.892
#'#
#'# $`MC error of power`
#'#     W    LR    RS    GR
#'# 0.002 0.002 0.002 0.002
#'#
#'# $`global deviation`
#'#     W    LR    RS    GR
#'# 0.102 0.107 0.105 0.109
#'#
#'# $`local deviation`
#'#         I1-C2  I2-C1  I2-C2  I3-C1  I3-C2 I4-C1 I4-C2  I5-C1  I5-C2
#'# group1  0.002 -0.997 -0.993  0.006  0.012 1.002 1.007  1.006  1.508
#'# group2 -0.007 -1.005 -1.007 -0.006 -0.009 0.993 0.984 -0.006 -0.510
#'#
#'# $`person score distribution in group 1`
#'#
#'#     1     2     3     4     5     6     7     8     9
#'# 0.112 0.130 0.131 0.129 0.122 0.114 0.101 0.091 0.070
#'#
#'# $`person score distribution in group 2`
#'#
#'#     1     2     3     4     5     6     7     8     9
#'# 0.091 0.108 0.117 0.122 0.122 0.121 0.115 0.110 0.093
#'#
#'# $`degrees of freedom`
#'# [1] 9
#'#
#'# $`noncentrality parameter`
#'#      W     LR     RS     GR
#'# 18.003 19.024 18.596 19.403
#'#
#'# $call
#'# powerPCM(alpha = 0.05, n_total = 200, persons1 = rnorm(10^6),
#'#          persons2 = rnorm(10^6), local_dev = local_dev)
#' }

powerPCM <- function(alpha = 0.05, n_total, persons1 = rnorm(10^6), persons2 = rnorm(10^6), local_dev){

  subfunc <- function(stats) {
    e <- stats/(sum(t1)+sum(t2))
    nc <- e * n_total * ((sum(t1)+sum(t2)) / (nrow(y1)+nrow(y2)))
    f <- function(nc){pchisq(qchisq(1-alpha, df), df, nc)}
    beta <- f(nc)
    # se <- sqrt((2 * df + 4 * (stats-df)) * (1/((sum(t1)+sum(t2))))^2 * (numDeriv::grad(f, nc) * n_total * ((sum(t1)+sum(t2)) / (nrow(y1)+nrow(y2))))^2)
    se <- sqrt((2*df + 4*stats) * (1/((sum(t1)+sum(t2))))^2 *
                 (numDeriv::grad(f, nc) * n_total * ((sum(t1)+sum(t2)) / (nrow(y1)+nrow(y2))))^2)
    power <- 1 - beta
    return(list('power' = round(power, digits = 3),
                'MC error of power' = round(se, digits = 3),
                'global deviation' = round(e, digits = 3),
                'noncentrality parameter' = round(nc, digits = 3) ))
  }

  call<-match.call()
  # thetas for rmvordlogis() for group 1 and group 2
  # a list with numeric vector elements, with free item-category parameters
  # last the discrimination parameter set to 1.
  func <- function(local_dev) lapply(local_dev, FUN = function(x) {append(x,1)})
  local_dev <- lapply(local_dev, func)


  y1 <- ltm::rmvordlogis(n = length(persons1), thetas = local_dev[[1]],
                    IRT = TRUE, model = "gpcm", link = "logit",
                    z.vals = persons1) - 1
  y2 <- ltm::rmvordlogis(n = length(persons2), thetas = local_dev[[2]],
                    IRT = TRUE, model = "gpcm", link = "logit",
                    z.vals = persons2) - 1

  df <- length(unlist(local_dev[[1]])) - length((local_dev[[1]])) - 1

  r1 <- psychotools::pcmodel(y1)
  r2 <- psychotools::pcmodel(y2)
  r3 <- psychotools::pcmodel(rbind(y1,y2))

  t1 <- table(factor(rowSums(y1),levels=1:df))
  t2 <- table(factor(rowSums(y2),levels=1:df))

  vstats <- do.call(c, invar_test_obj(r1=r1, r2=r2, r3=r3, model = "PCM"))

  res <- lapply(vstats, subfunc)

  results <- list(  'power' = unlist(do.call(cbind, res)[1,]),
                    'MC error of power' = unlist(do.call(cbind, res)[2,]),
                    'global deviation' = unlist(do.call(cbind, res)[3,]),
                    'local deviation' = round(rbind("group1"=r1$coefficients, "group2"= r2$coefficients), digits = 3),
                    'person score distribution in group 1' = round(t1/sum(t1), digits=3),
                    'person score distribution in group 2' = round(t2/sum(t2), digits=3),
                    'degrees of freedom' = df,
                    'noncentrality parameter' = unlist(do.call(cbind, res)[4,]),
                    "call" = call)
  return(results)
}
