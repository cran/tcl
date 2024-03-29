######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# sa_sizeRM
#
# Part of R/tlc - Testing in Conditional Likelihood Context package
#
# This file contains a routine that returns sample size for W,LRT, RS, GR
# given probabilities of errors of first and second kinds alpha and beta,
# and a given deviation from the hypothesis to be tested,
# i.e. that the difference of two group parameters equals zero,
# where one assumes all items presented to all persons.
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2021, Last Modified 25/10/2021
######################################################################
#' Sample size planning for tests of invariance of item parameters between two groups of persons in binary Rasch model
#'
#' Returns sample size for Wald (W), likelihood ratio (LR), Rao score (RS)
#' and gradient (GR) test given probabilities of errors of first and second kinds \eqn{\alpha} and
#' \eqn{\beta} as well as a deviation from the hypothesis to be tested. The hypothesis to be
#' tested assumes equal item parameters between two predetermined groups of persons. The alternative assumes
#' that at least one parameter differs between the two groups.
#'
#' In general, the sample size is determined from the assumption that the approximate distributions of
#' the four test statistics are from the family of noncentral \eqn{\chi^2} distributions with \eqn{df}
#' equal to the number of items minus 1, and noncentrality parameter \eqn{\lambda}. The latter is,
#' inter alia, a function of the sample size. Hence, the sample size can be determined from the condition
#' \eqn{\lambda = \lambda_0}, where \eqn{\lambda_0} is a predetermined constant  which depends on the probabilities of
#' the errors of the first and second kinds \eqn{\alpha} and \eqn{\beta}
#' (or power). More details about the distributions of the test statistics and the relationship between \eqn{\lambda},
#' power, and sample size can be found in Draxler and Alexandrowicz (2015).
#'
#'In particular, the determination of \eqn{\lambda} and the sample size, respectively, is based on a simple
#'Monte Carlo approach. As regards the concept of sample size a distinction between informative and total
#'sample size has to be made. In the conditional maximum likelihood context, the responses of persons
#'with minimum or maximum person score are completely uninformative. They do not contribute to the value o
#'f the test statistic. Thus, the informative sample size does not include these persons. The total sample
#'size is composed of all persons. The Monte Carlo approach used in the present problem to determine \eqn{\lambda}
#'and informative (and total) sample size can briefly be described as follows. Data (responses of a large number
#'of persons to a number of items) are generated given a user-specified scenario of a deviation from the hypothesis
#'to be tested. The hypothesis to be tested assumes equal item parameters between the two groups of persons.
#'A scenario of a deviation is given by a choice of the item parameters and the person parameters (to be drawn
#'randomly from a specified distribution) for each of the two groups. Such a scenario may be called local
#'deviation since deviations can be specified locally for each item. The relative group sizes are determined by
#'the choice of the number of person parameters for each of the two groups. For instance, by default \eqn{10^6} person
#'parameters are selected randomly for each group. In this case, it is implicitly assumed that the two groups of
#'persons are of equal size. The user can specify the relative groups sizes by choosing the lengths of the
#'arguments \code{persons1} and \code{persons2} appropriately. Note that the relative group sizes do have an impact on power
#'and sample size of the tests. The next step is to compute a test statistic \eqn{T} (Wald, LR, score, or gradient)
#'from the simulated data. The observed value \eqn{t} of the test statistic is then divided by the informative
#'sample size \eqn{n_{infsim}} observed in the simulated data. This yields the so-called global deviation
#'\eqn{e = t / n_{infsim}}, i.e., the chosen scenario of a deviation from the hypothesis to be tested being
#'represented by a single number. Let the informative sample size sought be denoted by \eqn{n_{inf}} (thus, this is
#'not the informative sample size observed in the sim. data). The noncentrality parameter \eqn{\lambda} can
#'be expressed by the product \eqn{n_{inf} * e}. Then, it follows from the condition \eqn{\lambda = \lambda_0} that
#'
#'\deqn{n_{inf} * e = \lambda_0}
#'
#'and
#'
#'\deqn{n_{inf} = \lambda_0 / e.}
#'
#'Note that the sample of size \eqn{n_{inf}} is assumed to be composed only of persons with informative person scores in both groups,
#'where the relative frequency distribution of these informative scores in each of both groups is considered to be equal
#'to the observed relative frequency distribution of informative scores in each of both groups in the simulated data. Note also that the
#'relative sizes of the two person groups are assumed to be equal to the
#'relative sizes of the two groups in the simulated data. By default, the two groups are equal-sized in the simulated
#'data, i.e., one yields \eqn{n_{inf} / 2} persons (with informative scores) in each of the two groups. The total
#'sample size \eqn{n_{total}} is obtained from the relation
#'\eqn{n_{inf} = n_{total} * pr}, where \eqn{pr} is the proportion or relative frequency of persons observed
#'in the simulated data with a minimum or maximum score. Basing the tests given a level \eqn{\alpha} on an informative
#'sample of size \eqn{n_{inf}} the probability of rejecting the hypothesis to be tested will be at least
#'\eqn{1 - \beta} if the true global deviation \eqn{\ge e}.
#'
#'Note that in this approach the data have to be generated only once. There are no replications needed. Thus, the
#'procedure is computationally not very time-consuming.
#'
#'Since \eqn{e} is determined from the value of the test statistic observed in the simulated data it has to be
#'treated as a realization of a random variable \eqn{E}. Consequently, \eqn{n_{inf}} is also a realization of a
#'random variable \eqn{N_{inf}}. Thus, the (realized) value \eqn{n_{inf}} need not be equal to the exact value of
#'the informative sample size that follows from the user-specified (predetermined) \eqn{\alpha}, \eqn{\beta}, and
#'scenario of a deviation from the hypothesis to be tested, i.e., the selected item parameters used for the
#'simulation of the data. If the CML estimates of these parameters computed from the simulated data are close
#'to the predetermined parameters \eqn{n_{inf}} will be close to the exact value. This will generally be the case
#'if the number of person parameters used for simulating the data, i.e., the lengths of the vectors \code{persons1}
#'and \code{persons2}, is large, e.g., \eqn{10^5} or even \eqn{10^6} persons. In such cases, the possible random
#'error of the computation procedure of \eqn{n_{inf}} based on the sim. data may not be of practical relevance any
#'more. That is why a large number (of persons for the simulation process) is generally recommended.
#'
#'For theoretical reasons, the random error involved in computing \eqn{n_{inf}} can be pretty well approximated.
#'A suitable approach is the well-known delta method. Basically, it is a Taylor polynomial of first order, i.e.,
#'a linear approximation of a function. According to it the variance of a function of a random variable can be
#'linearly approximated by multiplying the variance of this random variable with the square of the first
#'derivative of the respective function. In the present problem, the variance of the test statistic \eqn{T} is
#'(approximately) given by the variance of a noncentral \eqn{\chi^2} distribution.
#'Thus, \eqn{Var(T) = 2 (df + 2 \lambda)}, with \eqn{df} equal to the number of items minus 1 and
#'\eqn{\lambda = t}. Since the global deviation \eqn{e = (1 / n_{infsim}) * t} it
#'follows for the variance of the corresponding random variable \eqn{E} that \eqn{Var(E) = (1 / n_{infsim})^2 * Var(T)}.
#'Since \eqn{n_{inf} = f(e) = \lambda_0 / e} one obtains by the delta method (for the variance of the
#'corresponding random variable \eqn{N_{inf}})
#'
#'\deqn{Var(N_{inf}) = Var(E) * (f'(e))^2,}
#'
#'where \eqn{f'(e) = - \lambda_0 / e^2} is the derivative of \eqn{f(e)}. The square root of
#'\eqn{Var(N_{inf})} is then used to quantify the random error of the suggested Monte Carlo
#'computation procedure. It is called Monte Carlo error of informative sample size.
#'
#' @param alpha Probability of the error of first kind.
#' @param beta Probability of the error of second kind.
#' @param persons1 A vector of person parameters for group 1 (drawn from a specified distribution). By default
#' \eqn{10^6} parameters are drawn at random from the standard normal distribution. The larger this
#' number the more accurate are the computations. See Details.
#' @param persons2 A vector of person parameters for group 2 (drawn from a specified distribution). By default
#' \eqn{10^6} parameters are drawn at random from the standard normal distribution. The larger this
#' number the more accurate are the computations. See Details.
#' @param local_dev A list consisting of two vectors containing item parameters for the two person groups
#' representing a deviation from the hypothesis to be tested locally per item.
#'
#'@return A list of results.
#'  \item{informative sample size}{Informative sample size for each test omitting persons with min. and max score.}
#'  \item{MC error of sample size}{Monte Carlo error of informative sample size for each test.}
#'  \item{global deviation}{Global deviation computed from simulated data. See Details.}
#'  \item{local deviation}{CML estimates of free item parameters in both groups obtained from the simulated data.
#'  First item parameter set 0 in both groups.}
#'  \item{person score distribution in group 1}{Relative frequencies of person scores in group 1 observed in simulated data.
#'  Uninformative scores, i.e., minimum and maximum score, are omitted.
#'  Note that the person score distribution does also have an influence on the sample size.}
#'  \item{person score distribution in group 2}{Relative frequencies of person scores in group 2 observed in simulated data.
#'  Uninformative scores, i.e., minimum and maximum score, are omitted.
#'  Note that the person score distribution does also have an influence on the sample size.}
#'  \item{degrees of freedom}{Degrees of freedom \eqn{df}.}
#'  \item{noncentrality parameter}{Noncentrality parameter \eqn{\lambda} of \eqn{\chi^2} distribution from which sample size is determined.}
#'  \item{total sample size in group 1}{Total sample size in group 1 for each test. See Details.}
#'  \item{total sample size in group 1}{Total sample size in group 2 for each test. See Details.}
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
#'  }
#' @keywords sample_size_planning
#' @export
#' @seealso \code{\link{powerRM}}, and \code{\link{post_hocRM}}.
#' @examples
#' \dontrun{
#'##### Sample size of Rasch Model #####
#'
#' res <-  sa_sizeRM(local_dev = list( c(0, -0.5, 0, 0.5, 1) , c(0, 0.5, 0, -0.5, 1)))
#'
#'# > res
#'# $`informative sample size`
#'#   W  LR  RS  GR
#'# 159 153 155 151
#'#
#'# $`MC error of sample size`
#'#     W    LR    RS    GR
#'# 0.721 0.682 0.695 0.670
#'#
#'# $`global deviation`
#'#     W    LR    RS    GR
#'# 0.117 0.122 0.120 0.123
#'#
#'# $`local deviation`
#'#         Item2  Item3  Item4 Item5
#'# group1 -0.502 -0.005  0.497 1.001
#'# group2  0.495 -0.006 -0.501 0.994
#'#
#'# $`person score distribution in group 1`
#'#
#'#     1     2     3     4
#'# 0.249 0.295 0.268 0.188
#'#
#'# $`person score distribution in group 2`
#'#
#'#     1     2     3     4
#'# 0.249 0.295 0.270 0.187
#'#
#'# $`degrees of freedom`
#'# [1] 4
#'#
#'# $`noncentrality parameter`
#'# [1] 18.572
#'#
#'# $`total sample size in group 1`
#'#  W LR RS GR
#'# 97 93 94 92
#'#
#'# $`total sample size in group 2`
#'#  W LR RS GR
#'# 97 93 94 92
#'#
#'# $call
#'# sa_sizeRM(local_dev = list(c(0, -0.5, 0, 0.5, 1),
#'#                            c(0, 0.5, 0, -0.5, 1)))
#' }

sa_sizeRM <- function(alpha = 0.05, beta = 0.05, persons1 = rnorm(10^6),
                      persons2 = rnorm(10^6), local_dev){

  subfunc <- function(stats) {
    e <- stats /(sum(t1)+sum(t2))
    n <- ceiling(lambda0 / e)
    se <- sqrt((2*df + 4*stats) * lambda0^2 * (sum(t1) + sum(t2))^2 * (stats)^-4)
    n1 <- ceiling((n/2)/(sum(t1)/nrow(y1)))
    n2 <- ceiling((n/2)/(sum(t2)/nrow(y2)))
    return(list('informative sample size' = n,
                'MC error of sample size' = round(se, digits = 3),
                'global deviation' = round(e, digits = 3),
                'total sample size in group 1' = n1,
                'total sample size in group 2' = n2 ))
  }
  call<-match.call()

  y1 <- eRm::sim.rasch(persons = persons1, items = local_dev[[1]])
  y2 <- eRm::sim.rasch(persons = persons2, items = local_dev[[2]])

  func <- function (x) {beta - pchisq(qchisq(1-alpha,ncol(y1)-1), ncol(y1)-1, ncp = x)}
  lambda0 <- uniroot(f = func, interval = c(0,10^3), tol=.Machine$double.eps^0.5)$root

  r1 <- psychotools::raschmodel(y1)
  r2 <- psychotools::raschmodel(y2)
  r3 <- psychotools::raschmodel(rbind(y1,y2))

  df <- ncol(y1) - 1

  t1 <- table(factor(rowSums(y1),levels=1:(ncol(y1)-1)))
  t2 <- table(factor(rowSums(y2),levels=1:(ncol(y2)-1)))

  vstats <- do.call(c, invar_test_obj(r1=r1, r2=r2, r3=r3, model = "RM"))

  res <- lapply(vstats, subfunc)

  results <- list(  'informative sample size' = unlist(do.call(cbind, res)[1,]),
                    'MC error of sample size' = unlist(do.call(cbind, res)[2,]),
                    'global deviation' = unlist(do.call(cbind, res)[3,]),
                    'local deviation' = round(rbind("group1"=r1$coefficients, "group2"= r2$coefficients), digits = 3),
                    'person score distribution in group 1' = round(t1/sum(t1), digits=3),
                    'person score distribution in group 2' = round(t2/sum(t2), digits=3),
                    'degrees of freedom' = df,
                    'noncentrality parameter' = round(lambda0, digits = 3),
                    'total sample size in group 1' = unlist(do.call(cbind, res)[4,]),
                    'total sample size in group 2' = unlist(do.call(cbind, res)[5,]),
                    "call" = call)
  return(results)
}
