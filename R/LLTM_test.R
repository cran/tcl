######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# LLTM_test
#
# Part of R/tlc - Testing in Conditional Likelihood Context package
#
# This file contains a routine that computes Wald (W), Score (RS),
# likelihood ratio (LR) and gradient (GR) test
# for hypothisized linear restriction of item parameter space in Rasch model
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2019, Last Modified 02/05/2023
######################################################################
#' Testing linear restrictions on parameter space of item parameters of RM.
#'
#' Computes gradient (GR), likelihood ratio (LR), Rao score (RS) and Wald (W) test statistics for
#'   hypotheses defined by linear restrictions on parameter space of the item parameters of RM.
#'
#' The RM item parameters are assumed to be linear in the LLTM parameters.
#'  The coefficients of the linear functions are specified by a design matrix W. In this context,
#'  the LLTM is considered as a more parsimonious model than the RM. The LLTM parameters can be
#'  interpreted as the difficulties of certain cognitive operations needed to respond correctly
#'  to psychological test items. The item parameters of the RM are assumed to be linear combinations
#'  of these cognitive operations. These linear combinations are defined in the design matrix W.
#'
#' @param X Data matrix.
#' @param W Design matrix of LLTM.
#' @return A list of test statistics, degrees of freedom, and p-values.
#'  \item{test}{A numeric vector of gradient (GR), likelihood ratio (LR), Rao score (RS), and Wald test statistics.}
#'  \item{df}{Degrees of freedom.}
#'  \item{pvalue}{A vector of corresponding p-values.}
#'  \item{call}{The matched call.}
#' @references{
#'  Fischer, G. H. (1995). The Linear Logistic Test Model. In G. H. Fischer & I. W. Molenaar (Eds.),
#'  Rasch models: Foundations, Recent Developments, and Applications (pp. 131-155). New York: Springer.
#'
#'  Fischer, G. H. (1983). Logistic Latent Trait Models with Linear Constraints. Psychometrika, 48(1), 3-26.
#'  }
#' @keywords htest
#' @export
#' @seealso \code{\link{change_test}}, and \code{\link{invar_test}}.
#' @examples
#' \dontrun{
#' # Numerical example assuming no deviation from linear restriction
#'
#' # design matrix W defining linear restriction
#' W <- rbind(c(1,0), c(0,1), c(1,1), c(2,1))
#'
#'# assumed eta parameters of LLTM for data generation
#' eta <- c(-0.5, 1)
#'
#' # assumed vector of item parameters of RM
#' b <- colSums(eta * t(W))
#'
#' y <- eRm::sim.rasch(persons = rnorm(400), items = b - b[1])  # sum0 = FALSE
#'
#' res <- LLTM_test(X = y, W = W )
#'
#' res$test # test statistics
#' res$df # degrees of freedoms
#' res$pvalue # p-values
#'
#'}

LLTM_test <- function(X, W) {
  # X = observed data matrix
  # W design matrix

  call<-match.call()

  if (missing(W)) stop('Design matrix W is missing')
  else W <- as.matrix(W)

  model <- "RM"

  ###################################################
  #               main programm
  ###################################################

  y <- X

  r.RM <- eRm::RM( y, se = FALSE, sum0 = FALSE)                 # general model
  r.LLTM  <- eRm::LLTM( y,  W = W,  se = FALSE, sum0 = FALSE)   # specific model nested in RM

  #####################################################
  ########    compute  GR, LR, RS and Wald test      ##
  #####################################################

  # restricted item parameter estimates of RM (linear restriction imposed by W)
  e <- r.LLTM$etapar  # LLTM estimates
  eta.rest <- (colSums(e * t(W)) - colSums(e * t(W))[1])[2:nrow(W)] * -1

  eta.unrest <- r.RM$etapar  # unrestricted CML estimates of item parameters of RM

  ## score function and Hessian matrix evaluated at unrestricted estimates of item paramaters in RM
  unrest.y <- eRm_cml(X = y, eta = eta.unrest, model = "RM")

  ## score function and Hessian matrix evaluated at restricted estimates of item paramaters in RM
  rest.y <- eRm_cml(X = y, eta = eta.rest, model = "RM")

  # Fisher informaion matrix evaluated at unrestricted estimates of item paramaters in RM
  i <-  unrest.y$hessian

  # difference of eta parameters
  delta <-  eta.unrest - eta.rest

  Wald <- sum(colSums( delta * i) * delta) # classical Wald test

  GR <- sum(rest.y$scorefun * eta.unrest)  # gradient test statistic

  LR <- -2 * (r.LLTM$loglik - r.RM$loglik)  # likelihood ratio test statistic

  RS <- sum (colSums((rest.y$scorefun * solve(rest.y$hessian))) * rest.y$scorefun)

  test.stats <- c( GR, LR, RS, Wald)
  names(test.stats) <- c("GR", "LR", "RS", "W")

  df <- length(r.RM$etapar) - length(r.LLTM$etapar)
  pvalue <- 1 - (sapply(test.stats, stats::pchisq, df = df))

  pvalue <- pvalr(pvalue, digits = 3) # added 02.05.2023 AK

  res.list <- list("test" = round(test.stats, digits = 3),
                   "df" = df,
                   "pvalue" = pvalue,
                   "call" = call)

  return(res.list)
}
