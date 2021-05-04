######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# LLTM_test_obj
#
# Part of R/tlc - Testing in Conditional Likelihood Context package
#
# This file contains a routine that computes Wald (W), Score (RS),
# likelihood ratio (LR) and gradient (GR) test
# for hypothisized linear restriction of item parameter space in Rasch model
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2019, Last Modified 15/03/2021
######################################################################

LLTM_test_obj <- function(X, W, r.RM, r.LLTM) {
  # r.RM = general model
  # r.LLTM = specific model nested in RM
  # X = observed data matrix
  # W design matrix

  # call<-match.call()
  #
  if (missing(W)) stop('Design matrix W is missing')
  else W <- as.matrix(W)

  model <- "RM"

  ###################################################
  #               main programm
  ###################################################

  y <- X

  # r.RM <- eRm::RM( y, se = FALSE, sum0 = FALSE)                 # general model
  # r.LLTM  <- eRm::LLTM( y,  W = W,  se = FALSE, sum0 = FALSE)   # specific model nested in RM

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

  # test.stats <- c( GR, LR, RS, Wald)
  # names(test.stats) <- c("GR", "LR", "RS", "W")
  #
  # df <- length(r.RM$etapar) - length(r.LLTM$etapar)
  # pvalue <- 1 - (sapply(test.stats, stats::pchisq, df = df))

  # res.list <- list("test" = round(test.stats, digits = 3),
  #                  "df" = df,
  #                  "pvalue" = round(pvalue, digits = 3),
  #                  "call" = call)
  #
  # return(res.list)

  return(list("W"=Wald, "LR"=LR, "RS"=RS, "GR"=GR))
}
