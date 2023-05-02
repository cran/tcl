######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# invar_test_obj, invar_test_full, invar_test_RS
# tcl_fitobj
#
# Part of R/tlc - Testing in Conditional Likelihood Context package
#
# This file contains a routine that computes Wald (W), Score (RS),
# likelihood ratio (LR) and gradient (GR) test
# for hypothesis that the difference of two group parameters equals zero,
# where one assumes all items presented to all persons.
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2021, Last Modified 02/05/2023
######################################################################

invar_test_obj <- function(r1, r2, r3, model = "RM") {
  # r1, r2, r3 of fitted psychotools object

  y1 <- r1$data   # data group 1
  y2 <- r2$data   # data group 2

  ## CML estimates of eta parameters
  eta.rest     <- r3$coefficients     # CML eta estimates from total sample
  eta.unrest.1 <- r1$coefficients # CML eta estimates group 1
  eta.unrest.2 <- r2$coefficients # CML eta estimates group 2

  # score function and Hessian matrix evaluated at unrestricted eta estimates
  unrest.1 <- eRm_cml(X = y1, eta = eta.unrest.1, model = model)
  unrest.2 <- eRm_cml(X = y2, eta = eta.unrest.2, model = model)

  # score function and Hessian matrix evaluated at restricted eta estimates
  rest.1 <- eRm_cml(X = y1, eta = eta.rest, model = model)
  rest.2 <- eRm_cml(X = y2, eta = eta.rest, model = model)

  # Fisher informaion matrix
  i <- solve( solve(unrest.1$hessian ) + solve(unrest.2$hessian ))

  # difference of eta parameters in groups
  delta <- eta.unrest.1 - eta.unrest.2

  Wald <- sum( colSums(delta * i) * delta) # classical Wald test

  LR <- -2*(r3$loglik-r1$loglik-r2$loglik)

  RS <- sum (colSums( (rest.1$scorefun * solve(rest.1$hessian) )) * rest.1$scorefun) +
    sum (colSums( (rest.2$scorefun * solve(rest.2$hessian) )) * rest.2$scorefun)

  GR <- sum( (rest.1$scorefun * delta) )  # gradient test

  return(list("W"=Wald, "LR"=LR, "RS"=RS, "GR"=GR))

}


# subfunction switch for CML estimates
tcl_fitobj <- function(X, model) {

  model_switch <- function(option) {
  switch(option,
         "RM" =  eRm::RM( X, se = FALSE, sum0 = FALSE),
         "PCM" = eRm::PCM( X, se = FALSE, sum0 = FALSE),
         "RSM" = eRm::RSM( X, se = FALSE, sum0 = FALSE),
         stop("Invalid `option` value") )}

  return(obj = model_switch(option=model))
}


invar_test_full <- function(y, y1, y2, model, X, splitcr, del_pos) {
  ## CML estimates of eta parameters
  r  <- tcl_fitobj(X = y, model = model)
  r1 <- tcl_fitobj(X = y1, model = model)
  r2 <- tcl_fitobj(X = y2, model = model)

  eta.rest     <- r$etapar      # CML eta estimates from total sample
  eta.unrest.1 <- r1$etapar     # CML eta estimates group 1
  eta.unrest.2 <- r2$etapar     # CML eta estimates group 2

  if (model == "RM")  {
    df <- ncol(y1) - 1
  } else if (model == "PCM" || model == "RSM")  {
    df <- length(eta.rest)
  }


  ####################################
  ####      test statistics       ####
  ####################################

  # score function and Hessian matrix evaluated at unrestricted eta estimates
  unrest.1 <- eRm_cml(X = y1, eta = eta.unrest.1, model = model)
  unrest.2 <- eRm_cml(X = y2, eta = eta.unrest.2, model = model)

  # score function and Hessian matrix evaluated at restricted eta estimates
  rest.1 <- eRm_cml(X = y1, eta = eta.rest, model = model)
  rest.2 <- eRm_cml(X = y2, eta = eta.rest, model = model)

  # Fisher informaion matrix
  i <- solve( solve(unrest.1$hessian ) + solve(unrest.2$hessian ))

  # difference of eta parameters in groups
  delta <- eta.unrest.1 - eta.unrest.2

  Wald <- sum( colSums(delta * i) * delta) # classical Wald test

  GR <- sum( (rest.1$scorefun * delta) )  # gradient test

  LR <- -2*(r$loglik - r1$loglik - r2$loglik)

  # RS <- sum (colSums( (rest.1$scorefun * solve(rest.1$hessian) )) * rest.1$scorefun) +
  #   sum (colSums( (rest.2$scorefun * solve(rest.2$hessian) )) * rest.2$scorefun)

  RSobj <- RStest(X=X, splitcr = splitcr, model=model,del_pos=del_pos)
  RS <- RSobj$RS
  pvalue_RS <- RSobj$pvalue


  test.stats <- c( GR, LR, RS, Wald)
  names(test.stats) <- c("GR", "LR", "RS", "W")

  pvalue <- 1 - (sapply(test.stats, stats::pchisq, df = df))
  pvalue[3] <- pvalue_RS # if different df in full model

  pvalue <- pvalr(pvalue, digits = 3) # added AK 02.05.2023

  df_vec <- c(df,df,RSobj$df,df)
  names(df_vec) <- c("GR", "LR", "RS", "W")

  res.list <- list("test" = round(test.stats, digits = 3),
                   "df" = df_vec,
                   "pvalue" = pvalue)

  return(res.list)
}


RStest <- function(X,
                   splitcr ,
                   model,
                   ctr = c("psychotools", "eRm"),
                   eta_rest,
                   del_pos) {
  call <- match.call()
  if (length(ctr)==2) ctr <- ctr[1] # choose default = psychotools


  if (missing(del_pos)) del_pos <- tcl_datcheck_full(X,model)$del_pos
  if (any(!is.na(del_pos))) X <- X[,-del_pos]

  e <- tcl_splitcr_RStest(X = X, splitcr = splitcr, model = model)

  y  <- X
  y1 <- e$X_list[[1]]
  y2 <- e$X_list[[2]]

  if (model == "RM")  {

    if (missing(eta_rest)) {
      if (ctr=="psychotools") {
        r_rest <- psychotools::raschmodel(y, hessian = FALSE)
        eta_rest <-  c(0,-r_rest$coefficients)
      } else if (ctr=="eRm") {
        r_rest <- eRm::RM(y, se = FALSE, sum0 = FALSE)
        eta_rest <- r_rest$betapar
      }
    }# end missing

    df <- ncol(y) - 1

  } else if (model == "PCM" || model == "RSM")  {

    if (missing(eta_rest)) {
      if (ctr=="psychotools") {
        r_rest <- psychotools::pcmodel(y, hessian = FALSE)
        eta_rest <-  c(0,-r_rest$coefficients)
      } else if (ctr=="eRm") {
        r_rest <- eRm::PCM(X=y, se = FALSE, sum0 = FALSE)
        eta_rest <- r_rest$betapar
      }
    } # end missing

    df <- length(eta_rest) -1 # betapar
  } # end if

  # score function and Hessian matrix evaluated at restricted eta estimates
  rest_1 <- tcl_cml(X = y1, X0=y, eta = eta_rest) #, model = model)
  rest_2 <- tcl_cml(X = y2, X0=y, eta = eta_rest) #, model = model)


  RS <- sum (colSums( (rest_1$scorefun * solve(rest_1$hessian) )) * rest_1$scorefun) +
    sum (colSums( (rest_2$scorefun * solve(rest_2$hessian) )) * rest_2$scorefun)

  pvalue <- 1 - stats::pchisq(q=RS,df=df)


  res.list <- list( "RS" = round(RS, digits = 3),
                    "df" = df,
                    "pvalue" = round(pvalue, digits = 3),
                    "del_pos" = del_pos,
                    # "eta_rest" = round(eta_rest, digits = 3),
                    "call" = call)

  return(res.list)
}

#############################################
# depreciated functions from previous version
#############################################
# invar_test_full_old <- function(y, y1, y2, model) {
#   ## CML estimates of eta parameters
#   r  <- tcl_fitobj(X = y, model = model)
#   r1 <- tcl_fitobj(X = y1, model = model)
#   r2 <- tcl_fitobj(X = y2, model = model)
#
#   eta.rest     <- r$etapar      # CML eta estimates from total sample
#   eta.unrest.1 <- r1$etapar     # CML eta estimates group 1
#   eta.unrest.2 <- r2$etapar     # CML eta estimates group 2
#
#   if (model == "RM")  {
#     df <- ncol(y1) - 1
#   } else if (model == "PCM" || model == "RSM")  {
#     df <- length(eta.rest)
#   }
#
#
#   ####################################
#   ####      test statistics       ####
#   ####################################
#
#   # score function and Hessian matrix evaluated at unrestricted eta estimates
#   unrest.1 <- eRm_cml(X = y1, eta = eta.unrest.1, model = model)
#   unrest.2 <- eRm_cml(X = y2, eta = eta.unrest.2, model = model)
#
#   # score function and Hessian matrix evaluated at restricted eta estimates
#   rest.1 <- eRm_cml(X = y1, eta = eta.rest, model = model)
#   rest.2 <- eRm_cml(X = y2, eta = eta.rest, model = model)
#
#   # Fisher informaion matrix
#   i <- solve( solve(unrest.1$hessian ) + solve(unrest.2$hessian ))
#
#   # difference of eta parameters in groups
#   delta <- eta.unrest.1 - eta.unrest.2
#
#   Wald <- sum( colSums(delta * i) * delta) # classical Wald test
#
#   GR <- sum( (rest.1$scorefun * delta) )  # gradient test
#
#   LR <- -2*(r$loglik - r1$loglik - r2$loglik)
#
#   RS <- sum (colSums( (rest.1$scorefun * solve(rest.1$hessian) )) * rest.1$scorefun) +
#     sum (colSums( (rest.2$scorefun * solve(rest.2$hessian) )) * rest.2$scorefun)
#
#   test.stats <- c(GR, LR, RS, Wald)
#   names(test.stats) <- c("GR", "LR", "RS", "W")
#
#   pvalue <- 1 - (sapply(test.stats, stats::pchisq, df = df))
#
#   res.list <- list(test = round(test.stats, digits = 3), df = df,
#                    pvalue = round(pvalue, digits = 3))
#
#   return(res.list)
# }
#
# invar_test_RS <- function(y, y1, y2, model) {
#   ## CML estimates of eta parameters
#   r  <- tcl_fitobj(X = y, model = model)
#   # r1 <- tcl_fitobj(X = y1, model = model)
#   # r2 <- tcl_fitobj(X = y2, model = model)
#
#   eta.rest     <- r$etapar      # CML eta estimates from total sample
#   # eta.unrest.1 <- r1$etapar     # CML eta estimates group 1
#   # eta.unrest.2 <- r2$etapar     # CML eta estimates group 2
#
#   if (model == "RM")  {
#     df <- ncol(y) - 1
#   } else if (model == "PCM" || model == "RSM")  {
#     df <- length(eta.rest)
#   }
#
#
#   ####################################
#   ####      test statistics       ####
#   ####################################
#
#   # score function and Hessian matrix evaluated at unrestricted eta estimates
#   # unrest.1 <- eRm_cml(X = y1, eta = eta.unrest.1, model = model)
#   # unrest.2 <- eRm_cml(X = y2, eta = eta.unrest.2, model = model)
#
#   # score function and Hessian matrix evaluated at restricted eta estimates
#   rest.1 <- eRm_cml(X = y1, eta = eta.rest, model = model)
#   rest.2 <- eRm_cml(X = y2, eta = eta.rest, model = model)
#
#   # Fisher informaion matrix
#   # i <- solve( solve(unrest.1$hessian ) + solve(unrest.2$hessian ))
#
#   # difference of eta parameters in groups
#   # delta <- eta.unrest.1 - eta.unrest.2
#
#   # Wald <- sum( colSums(delta * i) * delta) # classical Wald test
#   #
#   # GR <- sum( (rest.1$scorefun * delta) )  # gradient test
#   #
#   # LR <- -2*(r$loglik - r1$loglik - r2$loglik)
#
#   RS <- sum (colSums( (rest.1$scorefun * solve(rest.1$hessian) )) * rest.1$scorefun) +
#     sum (colSums( (rest.2$scorefun * solve(rest.2$hessian) )) * rest.2$scorefun)
#
#
#   # test.stats <- c( GR, LR, RS, Wald)
#   test.stats <- c( NA, NA, RS, NA)
#   names(test.stats) <- c("GR", "LR", "RS", "W")
#
#   pvalue <- 1 - (sapply(test.stats, stats::pchisq, df = df))
#
#
#   res.list <- list("test" = round(test.stats, digits = 3),
#                    "df" = df,
#                    "pvalue" = round(pvalue, digits = 3))
#
#   return(res.list)
# }
