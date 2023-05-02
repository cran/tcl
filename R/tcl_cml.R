######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# tcl_cml
#
# Part of R/tlc - Testing in Conditional Likelihood Context package
#
# This file contains a routine related to computation of
# score funtion and Hessian matrix.
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2022, Last Modified 19/02/2022
######################################################################
#' @importFrom numDeriv hessian jacobian

tcl_cml <- function(X, X0, eta) {
  # beta_cml <- function(X, eta, model) {
  # X = observed data matrix
  # eta - item easiness parameter
  # model = RM, PCM, RSM, LLTM

  # ------------------ likelihood function -----------------
  # R code extracted from file functions_shoot.R
  # unpublished R package cmlPCMlasso_0.1.0 (Lasso Penalisierungen fürs PCM im CML Kontext)
  # Author: Can Gürer, UMIT Tirol
  # functions used to build total log.likelihood cltot:
  #     - cl_person_PCM(),
  #     - make_min0(),
  #     - make_catmax
  # Y Matrix (n x I) of responses.
  # y Vector of responses (lowest categories need to be 0).
  # sigma_p Vector of person specific parameters (sum parametrization).
  # catmax: vector with max categories per item (starting with cat 0)
  # no missings allowed

  make_min0 <- function(Y){
    mins <- apply(Y,2,function(x)min(x,na.rm=TRUE))
    return(t(t(Y) - mins))
  }

  # extract vector of max categories (after making min categories zero)
  make_catmax <- function(Y) {
    Y <- make_min0(Y)
    return(apply(Y,2,function(x)max(x,na.rm=TRUE)))
  }


  cl_person_PCM <- function(y, sigma_p, catmax){

    sigmas_p_list <- split(sigma_p,rep(1:length(catmax),catmax))
    non_missing <- !is.na(y)
    y <- y[non_missing]
    catmax <- catmax[non_missing]
    sigmas_p_list <- sigmas_p_list[non_missing]

    if((sum(y)==0)) {return(0)}
    if(sum(y)==sum(catmax)) {return(0)}

    esf_rp <- psychotools::elementary_symmetric_functions(lapply(sigmas_p_list,'-'),order = 0)

    positions <- (c(0,cumsum(catmax)[-length(catmax)]) + y)[y!=0] # exclude 0 answers
    selected_sigmas <- unlist(sigmas_p_list)[positions]

    cl <- sum(selected_sigmas) - log(esf_rp[[1]][sum(y)+1])
    return(cl)
  }

  cml <-function(X, sigma_p, catmax) {

    cltot <- sum(apply(X, 1, cl_person_PCM,sigma_p = sigma_p, catmax = catmax))
    return(cltot)
  }
  # ------------------ likelihood function ------------------



  ######################################################################
  # calculate score function and Hessian matrix using numDeriv package
  ######################################################################
  # catmax <-  make_catmax(Y=X)
  catmax <-  make_catmax(Y=X0) # restricted model
  sigma_p <- eta

  loglik <- cml(X=X,sigma_p=sigma_p,catmax=catmax) # conditional log likelihood evaluated at eta

  scorefun <- as.vector(numDeriv::jacobian(cml, x=sigma_p, X=X,catmax = catmax))
  r.hessian <- -numDeriv::hessian(cml, x=sigma_p, X=X,catmax = catmax)  # hessian * (-1) AK

  par.list <-      list( "eta" = eta,
                         "catmax" = catmax,
                         "call" = call)

  return( list( "scorefun" = scorefun,
                "hessian" = r.hessian,
                "loglik" = loglik,
                "par.list" = par.list,
                "call" = call))
}
