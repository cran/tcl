######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# invar_test
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
#' Test of invariance of item parameters between two groups.
#'
#' Computes gradient (GR), likelihood ratio (LR), Rao score (RS) and Wald (W) test statistics
#'   for hypothesis of equality of item parameters between two groups of persons against a two-sided
#'  alternative that at least one item parameter differs between the two groups.
#'
#'  Note that items are excluded for the computation of GR,LR, and W due to inappropriate
#'  response patterns within subgroups and for computation of RS due to inappropriate
#'  response patterns in the total data. If the model is identified from the total data but not from one
#'  or both subgroups only RS will be computed. If the model is not identified from the total data,
#'  no test statistic is computable.
#'
#' @param X Data matrix.
#' @param splitcr Split criterion which is either "mean", "median" or a numeric vector x.
#' \describe{
#'  \item{"mean"}{Corresponds to division of the sample according to the mean of the person score.}
#'  \item{"median"}{Corresponds to division of the sample according to the median of the person score.}
#'   \item{x}{Has length equal to number of persons and contains zeros and ones. It indicates group membership for every person.}
#'   }
#'
#' @param model RM, PCM, RSM
#' @return A list of test statistics, degrees of freedom, and p-values.
#'  \item{test}{A numeric vector of gradient (GR), likelihood ratio (LR), Rao score (RS), and Wald test statistics.}
#'  \item{df}{A numeric vector of corresponding degrees of freedom.}
#'  \item{pvalue}{A vector of corresponding p-values.}
#'  \item{deleted_items}{A list with numeric vectors of item numbers that were excluded before computing corresponding test statistics.}
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
#' @keywords htest
#' @export
#' @seealso \code{\link{change_test}}, and \code{\link{LLTM_test}}.
#' @examples
#' \dontrun{
#'##### Rasch Model #####
#'y <- eRm::sim.rasch(persons = rnorm(400), c(0,-3,-2,-1,0,1,2,3))
#'x <- c(rep(1,200),rep(0,200))
#'
#'res <- invar_test(y, splitcr = x, model = "RM")
#'
#'res$test # test statistics
#'res$df # degrees of freedoms
#'res$pvalue # p-values
#'res$deleted_items # excluded items
#'
#'$test
#'   GR    LR    RS     W
#'14.492 14.083 13.678 12.972
#'
#'$df
#'GR LR RS  W
#' 7  7  7  7
#'
#'$pvalue
#'   GR    LR    RS     W
#'"0.043" "0.050" "0.057" "0.073"
#'
#'$deleted_items
#'  $deleted_items$GR
#'  [1] "none"
#'
#'  $deleted_items$LR
#'  [1] "none"
#'
#'  $deleted_items$RS
#'  [1] "none"
#'
#'  $deleted_items$W
#'  [1] "none"
#'
#'
#'$call
#'invar_test(X = y, splitcr = x, model = "RM")
#'
#'}

invar_test <- function(X, splitcr = "median", model = "RM"){
  # X = observed data matrix comprised of k columns
  # splitcr... splitting criterion for 2 covariate groups.
  #   "median"  corresponds to a median person score split,
  #   "mean" corresponds to the mean person score split.
  #    vector of length n containing zeros or ones only for sample split
  #        (group 1 = '1', group 2 = '0')
  # model = RM, PCM, RSM

  call <- match.call()

  #---------------------------------------------------------------------
  # check of data matrix X for
  #    - inappropriate response patterns within subgroups
  #    - ill-conditioned data matrix (restricted and unrestricted models)
  #---------------------------------------------------------------------

  Xcheck<- tcl_datcheck_full(X,model)
  if (Xcheck$Xcheck == "none") return(list("test" = NA, "df" = NA,"pvalue" = NA,
                                           "deleted_items"=NA, "call" = call))
  del_pos_full <- Xcheck$del_pos  # check full model
  # if (any(!is.na(del_pos_full))) X <- X[,-del_pos_full]

  e <- tcl_splitcr(X = X, splitcr = splitcr, model = model)
  Xlist_check <- tcl_datcheck(X=X, Xlist = e$X.list, model = model)

  y  <- e$X.el
  y1 <- e$X.list[[1]]
  y2 <- e$X.list[[2]]

  ###################################################
  #               main programm
  ###################################################

  test_switch <- function(option) {
    switch(option,
           "full" = invar_test_full(y=y, y1=y1, y2=y2, model=model, X=X, splitcr = splitcr, del_pos=del_pos_full),
           # "RS" = invar_test_RS(y, y1, y2, model),
           "RS" =  RStest(X=X, splitcr = splitcr, model=model,del_pos=del_pos_full), # addded AK 20-02-2022
           "none" = list("test" = NA, "df" = NA,"pvalue" = NA,"deleted_items"=NA, "call" = call),
           stop("Invalid `option` value")
    )
  }

  res.list <- test_switch(option = Xlist_check)

  if(Xlist_check=="full") {
    if(is.na(e$del_pos)) {
      e$del_pos <- "none"
    } else {
      e$del_pos <-  paste0("I", e$del_pos)
    }

    if(is.na(del_pos_full))  {
      del_pos_full <- "none"
    } else {
      del_pos_full <- paste0("I", del_pos_full)
    }

    res.list$deleted_items <- list("GR"=e$del_pos, "LR"=e$del_pos, "RS"=del_pos_full,"W"=e$del_pos) # added AK 20-02-2022
  }

  if(Xlist_check=="RS") {
    test.stats <- c( NA, NA, res.list$RS, NA)
    names(test.stats) <- c("GR", "LR", "RS", "W")

    df_vec <- c(NA,NA,res.list$df,NA)
    names(df_vec) <- c("GR", "LR", "RS", "W")

    pvalue <- c(NA,NA,res.list$pvalue,NA)
    names(pvalue) <- c("GR", "LR", "RS", "W")

    if(is.na(del_pos_full))  {
      del_pos_full <- "none"
    } else {
      del_pos_full <- paste0("I", del_pos_full)
    }
    pvalue <- pvalr(pvalue, digits = 3)

    res.list <- list("test" = round(test.stats, digits = 3),
                     "df" = df_vec,
                     "pvalue" = pvalue,
                     "deleted_items" = list(  "GR"=NA,"LR"=NA,"RS"=del_pos_full,"W"=NA)) # addded AK 20-02-2022
  }

  res.list$call <- call
  return(res.list)
}
