######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# tcl_datcheck
#
# Part of R/tlc - Testing in Conditional Likelihood Context package
#
# function to check for ill-conditioned data in the RM
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2021, Last Modified 25/04/2021
######################################################################


tcl_datcheck <- function(X, Xlist, model = "RM") {

  # groupvec <- 1
  # mpoints <- 1
  # W <- NA
  # XWcheck <- suppressWarnings( datcheck ( X = X , W, mpoints, groupvec, model)  )   # inital check of X and W
  # X1Wcheck <- suppressWarnings( datcheck ( X = Xlist[[1]] , W, mpoints, groupvec, model) )
  # X2Wcheck <- suppressWarnings( datcheck ( X = Xlist[[2]] , W, mpoints, groupvec, model) )

  if (model == "RM") {

    func <- function(X) "Rm" %in% class(try(eRm::RM(X), silent = TRUE))

    #### Check for ill-conditioned data matrix X  in full model ######
    XWcheck <- !func(X)

        # if (XWcheck$ill_conditioned == TRUE ){
      if (XWcheck){
            warning(paste0(
              "\n",
              "\n",
              prettyPaste("Estimation stopped due to ill-conditioned data matrix X! Suspicious items in full model):"),
              "\n",
              "\n"
              # paste("No Estimation not possibke in full model!", collapse=" ")
              ),
              call. = FALSE, immediate.=TRUE)
            return ("none")
       } # end if

    #### Check for ill-conditioned data matrix X  in subgroup models ######
    X1Wcheck <- !func(X = Xlist[[1]])
    X2Wcheck <- !func(X = Xlist[[2]])

       if (X1Wcheck  || X2Wcheck){
          warning(paste0(
            "\n",
            "\n",
            prettyPaste("Estimation stopped due to ill-conditioned data matrix X! Suspicious items in group 1 and/or 2):"),
            "\n",
            "\n",
            paste("Estimation in full model only!", collapse=" ") ),
            call. = FALSE, immediate.=TRUE)
          return("RS")
        } # end if
  } # end if

  return("full")
}
