######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
#  get_tcl_arg, get_eRm_arg, pvalr, n_info
#
# Part of R/tlc - Testing in Conditional Likelihood Context package
#
# This file contains routines that function as helper functions
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2025, Last Modified 10/02/2025
######################################################################

# This file contains a routine that is adopted from orignal function
# "pvalr"
# original URL:https://stackoverflow.com/questions/23018256/printing-p-values-with-0-001

pvalr <- function(pvals,
                  sig.limit = .001,
                  digits = 3) {

  roundr <- function(x, digits = 1) {
    res <- sprintf(paste0('%.', digits, 'f'), x)
    zzz <- paste0('0.', paste(rep('0', digits), collapse = ''))
    res[res == paste0('-', zzz)] <- zzz
    res
  }

  func <- function(x, sig.limit) {
    if(!is.na(x)) {
      if (x < sig.limit) return(sprintf('< %s', format(sig.limit)))

      if (x > .1) {
        return(roundr(x, digits = 2))
      } else {
        return(roundr(x, digits = digits))
      }
    } else {
      return(NA)
    } # end if

  }
  sapply(pvals, FUN = func, sig.limit = sig.limit)

  # sapply(pvals, function(x, sig.limit) {
  #   if (x < sig.limit) return(sprintf('< %s', format(sig.limit)))
  #   if (x > .1)
  #     return(roundr(x, digits = 2)) else
  #       return(roundr(x, digits = digits))
  # }, sig.limit = sig.limit)
}


# Helper function to extract a specific argument from an object of class 'tcl_sa_size'.
# obj An object of class 'tcl_sa_size'.
# arg A character string specifying the name of the element to extract (default is "df").
# res = The value of the specified element in the object.
get_tcl_arg <- function(obj, arg = "df") {
    stopifnot("obj must be of class 'tcl_sa_size'" = inherits(obj, "tcl_sa_size"))
    # Extract arg using the specified argument
    res <- obj[[arg]]
  return(res)
}

# This function below is provided courtesy of Rainer W. Alexandrowicz.
n_info = function(obj) {
  stopifnot(inherits(obj,"LR"))

  G = length(obj$fitobj)    # no. of split groups
  N = rep(NA,G)             # informative sample sizes

  for (g in 1:G) {
    d = obj$fitobj[[g]]$X
    m = sum(apply(obj$X.list[[g]],2,max)) # na.rm?
    f = table(rowSums(d))     # score frequencies
    k = dim(d)[2]             # no. items
    s = names(f)              # scores
    x = which(s %in% c(0,m))  # exclude zero and perfect scores
    N[g] = sum(f[-x])         # informative sample size
  } # end for

  return(N)
} # end fun n_info

# informative sample size when binary data matrix is given
# n_info_data01 <- function(data) {
#   score_freq = table(rowSums(data))
#   max_score = as.character(dim(data)[2])
#   score_freq[names(score_freq) %in% c("0",max_score)] = 0 # no information if zero or total response
#   n = sum(score_freq)
#   return(n)
# }

n_info_data <- function(data) {
  X <- as.matrix(data)
  total_score <- rowSums(X, na.rm = TRUE)
  max_score <- sum(apply(X, 2, function(col) max(col, na.rm = TRUE)))
  score_freq <- table(total_score)
  keep_scores <- setdiff(as.numeric(names(score_freq)), c(0, max_score))
  N <- sum(score_freq[as.character(keep_scores)])
  return(N)
}
