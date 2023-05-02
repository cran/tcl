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
