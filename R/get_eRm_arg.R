#' @title Extract Arguments from an eRm Object
#'
#' @description This function extracts specific arguments from an object of class `'LR'` from the `eRm` package.
#' Depending on the selected argument, it retrieves degrees of freedom (`df`), local deviations (`local_dev`), or
#' informative sample size (`n_info`).
#'
#' @param obj An object of class `'LR'`, typically created using functions from the `eRm` package.
#' @param arg A character string specifying the argument to extract. Options are:
#'   \itemize{
#'     \item `"df"` (default): Extracts the degrees of freedom.
#'     \item `"local_dev"`: Extracts item parameters for the two person groups from the model. If more than two split groups are available, only the first two are selected.
#'     \item `"n_info"`: Computes and returns the informative sample size using `n_info()`.
#'   }
#'
#' @return The extracted argument value:
#'   \itemize{
#'     \item A numeric value if `arg = "df"`.
#'     \item A list containing local deviation parameters if `arg = "local_dev"`.
#'     \item A computed sample size if `arg = "n_info"`.
#'   }
#'
#' @details If multiple argument values are provided, `"df"` is selected by default. If an invalid `arg` is provided,
#' the function throws an error.
#'
#' @note If `obj` contains more than two split groups, only the first two will be selected for `"local_dev"`, with a message notifying the user.
#'
#' @examples
#' \dontrun{
#'   # Example usage with an LR object
#'   dat = eRm::sim.rasch(1000,10)
#'   mod = eRm::RM(dat)
#'
#'   obj <- eRm::LRtest(mod) # Create an LR object
#'   get_eRm_arg(obj, "df")      # Extract degrees of freedom
#'   get_eRm_arg(obj, "local_dev")  # Extract local deviations
#'   get_eRm_arg(obj, "n_info")  # Extract informative sample size
#' }
#'
#' @export
get_eRm_arg <- function(obj, arg = c("df","local_dev, n_info")) {
stopifnot("obj must be eRm of class 'LR'" = inherits(obj, "LR"))

if (length(arg) > 1) arg = "df"
# Extract argument based on specified key
if (arg == "df") {
  res <- obj[[arg]]
  message("Degrees of freedom (df) extracted from obj: ", res)
} else if (arg == "local_dev") {
  G = length(obj$fitobj)    # no. of split groups
  if (G > 2) cat("Only the first 2 of", G, "groups selected\n")
  res <- list(c(0, obj$fitobj[[1]]$etapar), c(0, obj$fitobj[[2]]$etapar))
  names(res) <- names(obj$fitobj)
} else if (arg == "n_info") {
  res <- n_info(obj)
  names(res) <- names(obj$fitobj)
  message("Informative sample size extracted from obj: ",
          paste(names(res), res, sep = " = ", collapse = ", "))
} else {
  stop("Invalid argument: ", arg)
}

return(res)
} # end func
