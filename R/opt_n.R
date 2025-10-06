##############################################################################
# UMIT Tirol -  Private University for Health Sciences and Health Technology
#   Institute of Psychology, Statistics and Psychometrics Working Group
#
# opt_n
#
# Part of R/tlc - Testing in Conditional Likelihood Context package
#
# This file contains a routine that
# calculates the informative sample size required to detect an effect of
# interest at given type I (alpha) and type II (beta) error levels using
# Wald (W), likelihood ratio (LR), Rao score (RS), or gradient (GR) test,
# either from a previous \code{invar_test()} result or from a specified effect size
# and degrees of freedom.
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2025, Last Modified 16/09/2025
######################################################################
#' Computes the optimal sample size for item parameter invariance tests.
#'
#'Computes the informative sample size given an effect of interest and type I and II
#' error probabilities (alpha and beta) for Wald (W), likelihood ratio (LR),
#' Rao score (RS), and gradient (GR) test.
#' The routine supports two modes:
#'  Either provide the return object of a previous call to \code{invar_test()}
#'  or provide the effect size of interest along with the degrees of freedom.
#'
#' @param invar_obj Return object of a previous call to \code{invar_test()}. Default is \code{NULL}.
#' If missing, values for \code{effect} and \code{df} need to be set manually.
#' @param effect       Numeric value representing the effect size. A real number
#' between 0 and 1, interpreted as a proportion of pseudo-variance between persons with
#' different covariate values (but the same person parameter). Default is \code{NULL}.
#' @param df        Degrees of freedom of the test. Default is \code{NULL}.
#' @param alpha     Type I error probability. Default is 0.05.
#' @param beta      Type II error probability. Default is 0.05.
#' @param n_range   A numeric vector specifying the sample sizes to be evaluated. Default is \code{10:10000}).
#'
#' @returns
#' A list of two elements:
#' \item{opt_n}{The required sample sizes for the four tests.}
#' \item{real_pow}{The realized power, as the sample sizes are rounded to the next integer.}
#' \item{call}{The matched call.}
#'
#' @details
#' The informative sample size is the number of observations realizing a score
#' greater than zero and less than the maximum possible score, as these two
#' values are not informative for the tests.
#'
#' Providing the return object of a previous call to \code{invar_test()} allows
#' using the results of a pilot study to obtain an empirical estimate of
#' parameter differences between the groups.
#'
#' The default search range of \code{10:10000} should suffice for most applications.
#' However, if the maximum is reached, a warning is given.
#'
#' If \code{effect} and \code{df} are provided, the sample sizes of all four tests will be
#' equal due to their asymptotic equivalence. If an \code{invar_obj} is provided,
#' the sample sizes will usually differ slightly.
#'
#' Note: The \code{invar_test()} function currently only supports a two-group split.
#'
#' @references
#' Draxler, C., & Kurz, A. (2025). Testing measurement invariance in a conditional likelihood framework by considering
#' multiple covariates simultaneously. \emph{Behavior Research Methods}, 57(1), 50.
#'
#' @seealso \code{\link{invar_test}}, \code{\link{p_curve}}, \code{\link{p_ncurve}}
#'
#' @examples
#' \dontrun{
#' # --- a priori mode:
#'
#'   opt_n(effect=0.3,df=20)       # n=102
#'
#'   opt_n(effect=0.001,df=300)    # Warning!
#'   opt_n(effect=0.001,df=300,n_range=1000:100000) # Warning disappears, n=91087
#'
#' # --- pilot sample mode:
#'
#'   library(eRm)
#'   opt_n(invar_test(raschdat1))
#'
#' # --- typical problem: items eliminated
#'
#'   ex2 = invar_test(pcmdat,model="PCM")
#'
#' # The following items were excluded for the computation of GR,LR, and W
#' # due to inappropriate response patterns within subgroups:
#' I2 I4 I1 I5
#'
#' > opt_n(ex2)
#'
#' # Parameters:
#' # alpha = 0.05
#' # beta = 0.05
#' # power = 0.95
#' # df = 7 7 19 7                 # note the different df!
#'
#' # Observed effects
#' #     GR     LR     RS      W
#' # 0.1295 0.1284 0.8462 0.1226   # note the effect differences!
#'
#' # Optimal Sample Size
#' #  GR  LR  RS   W
#' # 169 170  36 178               # note the different sample sizes!
#' #
#' # Realized Power
#' #    GR    LR    RS     W
#' # 0.950 0.950 0.952 0.950
#'}
#' @keywords sample_size_planning
#' @export
# ------------------------------------------------------------------------------
# This function below is provided courtesy of Rainer W. Alexandrowicz.
opt_n = function(invar_obj = NULL,
                 effect = NULL,
                 df = NULL,
                 alpha = 0.05,
                 beta = 0.05,
                 n_range = 10:10000) {

        call <- match.call()
        options("warn" = 1)

        if (all(is.null(invar_obj),
                is.null(effect),
                is.null(df))) stop("Pls. provide either invar_obj or effect+df.")

        if (is.null(effect)) {
            stopifnot(inherits(invar_obj,"tcl"))
            tst = names(invar_obj$test)
            effect = invar_obj$effect
            df  = invar_obj$df
            lab = "Observed effects"
        } else {
            # tst = c("GR","LR","RS","W")
            tst = c("W","LR","RS","LR")
            effect = rep(effect,4)
            df  = rep(df, 4)
            lab = "Planned effect"
        } # end if/else

        names(df)  = tst
        names(effect) = tst
        P = rep(NA,4)
        N = rep(NA,4)
        names(P) = tst
        names(N) = tst
        for (i in 1:4) {
             q = qchisq(1 - alpha, df[i])
             pow = 1 - pchisq(q, df[i], effect[i]*n_range)
             dif = abs(pow - (1 - beta))
             pos = which(dif == min(dif))
             P[i] = pow[pos]
             N[i] = n_range[pos]
        } # end for

        if (any(N == max(n_range))) warning("Maximum search range reached. Results might be inaccurate. Try extending the n_range parameter.",call. = FALSE)

        cat("\nParameters:\n")
        cat("alpha =",alpha,"\n")
        cat("beta =",beta,"\n")
        cat("power =",1 - beta,"\n")
        cat("df =",df,"\n\n")
        cat(lab,"\n")
        print(round(effect, 4))
        cat("\n\nOptimal Sample Size\n")
        print(N)
        cat("\n\nRealized Power\n")
        print(round(P, 3))

        return(invisible(list(opt_n = N, real_pow = P, call = call)))
} # end fun opt_n
