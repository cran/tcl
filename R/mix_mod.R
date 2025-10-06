##############################################################################
# UMIT Tirol -  Private University for Health Sciences and Health Technology
#   Institute of Psychology, Statistics and Psychometrics Working Group
#
# mix_mod
#
# Part of R/tlc - Testing in Conditional Likelihood Context package
#
# This file contains a routine that computes
# mixed logit model code allowing for missing responses except for covariates
# code depends on package psychotools, eRm, and numDeriv
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2025, Last Modified 16/09/2025
######################################################################
#' Mixed model considering the effects of multiple covariates  .
#'
#' Estimates and tests linear effects of multiple covariates on item parameters of the Rasch model (RM) simultaneously.
#'
#'The underlying model is a mixed-effects logit model with random person effects and ﬁxed item and covariates effects, i.e.,
#'
#' \deqn{
#' \log \frac{P(Y_{ij} = 1)}{1 - P(Y_{ij} = 1)}
#'   = \tau_i + \alpha_j + \sum_{p = 1}^q x_{ip} \delta_{jp}, \quad
#'   i = 1, \dots, n, \; j = 1, \dots, k, \; p = 1, \dots, q,
#' }
#' where \eqn{Y_{ij} \in \{0, 1\}}, \eqn{\tau_i} is a person parameter, \eqn{\alpha_j} is a baseline effect of item \eqn{j},
#' \eqn{\delta_{jp}} is an effect of covariate \eqn{p} on item \eqn{j}, and \eqn{x_{ip}} is a covariate
#' value observed for person \eqn{i} and covariate \eqn{p}. For identifiability, \eqn{\alpha_1 = 0}, \eqn{\delta_{1p} = 0 \;\forall p}.
#' Setting all \eqn{\delta} parameters (\eqn{\forall j, p}) to 0 yields the RM as a special case (with the \eqn{\alpha}s as the
#' item parameters of the RM).
#'
#' The \eqn{\alpha} and \eqn{\delta} parameters are estimated using a conditional maximum likelihood (CML)
#' approach and four different tests based on the conditional likelihood and derived from asymptotic theory
#' are provided, i.e., likelihood ratio (LR), Rao score (RS), Wald (W), and gradient (GR) test.
#' The hypothesis of interest is \eqn{\delta_{jp} = 0 \;\forall j, p} against the alternative that at least one \eqn{\delta}
#' parameter is not equal to \eqn{0}.
#' Furthermore, \eqn{Z} test statistics (i.e., standard normal distribution when
#' the true effect of a covariate on an item is \eqn{0}) are computed for each item and covariate separately.
#'
#' @param data Data matrix consisting of binary responses, i.e., 0s and 1s. Missing responses are NAs.
#'
#' @param Xcov Covariate matrix. Persons in rows and covariates in columns, e.g., age, gender, drug dosage, etc. In case of
#'one covariate Xcov must be a one-column matrix.
#'
#' @return A list of class \code{tcl} with the following components:
#'   \item{CML_estimates}{Conditional maximum likelihood (CML) estimates of item (\eqn{\alpha} and \eqn{\delta})
#'  parameters (easiness, attractiveness). The effects of the first item are set to \eqn{0}
#'  for identifiability.}
#'   \item{SE}{Standard errors of CML estimates.}
#'   \item{Z_statistics}{Z test statistics for each single parameter (\eqn{\alpha} and \eqn{\delta}), i.e.,
#'   testing the hypothesis that the true value of the respective parameter is \eqn{0}
#'   against the alternative of \eqn{\neq 0}.}
#'   \item{pvalue}{A matrix of two-sided p-values for the \eqn{Z} tests.}
#'  \item{loglik}{Conditional log-likelihood.}
#'   \item{tests}{A table summarizing the results of four tests (W, LR, RS, GR) of the hypothesis that the effects of all covariates on all items are all 0.
#'   The table contains the test statistic (\code{stat}), degrees of freedom (\code{df}), and two-sided p-value (\code{pvalue}) for each test.}
#'   \item{information_criteria}{AIC, BIC}
#'   \item{call}{The matched call.}
#'
#' @references{
#' Draxler, C., & Kurz, A. (2025). Testing measurement invariance in a conditional likelihood framework by considering
#' multiple covariates simultaneously. Behavior Research Methods, 57(1), 50.
#' }
#' @keywords htest
#' @export
#' @seealso \code{\link{invar_test}}, \code{\link{change_test}}, and \code{\link{LLTM_test}}.
#'
#' @examples
#' \dontrun{
#'##### Rasch Model #####
#' dat <- eRm::raschdat3
#' x1 <- c(rep(0,250), rep(1,250))
#' x2 <- runif(500,min = 0, max = 1)
#' X <- cbind(x1,x2)
#'
#' res <- mix_mod(data = dat, Xcov = X)
#' # $CML_estimates
#' #      1      2      3      4      5      6
#' # base 0 -0.596 -1.152 -1.804 -1.846 -2.353
#' # 1    0 -0.380 -0.403  0.072 -0.121 -0.452
#' # 2    0  0.814  0.780  0.612 -0.277  0.069
#' #
#' # $SE
#' #       1     2     3     4     5     6
#' # base NA 0.356 0.347 0.349 0.353 0.369
#' # 1    NA 0.314 0.303 0.301 0.307 0.320
#' # 2    NA 0.564 0.545 0.541 0.548 0.572
#' #
#' # $Z_statistics
#' #       1      2      3      4      5      6
#' # base NA -1.675 -3.320 -5.175 -5.223 -6.377
#' # 1    NA -1.210 -1.330  0.240 -0.396 -1.413
#' # 2    NA  1.443  1.431  1.132 -0.505  0.121
#' #
#' # $pvalue
#' #       1     2     3     4     5     6
#' # base NA 0.094 0.001 0.000 0.000 0.000
#' # 1    NA 0.226 0.183 0.810 0.692 0.158
#' # 2    NA 0.149 0.153 0.258 0.613 0.904
#' #
#' # $loglik
#' # [1] -993.8575
#' #
#' # $tests
#' #      stat df pvalue
#' # W  14.339 10  0.158
#' # LR 14.462 10  0.153
#' # RS 14.507 10  0.151
#' # GR 14.499 10  0.151
#' #
#' # $information_criteria
#' #           AIC      BIC
#' # [1,] 2007.715 2049.432
#' #
#' # $call
#' # mix_mod(data = dat, Xcov = X)
#' #
#' # attr(,"class")
#' # [1] "tcl"
#'}


mix_mod <- function(data, Xcov) {
  call <- match.call()

  y <- as.matrix(data) # must be a matrix
  x <- as.matrix(Xcov) # must be a matrix

  loglike <- function(b) {
    h <- function(i)
      b[1:k] + colSums(t(array(b[(k + 1):(k * (p + 1))], dim = c(k, p))) * x[i, ])
    l <- lapply(1:n, h)
    m <- t(simplify2array(l, higher = TRUE))

    z <- y
    z[y == 0] <- 1
    g <- apply(
      X = -m * z,
      MARGIN = 1,
      FUN = psychotools::elementary_symmetric_functions
    )
    r <- rowSums(y, na.rm = T)

    h1 <- function(j)
      g[[j]]$`0`[r[j] + 1]
    l1 <- lapply(1:n, h1)
    c <- unlist(l1)

    logl <- sum(rowSums(y * m, na.rm = T) - log(c))

    return(-logl)
  }

  score <- function(b) {
    h <- function(i)
      b[1:k] + colSums(t(array(b[(k + 1):(k * (p + 1))], dim = c(k, p))) * x[i, ])
    l <- lapply(1:n, h)
    m <- t(simplify2array(l, higher = TRUE))

    z <- y
    z[y == 0] <- 1
    g <- apply(
      X = -m * z,
      MARGIN = 1,
      FUN = psychotools::elementary_symmetric_functions,
      order = 1
    )
    r <- rowSums(y, na.rm = T)

    h1 <- function(j)
      g[[j]]$`0`[r[j] + 1]
    l1 <- lapply(1:n, h1)
    c <- unlist(l1)

    h2 <- function(j)
      g[[j]]$`1`[r[j] + 1, ]
    l2 <- lapply(1:n, h2)
    cp <- t(simplify2array(l2, higher = TRUE)) / c

    h3 <- function(j)
      colSums(cp * x[, j], na.rm = T)
    l3 <- lapply(1:ncol(x), h3)
    e <- c(colSums(cp, na.rm = T), unlist(l3))

    h4 <- function(j)
      colSums(y * x[, j], na.rm = T)
    l4 <- lapply(1:ncol(x), h4)
    o <- c(colSums(y, na.rm = T), unlist(l4))

    s <- o - e

    return(-s)
  }



  #info <- function(b){

  #h <- function(i) b[1:k] + colSums(t(array(b[(k+1):(k*(p+1))], dim = c(k, p))) * x[i,])
  #l <- lapply(1:n, h)
  #m <- t(simplify2array(l, higher = TRUE))

  #g <- apply(X = -m, MARGIN = 1, FUN = psychotools::elementary_symmetric_functions, order = 2)
  #r <- rowSums(y)



  #h1 <- function(i){
  #cp1 <- g[[i]]$'1' / g[[i]]$`0`
  #cp2 <- g[[i]]$'2' / g[[i]]$`0`

  #h2 <- function(j) cp2[,,j][r[i]+1,]
  #c2 <- simplify2array(lapply(1:k, h2))

  #Ip <- c2 - cp1[r[i]+1,] %o% cp1[r[i]+1,]

  #return(Ip)

  #}


  #list <- lapply(1:n, h1)
  #I <- Matrix::bdiag(list)

  #return(I)

  #}


  # Jacobian for each person

  #Jp <- function(i){

  #h <- function(){
  #a <- array(rep(0, k), dim = c(k, k))
  #diag(a) <- rep(1, k)
  #return(a)
  #}

  #a <- replicate(p, h())
  #a1 <- t(apply(X = a, FUN = cbind, MARGIN = 2))

  #a2 <- as.vector(t(array(x[i,], dim = c(p, k))))

  #a3 <- array(rep(0, k), dim = c(k, k))
  #diag(a3) <- rep(1, k)

  #a4 <- cbind(a3, t(t(a1) * a2))

  #return(a4)

  #}



  n <- nrow(x)
  k <- ncol(y)
  p <- ncol(x)


  #a3 <- simplify2array(lapply(1:n, Jp), higher = TRUE)
  #J <- apply(X = a3, FUN = cbind, MARGIN = 2)
  #J <- J[,-seq(from = 1, to = (p+1)*k, by = k)]  # total Jacobian matrix (for all persons)

  res <- optim(par = rep(0, (p + 1) * k),
               fn = loglike,
               gr = score,
               method = 'BFGS',
               hessian = T)  # unrestricted optimization of loglike
  b <- res$par  # unrestricted CML estimates (not normalized)
  mb <- t(array(b, dim = c(k, p + 1)))
  mb <- mb - mb[, 1]
  rownames(mb) <- c('base', 1:p)
  colnames(mb) <- 1:k
  # normalized estimates, parameters referring to first time point (or item) set to 0
  # first row refers to baseline parameters, second row to slope parameters of
  # first covariate, third row to slope parameters of second covariate, etc.

  urlogl <- res$value * -1  # unrestricted maximum of loglike
  #Iur <- t(J) %*% info(b = b) %*% J  # information matrix evaluated at unrestricted estimates
  Iur <- res$hessian[-seq(from = 1,
                          to = (p + 1) * k,
                          by = k), -seq(from = 1,
                                        to = (p + 1) * k,
                                        by = k)]  # information matrix evaluated at unrestricted estimates
  C <- solve(Iur)  # asymptotic covariance matrix of estimtes
  SE <- cbind(rep(NA, p + 1), t(array(sqrt(diag(C)), dim = c(k - 1, p + 1))))
  rownames(SE) <- c('base', 1:p)
  colnames(SE) <- 1:k
  # standard errors of all estimates

  z <- mb / SE  # signed Wald (z) tests of hypothesis that each parameter equals 0
  pv <- 1 - pchisq(z^2, 1)  # corresponding two-sided p values
  rownames(z) <- c('base', 1:p)
  colnames(z) <- 1:k
  rownames(pv) <- c('base', 1:p)
  colnames(pv) <- 1:k

  resr <- eRm::RM(y, sum0 = FALSE)  # restricted optimization of loglike
  br <- c(resr$betapar, rep(0, p * k))  # restricted CML estimates
  rlogl <- resr$loglik  # restricted maximum of loglike


  # likelihood ratio test

  df <- p * (k - 1)

  LR <- -2 * (rlogl - urlogl)
  plr <- 1 - pchisq(LR, df)


  # Gradient test

  G <- t(-score(b = br)) %*% b
  pg <- 1 - pchisq(G, df)


  # Rao score test

  #Ir <- t(J) %*% info(b = br) %*% J  # information matrix evaluated at restricted estimates
  Ir <- numDeriv::hessian(loglike, br)  # information matrix evaluated at restricted estimates
  RS <- t(score(b = br)[-seq(from = 1,
                             to = (p + 1) * k,
                             by = k)]) %*% solve(Ir[-seq(from = 1,
                                                         to = (p + 1) * k,
                                                         by = k), -seq(from = 1,
                                                                       to = (p + 1) * k,
                                                                       by = k)]) %*% score(b = br)[-seq(from = 1,
                                                                                                        to = (p + 1) * k,
                                                                                                        by = k)]
  prs <- 1 - pchisq(RS, df)

  # Wald test

  Cov <- C[-(1:(k - 1)), -(1:(k - 1))]  # asympt. covariance matrix of free slope parameters
  bf <- t(array(b[-(1:k)], dim = c(k, p)))
  W <- t(as.vector(t((bf - bf[, 1])[, -1]))) %*% solve(Cov) %*% as.vector(t((bf - bf[, 1])[, -1]))
  pw <- 1 - pchisq(W, df)

  t <- rbind(
    c(round(W, 3), df, round(pw, 3)),
    c(round(LR, 3), df, round(plr, 3)),
    c(round(RS, 3), df, round(prs, 3)),
    c(round(G, 3), df, round(pg, 3))
  )
  # colnames(t) <- c('test_statistic', 'df', 'p value')
  colnames(t) <- c('stat', 'df', 'pvalue')
  rownames(t) <- c('W',
                   'LR',
                   'RS',
                   'GR')
  # rownames(t) <- c('likelihood ratio test',
  #                  'gradient test',
  #                  'Rao score test',
  #                  'Wald test')

  # test.stats <- c( G, LR, RS, W)
  # names(test.stats) <- c("GR", "LR", "RS", "W")
  #
  # df_vec <- c(df,df,df,df)
  # names(df_vec) <- c("GR", "LR", "RS", "W")

  # pvalue <- c(pg, plr, prs, pw )
  # names( pvalue) <- c("GR", "LR", "RS", "W")
  # pvalue <- pvalr(pvalue, digits = 3)

  AIC <- -2 * urlogl + 2 * p * (k - 1)
  row <- rowSums(y, na.rm = T)
  ni <- length(row[(rowSums(y, na.rm = T) > 0 &
                      rowSums(y, na.rm = T) < k)])
  BIC <- -2 * urlogl + log(ni) * p * (k - 1)
  t1 <- array(c(AIC, BIC), dim = c(1, 2))
  colnames(t1) <- c('AIC', 'BIC')

  res.list <- list(
    'CML_estimates' = round(mb, 3),
    'SE' = round(SE, 3), # 'standard_errors' = round(SE, 3),
    'Z_statistics' = round(z, 3),
    'pvalue' = round(pv, 3), # 'two_sided_pvalues' = round(pv, 3),
    "loglik" = urlogl, # 'log likelihood'
    'tests' = t,
    'information_criteria' = t1, # 'information_criteria'
    "call" = call
  )

  res.list <- structure(
    res.list,
    class = "tcl"
  )

  return(res.list)
}
