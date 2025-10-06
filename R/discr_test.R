##############################################################################
# UMIT Tirol -  Private University for Health Sciences and Health Technology
#   Institute of Psychology, Statistics and Psychometrics Working Group
#
# mod_test
#
# Part of R/tlc - Testing in Conditional Likelihood Context package
#
# This file contains a routine that computes Wald (W), Score (RS),
# likelihood ratio (LR) and gradient (GR) test
# for hypothesis that the difference of item discrimination parameters equals zero,
# where one assumes all items presented to all persons.
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2025, Last Modified 17/09/2025
##############################################################################
#' Testing item discriminations
#'
#' Computes Wald (W), likelihood ratio (LR), Rao score (RS) and gradient (GR) test of
#'  hypothesis of equal item discriminations against the
#'   alternative that at least one item discriminates differently (only for binary data).
#'
#'The tests are based on the following model suggested in Draxler, Kurz, Guerer, and Nolte (2024)
#'
#'\deqn{ \text{logit} \big( E(Y) \big ) = \tau + \alpha + \delta (r - 1), }
#'
#' where \eqn{E(Y)} is the expected value of a binary response (of a person to an item),
#' \eqn{r = 1, \dots, k - 1} is the person score, i.e., number of correct responses of that person
#' when responding to \eqn{k} items, \eqn{\tau} is the respective person parameter and \eqn{\alpha} and
#' \eqn{\delta} are two parameters referring to the respective item. The parameter \eqn{\alpha}
#' represents a baseline, i.e., the easiness or attractiveness of the respective item in person score
#' group \eqn{r = 1}. The parameter \eqn{\delta} denotes the constant change of the attractiveness of that
#' item between successive person score groups. Thus, the model assumes a linear effect of the person
#' score \eqn{r} on the logit of the probability of a correct response.
#'
#' The four test statistics are derived from a conditional likelihood function in which the
#' \eqn{\tau} parameters are eliminated by conditioning on the observed person scores.
#' The hypothesis to be tested is formally given by setting all \eqn{\delta} parameters equal to \eqn{0}.
#' The alternative assumes that at least one \eqn{\delta} parameter is not equal to \eqn{0}.
#'
#'
#' @param data Data matrix.
#' @return A list of class \code{tcl} of test statistics, degrees of freedom, and p-values.
#'  \item{test}{A numeric vector of Wald (W), likelihood ratio (LR), Rao score (RS), and gradient (GR)statistics.}
#'  \item{df}{A numeric vector of corresponding degrees of freedom.}
#'  \item{pvalue}{A vector of corresponding p-values.}
#'  \item{data}{Data matrix.}
#'  \item{call}{The matched call.}
#'
#' @references{
#' Draxler, C., Kurz, A., Guerer, C., & Nolte, J. P. (2024). An Improved Inferential Procedure to Evaluate Item Discriminations
#' in a Conditional Maximum Likelihood Framework. \emph{Journal of Educational and Behavioral Statistics}, 49(3), 403-430.
#'
#'  }
#' @keywords htest
#' @export
#' @seealso \code{\link{invar_test}}, \code{\link{change_test}}, and \code{\link{LLTM_test}}.
#' @examples
#' \dontrun{
#'##### Dataset PISA Mathematics data.pisaMath {sirt} #####
#'
#' library(sirt)
#' data(data.pisaMath)
#' y <- data.pisaMath$data[, grep(names(data.pisaMath$data), pattern = "M" )]
#'
#' res <- discr_test(data = y)
#' # $test
#' #      W     LR     RS     GR
#' # 72.470 73.032 76.725 73.430
#' #
#' # $df
#' # W LR RS  GR
#' # 10 10 10 10
#' #
#' # $pvalue
#' #       W        LR        RS         GR
#' # "< 0.001" "< 0.001" "< 0.001" "< 0.001"
#' #
#' # $call
#' # discr_test(X = y)
#'
#'}

discr_test <- function(data) {
  # data = observed data matrix comprised of k columns

  call <- match.call()

  # general likelihood function for an arbitrary number of items k
  loglike <- function(b){
    # k <- ncol(y)      # will be computed in calling function
    # r <- rowSums(y)
    # y1 <- y[r > 0 & r < k,]
    # r1 <- rowSums(y1)

    m <- b[1:k] + b[(k+1):(2*k)] %o% (0:(k-2))
    g <- apply(X = m, MARGIN = 2, FUN = psychotools::elementary_symmetric_functions)
    i <- seq(from = 2, to = (k-1)*(k+1), by = k+1+1)
    ft <- table(factor(r1, levels = 1:(k-1)))
    c <- sum(ft * log(unlist(g)[i]))

    m1 <- array(-b[1:k], dim = c(k, nrow(y1)))
    m2 <- array(-b[(k+1):(2*k)], dim = c(k, nrow(y1)))
    m3 <- t(m2) * (r1 - 1)
    logl <- sum(rowSums(y1 * (t(m1) + m3))) - c

    return(-logl)
  }

  # score function, gradient
  score <- function(b){
    # k <- ncol(y)     # will be computed in calling function
    # r <- rowSums(y)
    ft <- table(factor(r, levels = 0:k))

    m <- b[1:k] + b[(k+1):(2*k)] %o% (0:(k-2))
    g <- apply(X = m, MARGIN = 2, FUN = psychotools::elementary_symmetric_functions, order = 1)

    f <- function(i) as.vector(ft) * (g[[i]]$'1' / g[[i]]$'0')
    a <- lapply(X = 1:(k-1), FUN = f)

    f1 <- function(j) a[[j]][j+1,]
    a1 <- sapply(X = 1:(k-1), FUN = f1)
    score1 <- colSums(y) - ft[k+1] - rowSums(a1)
    y1 <- y * (r-1)
    score2 <- colSums(y1[r > 1 & r < k,]) - rowSums(a1[,2:(k-1)] * (rep(1, k) %o% (1:(k-2))))
    score <- c(score1, score2)

    return(score)
  }

  # information function, block diagonal matrix, each block referring to one of k-1 person setClass("Class", slots = c(name = "type"))
  info <- function(b){

    # k <- ncol(y) # will be computed in calling function
    m <- b[1:k] + b[(k+1):(2*k)] %o% (0:(k-2))
    e <- apply(X = m, MARGIN = 2, FUN = psychotools::elementary_symmetric_functions, order = 2)

    f <- function(i){
      g <- e[[i]]$`0`
      g1 <- e[[i]]$`1`
      cp <- g1[2:k,] / g[2:k]
      return(cp)
    }

    f1 <- function(i){
      g <- e[[i]]$`0`
      g2 <- e[[i]]$`2`
      cp2 <- g2[2:k,,] / g[2:k]
      return(cp2)
    }

    cp <- lapply(X = 1:(k-1), FUN = f)
    cp2 <- lapply(X = 1:(k-1), FUN = f1)

    r <- rowSums(y)
    ft <- table(factor(r, levels = 1:(k-1)))

    f2 <- function(j) (ft[j] * (cp2[[j]][j,,] - cp[[j]][j,] %o% cp[[j]][j,]))[1:(k-1),1:(k-1)]

    l <- lapply(X = 1:(k-1), FUN = f2)
    I <- Matrix::bdiag(l)

    return(I)
  }

  # auxiliary functions
  f <- function(i,k){
    h1 <- array(rep(0, k-1), dim = c((k-1), (k-1)))
    diag(h1) <- rep(1, k-1)
    h2 <- h1 * array(i, dim = c((k-1), (k-1)))
    return(h2)
  }

  f1 <- function(k){
    h1 <- array(rep(0, k-1), dim = c((k-1), (k-1)))
    diag(h1) <- rep(1, k-1)
    return(h1)
  }

  f2 <- function(j,l1,l2) cbind(l1[,,j], l2[[j]])


  y <- data

  ###################################################
  #               main programm
  ###################################################
  k <- ncol(y)
  r <- rowSums(y)  # adapted, compute here for loglik(), score() and info()
  y1 <- y[r > 0 & r < k,]
  r1 <- rowSums(y1)

  df <- ncol(y) - 1

  # test statistics Wald, LR, RS, GT and LRmml

  # unrestricted and restricted maximization of loglike

  u <- optim(par = rep(0, 2*k), fn = loglike, gr = score, method = 'BFGS')
  re <- eRm::RM(y)
  # re <- psychotools::raschmodel(y)

  # score test statistic

  respar <- c(-re$betapar, rep(0, k)) # eRM

  s <- score(b = respar)[c(1:(k-1),(k+1):(2*k-1))]  # vector-valued score function for free parameters

  # r <- rowSums(y)
  ft <- table(factor(r, levels = 1:(k-1)))


  l1 <- replicate(k-1, f1(k=k)) # adapted AK
  l2 <- lapply(X = 0:(k-2), f,k=k) # adapted AK
  l3 <- lapply(X = 1:(k-1), f2, l1=l1, l2=l2) # adapted AK
  l4 <- simplify2array(x = l3, higher = TRUE)

  W <- apply(X = l4, MARGIN = 2, rbind) # Jacobian matrix, linear transformation from (k-1)*(k-1)-dimensional parameter space to 2*(k-1)-dim. space, k > 2

  # orginal (k-1)*(k-1)-dim. space refers to model assuming person score group specific item parameters, i.e.,
  # k-1 free item parameters for each of k-1 person scores. These are linearly restricted by allowing each item parameter to differ only by a constant
  # between successive person score groups

  # W is a matrix with (k-1)*(k-1) rows, each row referring to one of the free person score group specific item parameters,
  # and 2*(k-1) columns, each column referring to one of the free parameters of the linearly restriced model
  # the first k-1 item parameters per score group are treated as free parameters

  Ir <- t(W) %*% (info(b = respar)) %*% W
  # information matrix evaluated at restricted estimates

  # RS <- t(s) %*% solve(Ir) %*% s
  RS <- as.numeric(t(s) %*% solve(Ir) %*% s)


  # Wald test statistic

  uru <- u$par
  ur <- uru[(k+1):(2*k-1)] - uru[2*k]

  Iu <- t(W) %*% (info(b = uru)) %*% W
  # information matrix evaluated at unrestricted estimates

  Cov <- solve(Iu)[k:(2*k-2),k:(2*k-2)]

  # Wald <- t(ur) %*% solve(Cov) %*% ur
  Wald <- as.numeric( t(ur) %*% solve(Cov) %*% ur )


  # Likelihood ratio test statistic

  LR <- -2 * (re$loglik - -u$value)


  # Gradient test

  sg <- score(b = respar)
  GT <- sg %*% -uru

  test.stats <- c( GT, LR, RS, Wald)
  names(test.stats) <- c("GR", "LR", "RS", "W")

  pvalue <- 1 - (sapply(test.stats, stats::pchisq, df = df))
  pvalue <- pvalr(pvalue, digits = 3)

  df_vec <- c(df,df,df,df)
  names(df_vec) <- c("GR", "LR", "RS", "W")

  res.list <- list("test" = round(test.stats, digits = 3),
                   "df" = df_vec,
                   "pvalue" = pvalue,
                   "data" = y)

  res.list$call <- call

  # Define test order
  test_order <- c("W", "LR", "RS", "GR")

  # Reorder selected elements
  for (nm in c("test", "df", "pvalue")) {
    res.list[[nm]] <- res.list[[nm]][test_order]
  }

  res.list <- structure(
    res.list,
    class = "tcl"
  )

  return(res.list)
}
