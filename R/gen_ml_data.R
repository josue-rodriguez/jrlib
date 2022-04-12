#' gen_multilevel_data
#'
#' @param n A numeric specifying the number of units.
#' @param n_j A numeric specifying the number of observations per unit.
#' @param B fixed effects (intercept + slope)
#' @param tau Random effect SDs
#' @param r correlation for random effects
#' @param sigma residual SD
#' @param type categorical or continuous
#'
#' @importMethodsFrom Matrix t
#' @importFrom Matrix sparseMatrix KhatriRao
#' @export
gen_multilevel_data <- function(n = 5, n_j = 5,
                                B = c(1, 1),
                                tau = c(1, 1),
                                sigma = 1,
                                r = 0.3,
                                type = "categorical") {
  N <- n * n_j

  # design matrix + betas
  if (type == "categorical") {
    X <- cbind(1, sample(0:1, size = N, replace = TRUE))
  } else if (type == "continuous") {
    X <- cbind(1, rnorm(N))
  } else {
    warning("Invalid type selected. Using categorical predictor.")
    X <- cbind(1, sample(0:1, size = N, replace = TRUE))
  }


  # cor matrix for random effects
  R <- rbind(
    c(1, r),
    c(r, 1)
  )

  # vcov matrix for random effects
  tau <- diag(tau)
  Tau <- tau %*% R %*% tau

  # generate random effects and stack them for each regression
  u_mat <- mvtnorm::rmvnorm(n, rep(0, 2), Tau)

  # -- random effects
  u_list <- lapply(1:n, function(i) cbind(u_mat[i, 1:2]))
  u <- do.call(rbind, u_list)

  # generate factor levels + Z matrix
  # -- taken from lme4 documentation
  g <- gl(n, n_j)
  J <- as(g, Class = "sparseMatrix")
  t2 <- selectMethod("t", signature = "dgCMatrix")
  Ji <- t2(J)
  Z <- t(KhatriRao(t(Ji), t(X)))


  # vcov matrix for residuals
  eps <- rnorm(N, 0, sigma)

  # generate data
  y <- X %*% B + Z %*% u + eps
  y <- as.numeric(y)

  df <- data.frame(y = y, x = X[, 2], id = g)
  return(df)
}


