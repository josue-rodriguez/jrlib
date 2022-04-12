#' gen_melsm_data
#'
#' @importMethodsFrom Matrix t
#' @importFrom Matrix sparseMatrix KhatriRao
#' @export


gen_melsm_data <- function(n, n_i) {

  x1 = rnorm(n*n_i)
  x2 = rep(sample(0:1, n, replace = TRUE))

  Xloc <- cbind(int = 1, x1, x2 = rep(x2, each = n_i))
  Xscale <- cbind(int = 1, x1, x2 = rep(x2, each = n_i))
  Xtau <- cbind(int = 1, x2)


  Omega <- rbind(
    c(1, 0.1),
    c(0.1, 1)
  )


  iota_l <- exp(Xtau %*% rbind(-0.1, 0.2))
  iota_s <- exp(Xtau %*% rbind(-0.2, 0.1))

  tau <- cbind(iota_l, iota_s)



  Sigma_i <- array(data = NA, dim = c(2, 2, n))

  Theta <- matrix(data = NA, ncol = 2, nrow = n)

  for (j in 1:n) {
    Sigma_i[,,j] <- diag(tau[j, ]) %*% Omega %*% diag(tau[j, ])

    Theta[j, ] <- MASS::mvrnorm(n = 1, mu = c(0, 0), Sigma = Sigma_i[,,j])

  }


  Bsigma <- rbind(-0.2, 0.05, -0.2)
  sigma_ij <- exp(Xscale %*% Bsigma + rep(Theta[, 2], each = n_i))

  Bloc <- rbind(10, 2, 3)
  mu <- Xloc %*% Bloc +  rep(Theta[, 1], each = n_i)


  y <- rnorm(n * n_i, mu, sigma_ij)

  df <- data.frame(
    x1,
    x2 = factor(rep(x2, each = n_i)),
    y,
    id = rep(1:n, each = n_i),
    sigma = sigma_ij,
    mu = mu
  )

  return(df)
}


df <- gen_melsm_data(200, 100)

hist(df$y)

library(brms)
library(cmdstanr)

brm_fit <- brm(
  bf(
    y ~ x1 + x2 + (1|c|gr(id, by = x2)),
    sigma ~ x1 + x2 + (1|c|gr(id, by = x2))
  ),
  data = df,
  chains = 4,
  iter = 5000,
  file = "inst/models/brm_fit1.rds",
  backend = "cmdstanr",
  cores = parallel::detectCores()
)


brm_fit2 <- update(brm_fit, iter = 5000, chains = 4, backend = "cmdstanr",
                   cores = parallel::detectCores())

summary(brm_fit2)
plot(brm_fit2)

# exp(rbind(c(1,0), c(1,1)) %*% rbind(-0.1, 0.2))
# exp(rbind(c(1,0), c(1,1)) %*% rbind(-0.2, 0.1))
