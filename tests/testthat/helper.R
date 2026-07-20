
#' One line call to generate `knockoff.statistics()`(`stat_*()`) input for tests
#'
#' @param n number of observations
#' @param p number of covariates
#' @param p_b number of binary covariates
#' @param p_nn number of non-null covariates
#' @param seed
#'
#' @returns list of all components needed to generate knockoff statistics
#' for regression, classification, and survival outcomes
#'
make_ko_stat_input <- function(n, p, p_b, p_nn, seed) {
  #TODO usethis::use_package("withr")
  set.seed(seed)
  X <- generate_X(
    n = n, p = p, p_b = p_b,
    cov_type = "cov_equi", rho = 0.2
  )
  lp  <- generate_lp(X, p_nn = p_nn, a = 1)
  X_k <- knockoff(X)

  list(
    X   = X,
    X_k = X_k,
    lp  = lp,
    y_g = lp + rnorm(n),
    y_b = factor(rbinom(n, size = 1, prob = exp(lp) / (1 + exp(lp)))),
    y_s = simulWeib(N = n, lambda0 = 0.01, rho = 1, lp = lp),
    params = list(n = n, p = p, p_b = p_b, p_nn = p_nn, seed = seed)
  )
}


.valid_score <- function(test, X) {
  expect_true(is.numeric(test))
  expect_length(test, ncol(X))
  expect_named(test, colnames(X))
  expect_true(all(is.finite(test)))
}
