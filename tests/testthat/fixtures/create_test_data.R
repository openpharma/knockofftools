# ---------------------------------------------------------------------------
# (re)Generate test fixture data for knockofftools
# ---------------------------------------------------------------------------
#
# Run this script from the package root to (re)create the .rds fixtures:
#
#   Rscript tests/testthat/fixtures/create_test_data.R
#
# PRE-GENERATED FIXTURES
#
#   Pros
#   ----
#   - pre-computing knockoff.statistics() to avoid overhead in testing
#   - Fully deterministic: no seed sensitivity or platform-dependent RNG
#     differences can cause spurious failures.
#
#   Cons
#   ----
#   - Opacity: someone reading the tests cannot immediately see how the data
#     was created. This script exists to keep that transparent.
#   - Staleness: if knockofftools changes its knockoff algorithm or data
#     generation functions, the fixtures may become outdated. Re-run this
#     script after upgrading knockofftools to keep them in sync.
#   - Binary files in version control: .rds files are not human-readable,
#     making diffs unhelpful. Keep the fixtures small to minimize this cost.
#
# ---------------------------------------------------------------------------

make_fixture <- function(n, p, p_b, p_nn, seed) {
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
  ) %>%
  c(list(
    W_g = knockoff.statistics(y = .$y_g, X = .$X, type = "regression",     M = 50),
    W_b = knockoff.statistics(y = .$y_b, X = .$X, type = "classification", M = 50),
    W_s = knockoff.statistics(y = .$y_s, X = .$X, type = "survival",       M = 50)
  )
  )
}

fixtures_dir <- file.path("tests", "testthat", "fixtures")
dir.create(fixtures_dir, showWarnings = FALSE, recursive = TRUE)

# Single fixture used by all tests
saveRDS(
  make_fixture(n = 100, p = 10, p_b = 2, p_nn = 5, seed = 1),
  file.path(fixtures_dir, "test_data.rds")
)

cat("Fixtures written to", fixtures_dir, "\n")
