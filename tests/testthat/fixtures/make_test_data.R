# ---------------------------------------------------------------------------
# (re)Generate test fixture data for knockofftoolsExtra
# ---------------------------------------------------------------------------
#
# Run this script from the package root to (re)create the .rds fixtures:
#
#   Rscript tests/testthat/fixtures/make_test_data.R
#
# PRE-GENERATED FIXTURES
#
#   Pros
#   ----
#   - Tests run much faster: knockoff generation (knockofftools::knockoff) and
#     covariate simulation (knockofftools::generate_X) are the dominant cost in
#     each test; pre-generating them eliminates that overhead entirely.
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
#     making diffs unhelpful. Keep the fixtures small to minimise this cost.
#
# ---------------------------------------------------------------------------

library(knockofftools)

make_fixture <- function(n, p, p_b, p_nn, seed) {
  set.seed(seed)
  X <- knockofftools::generate_X(
    n = n, p = p, p_b = p_b,
    cov_type = "cov_equi", rho = 0.2
  )
  lp  <- knockofftools::generate_lp(X, p_nn = p_nn, a = 1)
  X_k <- knockofftools::knockoff(X)
  list(
    X   = X,
    X_k = X_k,
    lp  = lp,
    y_g = lp + rnorm(n),
    y_b = factor(rbinom(n, size = 1, prob = exp(lp) / (1 + exp(lp)))),
    y_s = knockofftools::simulWeib(N = n, lambda0 = 0.01, rho = 1, lp = lp)
  )
}

fixtures_dir <- file.path("tests", "testthat", "fixtures")
dir.create(fixtures_dir, showWarnings = FALSE, recursive = TRUE)

# Single fixture used by all tests
saveRDS(
  make_fixture(n = 100, p = 10, p_b = 2, p_nn = 5, seed = 1),
  file.path(fixtures_dir, "test_data_default.rds")
)

cat("Fixtures written to", fixtures_dir, "\n")
