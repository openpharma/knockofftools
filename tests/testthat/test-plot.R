# Copyright 2023 Novartis Institutes for BioMedical Research Inc.
#
# Licensed under the MIT License (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# https://www.mit.edu/~amini/LICENSE.md
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

test_that("plot.variable.selections works as expected", {

  set.seed(1)

  # Simulate 8 Gaussian covariate predictors and 2 binary factors:
  X <- generate_X(n=100, p=10, p_b=2, cov_type="cov_equi", rho=0.2)

  # create linear predictor with first 5 beta-coefficients = 1 (all other zero)
  lp <- generate_lp(X, p_nn = 5, a=1)

  # Gaussian

  # Simulate response from a linear model y = lp + epsilon, where epsilon ~ N(0,1):
  y <- lp + rnorm(100)

  # Calculate M independent knockoff feature statistics:
  W <- knockoff.statistics(y=y, X=X, type="regression", M=10)

  S <- variable.selections(W, error.type = "pfer", level = 1)

  expect_no_warning(p <- plot(S))
  expect_s3_class(p, "ggplot")

  S_miss_class <- S
  class(S_miss_class) <- setdiff(class(S), "variable.selections")
  expect_error(plot.variable.selections(S_miss_class))


})


test_that("plot snapshot", {

  skip_if_not_installed("vdiffr")
  skip_if_not(utils::packageVersion("ggplot2") >= "3.5.0")

  # tbd: use fixture
  # fixture <- readRDS(test_path("fixtures", "test_data.rds"))
  # W <- fixture$W_g

  set.seed(1)
  # Simulate 8 Gaussian covariate predictors and 2 binary factors:
  X <- generate_X(n=100, p=10, p_b=2, cov_type="cov_equi", rho=0.2)
  # create linear predictor with first 5 beta-coefficients = 1 (all other zero)
  lp <- generate_lp(X, p_nn = 5, a=1)
  # Gaussian
  # Simulate response from a linear model y = lp + epsilon, where epsilon ~ N(0,1):
  y <- lp + rnorm(100)
  # Calculate M independent knockoff feature statistics:
  W <- knockoff.statistics(y=y, X=X, type="regression", M=10)
  S <- variable.selections(W, error.type = "pfer", level = 1)

  p <- plot(S)

  vdiffr::expect_doppelganger(
    title = "plot-snapshot",
    fig = p
  )

})
