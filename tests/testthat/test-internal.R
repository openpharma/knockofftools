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

test_that("check_if_continuous works", {

  set.seed(1)
  X <- generate_X(n=1000, p=10, p_b=2, cov_type="cov_equi", rho=0.2)
  expect_equal(check_if_continuous(X), NULL)

  X$X1 <- sample(1:31, nrow(X), replace=TRUE)
  expect_equal(check_if_continuous(X), NULL)

  X$X1 <- sample(1:30, nrow(X), replace=TRUE)
  expect_warning(check_if_continuous(X))

})

# RANGER ----
test_that("ranger_importance_scores: output structure for regression and classification", {

  skip_if_not_installed("ranger")
  d <- readRDS(test_path("fixtures", "test_data.rds"))

  # regression: named numeric vector aligned with columns of X
  result_reg <- ranger_importance_scores(d$X, d$y_g, type = "regression", num.trees = 50)
  .valid_score(result_reg, d$X)

  # classification: same structure with binary factor response
  result_cls <- ranger_importance_scores(d$X, d$y_b, type = "classification", num.trees = 50)
  .valid_score(result_cls, d$X)

})

test_that("ranger_importance_scores: output structure for survival", {

  skip_if_not_installed("ranger")
  skip_if_not_installed("survival")

  d <- readRDS(test_path("fixtures", "test_data.rds"))

  result_surv <- ranger_importance_scores(d$X, d$y_s, type = "survival", num.trees = 50)
  .valid_score(result_surv, d$X)

})

test_that("ranger_importance_scores: ... forwarding", {

  skip_if_not_installed("ranger")
  d <- readRDS(test_path("fixtures", "test_data.rds"))

  # importance = "none" is a valid ranger setting,
  # but "permutation" is protected with a warning
  expect_warning(
    result_no_imp <- ranger_importance_scores(
      d$X, d$y_g, type = "regression", num.trees = 50,
      importance = "none"
    )
  )

  # unknown args should be silently filtered out by the intersect() call;
  # result should still have the correct structure
  result_unknown <- ranger_importance_scores(
    d$X, d$y_g, type = "regression",
    num.trees = 50,
    not_a_ranger_arg = TRUE
  )
  .valid_score(result_unknown, d$X)

})


test_that("ranger_importance_scores: protected arguments in dots raise warning", {

  skip_if_not_installed("ranger")
  d <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_warning(
    ranger_importance_scores(d$X, d$y_g, mtry = 1),
    "ignored"
  )

  expect_warning(
    ranger_importance_scores(d$X, d$y_g, importance = "none"),
    "ignored"
  )

})


test_that("ranger_importance_scores: invalid type raises an error", {

  skip_if_not_installed("ranger")
  d <- readRDS(test_path("fixtures", "test_data.rds"))
  expect_error(
    ranger_importance_scores(d$X, d$y_g, type = "invalid"),
    "should be one of"
  )

})

