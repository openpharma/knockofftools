# Package index

## All functions

- [`eval_fdp()`](https://openpharma.github.io/knockofftools/reference/eval_fdp.md)
  : False discovery proportion (fdp) as function of selection and known
  negatives:
- [`eval_tpp()`](https://openpharma.github.io/knockofftools/reference/eval_tpp.md)
  : True positive proportion (tpp) as function of selection and known
  positives:
- [`generate_X()`](https://openpharma.github.io/knockofftools/reference/generate_X.md)
  : Simulate Gaussian and binary covariate predictors
- [`generate_simdata()`](https://openpharma.github.io/knockofftools/reference/generate_simdata.md)
  : Generate simulated data set for vignettes
- [`generate_y()`](https://openpharma.github.io/knockofftools/reference/generate_y.md)
  : Simulate Gaussian response from a sparse regression model
- [`glasso_adjacency_matrix()`](https://openpharma.github.io/knockofftools/reference/glasso_adjacency_matrix.md)
  : Estimate adjacency matrix using graphical LASSO (glasso)
- [`knockoff()`](https://openpharma.github.io/knockofftools/reference/knockoff.md)
  : Knockoffs for general covariate data frames
- [`knockoff.statistics()`](https://openpharma.github.io/knockofftools/reference/knockoff.statistics.md)
  : Knockoff (feature) statistics:
- [`knockoffs_mx()`](https://openpharma.github.io/knockofftools/reference/knockoffs_mx.md)
  : Gaussian MX-knockoffs for continuous variables
- [`knockoffs_seq()`](https://openpharma.github.io/knockofftools/reference/knockoffs_seq.md)
  : Sequential knockoffs for continuous and categorical variables
- [`knockoffs_sparse_seq()`](https://openpharma.github.io/knockofftools/reference/knockoffs_sparse_seq.md)
  : Sparse sequential knockoff generation algorithm
- [`multi_select()`](https://openpharma.github.io/knockofftools/reference/multi_select.md)
  : Select variables based on the heuristic multiple selection algorithm
  from Kormaksson et al. 'Sequential knockoffs for continuous and
  categorical predictors: With application to a large psoriatic
  arthritis clinical trial pool.' Statistics in Medicine. 2021;1–16.
- [`plot(`*`<variable.selections>`*`)`](https://openpharma.github.io/knockofftools/reference/plot.variable.selections.md)
  : Heatmap of multiple variable selections ordered by importance
- [`selections_control_FDR()`](https://openpharma.github.io/knockofftools/reference/selections_control_FDR.md)
  : Controls the false discovery rate (FDR) given knockoff W-statistics.
- [`selections_control_PFER()`](https://openpharma.github.io/knockofftools/reference/selections_control_PFER.md)
  : Controls the per-familywise error rate (PFER) given knockoff
  W-statistics.
- [`selections_control_kFWER()`](https://openpharma.github.io/knockofftools/reference/selections_control_kFWER.md)
  : Controls the k-familywise error rate (k-FWER) given a vector of
  knockoff W-statistics.
- [`sim_glmnet()`](https://openpharma.github.io/knockofftools/reference/sim_glmnet.md)
  : Simulate from glmnet penalized regression model
- [`sim_simple()`](https://openpharma.github.io/knockofftools/reference/sim_simple.md)
  : Simple knockoff generator based on least squares fit (continuous
  variables) or multinomial logistic regression (factor variables)
  respectively. If X is empty, knockoffs are sampled from the marginal
  distribution of y
- [`simdata`](https://openpharma.github.io/knockofftools/reference/simdata.md)
  : Simulated dataset for knockofftools
- [`simulWeib()`](https://openpharma.github.io/knockofftools/reference/simulWeib.md)
  : Function that simulates response from Cox model with Weibull
  baseline hazard:
- [`stat_glmnet()`](https://openpharma.github.io/knockofftools/reference/stat_glmnet.md)
  : Knockoff (feature) statistics: Absolute elastic-net coefficient
  differences between original and knockoff variables
- [`stat_predictive_causal_forest()`](https://openpharma.github.io/knockofftools/reference/stat_predictive_causal_forest.md)
  : Causal forest based knockoff (feature) statistics that captues the
  predictive strength: Difference from importance scores derived by
  causal forest
- [`stat_predictive_glmnet()`](https://openpharma.github.io/knockofftools/reference/stat_predictive_glmnet.md)
  : Knockoff (feature) statistics that captues the predictive strength:
  Absolute coefficient differences between treatment original variables
  interaction terms and treatment knockoff variables interaction terms
- [`stat_random_forest()`](https://openpharma.github.io/knockofftools/reference/stat_random_forest.md)
  : Knockoff (feature) statistics: Random forest
- [`variable.selections()`](https://openpharma.github.io/knockofftools/reference/variable.selections.md)
  : Knockoff variable selection: Select the variables by controlling a
  user-specified error rate
