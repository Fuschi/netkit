#' Build a Signed Partial Network from a Feature Matrix
#'
#' @param x Numeric matrix or data frame with samples in rows and features in columns.
#' @param transform Optional preprocessing before covariance estimation. Use `"none"`
#'   for raw values or `"rank"` to rank-transform each column first.
#' @param covariance_estimator Use `"sample"` for the sample covariance,
#'   `"shrinkage"` for shrinkage covariance, or `"auto"` to start from the sample
#'   covariance and regularize automatically when it is problematic.
#' @param shrinkage_args Optional named list or named vector of arguments forwarded to
#'   [corpcor::cov.shrink()]. Supported entries are `lambda`, `lambda.var`, `w`,
#'   and `verbose`.
#' @param lambda Partial-network penalty. A single numeric value fits a single
#'   graphical lasso model. `NULL` creates a default lambda path and selects a value
#'   by StARS. A numeric vector is treated as an explicit lambda path for StARS.
#' @param nlambda Number of lambda values when `lambda = NULL`.
#' @param lambda_min_ratio Ratio between the smallest and largest lambda when
#'   `lambda = NULL`.
#' @param rep_num Number of subsamples used by StARS when a lambda path is fitted.
#' @param stars_thresh StARS instability threshold.
#' @param subsample_ratio Optional subsampling ratio passed to [pulsar::pulsar()].
#' @param seed Optional random seed used by StARS.
#' @param glasso_args Optional named list or named vector of arguments forwarded to
#'   [glasso::glasso()]. Use `lambda` rather than `glasso_args$rho`.
#' @param output One of `"graph"`, `"adjacency"`, `"matrix"`, `"both"`, or `"details"`.
#' @param verbose Logical; if `TRUE`, report the main estimation steps.
#'
#' @return Depending on `output`, a signed graph, matrix, adjacency, or a detail list.
#' @description Reconstruct a signed partial-correlation network from a feature matrix.
#'   The function estimates a covariance matrix, fits a graphical lasso model, converts
#'   the precision matrix to signed partial correlations, and returns the sparse partial
#'   network directly from the glasso solution. When `lambda` is `NULL` or a numeric
#'   vector, the function selects the penalty through StARS before fitting the final model.
#'
#' @section Shrinkage Covariance Estimation:
#' When `covariance_estimator = "shrinkage"` or when `covariance_estimator = "auto"`
#' falls back to regularization, the function uses [corpcor::cov.shrink()] before
#' graphical lasso estimation. This stabilizes the covariance matrix in singular or
#' badly conditioned settings but can also bias the resulting partial network.
#'
#' @references
#' Jin, Y., Notredame, C., & Erb, I. (2023). \emph{A unifying and general account
#' of symmetric power transformations for compositional data}. arXiv:2212.00496.
#' \url{https://arxiv.org/abs/2212.00496}
#' @export
build_partial_net <- function(
    x,
    transform = c("none", "rank"),
    covariance_estimator = c("sample", "shrinkage", "auto"),
    shrinkage_args = list(),
    lambda = NULL,
    nlambda = 20,
    lambda_min_ratio = 0.01,
    rep_num = 20,
    stars_thresh = 0.1,
    subsample_ratio = NULL,
    seed = NULL,
    glasso_args = list(),
    output = c("graph", "adjacency", "matrix", "both", "details"),
    verbose = FALSE
) {
  transform <- match.arg(transform)
  covariance_estimator <- match.arg(covariance_estimator)
  output <- match.arg(output)

  .verbose_inform("Preparing input matrix.", verbose)
  x <- .validate_input_matrix(x)
  x_prepared <- .prepare_matrix(x, corr_method = if (identical(transform, "rank")) "spearman" else "pearson")

  shrinkage_args <- .coerce_named_list(shrinkage_args, "shrinkage_args")
  glasso_args <- .coerce_named_list(glasso_args, "glasso_args")

  selection <- NULL
  lambda_path <- NULL

  if (is.null(lambda)) {
    .verbose_inform("Estimating covariance and building the default lambda path.", verbose)
    cov_estimate <- .estimate_partial_covariance(
      x_prepared,
      covariance_estimator = covariance_estimator,
      shrinkage_args = shrinkage_args,
      warn = TRUE
    )
    lambda_path <- .default_lambda_path(
      cov_estimate$covariance,
      nlambda = nlambda,
      lambda_min_ratio = lambda_min_ratio
    )
    .verbose_inform("Running StARS over the default lambda path.", verbose)
    selection <- .run_partial_stars(
      x = x,
      lambda = lambda_path,
      transform = transform,
      covariance_estimator = covariance_estimator,
      shrinkage_args = shrinkage_args,
      glasso_args = glasso_args,
      rep_num = rep_num,
      stars_thresh = stars_thresh,
      subsample_ratio = subsample_ratio,
      seed = seed
    )
    lambda <- selection$selected_lambda
    .verbose_inform("Selected lambda = {signif(lambda, 6)}.", verbose)
  } else if (length(lambda) > 1) {
    lambda_path <- as.numeric(lambda)
    if (any(!is.finite(lambda_path)) || any(lambda_path < 0)) {
      cli::cli_abort("{.arg lambda} must contain only finite numeric values >= 0.")
    }
    .verbose_inform("Running StARS over the supplied lambda path.", verbose)
    selection <- .run_partial_stars(
      x = x,
      lambda = lambda_path,
      transform = transform,
      covariance_estimator = covariance_estimator,
      shrinkage_args = shrinkage_args,
      glasso_args = glasso_args,
      rep_num = rep_num,
      stars_thresh = stars_thresh,
      subsample_ratio = subsample_ratio,
      seed = seed
    )
    lambda <- selection$selected_lambda
    .verbose_inform("Selected lambda = {signif(lambda, 6)}.", verbose)
  } else {
    lambda <- as.numeric(lambda)
    if (!is.finite(lambda) || lambda < 0) cli::cli_abort("{.arg lambda} must be a single finite numeric value >= 0.")
  }

  .verbose_inform("Estimating covariance for the final graphical lasso fit.", verbose)
  cov_estimate <- .estimate_partial_covariance(
    x_prepared,
    covariance_estimator = covariance_estimator,
    shrinkage_args = shrinkage_args,
    warn = TRUE
  )

  .verbose_inform("Fitting the final graphical lasso model.", verbose)
  fit <- .fit_partial_glasso(cov_estimate$covariance, lambda, glasso_args = glasso_args)
  assoc <- fit$matrix
  adj <- assoc
  graph <- as_signed_graph(adj)

  details <- list(
    method = "partial",
    transform = transform,
    covariance_estimator = covariance_estimator,
    covariance_estimator_used = cov_estimate$estimator_used,
    association_matrix = assoc,
    matrix = assoc,
    adjacency = adj,
    graph = graph,
    lambda = lambda,
    lambda_path = lambda_path,
    selection = selection,
    glasso_fit = fit$fit,
    covariance_rank = cov_estimate$covariance_rank,
    covariance_rcond = cov_estimate$covariance_rcond,
    sample_covariance_rank = cov_estimate$sample_rank,
    sample_covariance_rcond = cov_estimate$sample_rcond,
    shrinkage_lambda = if (!is.null(cov_estimate$shrinkage_fit)) attr(cov_estimate$shrinkage_fit, "lambda") else NULL,
    shrinkage_lambda_estimated = if (!is.null(cov_estimate$shrinkage_fit)) attr(cov_estimate$shrinkage_fit, "lambda.estimated") else NULL,
    shrinkage_lambda_var = if (!is.null(cov_estimate$shrinkage_fit)) attr(cov_estimate$shrinkage_fit, "lambda.var") else NULL,
    shrinkage_lambda_var_estimated = if (!is.null(cov_estimate$shrinkage_fit)) attr(cov_estimate$shrinkage_fit, "lambda.var.estimated") else NULL,
    shrinkage_args = if (identical(cov_estimate$estimator_used, "shrinkage") || identical(cov_estimate$estimator_used, "shrinkage_fallback")) shrinkage_args else NULL,
    glasso_args = glasso_args
  )

  switch(
    output,
    graph = graph,
    adjacency = adj,
    matrix = assoc,
    both = list(adjacency = adj, graph = graph),
    details = details
  )
}
