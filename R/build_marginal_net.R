#' Build a Signed Marginal Network from a Feature Matrix
#'
#' @param x Numeric matrix or data frame with samples in rows and features in columns.
#' @param corr_method Correlation flavor used to reconstruct the signed association matrix.
#' @param covariance_estimator Use `"sample"` for ordinary
#'   sample correlation, `"shrinkage"` for shrinkage correlation, or `"auto"` to
#'   start from the sample estimator and regularize automatically when the
#'   covariance matrix is problematic.
#' @param shrinkage_args Optional named list or named vector of arguments forwarded to
#'   [corpcor::cor.shrink()]. Supported entries are `lambda`, `w`, and `verbose`.
#' @param threshold_method Thresholding rule applied after marginal reconstruction.
#'   One of `"none"`, `"absolute"`, `"density"`, or `"pvalue"`.
#' @param threshold_value Either a single numeric threshold used with `threshold_method`,
#'   or a named list/vector with entries `abs`, `dens`, and/or `p`. When a named
#'   value is supplied, the corresponding thresholded edge sets are intersected.
#' @param adjust Multiple-testing correction used when p-value thresholding is requested.
#' @param output One of `"graph"`, `"adjacency"`, `"matrix"`, `"both"`, or `"details"`.
#'
#' @return Depending on `output`, a signed graph, matrix, adjacency, or a detail list.
#' @description Reconstruct a signed marginal association matrix, optionally
#'   apply thresholding, and return it as a matrix and/or graph. When
#'   `covariance_estimator = "sample"`, a singular or badly conditioned covariance
#'   matrix triggers a warning. When `covariance_estimator = "auto"`, the function
#'   warns and switches to shrinkage correlation instead.
#'
#' @section Shrinkage Covariance Estimation:
#' When `covariance_estimator = "shrinkage"` or when `covariance_estimator = "auto"`
#' falls back to regularization, the function uses [corpcor::cor.shrink()] to
#' stabilize the marginal association estimate. This can be useful when the sample
#' covariance is singular or badly conditioned, but it also introduces bias by
#' shrinking the association structure toward a more regular target.
#'
#' @references
#' Jin, Y., Notredame, C., & Erb, I. (2023). \emph{A unifying and general account
#' of symmetric power transformations for compositional data}. arXiv:2212.00496.
#' \url{https://arxiv.org/abs/2212.00496}
#' @export
build_marginal_net <- function(
    x,
    corr_method = c("pearson", "spearman"),
    covariance_estimator = c("sample", "shrinkage", "auto"),
    shrinkage_args = list(),
    threshold_method = c("none", "absolute", "density", "pvalue"),
    threshold_value = NULL,
    adjust = c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY"),
    output = c("graph", "adjacency", "matrix", "both", "details")
) {
  corr_method <- match.arg(corr_method)
  covariance_estimator <- match.arg(covariance_estimator)
  threshold_method <- match.arg(threshold_method)
  adjust <- match.arg(adjust)
  output <- match.arg(output)

  x <- .validate_input_matrix(x)
  x_prepared <- .prepare_matrix(x, corr_method = corr_method)

  shrinkage_args <- .coerce_named_list(shrinkage_args, "shrinkage_args")
  threshold_is_named <- !is.null(threshold_value) && !is.null(names(threshold_value))
  compute_p <- identical(threshold_method, "pvalue") || (threshold_is_named && "p" %in% names(threshold_value))

  covariance_check <- .check_covariance_matrix(stats::cov(x_prepared, use = "pairwise.complete.obs"))
  covariance_estimator_used <- covariance_estimator
  shrinkage_fit <- NULL

  if (identical(covariance_estimator, "sample")) {
    if (covariance_check$problematic) {
      cli::cli_warn(c(
        "Covariance matrix is singular or poorly conditioned.",
        "i" = "rank = {covariance_check$rank}, reciprocal condition number = {signif(covariance_check$rcond, 3)}",
        "i" = "Use {.arg covariance_estimator = \"shrinkage\"} or {.arg covariance_estimator = \"auto\"}."
      ))
    }

    estimate <- .pairwise_correlation(
      x_prepared,
      corr_method = if (identical(corr_method, "spearman")) "pearson" else corr_method,
      compute_p = compute_p
    )
    assoc <- estimate$matrix
    pval_mat <- estimate$p
  } else if (identical(covariance_estimator, "auto") && !covariance_check$problematic) {
    estimate <- .pairwise_correlation(
      x_prepared,
      corr_method = if (identical(corr_method, "spearman")) "pearson" else corr_method,
      compute_p = compute_p
    )
    assoc <- estimate$matrix
    pval_mat <- estimate$p
  } else {
    if (identical(covariance_estimator, "auto") && covariance_check$problematic) {
      cli::cli_warn(c(
        "Covariance matrix is singular or poorly conditioned; regularizing before marginal reconstruction.",
        "i" = "rank = {covariance_check$rank}, reciprocal condition number = {signif(covariance_check$rcond, 3)}"
      ))
      covariance_estimator_used <- "shrinkage_fallback"
    }

    if (compute_p) {
      cli::cli_abort("P-value thresholding is only available when the sample covariance estimator is used.")
    }

    shrinkage_estimate <- .estimate_shrinkage_correlation(x_prepared, shrinkage_args)
    shrinkage_fit <- shrinkage_estimate$fit
    assoc <- shrinkage_estimate$matrix
    pval_mat <- NULL
  }

  pval_mat_adjusted <- .adjust_pvalue_matrix(pval_mat, adjust = adjust)

  if (threshold_is_named) {
    adj <- .apply_threshold_list(
      matrix = assoc,
      threshold_value = threshold_value,
      pval_mat = pval_mat,
      adjust = adjust
    )
  } else {
    adj <- threshold_signed(
      matrix = assoc,
      method = threshold_method,
      value = threshold_value,
      pval_mat = pval_mat,
      adjust = adjust
    )
  }
  graph <- as_signed_graph(adj)

  details <- list(
    method = "marginal",
    corr_method = corr_method,
    covariance_estimator = covariance_estimator,
    covariance_estimator_used = covariance_estimator_used,
    association_matrix = assoc,
    matrix = assoc,
    adjacency = adj,
    graph = graph,
    p_value = pval_mat,
    p_value_adjusted = pval_mat_adjusted,
    threshold_method = threshold_method,
    threshold_value = threshold_value,
    covariance_rank = covariance_check$rank,
    covariance_rcond = covariance_check$rcond,
    shrinkage_lambda = if (!is.null(shrinkage_fit)) attr(shrinkage_fit, "lambda") else NULL,
    shrinkage_lambda_estimated = if (!is.null(shrinkage_fit)) attr(shrinkage_fit, "lambda.estimated") else NULL,
    shrinkage_args = if (identical(covariance_estimator_used, "shrinkage") || identical(covariance_estimator_used, "shrinkage_fallback")) shrinkage_args else NULL
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
