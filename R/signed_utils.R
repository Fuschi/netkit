#' Convert a Signed Adjacency Matrix to an `igraph`
#'
#' @param adj A square symmetric numeric matrix.
#' @param weighted Logical; if `TRUE`, store edge weights in `weight`.
#'
#' @return An undirected `igraph` object.
#' @description Convert a signed adjacency matrix into an undirected `igraph`.
#' @export
as_signed_graph <- function(adj, weighted = TRUE) {
  .assert_square_symmetric(adj, "adj")

  if (is.null(colnames(adj))) {
    colnames(adj) <- seq_len(ncol(adj))
  }
  if (is.null(rownames(adj))) {
    rownames(adj) <- colnames(adj)
  }

  if (all(adj == 0 | is.na(adj))) {
    g <- igraph::make_empty_graph(n = ncol(adj), directed = FALSE)
    igraph::V(g)$name <- colnames(adj)
    return(g)
  }

  igraph::graph_from_adjacency_matrix(
    adjmatrix = adj,
    mode = "undirected",
    weighted = weighted,
    diag = FALSE
  )
}

#' Threshold a Signed Association Matrix
#'
#' @param matrix A square symmetric signed association matrix.
#' @param method One of `"none"`, `"absolute"`, `"density"`, or `"pvalue"`.
#' @param value Threshold value interpreted according to `method`.
#' @param pval_mat Optional symmetric p-value matrix used when `method = "pvalue"`.
#' @param adjust Multiple-testing correction passed to [stats::p.adjust()].
#'
#' @return A signed adjacency matrix with diagonal set to zero.
#' @description Apply sign-preserving sparsification to a signed association matrix.
#' @export
threshold_signed <- function(
    matrix,
    method = c("none", "absolute", "density", "pvalue"),
    value = NULL,
    pval_mat = NULL,
    adjust = c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY")
) {
  method <- match.arg(method)
  adjust <- match.arg(adjust)
  .assert_square_symmetric(matrix, "matrix")

  adj <- matrix
  diag(adj) <- 0

  if (identical(method, "none")) {
    adj[is.na(adj)] <- 0
    return(adj)
  }

  if (is.null(value) || !is.numeric(value) || length(value) != 1) cli::cli_abort("{.arg value} must be a single numeric threshold.")

  if (identical(method, "absolute")) {
    if (value < 0 || value > 1) cli::cli_abort("For {.val absolute}, {.arg value} must be in [0, 1].")
    adj[abs(adj) < value | is.na(adj)] <- 0
    return(adj)
  }

  if (identical(method, "density")) {
    if (value <= 0 || value > 1) cli::cli_abort("For {.val density}, {.arg value} must be in (0, 1].")

    upper <- upper.tri(adj)
    values <- abs(adj[upper])
    valid <- !is.na(values)
    if (!any(valid)) {
      adj[,] <- 0
      return(adj)
    }

    target_edges <- ceiling(value * sum(upper))
    sorted <- sort(values[valid], decreasing = TRUE)
    cutoff <- if (target_edges < length(sorted)) sorted[target_edges] else min(sorted)
    adj[abs(adj) < cutoff | is.na(adj)] <- 0
    return(adj)
  }

  if (is.null(pval_mat)) cli::cli_abort("{.arg pval_mat} is required when {.arg method = \"pvalue\"}.")
  .assert_square_symmetric(pval_mat, "pval_mat")
  if (!all(dim(matrix) == dim(pval_mat))) cli::cli_abort("{.arg matrix} and {.arg pval_mat} must have the same dimensions.")
  if (value <= 0 || value > 1) cli::cli_abort("For {.val pvalue}, {.arg value} must be in (0, 1].")

  upper <- upper.tri(pval_mat)
  pvals <- pval_mat[upper]
  keep <- rep(FALSE, length(pvals))
  valid <- !is.na(pvals)
  if (any(valid)) {
    adjusted <- stats::p.adjust(pvals[valid], method = adjust)
    keep[valid] <- adjusted <= value
  }

  mask <- matrix(FALSE, nrow(matrix), ncol(matrix))
  mask[upper] <- keep
  mask <- mask | t(mask)
  adj[!mask | is.na(adj)] <- 0
  adj
}

#' Validate a Feature Matrix
#'
#' @param x A candidate feature matrix or data frame.
#'
#' @return A numeric matrix with the original columns preserved.
#' @keywords internal
#' @noRd
.validate_input_matrix <- function(x) {
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.numeric(x)) cli::cli_abort("{.arg x} must be a numeric matrix or data.frame.")
  if (is.null(colnames(x))) cli::cli_abort("{.arg x} must have column names.")
  if (ncol(x) < 2) cli::cli_abort("{.arg x} must contain at least two columns.")

  x
}

#' Coerce a Forwarded Argument Container to a Named List
#'
#' @param x Candidate list or named vector.
#' @param arg Name of the calling argument.
#'
#' @return A named list.
#' @keywords internal
#' @noRd
.coerce_named_list <- function(x, arg) {
  if (is.list(x)) return(x)
  if (is.atomic(x) && !is.null(names(x))) return(as.list(x))

  cli::cli_abort("{.arg {arg}} must be a named list or named vector.")
}

#' Validate Names in a Forwarded Argument List
#'
#' @param x Named list to validate.
#' @param allowed Character vector of supported names.
#' @param arg Name of the calling argument.
#'
#' @return The original list, invisibly.
#' @keywords internal
#' @noRd
.validate_forward_args <- function(x, allowed, arg) {
  unknown <- setdiff(names(x), allowed)
  if (length(unknown) > 0) cli::cli_abort("{.arg {arg}} contains unsupported entries: {.val {unknown}}.")

  invisible(x)
}

#' Assert That a Matrix Is Square and Symmetric
#'
#' @param x Matrix to validate.
#' @param arg Name of the calling argument.
#'
#' @return Called for its side effects; errors on invalid input.
#' @keywords internal
#' @noRd
.assert_square_symmetric <- function(x, arg) {
  if (!is.matrix(x)) cli::cli_abort("{.arg {arg}} must be a matrix.")
  if (!is.numeric(x)) cli::cli_abort("{.arg {arg}} must be numeric.")
  if (nrow(x) != ncol(x)) cli::cli_abort("{.arg {arg}} must be square.")
  if (!isSymmetric(x)) cli::cli_abort("{.arg {arg}} must be symmetric.")
}

#' Prepare a Matrix for Marginal Reconstruction
#'
#' @param x Numeric feature matrix.
#' @param corr_method Correlation flavor used downstream.
#'
#' @return A matrix with optional rank transform applied and constant columns removed.
#' @keywords internal
#' @noRd
.prepare_matrix <- function(x, corr_method) {
  if (identical(corr_method, "spearman")) {
    x <- apply(x, 2, rank, ties.method = "average", na.last = "keep")
    x <- matrix(x, nrow = nrow(x), dimnames = dimnames(x))
  }

  sds <- apply(x, 2, stats::sd, na.rm = TRUE)
  keep <- !(is.na(sds) | sds == 0)
  if (sum(keep) < 2) cli::cli_abort("At least two columns must have non-zero variance.")

  x[, keep, drop = FALSE]
}

#' Compute Pairwise Correlations and Optional P-values
#'
#' @param x Prepared numeric feature matrix.
#' @param corr_method Correlation method passed to [stats::cor()] or [stats::cor.test()].
#' @param compute_p Logical; if `TRUE`, also compute pairwise p-values.
#'
#' @return A list with components `matrix` and `p`.
#' @keywords internal
#' @noRd
.pairwise_correlation <- function(x, corr_method, compute_p = FALSE) {
  vars <- colnames(x)
  cor_mat <- matrix(NA_real_, ncol(x), ncol(x), dimnames = list(vars, vars))
  p_mat <- if (compute_p) matrix(NA_real_, ncol(x), ncol(x), dimnames = list(vars, vars)) else NULL

  cor_sub <- suppressWarnings(
    stats::cor(x, method = corr_method, use = "pairwise.complete.obs")
  )
  diag(cor_sub) <- 1
  cor_mat[,] <- cor_sub

  if (!compute_p) return(list(matrix = cor_mat, p = NULL))

  pairs <- utils::combn(seq_len(ncol(x)), 2, simplify = FALSE)
  for (pair in pairs) {
    xi <- x[, pair[1]]
    xj <- x[, pair[2]]
    ok <- stats::complete.cases(xi, xj)
    if (sum(ok) < 3) next

    test <- tryCatch(
      stats::cor.test(xi[ok], xj[ok], method = corr_method),
      error = function(e) NULL
    )
    if (is.null(test)) next

    p_mat[pair[1], pair[2]] <- test$p.value
    p_mat[pair[2], pair[1]] <- test$p.value
  }

  diag(p_mat) <- 0
  list(matrix = cor_mat, p = p_mat)
}

#' Estimate a Shrinkage Correlation Matrix
#'
#' @param x Prepared numeric feature matrix.
#' @param shrinkage_args Named list forwarded to [corpcor::cor.shrink()].
#'
#' @return A list with `matrix` and `fit`.
#' @keywords internal
#' @noRd
.estimate_shrinkage_correlation <- function(x, shrinkage_args = list()) {
  .validate_forward_args(
    shrinkage_args,
    allowed = c("lambda", "w", "verbose"),
    arg = "shrinkage_args"
  )

  shrinkage_call <- utils::modifyList(
    list(x = x, verbose = FALSE),
    shrinkage_args
  )
  fit <- do.call(corpcor::cor.shrink, shrinkage_call)
  cor_mat <- unclass(fit)
  diag(cor_mat) <- 1
  dimnames(cor_mat) <- list(colnames(x), colnames(x))

  list(matrix = cor_mat, fit = fit)
}

#' Estimate a Shrinkage Covariance Matrix
#'
#' @param x Prepared numeric feature matrix.
#' @param shrinkage_args Named list forwarded to [corpcor::cov.shrink()].
#'
#' @return A list with `matrix` and `fit`.
#' @keywords internal
#' @noRd
.estimate_shrinkage_covariance <- function(x, shrinkage_args = list()) {
  .validate_forward_args(
    shrinkage_args,
    allowed = c("lambda", "lambda.var", "w", "verbose"),
    arg = "shrinkage_args"
  )

  shrinkage_call <- utils::modifyList(
    list(x = x, verbose = FALSE),
    shrinkage_args
  )
  fit <- do.call(corpcor::cov.shrink, shrinkage_call)
  cov_mat <- unclass(fit)
  dimnames(cov_mat) <- list(colnames(x), colnames(x))

  list(matrix = cov_mat, fit = fit)
}

#' Apply Multiple Threshold Criteria by Intersection
#'
#' @param matrix A square symmetric signed association matrix.
#' @param threshold_value Named list or named vector with entries `abs`, `dens`, and/or `p`.
#' @param pval_mat Optional p-value matrix used when `p` is present.
#' @param adjust Multiple-testing correction passed to [stats::p.adjust()].
#'
#' @return A thresholded signed adjacency matrix.
#' @keywords internal
#' @noRd
.apply_threshold_list <- function(matrix, threshold_value, pval_mat = NULL, adjust = "none") {
  threshold_value <- .coerce_named_list(threshold_value, "threshold_value")
  .validate_forward_args(threshold_value, c("abs", "dens", "p"), "threshold_value")

  masks <- list()
  if (!is.null(threshold_value$abs)) {
    masks$abs <- threshold_signed(matrix, method = "absolute", value = threshold_value$abs) != 0
  }
  if (!is.null(threshold_value$dens)) {
    masks$dens <- threshold_signed(matrix, method = "density", value = threshold_value$dens) != 0
  }
  if (!is.null(threshold_value$p)) {
    masks$p <- threshold_signed(
      matrix,
      method = "pvalue",
      value = threshold_value$p,
      pval_mat = pval_mat,
      adjust = adjust
    ) != 0
  }

  if (length(masks) == 0) cli::cli_abort("{.arg threshold_value} must contain at least one of {.val abs}, {.val dens}, or {.val p}.")

  keep <- Reduce(`&`, masks)
  adj <- matrix * keep
  diag(adj) <- 0
  adj[is.na(adj)] <- 0
  adj
}

#' Adjust a Symmetric P-value Matrix
#'
#' @param pval_mat Symmetric p-value matrix.
#' @param adjust Multiple-testing correction passed to [stats::p.adjust()].
#'
#' @return A symmetric adjusted p-value matrix.
#' @keywords internal
#' @noRd
.adjust_pvalue_matrix <- function(pval_mat, adjust = "none") {
  if (is.null(pval_mat)) return(NULL)

  .assert_square_symmetric(pval_mat, "pval_mat")
  adj_p <- matrix(NA_real_, nrow(pval_mat), ncol(pval_mat), dimnames = dimnames(pval_mat))
  upper <- upper.tri(pval_mat)
  values <- pval_mat[upper]
  valid <- !is.na(values)

  if (any(valid)) {
    adjusted <- stats::p.adjust(values[valid], method = adjust)
    upper_idx <- which(upper)[valid]
    lower_idx <- which(t(upper))[valid]
    adj_p[upper_idx] <- adjusted
    adj_p[lower_idx] <- adjusted
  }

  diag(adj_p) <- 0
  adj_p
}

#' Check Whether a Covariance Matrix Is Problematic for Inversion
#'
#' @param covariance Covariance matrix to diagnose.
#' @param tol Numerical tolerance used for the reciprocal condition number.
#'
#' @return A list with logical flag `problematic`, matrix rank, and reciprocal condition number.
#' @keywords internal
#' @noRd
.check_covariance_matrix <- function(covariance, tol = 1e-10) {
  rank <- qr(covariance)$rank
  reciprocal_condition <- tryCatch(
    rcond(covariance),
    error = function(e) 0
  )

  problematic <- any(!is.finite(covariance)) ||
    rank < ncol(covariance) ||
    !is.finite(reciprocal_condition) ||
    reciprocal_condition < tol

  list(
    problematic = problematic,
    rank = rank,
    rcond = reciprocal_condition
  )
}
