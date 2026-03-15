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

#' Emit a Progress Message When Verbose Mode Is Enabled
#'
#' @param message Message text.
#' @param verbose Logical flag controlling message emission.
#'
#' @return Called for side effects only.
#' @keywords internal
#' @noRd
.verbose_inform <- function(message, verbose = FALSE) {
  if (isTRUE(verbose)) cli::cli_inform(message)
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

#' Convert a Precision Matrix to Partial Correlations
#'
#' @param precision Precision matrix.
#' @param vars Variable names.
#'
#' @return A symmetric partial-correlation matrix.
#' @keywords internal
#' @noRd
.precision_to_partial <- function(precision, vars = colnames(precision)) {
  scale <- sqrt(diag(precision))
  partial <- -precision / outer(scale, scale)
  partial <- (partial + t(partial)) / 2
  diag(partial) <- 1
  dimnames(partial) <- list(vars, vars)
  partial
}

#' Estimate the Covariance Matrix Used by Partial Reconstruction
#'
#' @param x Prepared numeric feature matrix.
#' @param covariance_estimator One of `"sample"`, `"shrinkage"`, or `"auto"`.
#' @param shrinkage_args Named list forwarded to [corpcor::cov.shrink()].
#' @param warn Logical; if `TRUE`, emit warnings for automatic fallback.
#'
#' @return A list with covariance diagnostics and the covariance matrix used downstream.
#' @keywords internal
#' @noRd
.estimate_partial_covariance <- function(
    x,
    covariance_estimator = c("sample", "shrinkage", "auto"),
    shrinkage_args = list(),
    warn = TRUE
) {
  covariance_estimator <- match.arg(covariance_estimator)
  sample_covariance <- stats::cov(x, use = "pairwise.complete.obs")
  sample_check <- .check_covariance_matrix(sample_covariance)

  if (identical(covariance_estimator, "sample")) {
    if (sample_check$problematic) {
      cli::cli_abort(c(
        "Covariance matrix is singular or poorly conditioned.",
        "i" = "rank = {sample_check$rank}, reciprocal condition number = {signif(sample_check$rcond, 3)}",
        "i" = "Use {.arg covariance_estimator = \"shrinkage\"} or {.arg covariance_estimator = \"auto\"}."
      ))
    }

    return(list(
      covariance = sample_covariance,
      estimator_requested = covariance_estimator,
      estimator_used = "sample",
      sample_rank = sample_check$rank,
      sample_rcond = sample_check$rcond,
      covariance_rank = sample_check$rank,
      covariance_rcond = sample_check$rcond,
      shrinkage_fit = NULL
    ))
  }

  if (identical(covariance_estimator, "auto") && !sample_check$problematic) {
    return(list(
      covariance = sample_covariance,
      estimator_requested = covariance_estimator,
      estimator_used = "sample",
      sample_rank = sample_check$rank,
      sample_rcond = sample_check$rcond,
      covariance_rank = sample_check$rank,
      covariance_rcond = sample_check$rcond,
      shrinkage_fit = NULL
    ))
  }

  if (identical(covariance_estimator, "auto") && sample_check$problematic && warn) {
    cli::cli_warn(c(
      "Covariance matrix is singular or poorly conditioned; regularizing before partial reconstruction.",
      "i" = "rank = {sample_check$rank}, reciprocal condition number = {signif(sample_check$rcond, 3)}"
    ))
  }

  shrinkage_estimate <- .estimate_shrinkage_covariance(x, shrinkage_args)
  shrinkage_check <- .check_covariance_matrix(shrinkage_estimate$matrix)

  list(
    covariance = shrinkage_estimate$matrix,
    estimator_requested = covariance_estimator,
    estimator_used = if (identical(covariance_estimator, "auto")) "shrinkage_fallback" else "shrinkage",
    sample_rank = sample_check$rank,
    sample_rcond = sample_check$rcond,
    covariance_rank = shrinkage_check$rank,
    covariance_rcond = shrinkage_check$rcond,
    shrinkage_fit = shrinkage_estimate$fit
  )
}

#' Build a Default Lambda Path for Partial Reconstruction
#'
#' @param covariance Covariance matrix used by graphical lasso.
#' @param nlambda Number of lambda values.
#' @param lambda_min_ratio Ratio between the smallest and largest lambda.
#'
#' @return A decreasing numeric lambda path.
#' @keywords internal
#' @noRd
.default_lambda_path <- function(covariance, nlambda = 20, lambda_min_ratio = 0.01) {
  if (!is.numeric(nlambda) || length(nlambda) != 1 || nlambda < 2) {
    cli::cli_abort("{.arg nlambda} must be a single numeric value >= 2.")
  }
  if (!is.numeric(lambda_min_ratio) || length(lambda_min_ratio) != 1 || lambda_min_ratio <= 0 || lambda_min_ratio > 1) {
    cli::cli_abort("{.arg lambda_min_ratio} must be a single numeric value in (0, 1].")
  }

  upper <- upper.tri(covariance)
  lambda_max <- max(abs(covariance[upper]), na.rm = TRUE)
  if (!is.finite(lambda_max) || lambda_max <= 0) lambda_max <- 1

  exp(seq(log(lambda_max), log(lambda_max * lambda_min_ratio), length.out = nlambda))
}

#' Fit a Graphical Lasso Partial Model
#'
#' @param covariance Covariance matrix used for graphical lasso.
#' @param lambda Penalty value passed as `rho`.
#' @param glasso_args Named list forwarded to [glasso::glasso()].
#'
#' @return A list containing the glasso fit and the partial-correlation matrix.
#' @keywords internal
#' @noRd
.fit_partial_glasso <- function(covariance, lambda, glasso_args = list()) {
  glasso_args <- .coerce_named_list(glasso_args, "glasso_args")
  .validate_forward_args(
    glasso_args,
    allowed = c("nobs", "zero", "thr", "maxit", "approx", "penalize.diagonal", "start", "w.init", "wi.init", "trace"),
    arg = "glasso_args"
  )
  if ("rho" %in% names(glasso_args)) {
    cli::cli_abort("{.arg glasso_args} must not contain {.val rho}; use {.arg lambda} instead.")
  }
  if (!is.numeric(lambda) || length(lambda) != 1 || !is.finite(lambda) || lambda < 0) {
    cli::cli_abort("{.arg lambda} must be a single finite numeric value >= 0.")
  }

  fit <- do.call(glasso::glasso, utils::modifyList(list(s = covariance, rho = lambda), glasso_args))
  partial <- .precision_to_partial(fit$wi, colnames(covariance))
  diag(partial) <- 0

  list(fit = fit, matrix = partial)
}

#' Build a Glasso Support Path for Pulsar/StARS
#'
#' @param data Feature matrix.
#' @param lambda Lambda path.
#' @param transform One of `"none"` or `"rank"`.
#' @param covariance_estimator One of `"sample"`, `"shrinkage"`, or `"auto"`.
#' @param shrinkage_args Named list forwarded to [corpcor::cov.shrink()].
#' @param glasso_args Named list forwarded to [glasso::glasso()].
#'
#' @return A list with a `path` element containing support matrices.
#' @keywords internal
#' @noRd
.partial_glasso_path <- function(
    data,
    lambda,
    transform = c("none", "rank"),
    covariance_estimator = c("sample", "shrinkage", "auto"),
    shrinkage_args = list(),
    glasso_args = list()
) {
  transform <- match.arg(transform)
  covariance_estimator <- match.arg(covariance_estimator)
  data <- .validate_input_matrix(data)
  x_prepared <- .prepare_matrix(data, corr_method = if (identical(transform, "rank")) "spearman" else "pearson")
  cov_estimate <- .estimate_partial_covariance(
    x_prepared,
    covariance_estimator = covariance_estimator,
    shrinkage_args = shrinkage_args,
    warn = FALSE
  )

  path <- lapply(lambda, function(rho) {
    fit <- .fit_partial_glasso(cov_estimate$covariance, rho, glasso_args = glasso_args)
    storage.mode(fit$matrix) <- "numeric"
    (fit$matrix != 0) * 1
  })

  list(path = path)
}

#' Run StARS Selection for the Partial Reconstruction Path
#'
#' @param x Prepared feature matrix.
#' @param lambda Lambda path.
#' @param transform One of `"none"` or `"rank"`.
#' @param covariance_estimator One of `"sample"`, `"shrinkage"`, or `"auto"`.
#' @param shrinkage_args Named list forwarded to [corpcor::cov.shrink()].
#' @param glasso_args Named list forwarded to [glasso::glasso()].
#' @param rep_num Number of subsamples used by StARS.
#' @param stars_thresh StARS instability threshold.
#' @param subsample_ratio Optional subsampling ratio.
#' @param seed Optional random seed.
#'
#' @return A list with the pulsar fit, selected lambda, lambda path, and selected index.
#' @keywords internal
#' @noRd
.run_partial_stars <- function(
    x,
    lambda,
    transform = c("none", "rank"),
    covariance_estimator = c("sample", "shrinkage", "auto"),
    shrinkage_args = list(),
    glasso_args = list(),
    rep_num = 20,
    stars_thresh = 0.1,
    subsample_ratio = NULL,
    seed = NULL
) {
  transform <- match.arg(transform)
  covariance_estimator <- match.arg(covariance_estimator)
  pulsar_fit <- pulsar::pulsar(
    data = x,
    fun = .partial_glasso_path,
    fargs = list(
      lambda = lambda,
      transform = transform,
      covariance_estimator = covariance_estimator,
      shrinkage_args = shrinkage_args,
      glasso_args = glasso_args
    ),
    criterion = "stars",
    thresh = stars_thresh,
    subsample.ratio = subsample_ratio,
    rep.num = rep_num,
    seed = seed,
    refit = FALSE
  )

  selected_index <- pulsar_fit$stars$opt.index
  list(
    pulsar = pulsar_fit,
    stars = pulsar_fit$stars,
    lambda = lambda,
    selected_index = selected_index,
    selected_lambda = lambda[selected_index]
  )
}
