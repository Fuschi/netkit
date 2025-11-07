# =========================
# Correlation & P-values
# =========================

#' Calculate correlation matrix (and p-values on demand)
#'
#' @description
#' Compute correlation and optionally p-values from a numeric matrix.
#' If `what = "cor"`, uses base `cor()` for speed. If `what = "cor_p"`, uses `cor.test()`
#' for each pair. Zero-variance columns are excluded from computation and return NA.
#'
#' @param x Numeric matrix or data.frame (samples in rows, variables in columns).
#' @param method One of "pearson", "spearman", or "kendall".
#' @param what Either "cor" (default) or "cor_p" to compute both correlation and p-values.
#'
#' @return A list with:
#'   - `cor`: correlation matrix (always returned)
#'   - `p`: p-value matrix (only if `what = "cor_p"`)
#'
#' @export
correlate_matrix <- function(
    x,
    method = c("pearson", "spearman", "kendall"),
    what = c("cor", "cor_p")
) {
  method <- match.arg(method)
  what <- match.arg(what)
  
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.numeric(x)) stop("`x` must be a numeric matrix or data.frame.")
  if (is.null(colnames(x))) stop("`x` must have column names.")
  
  p <- ncol(x)
  vars <- colnames(x)
  
  # 1. Identify zero-variance columns
  sds <- apply(x, 2, stats::sd, na.rm = TRUE)
  zero_var <- is.na(sds) | sds == 0
  nonzero_idx <- which(!zero_var)
  
  # Initialize full matrices
  cor_mat <- matrix(NA_real_, p, p, dimnames = list(vars, vars))
  pval_mat <- if (what == "cor_p") matrix(NA_real_, p, p, dimnames = list(vars, vars)) else NULL
  
  if (length(nonzero_idx) < 2) {
    warning("Correlation skipped: fewer than 2 variables with non-zero variance.")
    diag(cor_mat) <- 1
    if (!is.null(pval_mat)) diag(pval_mat) <- 0
    return(list(cor = cor_mat, p = pval_mat))
  }
  
  # Subset valid columns
  x_use <- x[, nonzero_idx, drop = FALSE]
  v_use <- colnames(x_use)
  
  if (what == "cor") {
    cor_sub <- suppressWarnings(stats::cor(x_use, method = method, use = "pairwise.complete.obs"))
    diag(cor_sub) <- 1
    cor_mat[v_use, v_use] <- cor_sub
    return(list(cor = cor_mat, p = NULL))
  }
  
  # cor_p mode - use cor.test() for each pair
  pairs <- utils::combn(seq_along(v_use), 2, simplify = FALSE)
  
  for (idx in pairs) {
    vi <- v_use[idx[1]]; vj <- v_use[idx[2]]
    xi <- x_use[, idx[1]]; xj <- x_use[, idx[2]]
    
    ok <- stats::complete.cases(xi, xj)
    if (sum(ok) < 2) next
    
    ct <- tryCatch(
      stats::cor.test(xi[ok], xj[ok], method = method),
      error = function(e) NULL
    )
    if (is.null(ct)) next
    
    r <- ct$estimate[[1]]
    p <- ct$p.value
    cor_mat[vi, vj] <- cor_mat[vj, vi] <- r
    pval_mat[vi, vj] <- pval_mat[vj, vi] <- p
  }
  
  # Set diagonal to 1 (correlation) and 0 (p-value) for valid variables
  diag(cor_mat)[v_use] <- 1
  if (!is.null(pval_mat)) diag(pval_mat)[v_use] <- 0
  
  list(cor = cor_mat, p = pval_mat)
}



# =========================
# Thresholds
# =========================

#' @title Threshold correlation matrix by absolute value
#' @description
#' Build an adjacency matrix by applying a minimum absolute correlation cutoff.
#'
#' @param cor_mat Correlation matrix (numeric, symmetric).
#' @param min_cor Minimum absolute correlation to include as an edge. Must be in \[0, 1\].
#'
#' @return Numeric symmetric adjacency matrix. Weights retain the (signed) correlation
#' wherever |r| >= min_cor; zeros elsewhere. Diagonal is set to 0.
#' @export
threshold_absolute <- function(cor_mat, min_cor) {
  if (!is.matrix(cor_mat)) cli::cli_abort("{.arg cor_mat} must be a matrix.")
  if (!is.numeric(cor_mat)) cli::cli_abort("{.arg cor_mat} must be numeric.")
  if (!isSymmetric(cor_mat)) cli::cli_abort("{.arg cor_mat} must be symmetric.")
  if (!is.numeric(min_cor) || length(min_cor) != 1 || min_cor < 0 || min_cor > 1)
    cli::cli_abort("{.arg min_cor} must be a number in [0, 1].")
  
  # Apply threshold
  adj <- cor_mat
  adj[abs(cor_mat) < min_cor | is.na(cor_mat)] <- 0
  diag(adj) <- 0
  adj
}



#' @title Threshold correlation matrix to achieve target edge density
#' @description
#' Keep the largest absolute correlations until the desired undirected edge density is reached.
#'
#' @param cor_mat Correlation matrix (numeric, symmetric).
#' @param density Target density in (0, 1\]. For p nodes, the number of edges is approximately
#'   \code{density * p*(p - 1)/2}.
#'
#' @return Numeric symmetric adjacency matrix with retained correlations as weights.
#' @export
threshold_density <- function(cor_mat, density) {
  if (!is.matrix(cor_mat)) cli::cli_abort("{.arg cor_mat} must be a matrix.")
  if (!is.numeric(cor_mat)) cli::cli_abort("{.arg cor_mat} must be numeric.")
  if (!isSymmetric(cor_mat)) cli::cli_abort("{.arg cor_mat} must be symmetric.")
  if (!is.numeric(density) || length(density) != 1 || density <= 0 || density > 1)
    cli::cli_abort("{.arg density} must be a number in (0, 1].")
  
  p <- ncol(cor_mat)
  total_possible_edges <- p * (p - 1) / 2
  target_edges <- ceiling(density * total_possible_edges)
  
  ut <- upper.tri(cor_mat)
  abs_vals <- abs(cor_mat[ut])
  valid <- !is.na(abs_vals)
  
  if (sum(valid) == 0) {
    adj <- matrix(0, p, p, dimnames = dimnames(cor_mat))
    return(adj)
  }
  
  # Get the threshold value to keep the top target_edges values
  sorted <- sort(abs_vals[valid], decreasing = TRUE)
  threshold <- if (target_edges < length(sorted)) sorted[target_edges] else min(sorted)
  
  # Build adjacency matrix
  mask <- abs(cor_mat) >= threshold & !is.na(cor_mat)
  adj <- cor_mat * mask
  diag(adj) <- 0
  adj
}


#' @title Threshold by p-value significance
#' @description
#' Build a (0/1) adjacency mask from a p-value matrix after multiple testing correction.
#' This function identifies *which* edges are statistically significant;
#' weights can be carried from the original correlation matrix downstream.
#'
#' @param pval_mat Square matrix of p-values (numeric, symmetric).
#' @param alpha Significance level in (0, 1].
#' @param adjust P-value adjustment method: one of "none", "holm", "hochberg", "hommel",
#'   "bonferroni", "BH", "BY", "fdr".
#'
#' @return Numeric symmetric adjacency matrix with values 1 (significant) or 0.
#' @export
threshold_pvalue <- function(pval_mat, alpha = 0.05,
                             adjust = c("none", "holm", "hochberg", "hommel",
                                        "bonferroni", "BH", "BY", "fdr")) {
  adjust <- match.arg(adjust)
  
  if (!is.matrix(pval_mat)) cli::cli_abort("{.arg pval_mat} must be a matrix.")
  if (!is.numeric(pval_mat)) cli::cli_abort("{.arg pval_mat} must be numeric.")
  if (!isSymmetric(pval_mat)) cli::cli_abort("{.arg pval_mat} must be symmetric.")
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha > 1)
    cli::cli_abort("{.arg alpha} must be a single number in (0, 1].")
  
  p <- nrow(pval_mat)
  if (p != ncol(pval_mat)) cli::cli_abort("{.arg pval_mat} must be square.")
  
  if (is.null(colnames(pval_mat))) colnames(pval_mat) <- seq_len(p)
  if (is.null(rownames(pval_mat))) rownames(pval_mat) <- colnames(pval_mat)
  
  ut <- upper.tri(pval_mat)
  pv_ut <- pval_mat[ut]
  
  if (length(pv_ut) == 0 || all(is.na(pv_ut))) {
    adj <- matrix(0, p, p, dimnames = dimnames(pval_mat))
    return(adj)
  }
  
  pv_adj <- stats::p.adjust(pv_ut, method = adjust)
  sig_ut <- as.numeric(pv_adj <= alpha)
  
  adj <- matrix(0, p, p)
  adj[ut] <- sig_ut
  adj <- adj + t(adj)
  diag(adj) <- 0
  dimnames(adj) <- dimnames(pval_mat)
  adj
}


# =========================
# Graph conversion
# =========================

#' @title Convert adjacency matrix to igraph
#' @description
#' Convert a symmetric numeric adjacency matrix into an undirected \code{igraph} object.
#' If all edges are zero (or NA), returns an empty graph with named vertices.
#'
#' @param adj Numeric symmetric adjacency matrix.
#' @param weighted Logical; if \code{TRUE}, non-zero entries are used as \code{weight}.
#'
#' @return An \code{igraph} object.
#' @export
adjacency_to_graph <- function(adj, weighted = TRUE) {
  if (!is.matrix(adj)) cli::cli_abort("{.arg adj} must be a matrix.")
  if (!is.numeric(adj)) cli::cli_abort("{.arg adj} must be numeric.")
  if (nrow(adj) != ncol(adj)) cli::cli_abort("{.arg adj} must be square.")
  if (!isSymmetric(adj)) cli::cli_abort("{.arg adj} must be symmetric.")
  
  # Ensure row/col names
  if (is.null(colnames(adj))) colnames(adj) <- seq_len(ncol(adj))
  if (is.null(rownames(adj))) rownames(adj) <- colnames(adj)
  
  # Handle empty network
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


# =========================
# High-level wrapper
# =========================

#' @title Build a correlation network in one call
#'
#' @description
#' High-level wrapper to construct a correlation-based network by computing pairwise correlations,
#' applying a thresholding strategy, and converting the resulting matrix into a graph.
#' Internally calls \code{\link{correlate_matrix}}, \code{\link{threshold_absolute}},
#' \code{\link{threshold_density}}, or \code{\link{threshold_pvalue}}, and finally \code{\link{adjacency_to_graph}}.
#' Designed for easy integration in workflows or use in functions like \code{set_omic()}.
#'
#' @param x A numeric matrix or data frame with samples in rows and variables (features/taxa) in columns.
#' @param cor_method A character string specifying the correlation method:
#'   one of \code{"pearson"}, \code{"spearman"}, or \code{"kendall"}.
#' @param thresh_method Thresholding method to apply:
#'   one of \code{"absolute"}, \code{"density"}, or \code{"p-value"}.
#' @param thresh_value A numeric threshold value, interpreted depending on \code{thresh_method}:
#'   \itemize{
#'     \item \code{"absolute"}: minimum absolute correlation, \eqn{\in [0, 1]}
#'     \item \code{"density"}: target edge density, \eqn{\in (0, 1]}
#'     \item \code{"p-value"}: significance level (alpha), \eqn{\in (0, 1]}
#'   }
#' @param adjust Character string indicating the p-value adjustment method
#'   (only used when \code{thresh_method = "p-value"}).
#'   One of: \code{"none"}, \code{"holm"}, \code{"hochberg"}, \code{"hommel"},
#'   \code{"bonferroni"}, \code{"BH"}, \code{"BY"}, or \code{"fdr"}.
#' @param output A character string specifying the output type:
#'   one of \code{"graph"}, \code{"adjacency"}, or \code{"both"}.
#'
#' @return
#' Depending on \code{output}:
#' \itemize{
#'   \item If \code{"graph"}: returns an \code{igraph} object.
#'   \item If \code{"adjacency"}: returns a numeric adjacency matrix.
#'   \item If \code{"both"}: returns a list with elements \code{adjacency} and \code{graph}.
#' }
#'
#' @seealso
#' \code{\link{correlate_matrix}}, \code{\link{threshold_absolute}},
#' \code{\link{threshold_density}}, \code{\link{threshold_pvalue}},
#' \code{\link{adjacency_to_graph}}
#'
#' @export
build_corr_net <- function(x,
                           cor_method    = c("pearson", "spearman", "kendall"),
                           thresh_method = c("absolute", "density", "p-value"),
                           thresh_value,
                           adjust        = c("none", "holm", "hochberg", "hommel",
                                             "bonferroni", "BH", "BY", "fdr"),
                           output        = c("graph", "adjacency", "both")) {
  
  cor_method    <- match.arg(cor_method)
  thresh_method <- match.arg(thresh_method)
  output        <- match.arg(output)
  adjust        <- match.arg(adjust)
  
  # Input checks
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.numeric(x)) cli::cli_abort("{.arg x} must be a numeric matrix or data.frame.")
  if (is.null(colnames(x))) cli::cli_abort("{.arg x} must have column names (variables/features).")
  if (missing(thresh_value))
    cli::cli_abort("{.arg thresh_value} is required: absolute [0,1], density (0,1], p-value (0,1].")
  
  # 1) Correlation
  res <- correlate_matrix(x = x, method = cor_method,
                          what = if (thresh_method == "p-value") "cor_p" else "cor")
  cor_mat <- res$cor
  
  # 2) Thresholding
  adj <- switch(thresh_method,
                "absolute" = {
                  if (!is.numeric(thresh_value) || thresh_value < 0 || thresh_value > 1)
                    cli::cli_abort("For {.val absolute}, {.arg thresh_value} must be in [0,1].")
                  threshold_absolute(res$cor, min_cor = thresh_value)
                },
                "density" = {
                  if (!is.numeric(thresh_value) || thresh_value <= 0 || thresh_value > 1)
                    cli::cli_abort("For {.val density}, {.arg thresh_value} must be in (0,1].")
                  threshold_density(res$cor, density = thresh_value)
                },
                "p-value" = {
                  if (!is.numeric(thresh_value) || thresh_value <= 0 || thresh_value > 1)
                    cli::cli_abort("For {.val p-value}, {.arg thresh_value} must be in (0,1].")
                  adj_pval <- threshold_pvalue(res$p, alpha = thresh_value, adjust = adjust)
                  adj_mat <- res$cor * (adj_pval > 0)
                  diag(adj_mat) <- 0
                  adj_mat
                }
  )
  
  # 3) Output
  graph <- adjacency_to_graph(adj, weighted = TRUE)
  switch(output,
         "graph" = graph,
         "adjacency" = adj,
         "both" = list(adjacency = adj, graph = graph))
}



