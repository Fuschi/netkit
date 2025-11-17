#' Compute a Graph Layout Using Only Positive Edges
#'
#' @description
#' Computes a graph layout by considering only the positively weighted edges
#' of the input graph. This is useful when visualizing signed networks in which
#' positive interactions (e.g., cooperation or similarity) are the primary
#' structure of interest. The default layout algorithm is Fruchterman–Reingold,
#' but any `igraph` layout specification can be used.
#'
#' @param graph An `igraph` or `tbl_graph` object.
#' @param algorithm A layout specification passed to `igraph::layout_()`,
#'   such as `igraph::with_fr()` or `igraph::with_kk()`.
#' @param circular A logical value required for compatibility with `ggraph()`.
#'   It is ignored by this function.
#' @param ... Additional arguments passed by `ggraph()`; ignored.
#'
#' @return A `data.frame` with columns `x` and `y` (and possibly others),
#'   suitable for use as a layout object in `ggraph()`.
#'
#' @export
layout_signed <- function(
    graph,
    algorithm = igraph::with_fr(),
    circular = FALSE,
    ...
) {
  
  if (inherits(graph, "tbl_graph")) {
    ig <- tidygraph::as.igraph(graph)
  } else if (inherits(graph, "igraph")) {
    ig <- graph
  } else {
    cli::cli_abort("{.arg graph} must be an {.cls igraph} or {.cls tbl_graph} object.")
  }
  
  # Keep only positive arc if the network is weighted
  if (igraph::is_weighted(ig)) {
    positive_edges <- which(igraph::E(ig)$weight > 0)
    ig <- igraph::subgraph.edges(ig, eids = positive_edges, delete.vertices = FALSE)
  }
  
  # Get the layout with the choosed algorithm
  coords <- igraph::layout_(ig, algorithm)
  
  # Return a data.frame with x,y columns
  coords <- as.data.frame(coords)
  names(coords)[1:2] <- c("x", "y")
  
  coords
}

