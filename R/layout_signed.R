#' Compute Graph Layout Using Only Positive Edges
#'
#' @description
#' Generates a graph layout by considering only the positively weighted edges of an input graph.
#' This is useful when visualizing signed graphs where positive interactions (e.g., cooperation, similarity)
#' are of primary interest. The default layout algorithm is Fruchterman-Reingold, but any `igraph` layout
#' function can be passed.
#'
#' @param graph An `igraph` object representing a graph (weighted or unweighted).
#' @param layout A layout specification (default: `with_fr()`). Any layout function from `igraph`
#'        that can be passed to `layout_()` is accepted.
#'
#' @details
#' If the graph is weighted, the function retains only edges with strictly positive weights.
#' Edges with zero or negative weights are ignored when computing the layout. The returned layout
#' matrix can be used in standard igraph plotting functions.
#'
#' @return A numeric matrix of layout coordinates (2 columns for 2D layouts; 3 columns for 3D).
#'
#' @examples
#' library(igraph)
#' g <- make_ring(10)
#' E(g)$weight <- c(1, -1, 1, 1, -1, 1, 1, -1, 1, 1)
#'
#' coords <- layout_signed(g)
#' plot(g, layout = coords)
#'
#' @export
layout_signed <- function(graph, layout = igraph::with_fr()) {
  if (!inherits(graph, "igraph")) {
    cli::cli_abort("{.arg graph} must be an igraph object.")
  }
  
  if (igraph::is_weighted(graph)) {
    positive_edges <- which(igraph::E(graph)$weight > 0)
    graph <- igraph::subgraph.edges(graph, eids = positive_edges, delete.vertices = FALSE)
  }
  
  coords <- igraph::layout_(graph, layout)
  return(coords)
}
