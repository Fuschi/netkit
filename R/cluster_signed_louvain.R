#' Community Detection in Signed Weighted Graphs (Gómez Modularity)
#'
#' @description
#' Detects communities in signed weighted graphs by optimizing the signed modularity
#' proposed by Gómez, Jensen, and Arenas (2009). This function wraps an external
#' executable that implements the signed modularity Louvain-like optimization,
#' capable of handling both positive and negative edge weights in undirected networks.
#'
#' The signed modularity is defined as:
#' \deqn{Q = Q^{+} - Q^{-}}
#' where \eqn{Q^{+}} represents the modularity for positive edges and
#' \eqn{Q^{-}} for negative edges (treated as repulsive interactions).
#'
#' @details
#' The algorithm identifies dense subgraphs by maximizing the difference between the
#' modularity of positive and negative edge weights. It is particularly suitable for
#' correlation or association networks where links can represent both cooperative and
#' antagonistic relationships.
#'
#' Internally, the function:
#' \enumerate{
#'   \item Writes the graph to a temporary file in Pajek (\code{.net}) format.
#'   \item Invokes the appropriate executable for the host operating system.
#'   \item Parses the output to build a standard \code{communities} object from \pkg{igraph}.
#' }
#'
#' The required executables must be present in the \code{exec/} directory of the package:
#' \itemize{
#'   \item \code{gomez_modularity_opt_Linux.exe} — Linux version
#'   \item \code{gomez_modularity_opt_Mac.exe} — macOS version
#'   \item \code{gomez_modularity_opt_Windows.exe} — Windows version
#' }
#'
#' At runtime, they are automatically located using:
#' \code{system.file("exec", package = "netkit")}
#'
#' Each executable must support the following command-line interface:
#' \preformatted{
#'   <executable> none WS l 1 <resistance> <penalty> <input_file> <output_file>
#' }
#'
#' @param g An \code{igraph} object representing a weighted and undirected graph.
#'          Edge weights can be positive or negative.
#' @param resistance Numeric; simulates a uniform self-loop resistance applied to all nodes.
#'        Default is \code{0}, corresponding to no resistance.
#' @param penalty Non-negative numeric; controls the weight of the null model term.
#'        Higher values tend to produce larger communities.
#' @param add.names Logical; if \code{TRUE} and vertex names exist, they are
#'        assigned to the membership vector. Default is \code{FALSE}.
#'
#' @return
#' An object of class \code{communities} (from \pkg{igraph}), with components:
#' \itemize{
#'   \item \code{membership}: a vector assigning each vertex to a community.
#'   \item \code{modularity}: the optimized signed modularity score \eqn{Q = Q^{+} - Q^{-}}.
#'   \item \code{algorithm}: the label \code{"signed modularity louvain"}.
#' }
#'
#' @section Executables:
#' The external binaries are distributed under the package’s \code{exec/} directory.
#' Ensure they are executable (use \code{chmod +x} on Unix-like systems).
#' These files are not compiled by R; they must be pre-built for each OS.
#'
#' @note
#' Temporary files (\code{.net} and \code{_res.txt}) are written to \code{tempdir()}
#' and automatically deleted after execution.
#'
#' @references
#' Gómez, S., Jensen, P., & Arenas, A. (2009). Analysis of community structure in networks
#' of correlated data. \emph{Physical Review E}, 80(1), 016114.
#' \doi{10.1103/PhysRevE.80.016114}
#'
#' @seealso
#' \code{\link[igraph]{communities}}, \code{\link[igraph]{make_clusters}},
#' \code{\link[igraph]{cluster_spinglass}}
#'
#' @examples
#' library(igraph)
#' g <- make_ring(10)
#' E(g)$weight <- rnorm(ecount(g), mean = 0, sd = 1)
#'
#' # Run community detection
#' comm <- cluster_signed_louvain(g)
#'
#' # Inspect results
#' membership(comm)
#' modularity(comm)
#' plot(comm, g)
#'
#' @export
cluster_signed_louvain <- function(g, resistance = 0, penalty = 1, add.names = FALSE) {
  # ---- Checks
  if (!inherits(g, "igraph")) cli::cli_abort("{.arg g} must be an igraph object.")
  if (!igraph::is_weighted(g) || igraph::is_directed(g)) {
    cli::cli_abort("Graph must be {.strong weighted} and {.strong undirected}.")
  }
  if (!is.numeric(resistance)) cli::cli_abort("{.arg resistance} must be numeric.")
  if (!is.numeric(penalty) || penalty < 0) {
    cli::cli_abort("{.arg penalty} must be a numeric value >= 0.")
  }
  if (add.names && !igraph::is_named(g)) {
    cli::cli_abort("Cannot use {.arg add.names = TRUE} if graph has no vertex names.")
  }
  
  # ---- OS & Executable
  os <- Sys.info()[["sysname"]]
  exec_name <- switch(
    os,
    "Linux"   = "gomez_modularity_opt_Linux.exe",
    "Windows" = "gomez_modularity_opt_Windows.exe",
    "Darwin"  = "gomez_modularity_opt_Mac.exe",
    cli::cli_abort("Unsupported operating system: {.val {os}}")
  )
  bin_path <- file.path(system.file("exec", package = "netkit", mustWork = TRUE), exec_name)
  
  if (!file.exists(bin_path)) {
    cli::cli_abort("Executable not found at {.file {bin_path}}.")
  }
  
  # ---- Write graph in Pajek format
  tmp_base <- tempfile("signed_mod_")
  file_graph <- paste0(tmp_base, ".net")
  file_result <- paste0(tmp_base, "_res.txt")
  
  igraph::write_graph(g, file = file_graph, format = "pajek")
  
  # ---- Execute command
  cmd <- sprintf('"%s" none WS l 1 %s %s %s %s',
                 bin_path, resistance, penalty, file_graph, file_result)
  
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  # ---- Parse result
  if (!file.exists(file_result)) {
    cli::cli_abort("No result file produced. External executable may have failed.")
  }
  
  lines <- readLines(file_result, warn = FALSE)
  if (length(lines) < 8) {
    cli::cli_abort("Output from executable is incomplete or corrupted.")
  }
  
  modularity <- suppressWarnings(as.numeric(strsplit(lines[3], " ")[[1]][3]))
  comm_num   <- suppressWarnings(as.integer(strsplit(lines[6], " ")[[1]][4]))
  
  if (is.na(modularity) || is.na(comm_num)) {
    cli::cli_abort("Failed to extract modularity or number of communities.")
  }
  
  comm <- rep(NA_integer_, igraph::vcount(g))
  for (i in seq_len(comm_num)) {
    members <- suppressWarnings(as.integer(strsplit(lines[7 + i], " ")[[1]][-1]))
    if (any(is.na(members))) {
      cli::cli_abort("Invalid node indices in community {.val {i}}.")
    }
    comm[members] <- i
  }
  
  unlink(c(file_graph, file_result), force = TRUE)
  
  res <- igraph::make_clusters(
    graph = g,
    membership = comm,
    algorithm  = "signed modularity louvain",
    modularity = modularity
  )
  if (add.names) names(res$membership) <- igraph::V(g)$name
  return(res)
}
