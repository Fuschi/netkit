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
#' At runtime, the wrapper first looks in the installed package
#' \code{exec/} directory and then falls back to the repository-local
#' \code{exec/} directory when the function is sourced from the project.
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
#' On Unix-like systems the wrapper attempts to make the selected executable
#' runnable if execute permissions are missing. These files are not compiled by R;
#' they must be pre-built for each OS.
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

  bin_path <- .resolve_gomez_executable()
  tmp_base <- tempfile("signed_mod_")
  file_graph <- paste0(tmp_base, ".net")
  file_result <- paste0(tmp_base, "_res.txt")

  on.exit(unlink(c(file_graph, file_result), force = TRUE), add = TRUE)
  igraph::write_graph(g, file = file_graph, format = "pajek")

  cmd_args <- c("none", "WS", "l", "1", resistance, penalty, file_graph, file_result)
  exec_log <- tempfile("signed_mod_log_")
  on.exit(unlink(exec_log, force = TRUE), add = TRUE)
  status <- system2(bin_path, args = as.character(cmd_args), stdout = exec_log, stderr = exec_log)
  exec_output <- readLines(exec_log, warn = FALSE)

  if (!identical(status, 0L)) {
    cli::cli_abort(c(
      "Gomez executable failed with exit status {.val {status}}.",
      if (length(exec_output)) "i" = paste(exec_output, collapse = "\n")
    ))
  }
  if (!file.exists(file_result)) {
    cli::cli_abort(c(
      "No result file was produced by the Gomez executable.",
      if (length(exec_output)) "i" = paste(exec_output, collapse = "\n")
    ))
  }

  parsed <- .parse_gomez_output(readLines(file_result, warn = FALSE), igraph::vcount(g))
  res <- igraph::make_clusters(
    graph = g,
    membership = parsed$membership,
    algorithm  = "signed modularity louvain",
    modularity = parsed$modularity
  )
  if (add.names) names(res$membership) <- igraph::V(g)$name
  res
}

#' @keywords internal
#' @noRd
.resolve_gomez_executable <- function() {
  os <- Sys.info()[["sysname"]]
  exec_name <- switch(
    os,
    "Linux" = "gomez_modularity_opt_Linux.exe",
    "Windows" = "gomez_modularity_opt_Windows.exe",
    "Darwin" = "gomez_modularity_opt_Mac.exe",
    cli::cli_abort("Unsupported operating system: {.val {os}}")
  )

  pkg_exec_dir <- system.file("exec", package = "netkit")
  local_exec_dir <- file.path(getwd(), "exec")
  candidates <- unique(normalizePath(
    c(
      if (nzchar(pkg_exec_dir)) file.path(pkg_exec_dir, exec_name),
      file.path(local_exec_dir, exec_name)
    ),
    winslash = "/",
    mustWork = FALSE
  ))
  bin_path <- candidates[file.exists(candidates)][1]

  if (is.na(bin_path)) {
    cli::cli_abort(c(
      "Could not locate the Gomez executable for {.val {os}}.",
      "i" = paste("Checked:", paste(candidates, collapse = ", "))
    ))
  }
  if (.Platform$OS.type != "windows" && file.access(bin_path, mode = 1) != 0) {
    Sys.chmod(bin_path, mode = "0755")
  }
  if (.Platform$OS.type != "windows" && file.access(bin_path, mode = 1) != 0) {
    cli::cli_abort(
      "Executable exists but is not runnable: {.file {bin_path}}. Try {.code chmod +x {bin_path}}."
    )
  }

  bin_path
}

#' @keywords internal
#' @noRd
.parse_gomez_output <- function(lines, n_vertices) {
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  modularity_line <- grep("^Q\\s*=", lines, value = TRUE)
  list_line <- grep("^Number of lists\\s*:", lines, value = TRUE)
  community_lines <- grep("^[0-9]+:\\s+", lines, value = TRUE)

  if (!length(modularity_line) || !length(list_line) || !length(community_lines)) {
    cli::cli_abort("Output from the Gomez executable is incomplete or could not be parsed.")
  }

  modularity <- suppressWarnings(as.numeric(sub("^Q\\s*=\\s*", "", modularity_line[1])))
  comm_num <- suppressWarnings(as.integer(sub("^Number of lists\\s*:\\s*", "", list_line[1])))
  if (is.na(modularity) || is.na(comm_num)) {
    cli::cli_abort("Failed to extract modularity or the number of communities from the Gomez output.")
  }
  if (length(community_lines) < comm_num) {
    cli::cli_abort("The Gomez output reported fewer community lines than expected.")
  }

  membership <- rep(NA_integer_, n_vertices)
  for (i in seq_len(comm_num)) {
    parts <- strsplit(community_lines[i], ":", fixed = TRUE)[[1]]
    members <- scan(text = parts[2], what = integer(), quiet = TRUE)
    if (!length(members)) cli::cli_abort("Community {.val {i}} has no parsed members.")
    if (any(members < 1L | members > n_vertices)) {
      cli::cli_abort("Community {.val {i}} contains invalid vertex indices.")
    }
    membership[members] <- i
  }
  if (anyNA(membership)) {
    cli::cli_abort("Some vertices were not assigned to any community by the Gomez executable.")
  }

  list(membership = membership, modularity = modularity, communities = comm_num)
}
