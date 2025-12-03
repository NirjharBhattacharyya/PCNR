#' PCN analysis utilities
#' @name pcn-analysis
NULL

#' Compute node-level metrics for a PCN
#' @param pcn A PCN object
#' @export
pcn_metrics <- function(pcn) {
  stopifnot(is.list(pcn), igraph::is.igraph(pcn$graph))
  g <- pcn$graph

  closeness_vec <- tryCatch(
    igraph::closeness(g, normalized = TRUE),
    error = function(e) {
      sp <- igraph::distances(g)
      inv <- ifelse(is.finite(sp) & sp > 0, 1/sp, 0)
      rowSums(inv) / (igraph::vcount(g) - 1)
    }
  )

  data.frame(
    resno = igraph::V(g)$resno,
    resname = igraph::V(g)$resname,
    degree = igraph::degree(g, mode = "all"),
    betweenness = igraph::betweenness(g, directed = FALSE, normalized = TRUE),
    closeness = as.numeric(closeness_vec),
    clustering_coeff = igraph::transitivity(g, type = "local", isolates = "zero"),
    row.names = igraph::V(g)$name,
    stringsAsFactors = FALSE
  )
}

#' Global summary metrics for a PCN
#' @param pcn A PCN object
#' @export
pcn_summary <- function(pcn) {
  stopifnot(is.list(pcn), igraph::is.igraph(pcn$graph))
  g <- pcn$graph

  comps <- igraph::components(g)
  giant <- igraph::induced_subgraph(g, which(comps$membership == which.max(comps$csize)))

  avg_path_length <- NA_real_
  if (igraph::is_connected(giant)) {
    avg_path_length <- igraph::mean_distance(giant, directed = FALSE, unconnected = FALSE)
  }

  list(
    num_nodes = igraph::vcount(g),
    num_edges = igraph::ecount(g),
    density = igraph::edge_density(g, loops = FALSE),
    num_components = comps$no,
    largest_component_size = max(comps$csize),
    average_clustering = igraph::transitivity(g, type = "average"),
    avg_path_length_giant = avg_path_length
  )
}

#' Detect residue communities (modules) in a PCN
#' @param pcn A PCN object
#' @param method One of "louvain","fast_greedy","walktrap".
#' @export
pcn_communities <- function(pcn, method = c("louvain","fast_greedy","walktrap")) {
  stopifnot(is.list(pcn), igraph::is.igraph(pcn$graph))
  method <- match.arg(method)
  g <- pcn$graph

  comm <- switch(method,
    louvain     = igraph::cluster_louvain(g),
    fast_greedy = igraph::cluster_fast_greedy(g),
    walktrap    = igraph::cluster_walktrap(g)
  )

  igraph::V(g)$community <- comm$membership
  pcn$graph <- g
  attr(comm, "pcn") <- pcn
  comm
}

#' Compare two PCNs (edge-level)
#' @param pcn1 First PCN
#' @param pcn2 Second PCN
#' @export
pcn_compare <- function(pcn1, pcn2) {
  stopifnot(is.list(pcn1), igraph::is.igraph(pcn1$graph),
            is.list(pcn2), igraph::is.igraph(pcn2$graph))
  g1 <- pcn1$graph; g2 <- pcn2$graph

  common_names <- intersect(igraph::V(g1)$name, igraph::V(g2)$name)
  g1c <- igraph::induced_subgraph(g1, vids = common_names)
  g2c <- igraph::induced_subgraph(g2, vids = common_names)

  ekey <- function(g) {
    if (igraph::ecount(g) == 0) return(character(0))
    ed <- igraph::as_edgelist(g, names = TRUE)
    apply(ed, 1, function(x) paste(sort(x), collapse = "--"))
  }

  k1 <- unique(ekey(g1c))
  k2 <- unique(ekey(g2c))

  only1 <- setdiff(k1, k2)
  only2 <- setdiff(k2, k1)
  inter <- intersect(k1, k2)
  denom <- length(k1) + length(k2) - length(inter)
  jacc <- if (denom == 0) NA_real_ else length(inter) / denom

  list(
    common_nodes = common_names,
    only_in_pcn1_edges = only1,
    only_in_pcn2_edges = only2,
    jaccard_edge_similarity = jacc
  )
}

#' Export a PCN graph to GraphML
#' @param pcn PCN object
#' @param file Output path ending with .graphml
#' @export
pcn_export_graphml <- function(pcn, file) {
  stopifnot(is.list(pcn), igraph::is.igraph(pcn$graph))
  igraph::write_graph(pcn$graph, file = file, format = "graphml")
  invisible(file)
}
