#' PCNR: Protein Contact Networks in R
#'
#' Build residue-level protein contact networks (PCNs) from PDB IDs/files,
#' compute adjacency, analyze centralities/communities, and export graphs.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom stats dist
#' @importFrom graphics plot
#' @importFrom utils write.csv read.table
#' @importFrom igraph graph_from_adjacency_matrix V as_edgelist
#' @importFrom igraph degree betweenness closeness transitivity components
#' @importFrom igraph induced_subgraph distances write_graph edge_density
#' @importFrom igraph layout_with_kk layout_with_fr layout_in_circle
#' @importFrom igraph cluster_louvain cluster_fast_greedy cluster_walktrap
NULL
