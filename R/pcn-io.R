#' Write adjacency matrix to CSV
#' @param pcn PCN object
#' @param file output path
#' @export
write_pcn_adjacency <- function(pcn, file) {
  utils::write.csv(pcn$adjacency, file = file, row.names = TRUE)
}

#' Write edge list to CSV
#' @param pcn PCN object
#' @param file output path
#' @export
write_pcn_edges <- function(pcn, file) {
  utils::write.csv(pcn_edges(pcn), file = file, row.names = FALSE)
}
