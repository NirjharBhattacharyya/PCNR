#' Quick plot for a PCN igraph
#' @param g igraph object (e.g., pcn$graph)
#' @param layout character: one of 'kk', 'fr', or 'circle'
#' @export
plot_pcn <- function(g, layout = "kk") {
  stopifnot(igraph::is.igraph(g))
  L <- switch(layout,
              kk = igraph::layout_with_kk,
              fr = igraph::layout_with_fr,
              circle = igraph::layout_in_circle,
              igraph::layout_with_kk)
  graphics::plot(g, layout = L(g), vertex.size = 3, vertex.label = NA, edge.arrow.mode = 0)
}
