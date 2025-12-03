#' Build a protein contact network from coordinates
#'
#' @param coords data.frame with columns: resno, resname, x, y, z
#' @param cutoff numeric distance threshold (Å) to define a contact
#' @param directed logical; whether to return a directed graph (default FALSE)
#' @return list with elements: graph (igraph), adjacency (matrix), coords (data.frame)
#' @export
build_pcn <- function(coords, cutoff = 7, directed = FALSE) {
  stopifnot(is.data.frame(coords), all(c("resno", "resname", "x", "y", "z") %in% names(coords)))
  stopifnot(is.numeric(cutoff), length(cutoff) == 1L, cutoff > 0)

  if (nrow(coords) < 2) stop("Need >= 2 residues to form a network.")

  d <- as.matrix(stats::dist(coords[, c("x", "y", "z")], method = "euclidean"))
  diag(d) <- Inf

  adj <- (d <= cutoff) * 1L
  diag(adj) <- 0L
  rownames(adj) <- colnames(adj) <- paste0(coords$resname, coords$resno)

  g <- igraph::graph_from_adjacency_matrix(adj, mode = if (directed) "directed" else "undirected")
  igraph::V(g)$resno <- coords$resno
  igraph::V(g)$resname <- coords$resname

  list(graph = g, adjacency = adj, coords = coords)
}

#' Extract residue coordinates from a bio3d pdb object
#'
#' @param pdb a `bio3d` pdb object (from `bio3d::read.pdb`)
#' @param chain single-letter chain identifier (e.g., "A")
#' @param atom atom name to select (e.g., "CA", "CB")
#' @return data.frame with resno, resname, x, y, z
#' @export
coords_from_pdb <- function(pdb, chain = "A", atom = "CA") {
  stopifnot(!is.null(pdb$atom))
  at <- as.data.frame(pdb$atom)

  needed <- c("elety", "chain", "x", "y", "z", "resno", "resid")
  if (!all(needed %in% names(at))) {
    stop("Unexpected PDB atom table format. Columns missing: ", paste(setdiff(needed, names(at)), collapse = ", "))
  }

  at <- at[at$chain == chain & at$elety == atom, , drop = FALSE]
  if (!nrow(at)) stop("No atoms matched chain=", chain, " and atom=", atom)

  at <- at[!duplicated(at$resno), ]

  data.frame(
    resno = at$resno,
    resname = at$resid,
    x = at$x,
    y = at$y,
    z = at$z,
    stringsAsFactors = FALSE
  )
}

#' Build a PCN from a PDB ID (download via bio3d)
#'
#' @param pdb_id four-character PDB identifier
#' @param chain chain identifier
#' @param atom atom name (e.g., "CA")
#' @param cutoff numeric Å cutoff
#' @param path directory to store downloaded PDB
#' @return list as in `build_pcn()`
#' @export
pcn_from_pdb <- function(pdb_id, chain = "A", atom = "CA", cutoff = 7, path = tempdir()) {
  stopifnot(is.character(pdb_id), nchar(pdb_id) >= 4)
  fn <- bio3d::get.pdb(pdb_id, path = path, pdbext = "pdb", overwrite = FALSE)
  pdb <- bio3d::read.pdb(fn)
  coords <- coords_from_pdb(pdb, chain = chain, atom = atom)
  build_pcn(coords, cutoff = cutoff, directed = FALSE)
}

#' Build a PCN from a local PDB file
#' @param file path to a .pdb file
#' @inheritParams pcn_from_pdb
#' @export
pcn_from_file <- function(file, chain = "A", atom = "CA", cutoff = 7) {
  stopifnot(file.exists(file))
  pdb <- bio3d::read.pdb(file)
  coords <- coords_from_pdb(pdb, chain = chain, atom = atom)
  build_pcn(coords, cutoff = cutoff, directed = FALSE)
}

#' Get adjacency matrix from a PCN object
#' @param pcn object returned by `build_pcn()`
#' @export
pcn_to_adjacency <- function(pcn) {
  stopifnot(is.list(pcn), !is.null(pcn$adjacency))
  pcn$adjacency
}

#' Get edge list (data.frame) from a PCN object
#' @param pcn object returned by `build_pcn()`
#' @return data.frame with columns: from, to (undirected)
#' @export
pcn_edges <- function(pcn) {
  g <- pcn$graph
  ed <- igraph::as_edgelist(g, names = TRUE)
  data.frame(from = ed[,1], to = ed[,2], stringsAsFactors = FALSE)
}
