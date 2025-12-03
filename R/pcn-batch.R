#' Batch processing of PDB jobs from a TSV
#'
#' @description
#' Read a tab-delimited file with columns: `PDB`, `CHAIN`, `ATOM`, `CUTOFF` and
#' run PCN construction for each row. Returns a list of results.
#'
#' @param tsv_file path to a tab-separated file
#' @param download_path directory for downloaded PDBs (default tempdir())
#' @param progress logical; print progress messages
#' @return named list of PCN objects
#' @export
pcn_batch <- function(tsv_file, download_path = tempdir(), progress = TRUE) {
  stopifnot(file.exists(tsv_file))
  tab <- utils::read.table(tsv_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  req <- c("PDB", "CHAIN", "ATOM", "CUTOFF")
  if (!all(req %in% names(tab))) stop("Batch file must contain columns: ", paste(req, collapse = ", "))

  out <- vector("list", nrow(tab))
  names(out) <- paste0(tab$PDB, "_", tab$CHAIN, "_", tab$ATOM, "_", tab$CUTOFF)

  for (i in seq_len(nrow(tab))) {
    pdb_id <- tab$PDB[i]
    chain  <- tab$CHAIN[i]
    atom   <- tab$ATOM[i]
    cutoff <- as.numeric(tab$CUTOFF[i])

    if (progress) message(sprintf("[%d/%d] %s chain %s atom %s cutoff %.2f", i, nrow(tab), pdb_id, chain, atom, cutoff))
    out[[i]] <- tryCatch(
      pcn_from_pdb(pdb_id = pdb_id, chain = chain, atom = atom, cutoff = cutoff, path = download_path),
      error = function(e) { warning("Failed for ", pdb_id, ": ", e$message); NULL }
    )
  }
  out
}
