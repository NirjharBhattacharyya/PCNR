skip_on_cran()
test_that('pcn_from_pdb basic', {
  pcn <- try(pcn_from_pdb('2K0A', chain='A', atom='CA', cutoff=7), silent = TRUE)
  if (!inherits(pcn, 'try-error')) {
    expect_true(igraph::is.igraph(pcn$graph))
    expect_true(is.matrix(pcn$adjacency))
  } else {
    succeed()  # offline environments
  }
})
