# PCNR

PCNR builds protein contact networks (PCNs) from PDB structures, performs network analysis, and exports graphs.

## Install from this folder

```r
# setwd('path/to/PCNR')  # the folder with DESCRIPTION
install.packages(c("devtools","roxygen2","bio3d","igraph","testthat","knitr","rmarkdown"))
devtools::document()
devtools::install(build_vignettes = FALSE, dependencies = FALSE, upgrade = "never")
```

## Quick start

```r
library(PCNR)
pcn <- pcn_from_pdb("2K0A", chain="A", atom="CA", cutoff=7)
plot_pcn(pcn$graph)
head(pcn_metrics(pcn))
pcn_summary(pcn)
```
