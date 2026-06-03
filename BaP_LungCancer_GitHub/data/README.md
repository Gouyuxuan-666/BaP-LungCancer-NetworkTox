# Data Sources

All gene expression data are downloaded directly from NCBI GEO at runtime.

## GEO Datasets

| Accession | Platform | Samples | Citation |
|-----------|----------|---------|----------|
| GSE10072 | Affymetrix HG-U133 Plus 2.0 | Tumor + paired normal | — |
| GSE115002 | Affymetrix HG-U133 Plus 2.0 | Lung cancer | — |
| GSE19804 | Affymetrix HG-U133 Plus 2.0 | Tumor + normal | — |
| GSE32863 | Affymetrix HG-U133 Plus 2.0 | Lung cancer | — |
| GSE43458 | Affymetrix HG-U133 Plus 2.0 | Lung cancer | — |
| GSE68465 | Affymetrix HG-U133A | Lung cancer | — |

## BaP Target Databases

| Database | URL | Method |
|----------|-----|--------|
| ChEMBL | https://www.ebi.ac.uk/chembl/ | Bioactivity data |
| SEA | https://sea.bkslab.org/ | Chemical similarity |
| SwissTargetPrediction | http://www.swisstargetprediction.ch/ | Shape/electrostatic similarity |

## Protein Structures (Molecular Docking)

- UniProt: https://www.uniprot.org/
- RCSB PDB: https://www.rcsb.org/

## Download Instructions

All GEO datasets can be downloaded programmatically using the R package `GEOquery`:

```r
library(GEOquery)
gse_ids <- c("GSE10072", "GSE115002", "GSE19804", "GSE32863", "GSE43458", "GSE68465")
for (id in gse_ids) {
  getGEO(id, destdir = "data/")
}
```

Alternatively, download manually from https://www.ncbi.nlm.nih.gov/geo/ and place the series matrix files in this directory.
