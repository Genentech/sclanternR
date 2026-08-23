# Introduction

🚧 🛠️ 🏗️ 👷

This repo provides the R analysis pipeline for our forthcoming preprint on lineage tracing with long reads.

Starting from single-cell allele counts, the method constructs a PCA embedding and distance metric between the single cells, that can be used for distance-based tree-building methods such as neighbor joining.

# Installation

```{r}
devtools::install_github("Genentech/sclanternR")
```

# Usage

See the
[vignette](https://jackkamm.github.io/sclanternR-vignette.html)
for an example on how to use the package.

# Overview of main objects and methods



Key functions:

* `CellAlleleCounts_from_dataframe()`: Reads a table produced by
  `sclantern-nf`, and creates a `CellAlleleCounts` object that
  represents the UMI count of each allele at each locus within each
  cell. Internally it's a large sparse matrix with cells as columns
  and locus-alleles as rows. The object has methods for slicing by
  cells or by loci.
* `fitFiniteMixtureModel()`: From a `CellAlleleCounts`, fits a
  multinomial mixture model, where each cluster has multinomial allele
  frequencies at each locus, and the cluster posterior probabilties
  for each cell are computed by EM.
* `estimateMiCutoff()`: After fitting the mixture model, this function
  selects informative loci via a scaled mutual information (selecting
  a cutoff for this metric).
* `fitCellDistancePipeline()`: Combines elements above in a pipeline
  to compute PCA and cell distances, used for downstream phylogeny
  building (e.g. via neighbor-joining). First a mixture model is fit
  on the cells, and informative loci selected. A further dimension
  reduction on the clusters is performed for smoothing (this is
  helpful for clusters with few cells). Then, allele frequencies are
  estimated in each cell by shrinking their empirical frequencies
  towards their expectations under the mixture model. A cell-level PCA
  is then performed and the squared-Euclidean distances returned as
  the phylogenetic distance metric. Intermediate PCA and Mutual
  Information objects are stored for visual inspection, and the
  pipeline can be rerun from intermediate steps using manually
  selected cutoffs.
