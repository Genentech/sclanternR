# Introduction

🚧 🛠️ 🏗️ 👷

This repo provides the R analysis pipeline for our [preprint](https://www.biorxiv.org/content/10.64898/2026.08.25.747112v1) on lineage tracing with long reads.

Starting from single-cell allele counts, the method constructs a PCA embedding and distance metric between the single cells, that can be used for distance-based tree-building methods such as neighbor joining.

# Installation

```{r}
devtools::install_github("Genentech/sclanternR")
```

# Usage

See the
[vignette](https://jackkamm.github.io/sclanternR-vignette.html)
for an example on how to use the package.
