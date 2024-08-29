LineageNet
================

# The LineageNet

To construct a lineage from cranial neural crest cells to non-osteogenic
specialized cells, we used the jaw masticatory muscles-tendon-bone as a
model and performed snMultiome and spatial transcriptomics sequencing
over consecutive days. Subsequently, we employed single-cell integration
analysis frameworks such as Seurat, Monocle, and Scanpy, along with
Huawei’s deep learning platform, to complete the spatiotemporal
multi-omics modeling.

We have encapsulated the above methods into an R package called
LineageNet. LineageNet includes the following three components:

1.  [Data quality control, dimensionality reduction, clustering, and
    cell
    annotation](./vignettes/Data-preprocessing-and-cell-annotation.Rmd)

2.  [Trajectory construction](./vignettes/Trajectory-construction.Rmd)

3.  [Dynamic network analysis of developmental
    trajectories](./vignettes/Dynamic-network.Rmd).

# Installation

You can install the development version of LineageNet like so: install
from github:

devtools::install_github(“YuHongHuang-lab/LineageNet”)
