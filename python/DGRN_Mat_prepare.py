import sys

import numpy as np
import pandas as pd
from scipy import sparse
import scanpy as sc
import os

# scMultiome_mat_file_path = "data/scMultiome_mat.csv"
# pseudotime_lineages_path = "result/pseudotime_lineages.csv"
# DEG_paths = ["result/DEgenes_MAST_sp4_0to1.csv", "result/DEgenes_MAST_sp4_0to2.csv"]
# HVGenes_Peaks_path = "result/HVGenes_Peaks.csv"
# adata_out = "result/adata_scMultiome_SEACells.h5ad"

scMultiome_mat_file_path = sys.argv[1]
pseudotime_lineages_path = sys.argv[2]
DEG_paths = [sys.argv[3], sys.argv[4]]
HVGenes_Peaks_path = sys.argv[5]
adata_out = sys.argv[6]


#HVGenes_Peaks_path
HVGenes_Peaks = pd.read_csv(HVGenes_Peaks_path)['x']

#DEGenes_Peaks
DEGenes_Peaks = pd.Series()
for path_ in DEG_paths:
    DE_ = pd.read_csv(path_)['primerid']
    DEGenes_Peaks = pd.concat([DEGenes_Peaks, DE_], ignore_index=True)

HVGenes_Peaks = HVGenes_Peaks[HVGenes_Peaks.isin(DEGenes_Peaks)]

#anndata
scMultiome_mat = pd.read_csv(scMultiome_mat_file_path)
scMultiome_mat = scMultiome_mat.set_index("cellid")
scMultiome_mat = scMultiome_mat.loc[:, scMultiome_mat.columns.isin(HVGenes_Peaks)]
adata = sc.AnnData(scMultiome_mat)

#pseudotime_lineages
pseudotime_lineages = pd.read_csv(pseudotime_lineages_path,
                                  index_col=0)

pseudotime_lineages.columns = [f"{col_}_lineage" for col_ in pseudotime_lineages.columns]
adata.obs = pseudotime_lineages.loc[adata.obs_names, ]

#DEG
adata.uns['lineages'] = pseudotime_lineages.columns.tolist()
for lineage_,path_ in zip(adata.uns['lineages'], DEG_paths):
    DEG_ = pd.read_csv(path_, index_col='primerid')
    DEG_ = DEG_.loc[DEG_.index.isin(adata.var_names), :]
    adata.var.loc[DEG_.index, lineage_+'_logFC'] = DEG_['logFC']

# adata.X = np.log2(adata.X+1)

adata.X = adata.X.astype(np.float32)
adata.write_h5ad(adata_out)



