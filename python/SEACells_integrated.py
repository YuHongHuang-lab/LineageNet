import numpy as np
import pandas as pd
import scanpy as sc
import sys

from SEACells.core import summarize_by_SEACell

rna_ad_path = sys.argv[1]
atac_ad_path = sys.argv[2]
tmp = sys.argv[3]
result_path = sys.argv[4]

def log_transform(ad, ps=0.1):
    ad.X.data = np.log2(ad.X.data + ps) - np.log2(ps)

rna_ad = sc.read_h5ad(rna_ad_path)
atac_ad = sc.read_h5ad(atac_ad_path)

#RNA
rna_meta_ad = summarize_by_SEACell(rna_ad, SEACells_label='SEACell',summarize_layer='X')
rna_meta_ad.obs['celltype'] = rna_ad.obs.groupby('SEACell').apply(lambda x: pd.Series(x['celltype']).mode())

# sc.pp.normalize_per_cell(rna_meta_ad, counts_per_cell_after=10000)
# log_transform(rna_meta_ad)

rna_mate_mat = pd.DataFrame(rna_meta_ad.X.toarray(),
                            columns=rna_meta_ad.var_names,
                            index=rna_meta_ad.obs_names)
rna_mate_mat['cellid'] = rna_mate_mat.index.to_series()

rna_mate_mat.to_csv(result_path+"rna_mate_mat.csv")


#ATAC
atac_ad.obs['SEACell'] = rna_ad.obs.loc[atac_ad.obs_names, 'SEACell']
atac_meta_ad = summarize_by_SEACell(atac_ad, SEACells_label='SEACell',summarize_layer='X')
atac_meta_ad.X = atac_meta_ad.X.toarray()
# sc.pp.normalize_per_cell(rna_meta_ad, counts_per_cell_after=10000)
# log_transform(rna_meta_ad)


atac_mate_mat = pd.DataFrame(atac_meta_ad.X,
                             columns=atac_meta_ad.var_names,
                             index=atac_meta_ad.obs_names)
atac_mate_mat['cellid'] = atac_mate_mat.index.to_series()

atac_mate_mat.to_csv(result_path+"atac_mate_mat.csv")


scMultiome_mat = pd.merge(atac_mate_mat,
                       rna_mate_mat,
                       how="inner",
                       on="cellid")

scMultiome_mat = scMultiome_mat.set_index("cellid")
scMultiome_mat.to_csv(result_path+"scMultiome_mat.csv")






