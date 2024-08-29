import cellrank as cr
import scanpy as sc
from cellrank.kernels import RealTimeKernel
from moscot.problems.time import TemporalProblem
from collections import Counter
import sys

data_path = sys.argv[1]
data_out = sys.argv[2]

#E10_Mand_E11_Mand
adata = cr.datasets.read(data_path)
adata.obs["day"] = adata.obs["time"].astype(float).astype("category")
adata.obs["day_numerical"] = adata.obs["day"].astype(float)
#*prepare problem
tp = TemporalProblem(adata)
tp = tp.score_genes_for_marginals(
    gene_set_proliferation="mouse", gene_set_apoptosis="mouse"
)
tp = tp.prepare(time_key="day", joint_attr = "X_pca")
#*solve
tp = tp.solve(epsilon=1e-3, tau_a=0.95, scale_cost="mean", threshold=0.1)

ct_desc = tp.cell_transition(
    1, 2,
    source_groups="celltype",
    target_groups="celltype",
    forward=True,
    key_added="transitions_1_2"
)

ct_desc.to_csv(data_out)







