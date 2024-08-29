import sys

import scanpy as sc
import SEACells.genescores

rna_ad_path = sys.argv[1]
tmp = sys.argv[2]

rna_ad = sc.read_h5ad(rna_ad_path)

rna_ad.var_names_make_unique()
sc.pl.highest_expr_genes(rna_ad, n_top=20)
sc.pp.normalize_total(rna_ad, target_sum=1e4)
sc.pp.log1p(rna_ad)
sc.pp.highly_variable_genes(rna_ad, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.neighbors(rna_ad, n_neighbors=10, n_pcs=30)
sc.tl.leiden(
    rna_ad,
    resolution=0.9,
    random_state=0,
    n_iterations=2,
    directed=False,
)

#SEACells
#RNA
model = SEACells.core.SEACells(ad=rna_ad,
                               build_kernel_on='X_pca',
                               n_SEACells=100,
                               n_waypoint_eigs=15,
                               convergence_epsilon=1e-5)

model.construct_kernel_matrix()
M = model.kernel_matrix

model.initialize_archetypes()

model.fit(min_iter=10, max_iter=50)

print(f'Ran for {len(model.RSS_iters)} iterations')
for _ in range(5):
    model.step()
print(f'Ran for {len(model.RSS_iters)} iterations')

model.plot_convergence(save_as=tmp+"convergence.pdf")
labels,weights = model.get_soft_assignments()

SEACells.plot.plot_2D(rna_ad,
                      key='X_umap',
                      colour_metacells=False,
                      save_as=tmp+"Umap_SCACells.pdf")
SEACells.plot.plot_SEACell_sizes(rna_ad,
                                 bins=5,
                                 save_as=tmp+"SEACell_sizes.pdf")
rna_ad.write_h5ad(rna_ad_path)




