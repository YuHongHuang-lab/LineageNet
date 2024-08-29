.libPaths('~/R/x86_64-pc-linux-gnu-library/4.3/')
reticulate::use_condaenv("/mnt/data/home/tycloud/software/anaconda3/envs/scFates/")
Sys.setenv(RETICULATE_PYTHON = "/mnt/data/home/tycloud/software/anaconda3/envs/scFates/bin/python3")
library(ggplot2)
library(devtools)
library(tidyverse)
library(pheatmap)
my_pal_science6 = c("#A9CA70","#C5D6F0","#F18C54",
                    "#9B59B6","#F7CAC9","#76D7C4")


# usethis::use_package("Seurat", type = "depends")
# usethis::use_package("dplyr", type = "depends")
# usethis::use_package("ggtree", type = "depends")
# usethis::use_package("monocle3", type = "depends")
# usethis::use_package(c("tidyfst", "tidyr",
#                        "viridis", "rlang",
#                        "reticulate","ggpubr",
#                        "ggplot2", "aplot",
#                        "ape", "Matrix"), type = "depends")
# usethis::use_package("tidyfst", type = "depends")
# usethis::use_package("tidyr", type = "depends")
# usethis::use_package("viridis", type = "depends")
# usethis::use_package("rlang", type = "depends")
# usethis::use_package("reticulate", type = "depends")
# usethis::use_package("ggpubr", type = "depends")
# usethis::use_package("ggplot2", type = "depends")
# usethis::use_package("aplot", type = "depends")
# usethis::use_package("ape", type = "depends")
# usethis::use_package("Matrix", type = "depends")


devtools::load_all()
seob = readRDS("/mnt/data/home/tycloud/workspace/scATAC/data/E12_E16_WT_Mand_RPCA_scRNA.rds")
# bad_genes = readRDS(file = "/mnt/data/home/tycloud/workspace/scATAC/data/bad_genes.rds")
# usethis::use_data(bad_genes)
# usethis::use_data_raw("bad_genes")
# usethis::use_data("data/bad_genes.rda")
data("bad_genes")

seob = RunSeurat(seob = seob,
                 split_by="sample",
                 integrated="RPCA",
                 outfeatures = bad_genes,
                 coln_add = "all",
                 cellid=NULL,
                 nfeatures=3000,
                 min_cells = 3,
                 min_features = 200,
                 npcs=20,
                 resolution=0.15,
                 k_anchor=3,
                 k_weight = 200,
                 use_old_RNA = F,
                 use_old_Seurat = F)

DimPlot(seob, reduction = "umap", group.by = "sample",label = T)+
  scale_color_manual(values = colorRampPalette(c("#F7903D", "#4D85BD", "#59A95A"))(5))
DimPlot(seob, reduction = "umap", group.by = "celltype",label = T)+
  scale_color_manual(values = colorRampPalette(c(my_pal_science6))(7))

plotDot(object = seob,
        features=unique(c('Sox9', 'Col2a1',
                          'Pecam1','Cdh5',
                          'Krt14','Epcam',
                          'Plp1','Sox10',
                          'Mitf',
                          'Prrx1', 'Twist1',
                          'Myod1','Myog')),
        mode=1, ID=c('celltype'),
        out_ID = c('celltype'),
        my_pal = colorRampPalette(c(my_pal_science6))(7),
        my_mal_gradientn = c("#4C3399","#F5F2EF","#FFE714"),
        cluster_ID = F,
        cluster_feature = F,
        x_angle=45)
cds = Runmonocl3(seob = seob,
                 cluster="seurat_clusters",
                 ncenter = 500)
cds <- order_cells(cds)
plot_cells(cds_Chai_Mes,
           color_cells_by = "seurat_clusters",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)+
  scale_color_manual(values = colorRampPalette(c(my_pal_science6))(length(unique(seob_Chai@meta.data$seurat_clusters))))


# integrated
seob_E10_Mand_flow = readRDS(file = "/mnt/data/home/tycloud/workspace/scATAC/data/E10_Mand_WT_CNCC_Flow_tmp_scRNA_seq.rds")
seob_E11_Mand_flow = readRDS(file = "/mnt/data/home/tycloud/workspace/scATAC/data/E11_Mand_WT_CNCC_Flow_tmp_scRNA_seq.rds")

seob_E10_Mand_E11_Mand = RunSeurat(seob = list("E10_Mand"=seob_E10_Mand_flow,
                                               "E11_Mand"=seob_E11_Mand_flow),
                                   split_by="sample",
                                   integrated="RPCA",
                                   coln_add = "all",
                                   cellid_add = F,
                                   cellid=NULL,
                                   nfeatures=3000,
                                   outfeatures = bad_genes,
                                   npcs=20,
                                   min_batch_num = 100,
                                   resolution=0.2,
                                   k_anchor=1,
                                   k_score = 120,
                                   k_filter = 600,
                                   k_weight = 100,
                                   use_old_RNA = F,
                                   use_old_Seurat = F)
seob_E10_Mand_E11_Mand@meta.data$celltype_velocyto = str_c(seob_E10_Mand_E11_Mand@meta.data$sample,
                                                           seob_E10_Mand_E11_Mand@meta.data$celltype,
                                                           sep = "_")

saveRDS(seob_E10_Mand_E11_Mand, file = "/mnt/data/home/tycloud/workspace/scATAC/data/E10_Mand_E11_Mand_WT_CNCC_SCT_Waddington_OT_scRNA_seq.rds")

seob_E10_Mand_E11_Mand = readRDS(file = "/mnt/data/home/tycloud/workspace/scATAC/data/E10_Mand_E11_Mand_WT_CNCC_SCT_Waddington_OT_scRNA_seq.rds")

DimPlot(seob_E10_Mand_E11_Mand, reduction = "umap", group.by = "sample",label = T)
DimPlot(seob_E10_Mand_E11_Mand, reduction = "umap", group.by = "celltype_velocyto",label = T)

seurat2anndata(obj = seob_E10_Mand_E11_Mand,
               outFile = "/mnt/data/home/tycloud/workspace/scATAC/data/E10_Mand_E11_Mand_WT_CNCC_SCT_Waddington_OT_scRNA_seq.h5ad",
               slot = 'data', main_layer = 'RNA',
               transfer_layers = NULL,
               drop_single_values = FALSE)

python_path = "/mnt/data/home/tycloud/mambaforge/envs/cellrank3/bin/python"
Wot_path = "python/Wot.py"
input_file = "/mnt/data/home/tycloud/workspace/scATAC/data/E10_Mand_E11_Mand_WT_CNCC_SCT_Waddington_OT_scRNA_seq.h5ad"
out_file = "/mnt/data/home/tycloud/workspace/scATAC/data/Waddington_OT_E10_Mand_E11_Mand.csv"

execute_py(python_path = python_path,
           script_path = Wot_path,
           input_file=input_file,
           out_file=out_file)
Tab_Wot = read.csv(out_file)%>%column_to_rownames(var = "X")
pheatmap(Tab_Wot,
         cluster_rows = F,
         cluster_cols = F)

load(file = "/mnt/data/home/tycloud/workspace/scATAC/data/trajectory_DAG.RData")
plot_trajectory_DAG(Tab_l = Tab_list_Wot,
                    Tab_l2 = NULL,
                    node_df = cell_df,
                    color_by="color",
                    cellid_by="cellid",
                    scale_ = NULL,
                    cut_ = 0.0,
                    order_ = rev(cell_df$cellid),
                    node_size = 2,
                    x_size = 3,
                    y_size = 3,
                    width_top = 0.5,
                    alpha_min = 0.5,
                    width_min = 0.3,
                    width_max = 0.8,
                    edge_color = c("#136783"),
                    node_color = my_pal_science5)

#1
tmp = "/mnt/data/home/tycloud/workspace/algorithms_raw/tmp/"
rna_ad_path = "/mnt/data/home/tycloud/workspace/algorithms_raw/data/PMBC_T.h5ad"
python_path = "/mnt/data/home/tycloud/mambaforge/envs/seacells/bin/python"

execute_py(python_path = "/mnt/data/home/tycloud/mambaforge/envs/seacells/bin/python",
           script_path = "/mnt/data/home/tycloud/workspace/LineageNet/LineageNet/python/SEACells_proprecess.py",
           rna_ad_path= "/mnt/data/home/tycloud/workspace/algorithms_raw/data/PMBC_T.h5ad",
           tmp="/mnt/data/home/tycloud/workspace/algorithms_raw/tmp/")


#2
rna_ad_path = "/mnt/data/home/tycloud/workspace/algorithms_raw/data/PMBC_T.h5ad"
atac_ad_path = "/mnt/data/home/tycloud/workspace/algorithms_raw/data/PMBC_ATAC_T.h5ad"
tmp = "/mnt/data/home/tycloud/workspace/algorithms_raw/tmp/"
result_path = "/mnt/data/home/tycloud/workspace/algorithms_raw/data/"

execute_py(python_path = "/mnt/data/home/tycloud/mambaforge/envs/seacells/bin/python",
           script_path = "/mnt/data/home/tycloud/workspace/LineageNet/LineageNet/python/SEACells_integrated.py",
           rna_ad_path = "/mnt/data/home/tycloud/workspace/algorithms_raw/data/PMBC_T.h5ad",
           atac_ad_path = "/mnt/data/home/tycloud/workspace/algorithms_raw/data/PMBC_ATAC_T.h5ad",
           tmp = "/mnt/data/home/tycloud/workspace/algorithms_raw/tmp/",
           result_path = "/mnt/data/home/tycloud/workspace/algorithms_raw/data/")


seob_T_SEACells = fread(file = "/mnt/data/home/tycloud/workspace/algorithms_raw/data/rna_mate_mat.csv",
                        quote = F, data.table = F, header = T)%>%
  column_to_rownames(var = "V1")
seob_T_SEACells = seob_T_SEACells[, colnames(seob_T_SEACells)!="cellid"]
seob_T_SEACells = t(seob_T_SEACells)
seob_T_SEACells = CreateSeuratObject(counts = seob_T_SEACells)

seob_T_SEACells = RunSeurat(seob_T_SEACells,
                            split_by=NULL,
                            coln_add=NULL,
                            cellid_add=T,
                            genes_add=NULL,
                            integrated="no",
                            integrated_assay="RNA",
                            cellid=NULL,
                            min_cells=3,
                            min_features=200,
                            max_features=0,
                            MT_filter = 20,
                            nfeatures=3000,
                            npcs=30,
                            resolution=0.5,
                            k_anchor=5,
                            k_filter=200,
                            k_score=30,
                            outfeatures=bad_genes,
                            min_batch_num=NULL,
                            k_weight=100)

seob_T_SEACells@meta.data$seurat_clusters =
  as.character(seob_T_SEACells@meta.data$seurat_clusters)

pseudotime_ = Runslingshot(seob=seob_T_SEACells,
                           start='2',
                           end=c('0','1'),
                           reduction="umap",
                           cluster="seurat_clusters")

exp_data = fread(file = "/mnt/data/home/tycloud/workspace/algorithms_raw/data/scMultiome_mat.csv", quote = F,
                 data.table = F)%>%column_to_rownames(var = "cellid")
exp_data = t(exp_data)
MAST_DEG(exp_data,
         pseudotime_lineages = pseudotime_,
         split_num=4,
         log_=T,
         out_dir = "/mnt/data/home/tycloud/workspace/algorithms_raw/result")
