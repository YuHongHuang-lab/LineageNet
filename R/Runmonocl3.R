#' Title a Warpper function of monocle3
#'
#' @param seob
#' @param cell_metadata
#' @param useumap
#' @param cluster
#' @param use_ncenter
#' @param ncenter
#' @param num_dim
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' cds_E8_5_E9 = Runmonocl3(seob = seob_E8_5_E9,
#' cluster = "celltype_velocyto")
#' }

Runmonocl3 = function(seob, cell_metadata,
                      useumap="seurat",
                      cluster="seurat_clusters",
                      use_ncenter=TRUE,
                      ncenter=300,
                      num_dim=50){
  if (class(seob)=="Seurat"){
    data <- GetAssayData(seob, assay = 'RNA', slot = 'counts')
    cell_metadata <- seob@meta.data
  }

  gene_annotation <- data.frame(gene_short_name = rownames(data))
  rownames(gene_annotation) = gene_annotation$gene_short_name

  cds <- new_cell_data_set(data,
                           cell_metadata = cell_metadata[colnames(data), ],
                           gene_metadata = gene_annotation[rownames(data), , drop=FALSE])
  cds <- preprocess_cds(cds, num_dim = num_dim)
  cds <- reduce_dimension(cds, preprocess_method = "PCA")

  if (useumap=="seurat"){
    cds.embed <- cds@int_colData$reducedDims$UMAP
    int.embed <- Embeddings(seob, reduction = "umap")
    int.embed <- int.embed[rownames(cds.embed),]
    cds@int_colData$reducedDims$UMAP <- int.embed
    plot_cells(cds, reduction_method="UMAP", color_cells_by=cluster,
               cell_size = 0.8)
  }


  cds <- cluster_cells(cds)

  if (isTRUE(use_ncenter)){
    cds <- learn_graph(cds,learn_graph_control=list(ncenter=ncenter),
                       use_partition = FALSE)
  }else {
    cds <- learn_graph(cds)
  }

  p = plot_cells(cds, label_groups_by_cluster = FALSE,
                 label_leaves = FALSE,
                 label_branch_points = FALSE,
                 color_cells_by = cluster,
                 cell_size = 0.8)
  print(p)
  return(cds)
}
