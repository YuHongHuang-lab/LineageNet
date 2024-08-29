#' Title
#' A Wrapper of seurat integrated, reduce and cluster
#' @param seob seurat object
#' @param split_by
#' @param coln_add
#' @param cellid_add
#' @param genes_add
#' @param integrated
#' @param integrated_assay
#' @param cellid
#' @param min_cells
#' @param min_features
#' @param max_features
#' @param MT_filter
#' @param nfeatures
#' @param npcs
#' @param resolution
#' @param k_anchor
#' @param k_filter Adjusting the value of k.filter can change
#' the number of cells considered as potential anchors.
#' Higher k.filter values imply stricter screening criteria,
#' which may reduce the number of anchors
#' but improve the specificity of integration;
#' lower k.filter values, on the other hand,
#' may increase the number of anchors,
#' which helps to capture a larger population of cells
#' but may increase false integrations
#' @param k_score The k.score parameter defines
#' the number of nearest neighbors to be used
#' when evaluating the quality of each anchor candidate.
#' This parameter is mainly used to evaluate the stability of
#' each anchor point:
#' For each potential anchor point, FindIntegrationAnchors calculates its average similarity score with its k.score nearest neighbors.
#' This score is used as a measure of the consistency and reliability
#'  of the anchor point, with a high score
#'  indicating that its neighbors have a high degree of
#'  consistency across datasets.
#'
#' @param outfeatures bad genes
#' @param min_batch_num
#' @param k_weight
#' @param use_old_RNA
#' @param use_old_Seurat
#' @param use_old_embeddings
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#'data("bad_genes")

#' seob = RunSeurat(seob = seob,
#'                  split_by="sample",
#'                  integrated="RPCA",
#'                  outfeatures = bad_genes,
#'                  coln_add = "all",
#'                  cellid=NULL,
#'                  nfeatures=3000,
#'                  min_cells = 3,
#'                  min_features = 200,
#'                  npcs=20,
#'                  resolution=0.15,
#'                  k_anchor=3,
#'                  k_weight = 200,
#'                  use_old_RNA = F,
#'                  use_old_Seurat = F)
#' }
#'
#'
RunSeurat = function(seob, split_by=NULL,
                     coln_add=NULL,
                     cellid_add=T,
                     genes_add=NULL,
                     integrated=NULL,
                     integrated_assay="RNA",
                     cellid=NULL,
                     min_cells=0,
                     min_features=0,
                     max_features=0,
                     MT_filter = 0,
                     nfeatures=3000,
                     npcs=30,
                     resolution=0.3,
                     k_anchor=5,
                     k_filter=200,
                     k_score=30,
                     outfeatures=NULL,
                     min_batch_num=NULL,
                     k_weight=100,
                     use_old_RNA=F,
                     use_old_Seurat=F,
                     use_old_embeddings=F){
  if (isFALSE(use_old_Seurat)){
    if (class(seob)=="list"){
      colid = lapply(seob, function(x)colnames(x@meta.data))
      comcolid = Reduce("intersect", colid)
      seob = lapply(seob, function(x){x@meta.data=x@meta.data[, comcolid];x})
      seob <- mapply(function(x, name) {
        x@meta.data[, 'splitid'] = name;x
      }, seob, names(seob), SIMPLIFY = FALSE)
      split_by = 'splitid'
      if (cellid_add){
        seob = merge(x=seob[[1]],y=seob[-1],
                     add.cell.ids=names(seob))
      }else {
        seob = merge(x=seob[[1]],y=seob[-1])
      }

      DefaultAssay(seob) = "RNA"
    }
    DefaultAssay(seob) = "RNA"

    seob = CreateSeuratObject(counts = seob@assays$RNA@counts,
                              meta.data = seob@meta.data,
                              min.cells = min_cells,
                              min.features = min_features)

    if (MT_filter>0){
      seob[["percent.mt"]] <- PercentageFeatureSet(seob, pattern = "(^MT-)|(^mt-)")
      seob <- subset(seob, subset = percent.mt < MT_filter)
    }

    if (max_features>0){
      seob <- subset(seob, subset = nFeature_RNA  < max_features)
    }

    if (!is.null(cellid)){
      seob = subset(seob, cells=cellid)
    }

    if (!is.null(outfeatures)){
      outfeatures_ = grep("(Rik$)|(^Gm)|(^ENS)", rownames(seob@assays$RNA@counts), value = T)
      outfeatures = unique(c(outfeatures, outfeatures_))
      features = setdiff(rownames(seob@assays$RNA@counts), outfeatures)
      seob = subset(seob, features=features)
    }



    seob_old = seob

    if (integrated=="RPCA"){
      if (integrated_assay == "RNA"){
        if (is.null(coln_add)){
          seob = CreateSeuratObject(counts = seob@assays$RNA@counts,
                                    meta.data = seob@meta.data[, split_by, drop=FALSE])
        }else if (coln_add=="all"){
          seob = CreateSeuratObject(counts = seob@assays$RNA@counts,
                                    meta.data = seob@meta.data)
        }else {
          seob = CreateSeuratObject(counts = seob@assays$RNA@counts,
                                    meta.data = seob@meta.data[, c(split_by, coln_add), drop=FALSE])
        }
        if (!is.null(genes_add)){genes_add = intersect(rownames(seob), genes_add)}
        seob_list = SplitObject(seob,
                                split.by = split_by)

        if (!is.null(min_batch_num)){
          seob_list = lapply(seob_list, FUN = function(x){
            if (ncol(x)>min_batch_num){
              return(x)
            }else {
              return(NULL)
            }
          })
          seob_list = seob_list[!sapply(seob_list, is.null)]
        }


        seob_list <- lapply(X = seob_list, FUN = function(x) {
          x <- NormalizeData(x)
          x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeatures)
        })

        features <- SelectIntegrationFeatures(object.list = seob_list,
                                              nfeatures = nfeatures)
        if (!is.null(genes_add)){features = unique(features, genes_add)}


        seob_list <- lapply(X = seob_list, FUN = function(x) {
          x <- ScaleData(x, features = features, verbose = TRUE)
          x <- RunPCA(x, features = features, npcs = npcs, verbose = TRUE)
        })


        anchors <- FindIntegrationAnchors(object.list = seob_list,
                                          anchor.features = features,
                                          reduction = "rpca",
                                          k.anchor = k_anchor,
                                          k.filter = k_filter,
                                          k.score = k_score,
                                          dims = 1:npcs)
        seob <- IntegrateData(anchorset = anchors,
                              k.weight = k_weight)
        rm(anchors);gc()

        DefaultAssay(seob) <- "integrated"
        seob <- ScaleData(seob, verbose = TRUE)
        seob_old = subset(seob_old, cells=colnames(seob))

        if (isTRUE(use_old_RNA)){
          seob[["RNA"]] = seob_old[["RNA"]]
        }else {
          DefaultAssay(seob) = "RNA"
          seob <- NormalizeData(seob)
        }
        DefaultAssay(seob) <- "integrated"
      }else if (integrated_assay == "integrated"){
        if (is.null(coln_add)){
          seob = CreateSeuratObject(counts = seob@assays$integrated@data,
                                    meta.data = seob@meta.data[, split_by, drop=FALSE])
        }else if (coln_add=="all"){
          seob = CreateSeuratObject(counts = seob@assays$integrated@data,
                                    meta.data = seob@meta.data)
        }else {
          seob = CreateSeuratObject(counts = seob@assays$integrated@data,
                                    meta.data = seob@meta.data[, c(split_by, coln_add), drop=FALSE])
        }
        if (!is.null(genes_add)){genes_add = intersect(rownames(seob), genes_add)}
        seob_list = SplitObject(seob,
                                split.by = split_by)

        if (!is.null(min_batch_num)){
          seob_list = lapply(seob_list, FUN = function(x){
            if (ncol(x)>min_batch_num){
              return(x)
            }else {
              return(NULL)
            }
          })
          seob_list = seob_list[!sapply(seob_list, is.null)]
        }


        seob_list <- lapply(X = seob_list, FUN = function(x) {
          x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeatures)
        })

        features <- SelectIntegrationFeatures(object.list = seob_list,
                                              nfeatures = nfeatures)
        if (!is.null(genes_add)){features = unique(features, genes_add)}


        seob_list <- lapply(X = seob_list, FUN = function(x) {
          x <- ScaleData(x, features = features, verbose = TRUE)
          x <- RunPCA(x, features = features, npcs = npcs, verbose = TRUE)
        })


        anchors <- FindIntegrationAnchors(object.list = seob_list,
                                          anchor.features = features,
                                          reduction = "rpca",
                                          k.anchor = k_anchor,
                                          k.filter = k_filter,
                                          k.score = k_score,
                                          dims = 1:npcs)
        seob <- IntegrateData(anchorset = anchors,
                              k.weight = k_weight)
        rm(anchors);gc()

        DefaultAssay(seob) <- "integrated"
        seob <- ScaleData(seob, verbose = TRUE)
        seob_old = subset(seob_old, cells=colnames(seob))

        if (isTRUE(use_old_RNA)){
          seob[["RNA"]] = seob_old[["RNA"]]
        }else {
          DefaultAssay(seob) = "RNA"
          seob <- NormalizeData(seob)
        }
        DefaultAssay(seob) <- "integrated"
      }
    }else if (integrated=="SCT"){
      if (is.null(coln_add)){
        seob = CreateSeuratObject(counts = seob@assays$RNA@counts,
                                  meta.data = seob@meta.data[, split_by, drop=FALSE])
      }else if (coln_add=="all"){
        seob = CreateSeuratObject(counts = seob@assays$RNA@counts,
                                  meta.data = seob@meta.data)
      }else {
        seob = CreateSeuratObject(counts = seob@assays$RNA@counts,
                                  meta.data = seob@meta.data[, c(split_by, coln_add), drop=FALSE])
      }
      if (!is.null(genes_add)){genes_add = intersect(rownames(seob), genes_add)}
      DefaultAssay(seob) = "RNA"
      seob_list = SplitObject(seob,
                              split.by = split_by)

      if (!is.null(min_batch_num)){
        seob_list = lapply(seob_list, FUN = function(x){
          if (ncol(x)>min_batch_num){
            return(x)
          }else {
            return(NULL)
          }
        })
        seob_list = seob_list[!sapply(seob_list, is.null)]
      }

      seob_list <- lapply(X = seob_list, FUN = function(x) {
        x <- SCTransform(x,
                         variable.features.n = nfeatures)
      })

      features = SelectIntegrationFeatures(object.list = seob_list,
                                           nfeatures = nfeatures)
      if (!is.null(genes_add)){features = unique(features, genes_add)}
      seob_list = PrepSCTIntegration(object.list = seob_list,
                                     anchor.features = features)
      anchors = FindIntegrationAnchors(object.list = seob_list,
                                       normalization.method = "SCT",
                                       anchor.features = features,
                                       k.anchor = k_anchor,
                                       k.filter = k_filter,
                                       k.score = k_score,)
      seob = IntegrateData(anchorset = anchors,
                           normalization.method = "SCT",
                           dims = 1:npcs,
                           k.weight = k_weight)
      rm(anchors);gc()

      DefaultAssay(seob) <- "integrated"
      #seob <- ScaleData(seob, verbose = TRUE)
    }else if (integrated=="no"){
      if (is.null(coln_add)){
        seob = CreateSeuratObject(counts = seob@assays$RNA@counts,
                                  meta.data = seob@meta.data[, split_by, drop=FALSE])
      }else if (coln_add=="all"){
        seob = CreateSeuratObject(counts = seob@assays$RNA@counts,
                                  meta.data = seob@meta.data)
      }else {
        seob = CreateSeuratObject(counts = seob@assays$RNA@counts,
                                  meta.data = seob@meta.data[, c(split_by, coln_add), drop=FALSE])
      }

      DefaultAssay(seob) = "RNA"
      if (!is.null(genes_add)){genes_add = intersect(rownames(seob), genes_add)}
      seob <- NormalizeData(seob)
      seob <- FindVariableFeatures(seob, selection.method = "vst", nfeatures = nfeatures)
      seob@assays$RNA@var.features = unique(c(seob@assays$RNA@var.features, genes_add))
      seob <- ScaleData(seob)
    }
  }else {
    #seob_E12_Mand_Chai <- Seurat::UpdateSeuratObject(object = Mouse1stArchE12)
    if (!is.null(cellid)){
      seob = subset(seob, cells=cellid)
    }

    if (!is.null(outfeatures)){
      features = setdiff(rownames(seob@assays$RNA@counts), outfeatures)
      seob = subset(seob, features=features)
    }
  }

  tryCatch({
    if (isTRUE(use_old_Seurat))DefaultAssay(seob)="integrated"
  },
  error=function(m)DefaultAssay(seob)="RNA")
  if (isFALSE(use_old_embeddings)){
    seob <- RunPCA(seob, npcs = npcs, verbose = TRUE)
    seob <- RunUMAP(seob, reduction = "pca", dims = 1:npcs,
                    return.model = TRUE)
  }

  seob <- FindNeighbors(seob, reduction = "pca", dims = 1:npcs)
  seob <- FindClusters(seob, resolution = resolution)
  DefaultAssay(seob) = "RNA"
  if (integrated=="SCT"){
    seob = NormalizeData(seob)
  }
  rm(seob_old);gc()
  return(seob)
}
