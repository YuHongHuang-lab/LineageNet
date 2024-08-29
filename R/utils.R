#' Title seurat2anndata convert seurat to adata
#'
#' @param obj
#' @param outFile
#' @param slot
#' @param main_layer
#' @param transfer_layers
#' @param drop_single_values
#' @param uns
#' @param obsm_add
#' @param obsm_add_name
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' bm_integrated_E11_E14 = seurat2anndata(obj = bm_integrated_E11_E14,
#' outFile = "data/E12_E14_WT_MT_CNCC_scVelo_scRNA_seq_E11_E14.h5ad",
#' slot = 'counts', main_layer = 'RNA',
#' transfer_layers = c("spliced", "unspliced", "ambiguous"),
#' drop_single_values = FALSE)
#' }

seurat2anndata <- function(
    obj,
    src=NULL,
    outFile = NULL,
    slot = 'counts',
    main_layer = 'RNA',
    transfer_layers = c("spliced", "unspliced", "ambiguous"),
    drop_single_values = TRUE,
    uns = NULL,
    obsm_add = NULL,
    obsm_add_name = NULL
) {

  if (!is.null(src)){
    com_id = intersect(colnames(obj), colnames(src))
    com_genes = intersect(rownames(obj), rownames(src))

    obj = obj[com_genes, com_id]
    src = src[com_genes, com_id]

    for (layer in transfer_layers) {
      obj[[layer]] = CreateAssayObject(counts = src@assays[[layer]]@counts)
    }

    rm(src); gc()
  }

  if (compareVersion(as.character(obj@version), '3.0.0') < 0)
    obj <- Seurat::UpdateSeuratObject(object = obj)

  X <- Seurat::GetAssayData(object = obj, assay = main_layer, slot = slot)

  obs <- .regularise_df(obj@meta.data, drop_single_values = drop_single_values)

  var <- .regularise_df(Seurat::GetAssay(obj, assay = main_layer)@meta.features, drop_single_values = drop_single_values)

  obsm <- NULL
  reductions <- names(obj@reductions)
  if (length(reductions) > 0) {
    obsm <- sapply(
      reductions,
      function(name) as.matrix(Seurat::Embeddings(obj, reduction=name)),
      simplify = FALSE
    )
    names(obsm) <- paste0('X_', tolower(names(obj@reductions)))
  }

  if (!is.null(obsm_add)){
    obsm[[obsm_add_name]] = obsm_add
  }

  layers <- list()
  for (layer in transfer_layers) {
    mat <- Seurat::GetAssayData(object = obj, assay = layer, slot = slot)
    layers[[layer]] <- Matrix::t(mat)
  }

  anndata <- reticulate::import('anndata', convert = FALSE)

  adata <- anndata$AnnData(
    X = Matrix::t(X),
    obs = obs,
    var = var,
    obsm = obsm,
    layers = layers,
    uns = uns
  )

  if (!is.null(outFile))
    adata$write_h5ad(outFile)

  adata
}

.regularise_df <- function(df, drop_single_values = TRUE) {
  if (ncol(df) == 0) df[['name']] <- rownames(df)
  if (drop_single_values) {
    k_singular <- sapply(df, function(x) length(unique(x)) == 1)
    if (sum(k_singular) > 0)
      warning(paste('Dropping single category variables:'),
              paste(colnames(df)[k_singular], collapse=', '))
    df <- df[, !k_singular, drop=F]
    if (ncol(df) == 0) df[['name']] <- rownames(df)
  }
  return(df)
}

#' Title Inferring Trajectories with Slingshot
#'
#' @param seob
#' @param start
#' @param end
#' @param reduction
#' @param cluster
#'
#' @return
#' @export
#'
#' @examples
Runslingshot = function(seob,
                        start=NULL,
                        end=NULL,
                        reduction="umap",
                        cluster="seurat_clusters"){
  library(slingshot)

  rd = Embeddings(seob, reduction=reduction)
  cl = as.numeric(seob@meta.data[, cluster])
  names(cl) = rownames(seob@meta.data)
  sds <- slingshot(rd, cl,
                   start.clus = start,
                   end.clus = end)
  pseudotime_ <- slingPseudotime(sds)
  colnames(pseudotime_) <- apply(pseudotime_, 2,function(x) paste0(start,'to',cl[which.max(x)]))
  row.names(pseudotime_) <- row.names(rd)
  return(pseudotime_)
}

#' Title Differential gene expression analysis by MAST
#'
#' @param exp_data col:cellid, row:gene/peak id
#' @param split_num
#' @param pseudotime_ should exports by slingshot
#' @param log_
#'
#' @return
#' @export
#'
#' @examples
MAST_DEG = function(exp_data,
                    pseudotime_lineages,
                    split_num,
                    log_=T,
                    out_dir = "result"){
  library(MAST)
  exp_data = as.matrix(exp_data)
  if (isTRUE(log_)){
    exp_data = log2(exp_data+1)
  }

  for(p in colnames(pseudotime_lineages))
  {
    lin <- pseudotime_lineages[, p]
    names(lin) <- row.names(pseudotime_lineages)
    lin <- lin[!is.na(lin)]
    lin <- sort(lin, decreasing=F)

    # split K groups based on pseudo-time
    len_per_group <- round(length(lin)/split_num)
    cell_groups <- c()
    for(i in 1:(split_num-1)) cell_groups <- c(cell_groups, rep(i,len_per_group))
    cell_groups <- c(cell_groups, rep(split_num, length(lin)-length(cell_groups)))
    names(cell_groups) <- names(lin)

    # expression data: gene * cell
    expData_p <- exp_data[, names(lin)]
    fData <- data.frame(names=rownames(expData_p), primerid=rownames(expData_p))
    rownames(fData) <- rownames(expData_p)
    cData <- data.frame(cond=cell_groups, wellKey=colnames(expData_p))
    rownames(cData) <- colnames(expData_p)

    sca <- FromMatrix(as.matrix(expData_p), cData, fData)
    colData(sca)$cdr <- scale(colSums(assay(sca)>0)) # cell detection rate
    cond <- factor(colData(sca)$cond)
    cond <- relevel(cond, 1) # set the reference level of the factor to be group 1 cells
    colData(sca)$condition <- cond

    # DE genes between the average expression of the start and end points of a lineage
    zlmdata <- zlm(~ condition + cdr, sca)
    avgLogFC <- 0
    for(i in seq(2,split_num))
    {
      end_group <- paste0('condition', i)
      summaryCond <- summary(zlmdata, doLRT=end_group) # test group i ('condition i')

      #
      summaryDt <- summaryCond$datatable
      fcHurdle <- merge(summaryDt[contrast==end_group & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                        summaryDt[contrast==end_group & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
      fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
      fcHurdle <- merge(fcHurdle, as.data.table(mcols(sca)), by='primerid')
      setnames(fcHurdle,'coef','logFC')
      fcHurdle <- fcHurdle[,-7]

      fcHurdle[fdr>0.01, 'logFC'] <- 0 # set logFC to 0 if fdr>0.01
      avgLogFC <- avgLogFC + abs(fcHurdle$logFC)
    }
    avgLogFC <- avgLogFC/(split_num-1)
    avgLogFC <- as.data.frame(avgLogFC, row.names=fcHurdle$primerid)
    avgLogFC$primerid <- row.names(avgLogFC)
    logFCAll <- as.data.frame(logFC(zlmdata)$logFC)
    logFCAll$primerid <- row.names(logFCAll)
    logFCAll <- merge(logFCAll, avgLogFC, by='primerid')
    # Delete the line containing NA
    NA_row <- unique(which(is.na(logFCAll), arr.ind=T)[, 1])
    if(length(NA_row)>0) logFCAll <- logFCAll[-NA_row, ]
    # Delete the line containing 0
    logFCAll <- logFCAll[logFCAll$avgLogFC!=0, ]

    setorder(logFCAll, -avgLogFC)
    setnames(logFCAll, 'avgLogFC', 'logFC')

    # write to file
    fwrite(logFCAll, paste0(out_dir, '/DEgenes_MAST_sp', split_num, '_', p, '.csv'))
  }
}


#' Title run python file in R
#'
#' @param python_path python interpreter
#' @param script_path python script
#' @param ... sys.argv must be order
#'
#' @return
#' @export
#'
#' @examples
execute_py = function(python_path=NULL,
                      script_path=NULL,
                      ...){
  args = c(script_path,
           unname(unlist(list(...))))
  system2(python_path, args = args, wait = T)
}

