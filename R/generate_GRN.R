#' Title utilize Qiuyue Yuan, Zhana Duren., *Nat Biotechnol*(2024). Apr 12
#' and Wang et al., *Nat Commun*(2023). Dec 20;14(1):8459
#'  method to construct the TF-CRE-Target prior network
#'
#' @param peaks chr1-10000-11000
#' @param Genes
#' @param input the data from Qiuyue Yuan, Zhana Duren., *Nat Biotechnol*(2024). Apr 12
#' @param out
#' @param genome_ human & mouse
#' @param min_overlap the peak_a intersect with peak_b
#' @param Start_ali
#' @param Sparse
#' @param motifList FIMO and motifmatchr
#' @param method
#' @param using_scan
#' @param strings_width
#' @param ncores
#' @param p_thresh motifScan 1e-4
#' @param Gene_network the data from Wang et al., *Nat Commun*(2023). Dec 20;14(1):8459
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' generate_GRN(peaks = peaks,
#' Genes = Genes,
#' input = "/mnt/data/home/tycloud/workspace/paper_code/LINGER/data/data_bulk/",
#' min_overlap = 300,
#' Start_ali = F,
#' motifList = pfmList,
#' method = "FIMO",
#' ncores = 12,
#' p_thresh=1e-5,
#' using_scan=F,
#' Gene_network="/mnt/data/home/tycloud/workspace/algorithms_raw/data/network_human.csv",
#' network_out="/mnt/data/home/tycloud/workspace/algorithms_raw/data/multiome_prior_network.csv")
#' }


generate_GRN = function(peaks,
                        Genes,
                        input,
                        out="out",
                        genome_="hg38",
                        min_overlap=400,
                        Start_ali=F,
                        Sparse=F,
                        motifList=NULL,
                        method="FIMO",
                        using_scan = "union",
                        strings_width=1000,
                        ncores=12,
                        p_thresh=1e-5,
                        Gene_network=NULL){

  if (!dir.exists(out)){
    dir.create(out)
  }

  peaks = str_replace_all(peaks, ":", "-")
  if (genome_=="hg38"){

    peaks_ref = fread(file = str_c(input, "hg38_hg19_pair.bed"), data.table = T)%>%
      .[, .(peaks_hg19=paste(V4,V5,V6,sep = '-'),
            peaks_hg38=paste(V1,V2,V3,sep = '-'))]
    peaks_ref = peaks_ref[!duplicated(peaks_ref$peaks_hg38), ]
    peaks_ref = setkey(peaks_ref, peaks_hg38)


    peaks_ol_df_hg38_hg19 = peaks_overlap(peaks,
                                          peaks_ref$peaks_hg38,
                                          min_overlap=min_overlap,
                                          dt=T)

    if (isTRUE(Start_ali)){
      start_src = str_replace(peaks_ol_df_hg38_hg19$peak_src, ':', '-')
      start_src = as.numeric(str_split(start_src, '-', simplify = T)[,2])
      start_ref = str_replace(peaks_ol_df_hg38_hg19$peak_ref, ':', '-')
      start_ref = as.numeric(str_split(start_ref, '-', simplify = T)[,2])

      peaks_ol_df_hg38_hg19 = peaks_ol_df_hg38_hg19[start_src==start_ref, ]
    }

    peaks_ol_df_hg38_hg19 = peaks_ol_df_hg38_hg19%>%
      .[, peaks_hg19:=peaks_ref[str_replace(peaks_ol_df_hg38_hg19$peak_ref,':','-'), "peaks_hg19"]]%>%
      .[, peaks_hg19:=str_replace(peaks_hg19, '-', ':')]



    Motif_binding_LINGER = generate_Motifpeaks_LINGER(peaks_ol=peaks_ol_df_hg38_hg19,
                                                      input=input,
                                                      out=out)

    # Motif_binding_LINGER = readRDS(file = str_c(out, "/Motif_Peak.rds"))


    Motif_binding_LINGER = Motif_binding_LINGER[, colnames(Motif_binding_LINGER)%in%Genes]
    Motif_binding_LINGER = Motif_binding_LINGER[rowSums(Motif_binding_LINGER)>0, ]


    TF = colnames(Motif_binding_LINGER)
    TF_motifList <- unname(sapply(motifList, function(x) x@name))
    motifList <- motifList[TF_motifList %in% TF]

    if (method=="FIMO"){
      Motif_binding_scan = generate_Motifpeaks_scan(peaks = peaks,
                                                    genome_=genome_,
                                                    meme_bin_path="/mnt/data/home/tycloud/meme/bin",
                                                    motifList=motifList,
                                                    method="FIMO",
                                                    strings_width=strings_width,
                                                    p_thresh=p_thresh)

    }else if (method=="motifmatchr"){
      Motif_binding_scan = generate_Motifpeaks_scan(peaks = peaks,
                                                    genome_=genome_,
                                                    motifList=motifList,
                                                    method="motifmatchr",
                                                    strings_width=strings_width)
    }

    Motif_binding_LINGER = tidyfst::mat_df(as.matrix(Motif_binding_LINGER))
    Motif_binding_LINGER = as.data.table(Motif_binding_LINGER)
    Motif_binding_LINGER = setnames(Motif_binding_LINGER, c("row", "col"), c("peak", "motif_id"))
    Motif_binding_LINGER = Motif_binding_LINGER[Motif_binding_LINGER$value>0, ]

    Motif_binding_LINGER = paste(Motif_binding_LINGER$peak,
                                 Motif_binding_LINGER$motif_id,
                                 sep = "_")

    Motif_binding_scan = paste(Motif_binding_scan$seqnames,
                               Motif_binding_scan$motif_id,
                               sep = "_")

    if (using_scan=="union"){
      Motif_peak = unique(c(Motif_binding_LINGER, Motif_binding_scan))
    }else if (using_scan=="inter"){
      Motif_peak = intersect(Motif_binding_LINGER, Motif_binding_scan)
    }
    Motif_peak = data.table(peak=str_split(Motif_peak, "_", simplify = T)[, 1],
                            Gene=str_split(Motif_peak, "_", simplify = T)[, 2])


    peak_Target = fread(file = str_c(input, "RE_gene_corr_hg19.bed"), data.table = T)%>%
      .[, .(peaks_hg19=paste(V1,':',V2,'-',V3,sep = ''),
            Gene=V4,
            value=V5)]

    peak_Target = merge(peaks_ol_df_hg38_hg19, peak_Target,
                        by = "peaks_hg19")
    peak_Target = peak_Target[peak_Target$peak_src%in%Motif_peak$peak &
                                peak_Target$Gene%in%Genes, ]

    Gene_network = fread(Gene_network, data.table = T)
    Gene_network = Gene_network[Gene_network$from%in%Genes &
                                  Gene_network$to%in%Genes, c("from", "to")]

    Motif_peak = Motif_peak[, .(from=Gene,
                                to=peak)]
    peak_Target = peak_Target[, .(from=peak_src,
                                  to=Gene)]

    multiome_prior_network = rbind(Motif_peak,
                                   peak_Target,
                                   Gene_network)
    write.csv(multiome_prior_network, file = network_out,
              row.names = F, quote = F)
  }
}

generate_Motifpeaks_LINGER = function(peaks_ol,out,input, Sparse=F){
  library(igraph)
  #read motif-target
  files = dir(input)[str_detect(dir(input), "MotifTarget_Matrix_")]
  Motif_TF = fread(file = file.path(input, "Match2.txt"),data.table = F)
  Motif_Weight = fread(file = file.path(input, "motifWeight.txt"),data.table = F)

  all_row_names <- NULL
  all_col_names <- NULL
  Motif_binding <- NULL

  for (file_ in files){
    chr_ = str_extract(file_, "chr[0-9XY]+")
    print(str_c("generate Motif_peak for ", chr_))
    Motif_df = fread(file = file.path(input, file_),data.table = T)
    setnames(Motif_df, "V1", "peak")
    setkey(Motif_df, peak)

    Motif_df = Motif_df[Motif_df$peak%in%peaks_ol$peaks_hg19, ]
    Motif_id = colnames(Motif_df)[2:length(Motif_df)]
    Motif_df = transpose(Motif_df, make.names = "peak")

    Motif_df = Motif_df[, Motif :=Motif_id]
    Motif_df = setcolorder(Motif_df, c("Motif", setdiff(names(Motif_df), "Motif")))

    Motif_TF_ = Motif_TF[!duplicated(Motif_TF$Match21), ]
    rownames(Motif_TF_) = Motif_TF_$Match21

    rownames(Motif_Weight) = Motif_Weight$motifName

    com_Motif = Reduce("intersect", list(Motif_df$Motif,
                                         rownames(Motif_TF_),
                                         rownames(Motif_Weight)))

    Motif_TF_ = Motif_TF_[com_Motif, ]
    Motif_Weight = Motif_Weight[com_Motif, "motifWeight", drop=F]

    Motif_df = Motif_df[Motif_df$Motif %in% com_Motif, ]
    Motif_df = Motif_df[match(com_Motif, Motif_df$Motif), ]

    peaks_id = colnames(Motif_df)[2:length(Motif_df)]

    Motif_df = as(Motif_df[,2:length(Motif_df)], "sparseMatrix")
    rownames(Motif_df) = com_Motif


    Motif_binding_ <- diag(as.numeric(1.0 / (t(Motif_Weight) + 0.1)))%*%Motif_df
    Motif_binding_ = log1p(Motif_binding_)
    rownames(Motif_binding_) = rownames(Motif_df)

    Motif_binding_ = as.data.table(Motif_binding_)
    Motif_binding_ = Motif_binding_[, Motif :=com_Motif]
    Motif_binding_ = setcolorder(Motif_binding_, c("Motif", setdiff(names(Motif_binding_), "Motif")))

    if (isTRUE(Sparse)){
      all_row_names <- union(all_row_names, rownames(Motif_binding_))
      all_col_names <- union(all_col_names, colnames(Motif_binding_))

      if (is.null(Motif_binding)){
        Motif_binding <- Matrix(0,
                                nrow = length(all_row_names),
                                ncol = length(all_col_names),
                                dimnames = list(all_row_names, all_col_names))
      }else {
        Motif_binding_extend <- Matrix(0, nrow = length(all_row_names), ncol = length(all_col_names),
                                       dimnames = list(all_row_names, all_col_names))
        Motif_binding_extend[rownames(Motif_binding), colnames(Motif_binding)] =
          Motif_binding
        Motif_binding = Motif_binding_extend
      }

      Motif_binding[rownames(Motif_binding_), colnames(Motif_binding_)] =
        Motif_binding_
      rm(Motif_binding_, Motif_df, Motif_binding_extend); gc()
    }else {
      if (is.null(Motif_binding)){
        Motif_binding = Motif_binding_
      }else {
        Motif_binding <- merge(Motif_binding,
                               Motif_binding_,
                               by = "Motif",
                               all = TRUE,
                               suffixes = c("", ".y"))
      }
      Motif_binding[is.na(Motif_binding)] = 0
      rm(Motif_binding_, Motif_df); gc()
    }

  }

  #sum motif
  Motif_id = Motif_binding$Motif
  Motif_binding = as(Motif_binding[, 2:ncol(Motif_binding)], "sparseMatrix")
  rownames(Motif_binding) = Motif_id

  Motif_TF_ = graph_from_data_frame(Motif_TF[, c("Match22","Match21")], directed = F)
  Motif_TF_ = as_adjacency_matrix(Motif_TF_, sparse = TRUE)
  Motif_TF_ = Motif_TF_[unique(Motif_TF$Match22), rownames(Motif_binding)]

  Motif_binding = Motif_TF_%*%Motif_binding

  #sum peak
  Motif_binding = t(Motif_binding)

  peaks_ol_ = graph_from_data_frame(peaks_ol[, c("peak_src","peaks_hg19")], directed = F)
  peaks_ol_ = as_adjacency_matrix(peaks_ol_, sparse = TRUE)
  peaks_ol_ = peaks_ol_[unique(peaks_ol$peak_src), rownames(Motif_binding)]

  Motif_binding = peaks_ol_%*%Motif_binding

  print("Save Sparse Motif-Peak data")
  saveRDS(Motif_binding, file = str_c(out, "/Motif_Peak.rds"))
  return(Motif_binding)
}

#' Title
#'
#' @param peaks peaks should be chr1-8888-8888
#' @param meme_bin_path
#' @param genome_
#' @param motifList if method=="FIMO", the length of motifList should
#' be < 1500
#' @param pfmList
#' @param method Parallel computation of "motifmatchr"
#' has not been tested and is not recommended
#'
#' @return
#' @export
#'
#' @examples
generate_Motifpeaks_scan = function(peaks,
                                    meme_bin_path="/mnt/data/home/tycloud/meme/bin",
                                    genome_,
                                    motifList=NULL,
                                    pfmList=NULL,
                                    method="FIMO",
                                    strings_width=1000,
                                    ncores=12,
                                    p_thresh=1e-5){

  if (method=="FIMO"){
    options(meme_bin = meme_bin_path)
    library(memes)
    library(universalmotif)

    motifList_universal = list()
    for (i in seq_along(motifList)) {
      print(i)
      motifList_universal[[i]] <- convert_motifs(motifList[[i]], class = "universalmotif")
    }

    peaks = peaks_to_string(peaks, genome_)

    peak_list <- base::split(peaks, ceiling(seq_along(peaks)/strings_width))

    if (ncores==1){
      fimo_results <- lapply(peak_list,
                             runFimo,
                             motifs = motifList_universal,
                             skip_matched_sequence = T,
                             thresh=p_thresh)
    }else {
      fimo_results <- mclapply(peak_list, runFimo,
                               motifs = motifList_universal,
                               skip_matched_sequence = TRUE,
                               mc.cores = ncores,
                               thresh=p_thresh)
      fimo_results <- lapply(fimo_results, function(x)as.data.table(x))
      fimo_results = do.call("rbind", fimo_results)

      fimo_results$seqnames = str_replace(fimo_results$seqnames,
                                          '-',':')
    }
    return(fimo_results)
  }else if (method=="motifmatchr"){
    library(motifmatchr)
    peak_list <- base::split(peaks, ceiling(seq_along(peaks)/strings_width))

    peak_list = lapply(peak_list, function(x){
      GRanges(seqnames = str_split(x, "-", simplify = T)[, 1],
              ranges = IRanges(start = as.integer(str_split(x, "-", simplify = T)[, 2]),
                               end = as.integer(str_split(x, "-", simplify = T)[, 3])))
    })

    motifmatchr_results <- lapply(peak_list, function(x){
      matchMotifs(motifList, x, genome=genome_)
    })
  }


}


#' Title
#' peaks and peaks_ref should be chr1-8888-8888
#' @param peaks
#' @param peaks_ref
#'
#' @return
#' @export
#'
#' @examples
peaks_overlap = function(peaks,
                         peaks_ref,
                         min_overlap,
                         dt=T){
  library(GenomicRanges)
  peaks = str_replace_all(peaks, ":", "-")
  peaks_ref = str_replace_all(peaks_ref, ":", "-")
  peaks = GRanges(seqnames = str_split(peaks, "-", simplify = T)[, 1],
                  ranges = IRanges(start = as.integer(str_split(peaks, "-", simplify = T)[, 2]),
                                   end = as.integer(str_split(peaks, "-", simplify = T)[, 3])))

  peaks_ref = GRanges(seqnames = str_split(peaks_ref, "-", simplify = T)[, 1],
                      ranges = IRanges(start = as.integer(str_split(peaks_ref, "-", simplify = T)[, 2]),
                                       end = as.integer(str_split(peaks_ref, "-", simplify = T)[, 3])))

  peaks_ol = findOverlaps(peaks, peaks_ref, minoverlap = min_overlap)

  peaks = as.data.frame(peaks)
  peaks$peak_src = str_c(peaks$seqnames,":",peaks$start,"-",peaks$end)

  peaks_ref = as.data.frame(peaks_ref)
  peaks_ref$peak_ref = str_c(peaks_ref$seqnames,":",peaks_ref$start,"-",peaks_ref$end)

  peaks_ol_df = data.frame(peak_src=peaks[queryHits(peaks_ol), "peak_src"],
                           peak_ref=peaks_ref[subjectHits(peaks_ol), "peak_ref"])
  if (isTRUE(dt)){
    peaks_ol_df = as.data.table(peaks_ol_df)
  }
  return(peaks_ol_df)
}

#' Title
#' peaks should be chr1-8888-8888
#'
#' @param peaks
#'
#' @return
#' @export
#'
#' @examples
peaks_to_string = function(peaks, genome_){
  library(GenomicRanges)
  library(Biostrings)

  peaks = peaks[!str_detect(peaks, "KI|GL")]
  peaks_names = peaks
  peaks = GRanges(seqnames = str_split(peaks, "-", simplify = T)[, 1],
                  ranges = IRanges(start = as.integer(str_split(peaks, "-", simplify = T)[, 2]),
                                   end = as.integer(str_split(peaks, "-", simplify = T)[, 3])))
  if (genome_=="hg38"){
    library(BSgenome.Hsapiens.UCSC.hg38)
    peaks_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, peaks)
    peaks_seq <- DNAStringSet(peaks_seq)
    names(peaks_seq) = peaks_names
    return(peaks_seq)
  }

}



#' Title run dictys dynamic subsets_rna
#'
#' @param python_path
#' @param ... the param of subsets_rna
#'
#' @return
#' @export
#'
#' @examples
generate_pseudotime_windows = function(python_path=NULL,
                                       ...){
  args_list = list(...)
  fi_traj <- args_list$fi_traj
  fi_traj_cell_rna <- args_list$fi_traj_cell_rna
  fi_coord_rna <- args_list$fi_coord_rna
  fo_subsets <- args_list$fo_subsets
  fo_subset_locs <- args_list$fo_subset_locs
  diro_subsets <- args_list$diro_subsets
  fo_subset_edges <- args_list$fo_subset_edges
  ncell <- as.character(args_list$ncell)
  noverlap <- as.character(args_list$noverlap)
  dmax <- as.character(args_list$dmax)

  args = c("-m", "dictys", "dynamic", "subsets_rna",
           fi_traj,
           fi_traj_cell_rna,
           fi_coord_rna,
           fo_subsets,
           fo_subset_locs,
           diro_subsets,
           fo_subset_edges,
           ncell,
           noverlap,
           dmax)
  system2(python_path, args = args, wait = T)
}


generate_windows_data = function(window_dir,
                                 seob,
                                 seob_ATAC,
                                 logFC_RNA=0.2,
                                 logFC_ATAC=0.2,
                                 de_pct_RNA=0.2,
                                 de_pct_ATAC=0.2){
  data = rhdf5::h5read(paste(window_dir, "subset_locs.h5", sep = "/"), "edges")
  data = data.frame(window_id=1:length(data),
                    order_in_branch=data)
  write.csv(data, file = paste(window_dir,"window_order_in_branch.csv",sep = '/'),
            row.names = F, quote = T)

  files = dir(window_dir)[str_detect(dir(window_dir), "Subset")]
  files <- files[order(as.numeric(gsub(".*?(\\d+)$", "\\1", files)))]

  cellid_df = data.table()
  DEG_df = data.table()
  DEpeak_df = data.table()
  for (id in seq(length(files))){
    print(paste("window: ", id))
    branch_id = data[id]
    file_tmp = files[id]
    file_path = file.path(window_dir, file_tmp)
    cellid = fread(paste(file_path, "names_rna.txt", sep = "/"), data.table = T)%>%
      .[, .(cellid=V1)]%>%.[,"windowid":=id]%>%.[,"branch_id":=branch_id]
    cellid_df = rbind(cellid_df, cellid)

    seob@meta.data$window = ifelse(rownames(seob@meta.data)%in%cellid$cellid,'1','0')
    seob_ATAC@meta.data$window = ifelse(rownames(seob_ATAC@meta.data)%in%cellid$cellid,'1','0')

    Idents(seob_T) = seob_T@meta.data$window
    DEG = FindAllMarkers(seob_T,
                         Mandn.pct = 0.2,
                         only.pos=TRUE,
                         logfc.threshold = logFC_RNA)
    DEG$de_pct = (DEG$pct.1-DEG$pct.2)/DEG$pct.1
    DEG = DEG[DEG$de_pct>de_pct_RNA, ]
    DEG = DEG[DEG$cluster=="1", "gene"]
    DEG = data.table(gene=DEG, window=id)
    DEG_df = rbind(DEG_df, DEG)
    write.table(DEG$gene, file = paste(file_path, "Genes.txt", sep = "/"),
                quote = F)


    Idents(seob_ATAC) = seob_ATAC@meta.data$window
    DEpeak = FindAllMarkers(seob_ATAC,
                            Mandn.pct = 0.2,
                            only.pos=TRUE,
                            logfc.threshold = logFC_ATAC)
    DEpeak$de_pct = (DEpeak$pct.1-DEpeak$pct.2)/DEpeak$pct.1
    DEpeak = DEpeak[DEpeak$de_pct>de_pct_ATAC, ]
    DEpeak = DEpeak[DEpeak$cluster=="1", "gene"]
    DEpeak = data.table(gene=DEpeak, window=id)
    DEpeak_df = rbind(DEpeak_df, DEpeak)

    write.table(DEpeak$gene, file = paste(file_path, "peaks.txt", sep = "/"),
                quote = F)

  }

}



