#' use ggplot2 to plot heatmap and bubble, can cluster row or col, add row and col metadata
#'
#' @param object a seurat object or data.frame,
#' if object is data.frame,it should be a two colmuns var(not a mat),
#' colinfo and rowinfo all should be in
#' @param features
#' @param model #' 1 metadata;
#' 2 3 metadata;
#' 3 Don't known how many data
#' @param my_pal
#' @param my_pal_10
#' @param my_pal_2
#' @param my_mal_gradientn
#' @param layer
#' @param ID '_' can't in ID
#' @param r_ID
#' @param out_ID only used in model 1 or 2
#' @param split_col
#' @param cluster_ID
#' @param cluster_feature must >2 for hcluster
#' @param flip_
#' @param plot_style
#' @param white_line_size
#' @param white_line_color
#' @param legend_position
#' @param legend_direction
#' @param y_size
#' @param x_size
#' @param x_linewidth
#' @param y_linewidth
#' @param x_angle
#' @param replace_
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' plotDot(object = seob,
#' features=unique(c('Sox9', 'Col2a1',
#'                   'Pecam1','Cdh5',
#'                  'Krt14','Epcam',
#'                   'Plp1','Sox10',
#'                   'Mitf',
#'                   'Prrx1', 'Twist1',
#'                   'Myod1','Myog')),
#' mode=1, ID=c('celltype'),
#' out_ID = c('celltype'),
#' my_pal = colorRampPalette(c(my_pal_science6))(7),
#' my_mal_gradientn = c("#4C3399","#F5F2EF","#FFE714"),
#' cluster_ID = F,
#' cluster_feature = F,
#' x_angle=45)
#' }


plotDot = function(object,
                   features,
                   model=1,
                   my_pal,
                   my_pal_10=c("#4292C6", "#EF3B2C"),
                   my_pal_2 = c("#F7903D", "#4D85BD", "#59A95A"),
                   my_mal_gradientn=NULL,
                   layer="data",
                   ID=NULL,
                   r_ID=NULL,
                   out_ID=NULL,
                   split_col=NULL,
                   cluster_ID=TRUE,
                   cluster_feature=TRUE,
                   flip_=F,
                   plot_style="bubble",
                   white_line_size=0.5,
                   white_line_color="white",
                   legend_position="rigth",
                   legend_direction="vertical",
                   y_size=10,
                   x_size=10,
                   x_linewidth=0.5,
                   y_linewidth=0.5,
                   x_angle=0,
                   replace_=F){

  if (model==3&class(object)!="Seurat"){
    cluster_ID = F
    cluster_feature = F
  }
  #source("/mnt/data/home/tycloud/workspace/scATAC/utils.R")

  my_pal_ = my_pal
  if (is.null(my_mal_gradientn)){
    my_mal_gradientn = viridis::viridis(20)
  }else {
    my_mal_gradientn = colorRampPalette(my_mal_gradientn)(20)
  }


  #dat = DotPlot(object,features = features, group.by = ID[1])$data

  if (replace_){
    if (class(object)=="Seurat"){
      for (name_ in ID){
        object@meta.data[,name_] =
          str_replace_all(object@meta.data[,name_],"_","-")
      }
    }
  }

  if (model==1){
    cells <- unlist(x = Seurat::CellsByIdentities(object = object, idents = NULL))
    data.features <- Seurat::FetchData(object = object, vars = features,
                                       cells = cells,
                                       layer = layer)
    data.features$id <- object[[ID[1], drop = TRUE]][cells, drop = TRUE]

    if (!is.factor(x = data.features$id)) {
      data.features$id <- factor(x = data.features$id)
    }
    id.levels <- levels(x = data.features$id)
    data.features$id <- as.vector(x = data.features$id)

    # "id" column must be in the end of data.features
    data_plot = avg_pct_func_(data.features)

    avg.exp.scaled <- sapply(X = unique(x = data_plot$features.plot),
                             FUN = function(x) {
                               data.use <- data_plot[data_plot$features.plot ==
                                                       x, "avg.exp"]

                               data.use <- scale(x = data.use)
                               data.use <- MinMax(data = data.use, min = -2.5,
                                                  max = 2.5)
                               return(data.use)
                             })
    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))

    data_plot$avg.exp.scaled <- avg.exp.scaled
    data_plot$features.plot <- factor(x = data_plot$features.plot,
                                      levels = features)
    data_plot$pct.exp[data_plot$pct.exp < 0] <- NA
    data_plot$pct.exp <- data_plot$pct.exp * 100



    mat <- data_plot %>% dplyr::select(-pct.exp, -avg.exp.scaled) %>%
      tidyr::pivot_wider(names_from = features.plot, values_from = avg.exp ) %>%
      data.frame()

    row.names(mat) <- mat$id
    mat <- mat[,-1]
    clust_col <- hclust(dist(mat %>% as.matrix()%>% t()))
    ddgram_col <- ape::as.phylo(clust_col)

    clust_row <- hclust(dist(mat %>% as.matrix() ))
    ddgram_row <- ape::as.phylo(clust_row)

    ggtree_plot_col <- ggtree::ggtree(ddgram_col,branch.length="none")+
      geom_tippoint(color="#F7903D", shape=16, size=3)
    ggtree_plot_row <- ggtree::ggtree(ddgram_row,branch.length="none")+
      geom_tippoint(color="#4D85BD", shape=16, size=3)+
      layout_dendrogram()

    if (isTRUE(cluster_ID)){
      data_plot$id <-  factor(as.character(data_plot$id) , levels = clust_row$labels[clust_row$order])
      level_row = clust_row$order
    }else {
      data_plot$id = factor(data_plot$id, levels = id.levels)
      level_row = 1:length(clust_row$order)
    }

    if (isTRUE(cluster_feature)){
      data_plot$features.plot<-  factor(as.character(data_plot$features.plot) , levels = clust_col$labels[clust_col$order])
    }

    if (length(my_pal)<length(unique(data_plot$id))){
      if (is.null(my_pal)){
        my_pal = scales::hue_pal()(length(unique(data_plot$id)))
      }else {
        my_pal = colorRampPalette(my_pal)(length(unique(data_plot$id)))
      }
    }

    group =
      ggplot(data_plot, aes(id,1,fill=id))+
      geom_tile()+
      scale_y_discrete(expand = c(0,0),position="right")+
      scale_x_discrete(expand=c(0,0)) +
      scale_fill_manual(values=my_pal[level_row])+
      theme_void()+
      theme(axis.text.y=element_text(color="black",size=10))+
      theme(  legend.key=element_blank(),
              legend.text = element_text(color="black",size=10),
              legend.spacing.x=unit(0.1,'cm'),
              legend.key.width=unit(0.5,'cm'),
              legend.key.height=unit(0.5,'cm'),
              legend.background=element_blank())

    if (plot_style=="bubble"){
      dotplot <- ggplot(data_plot,
                        aes(x=id,y = features.plot,size = pct.exp, color = avg.exp.scaled))+
        geom_point() +
        scale_size("% detected", range = c(0,6)) +
        scale_y_discrete(position = "right") +
        scale_x_discrete(position = "top") +
        scale_color_gradientn(colours = my_mal_gradientn,
                              guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"),
                              name = "Average\nexpression") +
        theme_bw() +
        xlab("") + ylab("") +
        theme(
          axis.text.x.top = element_text(size=12, angle=x_angle, hjust=0.5, color="black",face="bold"),
          axis.ticks.x = element_blank(),

          axis.text.y = element_text(size=12, color="skyblue",face="bold"),
          axis.ticks.y = element_blank(),
          axis.text.y.right = element_text(size=8, angle=0, hjust=0.5, color="black",face="bold"),
          axis.line = element_line(colour = 'grey30',size = 0.2),
          panel.grid = element_blank()
        )
    }else if(plot_style=="heatmap"){
      color_map = scales::col_numeric(palette = c(my_pal_[1],
                                                  my_pal_[2],
                                                  my_pal_[3]),
                                      domain = range(data_plot$avg.exp.scaled))
      data_plot$color <- color_map(data_plot$avg.exp.scaled)

      dotplot <- ggplot(data_plot, aes(x = features.plot, y = id)) +
        geom_tile(aes(fill = avg.exp.scaled), size = white_line_size, colour = white_line_color) +
        scale_fill_gradientn(colors = c(my_pal_[1], my_pal_[1],my_pal_[2], my_pal_[3])) +
        scale_y_discrete(position = "right") +
        scale_x_discrete(position = "top") +
        theme_bw() +
        xlab("") + ylab("") +
        theme(
          axis.text.x.top = element_text(size = x_size, angle = x_angle, hjust = 0.5, color = "black", face = "bold"),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = y_size, color = "skyblue", face = "bold"),
          axis.ticks.y = element_blank(),
          axis.text.y.right = element_text(size = y_size, angle = 0, hjust = 0, vjust = 0,
                                           color = "black", face = "bold"),
          axis.line = element_line(colour = 'grey30', size = 0.2),
          panel.grid = element_blank()
        )
    }

    if (isTRUE(cluster_ID)&isFALSE(cluster_feature)){
      dotplot%>%
        aplot::insert_top(ggtree_plot_row, height = 0.2)%>%
        aplot::insert_bottom(group, height = 0.03)
    }else if (isFALSE(cluster_ID)&isTRUE(cluster_feature)){
      dotplot%>%aplot::insert_left(ggtree_plot_col, width = 0.2)%>%
        aplot::insert_bottom(group, height = 0.03)
    }else if (isTRUE(cluster_ID)&isTRUE(cluster_feature)){
      dotplot%>%aplot::insert_left(ggtree_plot_col, width = 0.2)%>%
        aplot::insert_top(ggtree_plot_row, height = 0.2)%>%
        aplot::insert_bottom(group, height = 0.03)
    }else {
      if (plot_style=="bubble"){
        dotplot%>%aplot::insert_bottom(group, height = 0.03)
      }else {
        dotplot
      }
    }

  }else if (model==2){
    cells <- unlist(x = Seurat::CellsByIdentities(object = object, idents = NULL))
    data.features <- Seurat::FetchData(object = object, vars = features,
                                       cells = cells)
    data.features$id <- object[[ID[1], drop = TRUE]][cells, drop = TRUE]

    # "id" column must be in the end of data.features
    data_plot = avg_pct_func_(data.features)

    # row and col order
    mat <- data_plot %>% dplyr::select(-pct.exp) %>%
      tidyr::pivot_wider(names_from = features.plot, values_from = avg.exp ) %>%
      data.frame()

    row.names(mat) <- mat$id
    mat <- mat[,-1]
    clust_col <- hclust(dist(mat %>% as.matrix()%>% t()))
    ddgram_col <- ape::as.phylo(clust_col)

    clust_row <- hclust(dist(mat %>% as.matrix() ))
    ddgram_row <- ape::as.phylo(clust_row)

    ggtree_plot_col <- ggtree::ggtree(ddgram_col,branch.length="none")+
      geom_tippoint(color="#F7903D", shape=16, size=3)
    ggtree_plot_row <- ggtree::ggtree(ddgram_row,branch.length="none")+
      geom_tippoint(color="#4D85BD", shape=16, size=3)+
      layout_dendrogram()

    id_level <-  clust_row$labels[clust_row$order]
    features_level<-  clust_col$labels[clust_col$order]

    #mode1 2
    data.features <- Seurat::FetchData(object = object, vars = features,
                                       cells = cells)
    data.features$id_1 <- object[[ID[1], drop = TRUE]][cells, drop = TRUE]
    data.features$id_2 <- object[[ID[2], drop = TRUE]][cells, drop = TRUE]
    data.features$id_3 <- object[[ID[3], drop = TRUE]][cells, drop = TRUE]

    data.features = data.features%>%tidyr::unite("id", id_1:id_3,
                                                 sep = "_",
                                                 remove = FALSE)
    data.features = data.features[, 1:(ncol(data.features)-3)]
    data_plot = avg_pct_func_(data.features)

    avg.exp.scaled <- sapply(X = unique(x = data_plot$features.plot),
                             FUN = function(x) {
                               data.use <- data_plot[data_plot$features.plot ==
                                                       x, "avg.exp"]

                               data.use <- scale(x = data.use)
                               data.use <- MinMax(data = data.use, min = -2.5,
                                                  max = 2.5)
                               return(data.use)
                             })
    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))

    data_plot$avg.exp.scaled <- avg.exp.scaled
    data_plot$features.plot <- factor(x = data_plot$features.plot,
                                      levels = features)
    data_plot$pct.exp[data_plot$pct.exp < 0] <- NA
    data_plot$pct.exp <- data_plot$pct.exp * 100

    data_plot = data_plot%>%tidyr::separate(col = id,
                                            sep = "_",
                                            into = out_ID,
                                            remove = FALSE)

    data_plot[, out_ID[1]] <-  factor(data_plot[, out_ID[1]], levels = id_level)
    data_plot$features.plot<-  factor(as.character(data_plot$features.plot) , levels = features_level)

    data_plot = data_plot[order(data_plot[, out_ID[1]],
                                data_plot[, out_ID[2]],
                                data_plot[, out_ID[3]]), ]
    data_plot$id = factor(data_plot$id, levels = unique(data_plot$id))

    if (any(is.na(as.numeric(clust_row$labels[clust_row$order])))){
      level_row = clust_row$order
    }else {
      level_row = as.numeric(clust_row$labels[clust_row$order])+1
    }

    if (length(my_pal)<length(unique(data_plot$id))){
      my_pal = scales::hue_pal()(length(unique(data_plot$id)))
    }

    x = out_ID[1]
    group1 =
      ggplot(data_plot, aes(x = id, y = 1,
                            fill=!!as.symbol(x)))+
      geom_tile()+
      scale_y_discrete(expand = c(0,0),position="right")+
      scale_x_discrete(expand=c(0,0)) +
      scale_fill_manual(values=my_pal[level_row])+
      theme_void()+
      theme(axis.text.y=element_text(color="black",size=10))+
      theme(  legend.key=element_blank(),
              legend.text = element_text(color="black",size=10),
              legend.spacing.x=unit(0.1,'cm'),
              legend.key.width=unit(0.5,'cm'),
              legend.key.height=unit(0.5,'cm'),
              legend.background=element_blank())

    x = out_ID[2]
    group2 =
      ggplot(data_plot, aes(x = id, y = 1,
                            fill=!!as.symbol(x)))+
      geom_tile()+
      scale_y_discrete(expand = c(0,0),position="right")+
      scale_x_discrete(expand=c(0,0)) +
      scale_fill_manual(values=my_pal_10)+
      theme_void()+
      theme(axis.text.y=element_text(color="black",size=10))+
      theme(  legend.key=element_blank(),
              legend.text = element_text(color="black",size=10),
              legend.spacing.x=unit(0.1,'cm'),
              legend.key.width=unit(0.5,'cm'),
              legend.key.height=unit(0.5,'cm'),
              legend.background=element_blank())

    x = out_ID[3]
    group3 =
      ggplot(data_plot, aes(x = id, y = 1,
                            fill=!!as.symbol(x)))+
      geom_tile()+
      scale_y_discrete(expand = c(0,0),position="right")+
      scale_x_discrete(expand=c(0,0)) +
      scale_fill_manual(values=my_pal_2)+
      theme_void()+
      theme(axis.text.y=element_text(color="black",size=10))+
      theme(  legend.key=element_blank(),
              legend.text = element_text(color="black",size=10),
              legend.spacing.x=unit(0.1,'cm'),
              legend.key.width=unit(0.5,'cm'),
              legend.key.height=unit(0.5,'cm'),
              legend.background=element_blank())

    data_plot_list = split(data_plot, data_plot[, split_col])

    data_plot_list = sapply(1:length(data_plot_list),
                            function(x){
                              data_plot_tmp = data_plot_list[[x]]
                              data_plot_tmp = data_plot_tmp[order(data_plot_tmp$avg.exp, decreasing = TRUE), ]
                              avg.exp_c = as.character(unique(data_plot_tmp$avg.exp))
                              col = rev(grDevices::colorRampPalette(c('white', my_pal_10[x]),
                                                                    bias=1)(n=length(avg.exp_c)))
                              names(col) = avg.exp_c
                              data_plot_tmp$col = col[as.character(data_plot_tmp$avg.exp)]
                              return(data_plot_tmp)
                            },
                            simplify = FALSE)

    data_plot = do.call("rbind", data_plot_list)

    dotplot <- ggplot(data_plot,
                      aes(x=id,y = features.plot,size = pct.exp))+
      geom_point(colour = data_plot$col) +
      scale_size("% detected", range = c(0,6)) +
      scale_y_discrete(position = "right") +
      scale_x_discrete(position = "top") +
      scale_color_gradientn(colours = my_mal_gradientn,
                            guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"),
                            name = "Average\nexpression") +
      theme_bw() +
      xlab("") + ylab("") +
      theme(
        axis.text.x.top = element_text(size=12, color="black",face="bold"),
        axis.ticks.x = element_blank(),

        axis.text.y = element_text(size=12, color="skyblue",face="bold"),
        axis.ticks.y = element_blank(),
        axis.text.y.right = element_text(size=8, angle=0, hjust=0.5, color="black",face="bold"),
        axis.line = element_line(colour = 'grey30',size = 0.2),
        panel.grid = element_blank()
      ) + ggpubr::rotate_x_text(angle = 45, hjust = 0, vjust = 0)

    dotplot%>%aplot::insert_left(ggtree_plot_col, width = 0.2)%>%
      aplot::insert_bottom(group3, height = 0.03)%>%
      aplot::insert_bottom(group2, height = 0.03)%>%
      aplot::insert_bottom(group1, height = 0.03)
  }else if (model==3){

    if (class(object)=="Seurat"){
      cells <- unlist(x = Seurat::CellsByIdentities(object = object, idents = NULL))

      data.features <- Seurat::FetchData(object = object, vars = features,
                                         cells = cells)
      data.features$id <- object[[ID[1], drop = TRUE]][cells, drop = TRUE]

      # "id" column must be in the end of data.features
      data_plot = avg_pct_func_(data.features)

      # row and col order
      mat <- data_plot %>% dplyr::select(-pct.exp) %>%
        tidyr::pivot_wider(names_from = features.plot, values_from = avg.exp ) %>%
        data.frame()

      row.names(mat) <- mat$id
      mat <- mat[,-1]
    }else {
      data_plot = object
      colnames(data_plot)[which(!colnames(data_plot)%in%c(ID, r_ID))] = "value"
      row = ID[1]
      col = r_ID[2]
      param_list = list(df=rlang::expr(data_plot),
                        row=as.symbol(ID[1]),
                        col=as.symbol(r_ID[1]),
                        value=rlang::expr(value))
      mat = do.call(tidyfst::df_mat, args = param_list)
      mat[is.na(mat)] = 0
    }

    clust_col <- hclust(dist(mat %>% as.matrix()%>% t()))
    ddgram_col <- ape::as.phylo(clust_col)

    clust_row <- hclust(dist(mat %>% as.matrix() ))
    ddgram_row <- ape::as.phylo(clust_row)

    ggtree_plot_col <- ggtree::ggtree(ddgram_col,branch.length="none")+
      geom_tippoint(color="#F7903D", shape=16, size=3)
    ggtree_plot_row <- ggtree::ggtree(ddgram_row,branch.length="none")+
      geom_tippoint(color="#4D85BD", shape=16, size=3)+
      layout_dendrogram()

    if (cluster_ID){
      id_level <-  clust_row$labels[clust_row$order]
    }else if(class(object)=="Seurat"&class(object@meta.data[,ID[1]])=="factor") {
      id_level = levels(object@meta.data[,ID[1]])
    }else {
      id_level <-  clust_row$labels[clust_row$order]
    }

    features_level<-  clust_col$labels[clust_col$order]

    if (class(object)=="Seurat"){

      if (isTRUE(cluster_feature)){
        data_plot$features.plot<-  factor(as.character(data_plot$features.plot) , levels = clust_col$labels[clust_col$order])
      }

      for (name in ID){
        data.features[, name] = object[[name, drop = TRUE]][cells, drop = TRUE]
      }

      data.features = data.features%>%tidyr::unite("id", ID[1]:ID[length(ID)],
                                                   sep = "_",
                                                   remove = FALSE)
      data.features = data.features[, 1:(ncol(data.features)-length(ID))]
      data_plot = avg_pct_func_(data.features)

      avg.exp.scaled <- sapply(X = unique(x = data_plot$features.plot),
                               FUN = function(x) {
                                 data.use <- data_plot[data_plot$features.plot ==
                                                         x, "avg.exp"]

                                 data.use <- scale(x = data.use)
                                 data.use <- MinMax(data = data.use, min = -2.5,
                                                    max = 2.5)
                                 return(data.use)
                               })
      avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
      data_plot$avg.exp.scaled <- avg.exp.scaled
      data_plot$pct.exp[data_plot$pct.exp < 0] <- NA
      data_plot$pct.exp <- data_plot$pct.exp * 100
      data_plot = data_plot%>%tidyr::separate(col = id,
                                              sep = "_",
                                              into = out_ID,
                                              remove = FALSE)

      data_plot[, ID[1]] <-  factor(data_plot[, ID[1]], levels = id_level)

      order_rows = lapply(ID, function(col_name){
        data_plot[[col_name]]
      })


      data_plot = data_plot[do.call(order, order_rows), ]
      data_plot$id = factor(data_plot$id, levels = unique(data_plot$id))


    }else {

      if (isTRUE(cluster_feature)){
        data_plot$features.plot<-  factor(as.character(data_plot[, r_ID[1]]) , levels = features_level)
      }else {
        data_plot$features.plot<-  factor(as.character(data_plot[, r_ID[1]]) , levels = features)
        data_plot = data_plot[order(data_plot$features.plot), ]
        data_plot_2 = data_plot[!duplicated(data_plot[, r_ID[1]]), ]
      }

      order_rows = lapply(ID[2:length(ID)], function(col_name){
        data_plot[[col_name]]
      })
      data_plot = data_plot[do.call(order, order_rows), ]
      data_plot$id = factor(data_plot[, ID[1]], levels = unique(data_plot[, ID[1]]))
      data_plot = data_plot[order(data_plot$id), ]
      data_plot_3 = data_plot[!duplicated(data_plot[, ID[1]]), ]
    }

    if (cluster_feature){
      level_row = clust_row$order
    }else {
      if (class(object)=="Seurat"){
        level_row = 1:length(clust_row$order)
      }else {
        level_row = 1:length(unique(data_plot[, ID[1]]))
      }
    }

    #scale for each split_col
    if (!is.null(split_col)){
      data_plot_list = split(data_plot, data_plot[, split_col])

      data_plot_list = sapply(1:length(data_plot_list),
                              function(x){
                                data_plot_tmp = data_plot_list[[x]]
                                data_plot_tmp = data_plot_tmp[order(data_plot_tmp$avg.exp, decreasing = TRUE), ]
                                avg.exp_c = as.character(unique(data_plot_tmp$avg.exp))
                                col = rev(grDevices::colorRampPalette(c('white', my_pal_10[x]),
                                                                      bias=1)(n=length(avg.exp_c)))
                                names(col) = avg.exp_c
                                data_plot_tmp$col = col[as.character(data_plot_tmp$avg.exp)]
                                return(data_plot_tmp)
                              },
                              simplify = FALSE)

      data_plot = do.call("rbind", data_plot_list)
    }

    if (plot_style=="bubble"){
      dotplot <- ggplot(data_plot,
                        aes(x=id,y = features.plot,size = pct.exp,, color = avg.exp.scaled))+
        geom_point() +
        scale_size("% detected", range = c(0,6)) +
        scale_y_discrete(position = "right") +
        scale_x_discrete(position = "top") +
        scale_color_gradientn(colours = my_mal_gradientn,
                              guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"),
                              name = "Average\nexpression") +
        theme_bw() +
        xlab("") + ylab("") +
        theme(
          axis.text.x.top = element_text(size=x_size, color="black",face="bold"),
          axis.ticks.x = element_blank(),

          axis.text.y = element_text(size=y_size, color="skyblue",face="bold"),
          axis.ticks.y = element_blank(),
          axis.text.y.right = element_text(size=y_size, angle=0, hjust=0.5, color="black",face="bold"),
          axis.line = element_line(colour = 'gray30',size = 0.2),
          panel.grid = element_blank()
        ) + ggpubr::rotate_x_text(angle = 45, hjust = 0, vjust = 0)
    }else if (plot_style=="heatmap"){
      dotplot = ggplot(data_plot, aes(x = features.plot, y = id))+
        geom_tile(aes(fill = value),size=white_line_size,
                  colour = white_line_color)+
        scale_fill_gradient2(low = "#2f5688",
                             high = "#CC0000",
                             mid = "white")+
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.background = element_rect(fill = alpha("blue",0.05), color = NA),
              #axis.line = element_line(color = "black", linewidth = 0.4),
              #axis.line.x = element_line(color = "black", linewidth = 0.4)
              legend.key = element_blank(),
              legend.text = element_text(colour = "black", size = 10),
              legend.background = element_blank(),
              legend.position = legend_position,
              legend.direction = legend_direction)

    }
    if ((cluster_ID&plot_style=="heatmap")|
        (cluster_feature&plot_style=="bubble")){
      dotplot = dotplot%>%aplot::insert_left(ggtree_plot_col, width = 0.2)
    }

    if (class(object)=="Seurat"){
      ID_new = ID
      data_plot_3 = data_plot
    }else {
      ID_new = ID[2:length(ID)]
      r_ID_new = r_ID[2:length(r_ID)]
    }
    for (name in rev(ID_new)){
      elem_num = length(unique(data_plot[, name]))
      if (elem_num==2){
        my_pal_tmp = my_pal_10
      }else if (elem_num==3){
        my_pal_tmp = my_pal_2
      }else {
        if (length(my_pal)<length(unique(data_plot[, name]))){
          if (is.null(my_pal)){
            my_pal = scales::hue_pal()(length(unique(data_plot[, name])))
          }else {
            my_pal = colorRampPalette(my_pal)(length(unique(data_plot[, name])))
          }
        }
        my_pal_tmp = my_pal[level_row]
      }

      group =
        ggplot(data_plot_3, aes(x = id, y = 1,
                                fill=!!as.symbol(name)))+
        geom_tile()+
        scale_y_discrete(expand = c(0,0),position="right")+
        scale_x_discrete(expand=c(0,0)) +
        scale_fill_manual(values=my_pal_tmp)+
        theme_void()+
        theme(axis.text.y=element_text(color="black",size=y_size))+
        theme(  legend.key=element_blank(),
                legend.text = element_text(color="black",size=10),
                legend.spacing.x=unit(0.1,'cm'),
                legend.key.width=unit(0.5,'cm'),
                legend.key.height=unit(0.5,'cm'),
                legend.background=element_blank())

      if (plot_style=="bubble"){
        dotplot = dotplot%>%
          aplot::insert_bottom(group, height = 0.03)
      }else if (plot_style=="heatmap") {
        group = group+coord_flip()+theme(axis.text.y = element_blank())
        if (name==ID_new[1]){
          group = group+coord_flip()+theme(axis.text.y = element_text(size = y_size),
                                           axis.ticks.y = element_line(color = "black", linewidth = y_linewidth),
                                           axis.ticks.length.y = unit(0.5, "mm"))
        }
        dotplot = aplot::insert_left(dotplot, group, width = 0.01)
      }

    }

    if (plot_style=="heatmap"){
      for (name in r_ID_new){
        elem_num = length(unique(data_plot[, name]))
        if (elem_num==2){
          my_pal_tmp = my_pal_10
        }else if (elem_num==3){
          my_pal_tmp = my_pal_2
        }else {
          if (length(my_pal)<length(unique(data_plot_2[, name]))){
            my_pal_tmp = scales::hue_pal()(length(unique(data_plot_2[, name])))
          }
          my_pal_tmp = my_pal[1:elem_num]
        }

        group =
          ggplot(data_plot_2, aes(x = features.plot, y = 1,
                                  fill=!!as.symbol(name)))+
          geom_tile()+
          scale_y_discrete(expand = c(0,0),position="right")+
          scale_x_discrete(expand=c(0,0)) +
          scale_fill_manual(values=my_pal_tmp)+
          theme_void()+
          theme(axis.text.y=element_text(color="black",size=y_size))+
          theme(  legend.key=element_blank(),
                  legend.text = element_text(color="black",size=10),
                  legend.spacing.x=unit(0.1,'cm'),
                  legend.key.width=unit(0.5,'cm'),
                  legend.key.height=unit(0.5,'cm'),
                  legend.background=element_blank())

        if (name == rev(r_ID_new)[1]){
          group = group+theme(axis.text.x = element_text(angle = 90,
                                                         vjust = 0,
                                                         hjust = 1,
                                                         size = x_size),
                              axis.ticks.x = element_line(color = "black", linewidth = x_linewidth),
                              axis.ticks.length.x = ggplot2::unit(0.5, "mm"))


        }
        dotplot = dotplot%>%aplot::insert_bottom(group, height = 0.03)
      }
      print(dotplot)
    }
    print(dotplot)
  }
}

avg_pct_func_ = function(data.features){
  data_plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident,
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = Seurat::PercentAbove,
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data_plot) <- unique(x = data.features$id)

  data_plot <- lapply(X = names(x = data_plot), FUN = function(x) {
    data.use <- as.data.frame(x = data_plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data_plot <- do.call(what = "rbind", args = data_plot)
}

#' Title plot_trajectory_DAG Implemented via the ggraph package
#'
#' @param Tab_l the frist elemnet must be a single root,
#' The list element can't be all = 0,
#' the names of list must be splited by ":"
#' @param scale_ "src", "obj" or NULL
#' @param order_ the order in y axis: from low tp high
#' @param cut_ The minimum transformation probability that satisfies the drawing
#' @param trans_link The value of links within same time; eg: E9_p,E9_l
#' @param node_color has been set in node_info
#'
#' @return
#' @export
#'
#' @examples


plot_trajectory_DAG = function(Tab_l,
                               Tab_l2,
                               scale_=NULL,
                               order_=NULL,
                               node_df=NULL,
                               color_by=NULL,
                               cellid_by=NULL,
                               cut_=0.1,
                               trans_link=0.3,
                               width_min=0.1,
                               width_max=0.3,
                               width_2=0.2,
                               width_top=0.5,
                               alpha_min=0.1,
                               x_size=3,
                               y_size=3,
                               node_size=2,
                               node_color="gray40",
                               edge_color){
  library(ggraph)
  library(igraph)
  library(aplot)

  order_n = 1:length(order_)
  names(order_n) = order_

  trajectory_DAG = trajectory_DAG_(Tab_l,
                                   cut_=cut_,
                                   trans_link=trans_link)
  Tab_ts = trajectory_DAG[["Tab_ts"]]
  time_ = trajectory_DAG[["time_"]]
  node_visited = trajectory_DAG[["node_visited"]]

  stopifnot(length(order_)==length(node_visited))

  if (!is.null(Tab_l2)){
    trajectory_DAG2 = trajectory_DAG_(Tab_l2,
                                      cut_=cut_,
                                      trans_link=trans_link)
    Tab_ts2 = trajectory_DAG2[["Tab_ts"]]

    Tab_ts = merge(Tab_ts, Tab_ts2,
                   by=c("src","obj"),
                   all.x=T)
    Tab_ts$new_column <- ifelse(!is.na(Tab_ts$value.y), 1, 0)

    Tab_ts <- Tab_ts %>%
      mutate(emphasize = new_column,
             width = ifelse(emphasize == 1, width_2,
                            width_min + (width_max - width_min) * value.x),
             color = ifelse(emphasize == 1, edge_color[2],
                            alpha(edge_color[1],alpha_min + (1 - alpha_min) * value.x)))

    Tab_ts <- dplyr::select(Tab_ts,
                            src=src,
                            obj=obj,
                            value=value.x,
                            emphasize=emphasize,
                            width=width,
                            color=color)
  }else {
    Tab_ts <- Tab_ts %>%
      mutate(width = width_min + (width_max - width_min) * value,
             color = alpha(edge_color[1],alpha_min + (1 - alpha_min) * value))

    Tab_ts <- dplyr::select(Tab_ts,
                            src=src,
                            obj=obj,
                            value=value,
                            width=width,
                            color=color)
  }

  time_num = rev(1:(length(time_)))
  names(time_num) = time_

  node_ll_names = unique(c(Tab_ts$src, Tab_ts$obj))
  node_ll = 1:length(node_ll_names)
  names(node_ll) = node_ll_names
  Tab_ts$from = unname(node_ll[Tab_ts$src])
  Tab_ts$to = unname(node_ll[Tab_ts$obj])

  node_info = data.frame(id=unname(node_ll),
                         label=names(node_ll))
  node_info$id = as.numeric(node_info$id)

  if (!is.null(node_df)){
    colnames(node_df)[which(colnames(node_df)==cellid_by)] = "cellid"
    colnames(node_df)[which(colnames(node_df)==color_by)] = "color"
    node_info$label_old = sub(".*:", "", node_info$label)
    node_info <- node_info %>%
      dplyr::left_join(node_df,by = c("label_old"="cellid"))
  }

  coords_x = unname(order_n[str_split(node_info$label,":",simplify = T)[,2]])
  coords_y = unname(time_num[str_split(node_info$label,":",simplify = T)[,1]])

  Tab_ts <- Tab_ts %>%
    mutate(from = unname(node_ll[src]),
           to = unname(node_ll[obj])) %>%
    dplyr::select(from, to, value, width, color)


  #Tab_ts = Tab_ts[, c("from", "to", "value")]
  #Tab_ts$from = as.numeric(Tab_ts$from)
  #Tab_ts$to = as.numeric(Tab_ts$to)

  tree <- graph_from_data_frame(Tab_ts, vertices = node_info)

  p = ggraph(tree, layout = 'manual',
             x=coords_x,
             y=coords_y,
             circular = FALSE) +
    geom_edge_diagonal(aes(start_cap = label_rect(node1.name),
                           end_cap = label_rect(node2.name),
                           color=color,
                           width=width))+
    geom_node_point(aes(color=color),size=node_size)+
    scale_edge_alpha_identity()+
    scale_edge_color_identity()+
    scale_color_identity()+
    scale_y_reverse()+
    coord_flip()+theme_void()+
    theme(legend.position = "top")

  #add x_axis development time
  p = p+annotate("text", y = unique(unname(time_num)),
                 x = max(coords_x)+1,
                 label=c(names(time_num)),
                 hjust=1,
                 size=x_size)
  #add y_axis cell states
  p = p+annotate("text", x = sort(unique(coords_x)),
                 y = 0,
                 label=order_,
                 hjust=0.0,
                 size=y_size)
  #Add inter- and intra-timing plotting layers

  for (i in rev(seq_along(unique(coords_y)))) {
    p = p + annotate("rect", xmin = max(coords_x)+2, xmax = max(coords_x)+2+width_top,
                     ymin = unique(coords_y)[i] - 0.5, ymax = unique(coords_y)[i] + 0.5,
                     fill = ifelse(i %% 2 == 0, "gray60", "gray90"), alpha = 0.6)
  }
  print(p)
}

#' Title a wrap function to prepare transmit net data for plot_trajectory_DAG
#'
#' @param Tab_l samilar to plot_trajectory_DAG
#'
#' @return
#' @export
#'
#' @examples
trajectory_DAG_ = function(Tab_l,cut_,trans_link){
  Tab_ts = data.frame()
  node_visited = c()
  node_old = c()
  time_ = c()
  for (name in names(Tab_l)){
    Tab_t = Tab_l[[name]]
    Tab_t[Tab_t<cut_] = 0
    src_ = rownames(Tab_t)
    obj_f = colnames(Tab_t)[1]
    obj_l = colnames(Tab_t)[ncol(Tab_t)]
    Tab_t = rownames_to_column(Tab_t, var = "src")

    if (obj_f==obj_l){
      cols_expr=rlang::expr(!!rlang::sym(obj_f))
    }else {
      cols_expr= rlang::expr(!!rlang::sym(obj_f):!!rlang::sym(obj_l))
    }

    Tab_t = pivot_longer(data = Tab_t,
                         cols = !!cols_expr,
                         names_to = c("obj"),
                         values_to = "value")
    Tab_t = Tab_t[Tab_t$value>0, ]

    if (isTRUE(grepl(":", name))){
      time_p = str_split(name, ":", simplify = T)[1]
      time_p = str_c(time_p,"n",sep = "_")
      time_n = str_split(name, ":", simplify = T)[2]
      time_n = str_c(time_n,"p",sep = "_")
      Tab_t$src = str_c(time_p,Tab_t$src,sep = ":")
      Tab_t$obj = str_c(time_n,Tab_t$obj,sep = ":")

    }else {
      if(!is.null(node_visited)){
        src_=node_old
      }
      time_p = str_c(name, "p", sep = "_")
      time_n = str_c(name, "n", sep = "_")
      Tab_n = data.frame(src=str_c(time_p,src_,sep = ":"),
                         obj=str_c(time_n,src_,sep = ":"),
                         value=rep(trans_link,length(src_)))

      Tab_t$src = str_c(time_p,Tab_t$src,sep = ":")
      Tab_t$obj = str_c(time_n,Tab_t$obj,sep = ":")
      Tab_t = rbind(Tab_n, Tab_t)
    }

    time_ = unique(c(time_, setdiff(c(time_p, time_n),time_)))
    Tab_ts = rbind(Tab_ts, Tab_t)
    node_old = unique(sub(".*:", "", Tab_t$obj))
    node_visited  = unique(c(node_visited, node_old))
  }
  return(list("Tab_ts"=Tab_ts,
              "time_"=time_,
              "node_visited"=node_visited))
}
