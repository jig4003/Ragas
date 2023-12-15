#' Matrix Plot
#'
#' An internal function to make matrix plot based on Seurat object.
#'
#' @param object A Seurat object
#' @param markers.list A data.frame that contains marker data from \code{\link[Seurat]{FindAllMarkers}}
#' @param n.genes Number of genes to plot for each identity group
#' @param type The type of markers, either "all" or "conserved", which correspond results from  \code{\link[Seurat]{FindAllMarkers}}
#' and \code{\link[Seurat]{FindConservedMarkers}}, respectively. In the current version, only \code{FindAllMarkers} is supported (default: "all")
#' @param up.genes Whether to plot the up-regulated genes (default: TRUE)
#' @param down.genes Whether to plot the down-regulated genes (default: FALSE)
#' @param heatmap.cols Colors to plot the heatmap (optional). A character vector with colors for low, medium and high expression.
#' @param min.exp The minimum expression value to plot (default: -2)
#' @param max.exp The maximum expression value to plot (default: 2)
#' @param column.fontsize Font size for column texts (default: 12)
#' @param row.fontsize Font size for row texts (default: 12)
#' @param column.anno.cols Named character vector with colors for column annotations (default: NULL). Names should match Idents of object.
#' @param column.anno.name.fontsize Font size for column annotations, such as clusters (default: 15)
#' @param column.anno.name.rot Rotation of column annotations (default: 45)
#' @param column.anno.name.just Alignment of column annotations (default: "center")
#' @param legend.label.fontsize Font size for legend label (default: 13)
#' @param legend.title.fontsize Font size for legend title (15)
#' @param heatmap.width Width of the heatmap (default: 15)
#' @param heatmap.height Height of the heatmap (default: 10)
#'
#' @importFrom dplyr %>% summarise group_by
#' @importFrom reshape2 melt dcast
#' @importFrom scales hue_pal
#' @importFrom ComplexHeatmap HeatmapAnnotation anno_text Heatmap
#' @importFrom grid gpar unit
#' @importFrom circlize colorRamp2
#' @importFrom gplots colorpanel
#'
#' @return A \code{\link[ComplexHeatmap]{Heatmap-class}} object
#'

MatrixPlot <- function(object, markers.list, n.genes = 3, type = "all", up.genes = TRUE, down.genes = FALSE,
                       heatmap.cols = NULL,
                       min.exp=-2,max.exp=2,column.fontsize=12,row.fontsize=12,
                       column.anno.cols = NULL,
                       column.anno.name.fontsize=15,column.anno.name.rot=45,column.anno.name.just="center",
                       legend.label.fontsize=13,legend.title.fontsize=15,
                       heatmap.width=15,heatmap.height=10)
{

  if(type == "conserved")
  {
    #retrive top N genes
    top.markers.list <- list() ; cluster2remove = NULL
    for(i in 1:length(markers.list))
    {
      idx.lfc <- grep("_avg_log2FC",colnames(markers.list[[i]]))
      markers.list[[i]] <- cbind(markers.list[[i]], avg.log2FC.all = apply(markers.list[[i]][,idx.lfc],1,mean))
      #markers.sig <- markers.list[[i]][which(markers.list[[i]]$minimump_p_val < 0.05),]
      markers.sig <- markers.list[[i]]
      markers.sig <- markers.sig[order(markers.sig$avg.log2FC.all,decreasing = T),]
      if(up.genes==TRUE & down.genes ==FALSE)
      {
        markers.sig.1 <- markers.sig[which(markers.sig$avg.log2FC.all > 0),]
        if(nrow(markers.sig.1)!=0)
        {
          condition <- nrow(markers.sig.1) > n.genes
          top.markers.list[[i]] <- switch(2-condition,rownames(markers.sig.1)[1:n.genes],rownames(markers.sig.1))
        }else
        {
          cluster2remove <- c(cluster2remove,names(markers.list)[i])
        }

      }else if(up.genes==TRUE & down.genes ==TRUE)
      {
        stop("Error: Cannot set both up.genes and down.genes to TRUE")
      }else if(up.genes==FALSE & down.genes ==TRUE)
      {
        markers.sig.1 <- markers.sig[which(markers.sig$avg.log2FC.all < 0),]
        if(nrow(markers.sig.1)!=0)
        {
          markers.sig.1 <- markers.sig.1[order(abs(markers.sig.1$avg.log2FC.all),decreasing = TRUE),]
          condition <- nrow(markers.sig.1) > n.genes
          top.markers.list[[i]] <- switch(2-condition,rownames(markers.sig.1)[1:n.genes],rownames(markers.sig.1))
        }else
        {
          cluster2remove <- c(cluster2remove,names(markers.list)[i])
        }

      }
    }

    if(!is.null(cluster2remove))
    {
      idx.cluster2remove <- match(cluster2remove,names(markers.list))
      top.markers.list <- top.markers.list[-idx.cluster2remove]
      names(top.markers.list) <- names(markers.list)[-idx.cluster2remove]
    }else
    {
      names(top.markers.list) <- names(markers.list)
    }


    dat.assay <- GetAssayData(object,slot = "scale.data")
    seurat_clusters <- Idents(object)

    top.markers.avg.exp <- list()
    for(i in 1:length(top.markers.list))
    {

      exp.dat <- dat.assay[match(top.markers.list[[i]],rownames(dat.assay)),]
      if(is.null(nrow(exp.dat)))
      {
        exp.dat <- t(as.matrix(exp.dat))
        rownames(exp.dat) <- top.markers.list[[i]]
      }
      exp.dat.1 <- data.frame(seurat_clusters=rep(seurat_clusters,each=nrow(exp.dat)),
                              exp=melt(exp.dat))
      top.markers.avg.exp[[i]] <-  cbind(Cluster=names(top.markers.list)[i],
                                         as.data.frame(exp.dat.1 %>% group_by(seurat_clusters,exp.Var1) %>%
                                                         summarise(avg.exp = mean(exp.value),.groups = 'drop')))

    }
    names(top.markers.avg.exp) <- names(top.markers.list)

    top.markers.avg.exp.1 <- do.call("rbind",top.markers.avg.exp)
    if(!is.null(cluster2remove))
    {
      idx <- which(top.markers.avg.exp.1$seurat_clusters==gsub("cluster","",cluster2remove))
      top.markers.avg.exp.1 <- top.markers.avg.exp.1[-idx,]
    }
    top.markers.avg.exp.1$seurat_clusters <- droplevels(top.markers.avg.exp.1$seurat_clusters)

  }else if(type == "all")
  {
    markers.clusters <- levels(markers.list$cluster)
    #retrive top N genes
    top.markers.list <- list() ; cluster2remove = NULL
    for(i in 1:length(markers.clusters))
    {
      markers.sig <- markers.list[which(markers.list$cluster==markers.clusters[i]),]
      markers.sig <- markers.sig[order(markers.sig$avg_log2FC,decreasing = T),]

      if(up.genes==TRUE & down.genes ==FALSE)
      {
        markers.sig.1 <- markers.sig[which(markers.sig$avg_log2FC > 0),]
        if(nrow(markers.sig.1)!=0)
        {
          condition <- nrow(markers.sig.1) > n.genes
          top.markers.list[[i]] <- switch(2-condition,markers.sig.1$gene[1:n.genes],markers.sig.1$gene)
        }else
        {
          cluster2remove <- c(cluster2remove,markers.clusters[i])
        }

      }else if(up.genes==TRUE & down.genes ==TRUE)
      {
        stop("Error: Cannot set both up.genes and down.genes to TRUE")
      }else if(up.genes==FALSE & down.genes ==TRUE)
      {
        markers.sig.1 <- markers.sig[which(markers.sig$avg_log2FC < 0),]
        if(nrow(markers.sig.1)!=0)
        {
          markers.sig.1 <- markers.sig.1[order(abs(markers.sig.1$avg_log2FC),decreasing = TRUE),]
          condition <- nrow(markers.sig.1) > n.genes
          top.markers.list[[i]] <- switch(2-condition,markers.sig.1$gene[1:n.genes],markers.sig.1$gene)
        }else
        {
          cluster2remove <- c(cluster2remove,markers.clusters[i])
        }

      }
    }

    if(!is.null(cluster2remove))
    {
      idx.cluster2remove <- match(cluster2remove,markers.clusters)
      top.markers.list <- top.markers.list[-idx.cluster2remove]
      names(top.markers.list) <- paste("cluster",markers.clusters[-idx.cluster2remove],sep="")
    }else
    {
      names(top.markers.list) <- paste("cluster",markers.clusters,sep="")
    }

    dat.assay <- GetAssayData(object,slot = "scale.data")
    #seurat_clusters <- object$seurat_clusters
    seurat_clusters <- Idents(object)

    top.markers.avg.exp <- list()
    for(i in 1:length(top.markers.list))
    {

      exp.dat <- dat.assay[match(top.markers.list[[i]],rownames(dat.assay)),]
      if(is.null(nrow(exp.dat)))
      {
        exp.dat <- t(as.matrix(exp.dat))
        rownames(exp.dat) <- top.markers.list[[i]]
      }
      exp.dat.1 <- data.frame(seurat_clusters=rep(seurat_clusters,each=nrow(exp.dat)),
                              exp=melt(exp.dat))
      top.markers.avg.exp[[i]] <-  cbind(Cluster=names(top.markers.list)[i],
                                         as.data.frame(exp.dat.1 %>% group_by(seurat_clusters,exp.Var1) %>%
                                                         summarise(avg.exp = mean(exp.value),.groups = 'drop')))

    }
    names(top.markers.avg.exp) <- names(top.markers.list)

    top.markers.avg.exp.1 <- do.call("rbind",top.markers.avg.exp)
    if(!is.null(cluster2remove))
    {
      idx <- which(top.markers.avg.exp.1$seurat_clusters==gsub("cluster","",cluster2remove))
      top.markers.avg.exp.1 <- top.markers.avg.exp.1[-idx,]
    }
    top.markers.avg.exp.1$seurat_clusters <- droplevels(top.markers.avg.exp.1$seurat_clusters)

  }

  plot.dat <- dcast(top.markers.avg.exp.1,Cluster+exp.Var1 ~ seurat_clusters,value.var = "avg.exp")
  clusters <- plot.dat$Cluster
  if(sum(clusters %in% c('0','1')) > 0)
  {
    clusters <- as.numeric(gsub("cluster","",clusters))
    plot.dat <- plot.dat[order(clusters,decreasing = F),]
    clusters <- clusters[order(clusters, decreasing = F)]
  }else
  {
    clusters <- gsub("cluster","",clusters)
  }

  #order the markers by log fold change
  plot.dat.new <- NULL
  for(i in 1:length(top.markers.list))
  {
    tmp.dat <- plot.dat[which(plot.dat$Cluster==names(top.markers.list)[i]),]
    idx <- match(top.markers.list[[i]],tmp.dat$exp.Var1,)
    tmp.dat <- tmp.dat[idx,]
    plot.dat.new <- rbind(plot.dat.new,tmp.dat)
  }


  plot.dat.1 <- t(plot.dat.new[,-c(1:2)])
  colnames(plot.dat.1) <- plot.dat.new$exp.Var1


  dup.idx <- lapply(colnames(plot.dat.1),function(x){which(colnames(plot.dat.1)==x)})
  names(dup.idx) <- colnames(plot.dat.1)

  dup.times <- lapply(dup.idx,function(x){length(x)})
  dup.names <- list(names(dup.times),unlist(dup.times))
  dup.names <- unlist(mapply(FUN = function(x,y){paste(x ,c(0:(y-1)),sep=".")},dup.names[[1]],dup.names[[2]]))
  dup.names <- gsub(".0$","",dup.names) ; names(dup.names) <- NULL


  colnames(plot.dat.1)[unlist(dup.idx)] <- dup.names
  plot.dat.new$name <- colnames(plot.dat.1)


  h.clust <- hclust(dist(plot.dat.1 %>% as.matrix())) # hclust with distance matrix
  ddgram.row <- as.dendrogram(h.clust)


  idx <- unlist(lapply(paste("cluster",h.clust$labels[h.clust$order],sep=""),function(x){which(plot.dat.new$Cluster==x)}))
  plot.dat.new <- plot.dat.new[idx,]
  plot.dat.2 <- plot.dat.1[,plot.dat.new$name]


  clusters <- plot.dat.new$Cluster
  if(sum(clusters %in% c('0','1')) > 0)
  {
    clusters <- as.numeric(gsub("cluster","",clusters))
  }else
  {
    clusters <- gsub("cluster","",clusters)
  }
  anno.col <- data.frame(Clusters=factor(clusters,levels = unique(clusters)))

  if(is.null(column.anno.cols))
  {
    cols.clusters <- hue_pal()(length(levels(anno.col$Clusters)))
    names(cols.clusters) <- levels(anno.col$Clusters)
  }else
  {
    if(length(column.anno.cols) < length(levels(anno.col$Clusters)))
    {
      print('Length of column.anno.cols should match the number of clusters')
      stop()
    }else
    {
      cols.clusters <- column.anno.cols
      idx <- match(names(cols.clusters), anno.col$Clusters)
      if(sum(is.na(idx)) >=1)
      {
        stop("Names of \"column.anno.cols\" should match Idents")
      }else
      {
        cols.clusters <- cols.clusters[levels(anno.col$Clusters)]
      }
      #names(cols.clusters) <- levels(anno.col$Clusters)
    }
  }


  clus.name <- anno.col$Clusters
  idents.seurat.obj <- levels(Idents(object))
  idents.seurat.obj <- idents.seurat.obj[match(levels(clus.name),idents.seurat.obj)]

  clus.name <- factor(clus.name,labels = idents.seurat.obj)
  clust.size <-NULL
  for(clus in 1:length(levels(clus.name)))
  {
    idx <- which(clus.name==levels(clus.name)[clus])
    names(idx) <- c(1:length(idx))
    size <- rep(0,length(idx))
    size[floor(mean(as.numeric(names(idx))))] = column.anno.name.fontsize
    clust.size <- c(clust.size,size)
  }

  top.anno <- HeatmapAnnotation(foo = anno_text(clus.name, gp = gpar(fontsize = clust.size),
                                                rot=column.anno.name.rot,location = unit(0,"npc"),just=column.anno.name.just),
                                df=anno.col,
                                col = list(Clusters=cols.clusters),
                                annotation_height = unit(c(0.5, 0.5), "cm"),
                                show_annotation_name = F,
                                annotation_name_offset = unit(5, "mm"),
                                annotation_name_rot = c(0, 0),
                                gp = gpar(col = "white"),
                                gap = unit(5, "points"),
                                annotation_legend_param = list(labels_gp = gpar(fontsize=legend.label.fontsize),title_gp = gpar(fontsize=legend.title.fontsize,fontface="bold")))

  if(is.null(heatmap.cols))
  {
    heatmap.cols <- c('blue','white','red')
  }

  p <- Heatmap(plot.dat.2,cluster_rows = ddgram.row, cluster_columns = F,top_annotation = top.anno,show_column_names = T,name="average\nexpression",
               row_dend_side = "right", row_names_side = "left",rect_gp = gpar(col = "black"), height =  unit(0.5, "npc"),
               column_labels = plot.dat.new$exp.Var1,column_names_gp = gpar(fontsize=column.fontsize),row_names_gp = gpar(fontsize = row.fontsize),
               heatmap_width = unit(heatmap.width, "cm"),heatmap_height = unit(heatmap.height,"cm"),
               colorRamp2(breaks=seq(min.exp,max.exp,length.out=21),colors = colorpanel(21, low=heatmap.cols[1], mid=heatmap.cols[2], high=heatmap.cols[3])),
               heatmap_legend_param = list(at= seq(min.exp,max.exp),legend_height = unit(5,"cm"),
                                           labels_gp = gpar(fontsize=legend.label.fontsize),title_gp = gpar(fontsize=legend.title.fontsize,fontface="bold")))

  return(p)
}
