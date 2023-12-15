#' Summarized Heatmap for Seurat or post-integration (pi) object
#'
#' An internal function to plot heatmap that groups by identity class or additional meta data. Use its wrapper function \code{\link[Ragas]{RunSummarizedHeatmap}} for plotting.
#'
#' SummarizedHeatmap takes data from the data slot of the Seurat object for the given features and split the cells
#' and summarize their average expression per identity (or by identity and additional metadata as determined by the split.by argument)
#' and scaled across all columns.
#'
#' @param object A Seurat or pi object
#' @param features A character vector of features
#' @param assay Name of the assay to use (default: "RNA")
#' @param split.by Factor to split the columns by
#' @param heatmap.cols Colors for filling the heatmap
#' @param row.annotation,column.annotation Data.frame for row and column annotation
#' @param clust.row,clust.column Whether to cluster the rows/columns (default: TRUE)
#' @param row.annotation.cols,column.annotation.cols Colors for row/Column annotations; If set to NULL, \code{link[randomcoloR]{distinctColorPalette}} will be used to generate random colors
#' @param show.row.names,show.column.names Whether to show row/column names (default: TRUE)
#' @param min.exp,max.exp Minimum/maximum scaled expression value to plot (default: -2/2)
#' @param column.names.rotation Rotation for column names (default: 90)
#' @param row.fontsize,column.fontsize Row/column font size (default: 10/8)
#' @param legend.label.fontsize,legend.title.fontsize,annotation.name.fontsize Font sizes for legend label/title and annotation name (default: 10)
#' @param heatmap.width,heatmap.height Width and height of heatmap (default: 12/10)
#' @param random.col.seed Set seed to control row and column annotation colors (default: 1)
#' @param ... Extra parameters passed to \code{\link[ComplexHeatmap]{Heatmap}}
#'
#' @return A \code{\link[ComplexHeatmap]{Heatmap-class}} object
#'
#' @references Hao and Hao et al. (2021). “Integrated analysis of multimodal single-cell data.” Cell. doi:10.1016/j.cell.2021.04.048
#' @references Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics.
#'
#' @importFrom dplyr %>% group_by
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom ComplexHeatmap HeatmapAnnotation Heatmap
#' @importFrom grid gpar unit
#' @importFrom circlize colorRamp2
#' @importFrom gplots colorpanel
#'
#' @keywords internal



SummarizedHeatmap <- function(object, features, assay = 'RNA', split.by = NULL,
                              heatmap.cols = NULL,
                              clust.row = TRUE,clust.column = FALSE,
                              row.annotation = NULL,row.annotation.cols = NULL,
                              column.annotation = NULL,  column.annotation.cols = NULL,
                              show.row.names = TRUE, show.column.names = TRUE,
                              min.exp=-2,max.exp=2,column.names.rotation=45,
                              column.fontsize=12,row.fontsize=12,
                              legend.label.fontsize=13,legend.title.fontsize=15,
                              annotation.name.fontsize=8, heatmap.width=12,heatmap.height=14,random.col.seed = 1,
                              ...)
{

  # Check if features are unique
  if(length(unique(features)) != length(features))
  {
    print('Features are not unique!')
    stop()
  }

  all.features <- features
  idx <- match(features, rownames(object))
  if(sum(is.na(idx))>0)
  {
    print(paste("Warning: Features",paste(features[which(is.na(idx))],collapse = ", "),"not found! Omitting these features!!"))
    features <- features[which(!is.na(idx))]
    idx <- match(features,rownames(object))
  }

  dat <- GetAssayData(object, slot="data", assay = assay)
  if(length(idx) > 1)
  {
    dat <- as.data.frame(dat[idx,])
  }else if(length(idx) == 1)
  {
    dat <- t(as.data.frame(dat[idx,]))
    rownames(dat) <- features
  }else if(length(idx) == 0)
  {
    print(paste('Need at least 1 valid feature to plot heatmap!'))
    stop()
  }

  seurat.clusters <- Idents(object)

  if(is.null(split.by))
  {
    dat.1 <- suppressMessages(data.frame(seurat.clusters=rep(seurat.clusters,each=nrow(dat)),
                                         Gene = rep(rownames(dat),ncol(dat)),
                                         exp=reshape2::melt(dat)))

    dat.2 <- as.data.frame(dat.1 %>% group_by(seurat.clusters,Gene) %>%
                             dplyr::summarise(avg.exp = mean(exp.value),.groups = 'drop'))


    plot.dat <- reshape2::dcast(dat.2,Gene ~ seurat.clusters,value.var = "avg.exp")
    rownames(plot.dat) <- plot.dat$Gene
    plot.dat <- plot.dat[,-match("Gene",colnames(plot.dat))]
    plot.dat <- plot.dat[match(features,rownames(plot.dat)),]
    plot.dat <- as.matrix(plot.dat)

    idx <- which(rowSums(plot.dat)==0)
    if(length(idx)!=0)
    {
      r.names <- rownames(plot.dat)[-idx]
      plot.dat <- plot.dat[-idx,]
      if(!is.matrix(plot.dat))
      {
        plot.dat <- t(as.data.frame(plot.dat))
        rownames(plot.dat) <- r.names
      }
    }

    if(nrow(plot.dat) == 0)
    {
      stop('Features do not have any expression!')
    }

    plot.dat.z <- (plot.dat - apply(plot.dat,1,mean))/apply(plot.dat,1,sd)

    if(!is.null(column.annotation))
    {
      col.anno.names <- colnames(column.annotation)
      col.anno.idx <- match(colnames(plot.dat.z),column.annotation[,1])
      column.annotation <- data.frame(column.annotation[col.anno.idx,])
      colnames(column.annotation) <- col.anno.names
      column.annotation$Cluster <- factor(column.annotation$Cluster, levels = levels(Idents(object)))

      if(is.null(column.annotation.cols))
      {
        set.seed(random.col.seed)
        column.annotation.cols <- distinctColorPalette(length(unique(column.annotation[,1])))
        names(column.annotation.cols) <- unique(column.annotation[,1])
        column.annotation.cols <- list(column.annotation.cols)
        names(column.annotation.cols) <- colnames(column.annotation)[1]
      }

      top.anno <- HeatmapAnnotation(df=column.annotation,
                                    col = column.annotation.cols,
                                    annotation_height = unit(rep(0.5, ncol(column.annotation)), "cm"),
                                    show_annotation_name = T,
                                    annotation_name_offset = unit(1, "mm"),
                                    annotation_name_gp = gpar(fontface ="bold",fontsize=annotation.name.fontsize),
                                    annotation_name_rot = rep(0, ncol(column.annotation)),
                                    gp = gpar(col = "white"),
                                    annotation_legend_param = list(labels_gp = gpar(fontsize=legend.label.fontsize),
                                                                   title_gp = gpar(fontsize=legend.title.fontsize,fontface="bold")))
    }else if(is.null(column.annotation))
    {
      top.anno <- NULL
    }

    if(!is.null(row.annotation))
    {
      row.anno.names <- colnames(row.annotation)
      row.anno.idx <- match(rownames(plot.dat.z),all.features)
      row.annotation <- data.frame(row.annotation[row.anno.idx,])
      colnames(row.annotation) <- row.anno.names

      if(is.null(row.annotation.cols))
      {
        set.seed(random.col.seed)
        row.annotation.cols <- distinctColorPalette(length(unique(row.annotation[,1])))
        names(row.annotation.cols) <- unique(row.annotation[,1])
        row.annotation.cols <- list(row.annotation.cols)
        names(row.annotation.cols) <- colnames(row.annotation)[1]
      }

      left.anno <- HeatmapAnnotation(df = row.annotation,
                                     col = row.annotation.cols,
                                     annotation_height = unit(0.5, "cm"),
                                     show_annotation_name = F,
                                     annotation_name_offset = unit(1, "mm"),
                                     annotation_name_gp = gpar(fontsize=0),
                                     annotation_name_rot = 90,
                                     gp = gpar(col = "white"),which = "row",
                                     annotation_legend_param = list(labels_gp = gpar(fontsize=legend.label.fontsize),
                                                                    title = 'Row\nannotation',
                                                                    title_gp = gpar(fontsize=legend.title.fontsize,fontface="bold")))
    }else if(is.null(row.annotation))
    {
      left.anno <- NULL
    }

    if(is.null(heatmap.cols))
    {
      heatmap.cols <- c('blue','white','red')
    }

    p <- Heatmap(plot.dat.z, cluster_rows = clust.row, cluster_columns = clust.column,
                 left_annotation = left.anno,top_annotation = top.anno,
                 show_row_names = show.row.names, show_column_names = show.column.names, name="Average\nexpression",
                 row_dend_side = "right", row_names_side = "right",rect_gp = gpar(col = "black"),
                 column_names_gp = gpar(fontsize=column.fontsize),column_names_rot = column.names.rotation,
                 row_names_gp = gpar(fontsize = row.fontsize),
                 heatmap_width = unit(heatmap.width, "cm"),heatmap_height = unit(heatmap.height,"cm"),
                 colorRamp2(breaks=seq(min.exp,max.exp,length.out=21),colors = colorpanel(21, low=heatmap.cols[1], mid=heatmap.cols[2], high=heatmap.cols[3])),
                 heatmap_legend_param = list(at= seq(min.exp,max.exp),legend_height = unit(5,"cm"),
                                             labels_gp = gpar(fontsize=legend.label.fontsize),
                                             title_gp = gpar(fontsize=legend.title.fontsize,fontface="bold")),...)

  }else if(!is.null(split.by))
  {
    group.var <- split.by
    group <- object[[group.var]][,1]

    dat.1 <- suppressMessages(data.frame(seurat.clusters=rep(seurat.clusters,each=nrow(dat)),group = rep(group,each=nrow(dat)),
                                         Gene = rep(rownames(dat),ncol(dat)),
                                         exp=reshape2::melt(dat)))

    dat.2 <- as.data.frame(dat.1 %>% group_by(seurat.clusters,group,Gene) %>%
                             dplyr::summarise(avg.exp = mean(exp.value),.groups = 'drop'))

    plot.dat <- reshape2::dcast(dat.2,Gene ~ seurat.clusters+group,value.var = "avg.exp")
    rownames(plot.dat) <- plot.dat$Gene
    plot.dat <- plot.dat[,-match("Gene",colnames(plot.dat))]
    plot.dat <- plot.dat[match(features,rownames(plot.dat)),]
    plot.dat <- as.matrix(plot.dat)

    idx <- which(rowSums(plot.dat)==0)
    if(length(idx)!=0)
    {
      r.names <- rownames(plot.dat)[-idx]
      plot.dat <- plot.dat[-idx,]
      if(!is.matrix(plot.dat))
      {
        plot.dat <- t(as.data.frame(plot.dat))
        rownames(plot.dat) <- r.names
      }
    }


    if(nrow(plot.dat) == 0)
    {
      stop('Features do not have any expression!')
    }

    plot.dat.z <- (plot.dat - apply(plot.dat,1,mean))/apply(plot.dat,1,sd)


    if(!is.null(column.annotation))
    {
      if(ncol(column.annotation) == 2)
      {
        col.anno.names <- colnames(column.annotation)
        col.anno.idx <- match(colnames(plot.dat.z),apply(column.annotation,1,function(x){paste(x,collapse = '_')}))
        if(is.na(sum(col.anno.idx)))
        {
          col.anno.idx <- match(colnames(plot.dat.z),apply(column.annotation,1,function(x){paste(rev(x),collapse = '_')}))
        }

        column.annotation <- data.frame(column.annotation[col.anno.idx,])
        colnames(column.annotation) <- col.anno.names
        column.annotation$Cluster <- factor(column.annotation$Cluster, levels = levels(Idents(object)))
      }else if(ncol(column.annotation) >2)
      {
        temp <- column.annotation[,match(c('Cluster',split.by),colnames(column.annotation))]
        col.anno.names.orig <- colnames(column.annotation)
        col.anno.names <- colnames(temp)
        col.anno.idx <- match(colnames(plot.dat.z),apply(temp,1,function(x){paste(x,collapse = '_')}))
        if(is.na(sum(col.anno.idx)))
        {
          col.anno.idx <- match(colnames(plot.dat.z),apply(temp,1,function(x){paste(rev(x),collapse = '_')}))
        }

        column.annotation <- data.frame(column.annotation[col.anno.idx,])
        colnames(column.annotation) <- col.anno.names.orig
        column.annotation$Cluster <- factor(column.annotation$Cluster, levels = levels(Idents(object)))
      }



      if(is.null(column.annotation.cols))
      {
        set.seed(random.col.seed)
        column.annotation.cols <- list()
        for(i in 1:ncol(column.annotation))
        {
          cols <- distinctColorPalette(length(unique(column.annotation[,i])))
          names(cols) <- unique(column.annotation[,i])
          column.annotation.cols[[i]] <- cols
        }
        names(column.annotation.cols) <- colnames(column.annotation)
      }

      top.anno <- HeatmapAnnotation(df=column.annotation,
                                    col = column.annotation.cols,
                                    annotation_height = unit(rep(0.5, ncol(column.annotation)), "cm"),
                                    show_annotation_name = T,
                                    annotation_name_offset = unit(1, "mm"),
                                    annotation_name_gp = gpar(fontface ="bold",fontsize=annotation.name.fontsize),
                                    annotation_name_rot = rep(0, ncol(column.annotation)),
                                    gp = gpar(col = "white"),
                                    annotation_legend_param = list(labels_gp = gpar(fontsize=legend.label.fontsize),title_gp = gpar(fontsize=legend.title.fontsize,fontface="bold")))
    }else
    {
      top.anno <- NULL
    }


    if(!is.null(row.annotation))
    {
      row.anno.names <- colnames(row.annotation)
      row.anno.idx <- match(rownames(plot.dat.z),all.features)
      row.annotation <- data.frame(row.annotation[row.anno.idx,])
      colnames(row.annotation) <- row.anno.names

      if(is.null(row.annotation.cols))
      {
        set.seed(random.col.seed)
        row.annotation.cols <- distinctColorPalette(length(unique(row.annotation[,1])))
        names(row.annotation.cols) <- unique(row.annotation[,1])
        row.annotation.cols <- list(row.annotation.cols)
        names(row.annotation.cols) <- colnames(row.annotation)[1]
      }

      left.anno <- HeatmapAnnotation(df = row.annotation,
                                     col = row.annotation.cols,
                                     annotation_height = unit(0.5, "cm"),
                                     show_annotation_name = F,
                                     annotation_name_offset = unit(1, "mm"),
                                     annotation_name_gp = gpar(fontsize=0),
                                     annotation_name_rot = 90,
                                     gp = gpar(col = "white"),which = "row",
                                     annotation_legend_param = list(labels_gp = gpar(fontsize=legend.label.fontsize),
                                                                    title = 'Row\nannotation',
                                                                    title_gp = gpar(fontsize=legend.title.fontsize,fontface="bold")))
    }else if(is.null(row.annotation))
    {
      left.anno <- NULL
    }

    if(is.null(heatmap.cols))
    {
      heatmap.cols <- c('blue','white','red')
    }

    p <- Heatmap(plot.dat.z,cluster_rows = clust.row, cluster_columns = clust.column,
                 left_annotation = left.anno,top_annotation = top.anno,
                 show_row_names = show.row.names, show_column_names = show.column.names,name="Average\nexpression",
                 row_dend_side = "right", row_names_side = "right",rect_gp = gpar(col = "black"),
                 column_names_gp = gpar(fontsize=column.fontsize),column_names_rot = column.names.rotation,
                 row_names_gp = gpar(fontsize = row.fontsize),
                 heatmap_width = unit(heatmap.width, "cm"),heatmap_height = unit(heatmap.height,"cm"),
                 colorRamp2(breaks=seq(min.exp,max.exp,length.out=21),colors = colorpanel(21, low=heatmap.cols[1], mid=heatmap.cols[2], high=heatmap.cols[3])),
                 heatmap_legend_param = list(at= seq(min.exp,max.exp),legend_height = unit(5,"cm"),
                                             labels_gp = gpar(fontsize=legend.label.fontsize),
                                             title_gp = gpar(fontsize=legend.title.fontsize,fontface="bold")),...)

  }
  return(p)
}
