#' AnnotatedDotPlot for Seurat or post-integration (pi) object
#'
#' The AnnotatedDotPlot function enhances the original \code{\link[Seurat]{DotPlot}} function from Seurat by:
#' (1) adding feature annotations to the columns); (2) allowing clustering and visualzation of dendrograms of
#' identities on the rows. Use its wrapper function \code{\link[Ragas]{RunAnnotatedDotPlot}} for plotting.
#'
#' @param object A Seurat object
#' @param features Genes/features to plot
#' @param cols Colors to plot. Same as the "cols" argument of \code{\link[Seurat]{DotPlot}} (default: c("lightgrey", "blue"))
#' @param split.by Factor to split the groups by. Same as the "split.by" argument of \code{\link[Seurat]{DotPlot}}
#' @param group.by Factor to group the cells by. Same as the "group.by" argument of \code{\link[Seurat]{DotPlot}}
#' @param clust.row Whether to cluster the rows/identities (default: TRUE)
#' @param clust.column Whether to cluster the columns/features (default: FALSE)
#' @param column.annotation A data.frame contains features and their annotation groups
#' @param column.annotation.cols Colors for column annotation
#' @param column.fontsize Size of column text (default: 12)
#' @param row.fontsize Size of row text (default: 12)
#' @param column.fontface,row.fontface Fontface for column/row labels either "plain", "bold", "italic" or "bold.italic" (default: "plain")
#' @param legend.label.fontsize Size of the legend labels (default: 13)
#' @param legend.title.fontsize Size of the legend title (default: 15)
#' @param ... Extra parameters passed to DotPlot
#'
#' @return A \code{\link[ggplot2]{ggplot}} object
#'
#' @references Hao and Hao et al. (2021). “Integrated analysis of multimodal single-cell data.” Cell. doi:10.1016/j.cell.2021.04.048
#'
#' @importFrom ggplot2 theme element_blank element_text scale_y_discrete ggplot aes geom_tile scale_fill_manual
#' @importFrom patchwork plot_spacer plot_layout
#' @importFrom cowplot theme_nothing plot_grid get_legend
#' @importFrom ggtree ggtree
#' @importFrom dplyr %>%
#' @importFrom aplot xlim2
#' @importFrom randomcoloR distinctColorPalette
#
#' @keywords internal


AnnotatedDotPlot <- function(object,
                             features,
                             cols = c("lightgrey", "blue"),
                             split.by = NULL,
                             group.by = NULL,
                             clust.row = TRUE,
                             clust.column = FALSE,
                             column.annotation = NULL,
                             column.annotation.cols = NULL,
                             column.fontsize = 12,
                             row.fontsize = 12,
                             column.fontface = "plain",
                             row.fontface = "plain",
                             legend.label.fontsize=13,
                             legend.title.fontsize=15,...)
{

  if(!is.null(split.by) & clust.row == TRUE)
  {
    print('Cannot use split.by and cluster rows together!')
    print('Running without clustering rows...')
    clust.row <- FALSE
  }

  if(!is.null(group.by))
  {
    Idents(object) <- object@meta.data[,group.by]
  }

  p <- DotPlot(object = object, features = features, cols=cols, split.by = split.by, group.by = group.by, ...) + RotatedAxis() +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(size = column.fontsize, face = column.fontface, angle = 90, vjust = 0.5),
          axis.text.y = element_text(size = row.fontsize, face = row.fontface),
          legend.title = element_text(size = legend.title.fontsize),
          legend.text = element_text(size = legend.label.fontsize)) +
    scale_y_discrete(position = "right")

  if(is.null(split.by))
  {
    dat <- matrix(0,nrow = length(levels(Idents(object))), ncol = length(features))
    rownames(dat)<-  levels(Idents(object))
    colnames(dat) <- features

    for(i in 1:nrow(dat))
    {
      idx <- which(p$data$id==rownames(dat)[i])
      tmp <- p$data$pct.exp[idx] ; names(tmp) = p$data$features.plot[idx]
      dat[i,names(tmp)] <- p$data$pct.exp[idx]
    }
  }else
  {
    if(!is.factor(object@meta.data[,split.by]))
    {
      object@meta.data[,split.by] <- factor(object@meta.data[,split.by])
    }

    dat <- matrix(0,nrow = length(levels(Idents(object)))* length(levels(object@meta.data[,split.by])),
                  ncol = length(features))
    rownames(dat)<- paste(rep(levels(Idents(object)),each = length(levels(object@meta.data[,split.by]))),
                          levels(object@meta.data[,split.by]),sep='_')
    colnames(dat) <- features

    for(i in 1:nrow(dat))
    {
      idx <- which(p$data$id==rownames(dat)[i])
      tmp <- p$data$pct.exp[idx] ; names(tmp) = p$data$features.plot[idx]
      dat[i,names(tmp)] = p$data$pct.exp[idx]
    }
  }



  if(clust.column==TRUE & clust.row == TRUE)
  {
    v_clust <- hclust(dist(dat %>% as.matrix() %>% t())) # hclust with distance matrix
    ddgram_col <- as.dendrogram(v_clust)
    ggtree_plot_col <- ggtree(ddgram_col) + layout_dendrogram()

    h_clust <- hclust(dist(dat %>% as.matrix())) # hclust with distance matrix
    ddgram_row <- as.dendrogram(h_clust)
    ggtree_plot_row <- ggtree(ddgram_row)

    if(is.null(column.annotation))
    {
      tmp <- object
      Idents(tmp) <- factor(Idents(tmp),levels = h_clust$labels[h_clust$order])

      if(!is.null(group.by))
      {
        group.by <- NULL
      }

      p.dot <- DotPlot(tmp, features = v_clust$labels[v_clust$order],cols=cols, split.by = split.by, group.by = group.by,...) + RotatedAxis() +
        theme(axis.title = element_blank(),
              axis.text.x = element_text(size = column.fontsize, face = column.fontface, angle = 90, vjust = 0.5),
              axis.text.y = element_text(size = row.fontsize, face = row.fontface),
              legend.title = element_text(size = legend.title.fontsize),
              legend.text = element_text(size = legend.label.fontsize)) +
        scale_y_discrete(position = "right")

      final.plot <- plot_spacer() + ggtree_plot_col +
                    ggtree_plot_row + p.dot+
                    plot_layout(ncol = 2, widths = c(0.7, 4), heights = c(0.9, 4, 1))

    }else
    {
      tmp <- object
      Idents(tmp) <- factor(Idents(tmp),levels = h_clust$labels[h_clust$order])

      if(!is.null(group.by))
      {
        group.by <- NULL
      }

      p.dot <- DotPlot(tmp, features = v_clust$labels[v_clust$order],cols=cols, split.by = split.by, group.by = group.by,...) + RotatedAxis() +
        theme(axis.title = element_blank(),
              axis.text.x = element_text(size = column.fontsize, face = column.fontface, angle = 90, vjust = 0.5),
              axis.text.y = element_text(size = row.fontsize, face = row.fontface),
              legend.title = element_text(size = legend.title.fontsize),
              legend.text = element_text(size = legend.label.fontsize)) +
        scale_y_discrete(position = "right")

      if(is.null(column.annotation.cols))
      {
        set.seed(1)
        column.annotation.cols <- distinctColorPalette(length(unique(column.annotation$annotation)))
      }
      labels <- ggplot(column.annotation,aes(x = features, y = 1, fill = annotation))+
        geom_tile() + scale_fill_manual(values = column.annotation.cols) + theme_nothing() +
        xlim2(p.dot)

      legend <- plot_grid(get_legend(labels + theme(legend.position="bottom",
                                                    legend.title = element_text(size = legend.title.fontsize),
                                                    legend.text = element_text(size = legend.label.fontsize))))

      final.plot <- plot_spacer() + ggtree_plot_col +
                    plot_spacer() + labels +
                    ggtree_plot_row+ p.dot+
                    plot_spacer() +  legend +
                    plot_layout(ncol = 2, widths = c(0.7, 4), heights = c(0.9, 0.2, 4, 1))
    }
  }else if(clust.column == FALSE & clust.row == TRUE)
  {
    h_clust <- hclust(dist(dat %>% as.matrix())) # hclust with distance matrix
    ddgram_row <- as.dendrogram(h_clust)
    ggtree_plot_row <- ggtree(ddgram_row)

    if(is.null(column.annotation))
    {
      tmp <- object
      Idents(tmp) <- factor(Idents(tmp),levels = h_clust$labels[h_clust$order])

      if(!is.null(group.by))
      {
        group.by <- NULL
      }

      p.dot <- DotPlot(tmp, features = features,cols=cols, split.by = split.by,group.by = group.by,...) + RotatedAxis() +
        theme(axis.title = element_blank(),
              axis.text.x = element_text(size = column.fontsize, face = column.fontface, angle = 90, vjust = 0.5),
              axis.text.y = element_text(size = row.fontsize, face = row.fontface),
              legend.title = element_text(size = legend.title.fontsize),
              legend.text = element_text(size = legend.label.fontsize)) +
        scale_y_discrete(position = "right")

      final.plot <- ggtree_plot_row + p.dot+
                    plot_layout(ncol = 2, widths = c(0.7, 4))
    }else
    {
      tmp <- object
      Idents(tmp) <- factor(Idents(tmp),levels = h_clust$labels[h_clust$order])

      if(!is.null(group.by))
      {
        group.by <- NULL
      }

      p.dot <- DotPlot(tmp, features = features,cols=cols, split.by = split.by,group.by = group.by,...) + RotatedAxis() +
        theme(axis.title = element_blank(),
              axis.text.x = element_text(size = column.fontsize, face = column.fontface, angle = 90, vjust = 0.5),
              axis.text.y = element_text(size = row.fontsize, face = row.fontface),
              legend.title = element_text(size = legend.title.fontsize),
              legend.text = element_text(size = legend.label.fontsize)) +
        scale_y_discrete(position = "right")

      if(is.null(column.annotation.cols))
      {
        set.seed(1)
        column.annotation.cols <- distinctColorPalette(length(unique(column.annotation$annotation)))
      }
      labels <- ggplot(column.annotation,aes(x = features, y = 1, fill = annotation))+
        geom_tile() + scale_fill_manual(values = column.annotation.cols) + theme_nothing() +
        xlim2(p.dot)

      legend <- plot_grid(get_legend(labels + theme(legend.position="bottom",
                                                    legend.title = element_text(size = legend.title.fontsize),
                                                    legend.text = element_text(size = legend.label.fontsize))))

      final.plot <- plot_spacer() +  labels +
                    ggtree_plot_row +p.dot+
                    plot_spacer()  + legend +
                    plot_layout(ncol = 2, widths = c(0.7, 4), heights = c(0.2, 4, 1))
    }
  }else if(clust.column==TRUE & clust.row == FALSE)
  {
    v_clust <- hclust(dist(dat %>% as.matrix() %>% t())) # hclust with distance matrix
    ddgram_col <- as.dendrogram(v_clust)
    ggtree_plot_col <- ggtree(ddgram_col) + layout_dendrogram()

    if(is.null(column.annotation))
    {
      tmp <- object

      p.dot <- DotPlot(tmp, features = v_clust$labels[v_clust$order],cols=cols, split.by = split.by,group.by = group.by,...) + RotatedAxis() +
        theme(axis.title = element_blank(),
              axis.text.x = element_text(size = column.fontsize, face = column.fontface, angle = 90, vjust = 0.5),
              axis.text.y = element_text(size = row.fontsize, face = row.fontface),
              legend.title = element_text(size = legend.title.fontsize),
              legend.text = element_text(size = legend.label.fontsize)) +
        scale_y_discrete(position = "right")

      final.plot <- ggtree_plot_col +
                    p.dot+
                    plot_layout(ncol = 1, widths = c(0.7,4), heights = c(0.9, 4))
    }else
    {
      tmp <- object

      p.dot <- DotPlot(tmp, features = v_clust$labels[v_clust$order],cols=cols, split.by = split.by,group.by = group.by,...) + RotatedAxis() +
        theme(axis.title = element_blank(),
              axis.text.x = element_text(size = column.fontsize, face = column.fontface, angle = 90, vjust = 0.5),
              axis.text.y = element_text(size = row.fontsize, face = row.fontface),
              legend.title = element_text(size = legend.title.fontsize),
              legend.text = element_text(size = legend.label.fontsize)) +
        scale_y_discrete(position = "right")

      if(is.null(column.annotation.cols))
      {
        set.seed(1)
        column.annotation.cols <- distinctColorPalette(length(unique(column.annotation$annotation)))
      }
      labels <- ggplot(column.annotation,aes(x = features, y = 1, fill = annotation))+
        geom_tile() + scale_fill_manual(values = column.annotation.cols) + theme_nothing() +
        xlim2(p.dot)

      legend <- plot_grid(get_legend(labels + theme(legend.position="bottom",
                                                    legend.title = element_text(size = legend.title.fontsize),
                                                    legend.text = element_text(size = legend.label.fontsize))))

      final.plot <-  ggtree_plot_col +
                    labels +
                     p.dot +
                     legend +
                     plot_layout(ncol = 1, widths = c(0.7, 4), heights = c(0.9, 0.2, 4, 1))
    }
  }else if(clust.column == FALSE & clust.row == FALSE)
  {

    if(is.null(column.annotation))
    {
      tmp <- object
      p.dot <- DotPlot(tmp, features = features,cols=cols, split.by = split.by,group.by = group.by,...) + RotatedAxis() +
        theme(axis.title = element_blank(),
              axis.text.x = element_text(size = column.fontsize, face = column.fontface, angle = 90, vjust = 0.5),
              axis.text.y = element_text(size = row.fontsize, face = row.fontface),
              legend.title = element_text(size = legend.title.fontsize),
              legend.text = element_text(size = legend.label.fontsize)) +
        scale_y_discrete(position = "right")

      final.plot <- p.dot
    }else
    {
      tmp <- object
      p.dot <- DotPlot(tmp, features = features,cols=cols, split.by = split.by,group.by = group.by,...) + RotatedAxis() +
        theme(axis.title = element_blank(),
              axis.text.x = element_text(size = column.fontsize, face = column.fontface, angle = 90, vjust = 0.5),
              axis.text.y = element_text(size = row.fontsize, face = row.fontface),
              legend.title = element_text(size = legend.title.fontsize),
              legend.text = element_text(size = legend.label.fontsize)) +
        scale_y_discrete(position = "right")

      if(is.null(column.annotation.cols))
      {
        set.seed(1)
        column.annotation.cols <- distinctColorPalette(length(unique(column.annotation$annotation)))
      }
      labels <- ggplot(column.annotation,aes(x = features, y = 1, fill = annotation))+
        geom_tile() + scale_fill_manual(values = column.annotation.cols) + theme_nothing() +
        xlim2(p.dot)

      legend <- plot_grid(get_legend(labels + theme(legend.position="bottom",
                                                    legend.title = element_text(size = legend.title.fontsize),
                                                    legend.text = element_text(size = legend.label.fontsize))))

      final.plot <-  labels +
                     p.dot+
                     legend +
        plot_layout(ncol = 1,  heights = c(0.2,  4, 1))
    }
  }

  return(final.plot)
}
