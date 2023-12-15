#' Stacked Violin Plot
#'
#' Main function to plot stacked violin plot
#'
#' @param object A Seurat object
#' @param features Features to plot
#' @param feature.annotation Annotation for features
#' @param feature.annotation.cols Annotation color
#' @param assay Name of the assay to use
#' @param split.by Factor to split the columns by
#' @param color.by One of the five ways to color for the violin plot: "features", "clusters", "median.exp" (median expression), "mean.exp" (mean expression), "split.var" (the variable that splits the violin plot).
#' When \code{split.by} is set (not to NULL), \code{color.by} must set to "split.var"; and vice versa.
#' @param clust.row,clust.column Whether to clusters row and/or columns
#' @param add.points Whether to plot data
#' @param points.size Size for data points
#' @param column.fontsize,row.fontsize,row.title.fontsize,legend.fontsize,legend.title.fontsize,features.fontsize Font size for figure and legend
#' @param column.names.rotation Rotation for column names
#' @param axis.text.hjust,axis.text.vjust Horizontal/vertical justification (in [0, 1]) for axis text (default: 0.5)
#' @param colors The colors to fill the violin plot when \code{color.by} is set to "median.exp" or "mean.exp"
#' @param random.annotation.col.seed Seed to generate random colors for feature annotations
#'
#' @return A \code{\link[ggplot2]{ggplot}} object
#'
#' @references Hao and Hao et al. (2021). “Integrated analysis of multimodal single-cell data.” Cell. doi:10.1016/j.cell.2021.04.048
#'
#' @importFrom dplyr %>% group_by
#' @importFrom ggplot2 ggplot scale_fill_continuous scale_fill_manual aes geom_violin facet_grid vars theme element_blank element_text element_rect scale_y_continuous scale_y_discrete labs geom_tile guides guide_legend geom_jitter position_jitterdodge
#' @importFrom ggtree ggtree layout_dendrogram
#' @importFrom scales hue_pal
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom cowplot plot_grid get_legend theme_nothing
#' @importFrom patchwork plot_spacer plot_layout
#'
#' @keywords internal
#'
StackedVlnPlot <- function(object, features, feature.annotation = NULL,feature.annotation.cols = NULL,
                           assay = 'RNA', split.by = NULL, color.by = 'features', clust.row = FALSE,clust.column = TRUE,
                           add.points = FALSE,points.size = 0.1,
                           column.fontsize = 12, row.fontsize = 8, row.title.fontsize = 15,
                           legend.fontsize = 12, legend.title.fontsize = 15, features.fontsize = 12,
                           column.names.rotation = 0,
                           axis.text.hjust = 0.5, axis.text.vjust = 0.5,
                           colors = NULL,random.annotation.col.seed = 1)
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
    if(!is.null(feature.annotation))
    {
      feature.annotation <- feature.annotation[match(features, feature.annotation$features),]
    }
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
    dat.1 <- data.frame(seurat.clusters=rep(seurat.clusters,each=nrow(dat)),
                        Gene = rep(rownames(dat),ncol(dat)),
                        exp=suppressMessages(reshape2::melt(dat)))

    dat.2 <- as.data.frame(dat.1 %>% group_by(seurat.clusters,Gene) %>%
                             dplyr::mutate(median.exp = median(exp.value), mean.exp = mean(exp.value)))

  }else
  {
    split.dat <- object@meta.data[,split.by]
    dat.1 <- data.frame(seurat.clusters=rep(seurat.clusters,each=nrow(dat)),
                        Gene = rep(rownames(dat),ncol(dat)), Split = rep(split.dat, each = nrow(dat)),
                        exp=suppressMessages(reshape2::melt(dat)))
    dat.2 <- as.data.frame(dat.1 %>% group_by(seurat.clusters,Gene) %>%
                             dplyr::mutate(clus.median.exp = median(exp.value), clus.mean.exp = mean(exp.value)))
    dat.2 <- as.data.frame(dat.2 %>% group_by(seurat.clusters,Gene,Split) %>%
                             dplyr::mutate(median.exp = median(exp.value), mean.exp = mean(exp.value)))
  }

  if(color.by == 'features')
  {
    fill.var <- 'Gene'
  }else if(color.by == 'clusters')
  {
    fill.var <- 'seurat.clusters'
  }else if(color.by == 'median.exp')
  {
    fill.var <- 'median.exp'
  }else if(color.by == 'mean.exp')
  {
    fill.var <- 'mean.exp'
  }else if(color.by == 'split.var' & !is.null(split.by))
  {
    fill.var <- 'Split'
  }else if(color.by == 'split.var' & is.null(split.by))
  {
    stop('Error: Cannot color by split.var if split.by is set to NULL!')
  }else
  {
    stop('Error: color.by must be one of features, clusters, median.exp, mean.exp or split.var')
  }

  if(fill.var == 'median.exp' | fill.var == 'mean.exp')
  {
    if(is.null(colors))
    {
      color.values <- scale_fill_continuous(low = "white", high = "blue")
    }else
    {
      color.values <- scale_fill_continuous(low = colors[1], high = colors[2])
    }
    legend.position = 'bottom'
    legend.title = gsub('.','\n',color.by, fixed = TRUE)
  }else if(fill.var == 'Split')
  {
    if(is.null(colors))
    {
      color.values <- scale_fill_manual(values = hue_pal()(length(unique(dat.2$Split))))
    }else
    {
      if(length(colors)!=length(unique(dat.2$Split)))
      {
        stop('Error: number of colors provided does not match with number of levels in split.by')
      }else
      {
        color.values <- scale_fill_manual(values = colors)
      }
    }
    legend.position = 'bottom'
    legend.title = 'Group'

  }else
  {
    color.values <- NULL
    legend.position <- 'none'
    legend.title = NULL
  }

  if(isTRUE(add.points))
  {
    p.add.points <- geom_jitter(position = position_jitterdodge(seed = 1),size = points.size)

  }else
  {
    p.add.points <- NULL
  }

  if(clust.row == FALSE & clust.column == FALSE)
  {
    dat.2$Gene <- factor(dat.2$Gene, levels = features)
    if(is.null(feature.annotation))
    {
      p <- ggplot(dat.2, aes(seurat.clusters, exp.value, fill = get(fill.var))) +
        geom_violin(scale = "width", adjust = 1, trim = TRUE) + p.add.points +
        facet_grid(rows =vars(Gene), scales = 'free_y' ) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = column.fontsize, angle = column.names.rotation, hjust = axis.text.hjust, vjust = axis.text.vjust),
              axis.title.y = element_text(size = row.title.fontsize), axis.text.y = element_text(size = row.fontsize),
              legend.title = element_text(size = legend.title.fontsize, hjust = 0.5), legend.text = element_text(size = legend.fontsize),
              legend.position = legend.position,legend.title.align  = 0.5,
              strip.background = element_blank(),
              strip.text.y = element_text(angle = 0,size = features.fontsize, hjust = 0),
              panel.background = element_rect(fill = NA, color = 'black')) +
        scale_y_continuous(expand = c(0,0)) +
        labs(y = 'Expression', fill = legend.title) + color.values

      final.plot <- p
    }else
    {
      p <- ggplot(dat.2, aes(seurat.clusters, exp.value, fill = get(fill.var))) +
        geom_violin(scale = "width", adjust = 1, trim = TRUE) + p.add.points +
        facet_grid(rows =vars(Gene), scales = 'free_y' ) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = column.fontsize, angle = column.names.rotation, hjust = axis.text.hjust, vjust = axis.text.vjust),
              axis.title.y = element_text(size = row.title.fontsize), axis.text.y = element_text(size = row.fontsize),
              legend.title = element_text(size = legend.title.fontsize), legend.text = element_text(size = legend.fontsize),
              legend.position = legend.position,
              strip.background = element_blank(),
              strip.text.y = element_text(angle = 0,size = features.fontsize, hjust = 0),
              panel.background = element_rect(fill = NA, color = 'black')) +
        scale_y_continuous(expand = c(0,0)) +
        labs(y = 'Expression', fill = legend.title) + color.values

      if(fill.var == 'median.exp' | fill.var == 'mean.exp' | fill.var == 'Split')
      {
        legend.p <- plot_grid(get_legend(p + theme(legend.position="bottom",legend.title.align  = 0.5,
                                                   legend.title = element_text(size = legend.title.fontsize, hjust = 0.5),
                                                   legend.text = element_text(size = legend.fontsize))))
        p <- p + theme(legend.position = 'none')
      }

      feature.annotation$features <- factor(feature.annotation$features,levels = rev(levels(dat.2$Gene)))
      if(is.null(feature.annotation.cols))
      {
        set.seed(random.annotation.col.seed)
        feature.annotation.cols <- distinctColorPalette(length(unique(feature.annotation$annotation)))
      }

      labels <- ggplot(feature.annotation,aes(x = 1, y = features))+
        geom_tile(aes(fill = annotation)) +
        scale_fill_manual(values = feature.annotation.cols) +
        scale_y_discrete(expand = c(0,0)) +
        theme_nothing()

      legend <- plot_grid(get_legend(labels + theme(legend.position="bottom",
                                                    legend.title = element_text(size = legend.title.fontsize),
                                                    legend.text = element_text(size = legend.fontsize))+
                                       guides(fill=guide_legend(nrow = 3))))

      if(fill.var == 'median.exp' | fill.var == 'mean.exp' | fill.var == 'Split')
      {
        final.plot <- p +  labels + legend.p +
          plot_spacer() +
          legend + plot_spacer()+
          plot_layout(ncol = 2, widths = c(0.95, 0.05,0.15), heights = c(0.8, 0.05,0.15))
      }else
      {
        final.plot <- p +  labels  +
          legend + plot_spacer()+
          plot_layout(ncol = 2, widths = c(0.95, 0.05), heights = c(0.8, 0.2))
      }
    }
  }else if(clust.row == TRUE & clust.column == FALSE)
  {
    if(is.null(split.by))
    {
      dat <- reshape2::dcast(dat.2,Gene ~ seurat.clusters,value.var = "mean.exp", fun.aggregate = mean)
    }else
    {
      dat <- reshape2::dcast(dat.2,Gene ~ seurat.clusters,value.var = "clus.mean.exp", fun.aggregate = mean)
    }
    rownames(dat) <- dat$Gene
    dat <- dat[,-1]

    h_clust <- hclust(dist(dat %>% as.matrix())) # hclust with distance matrix
    ddgram_row <- as.dendrogram(h_clust)
    ggtree_plot_row <- ggtree(ddgram_row)

    dat.2$Gene <- factor(dat.2$Gene, levels = rev(h_clust$labels[h_clust$order]))
    if(is.null(feature.annotation))
    {
      p <- ggplot(dat.2, aes(seurat.clusters, exp.value, fill = get(fill.var))) +
        geom_violin(scale = "width", adjust = 1, trim = TRUE) + p.add.points +
        facet_grid(rows =vars(Gene), scales = 'free_y' ) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = column.fontsize, angle = column.names.rotation, hjust = axis.text.hjust, vjust = axis.text.vjust),
              axis.title.y = element_text(size = row.title.fontsize), axis.text.y = element_text(size = row.fontsize),
              legend.title = element_text(size = legend.title.fontsize, hjust = 0.5), legend.text = element_text(size = legend.fontsize),
              legend.position = legend.position,legend.title.align  = 0.5,
              strip.background = element_blank(),
              strip.text.y = element_text(angle = 0,size = features.fontsize, hjust = 0),
              panel.background = element_rect(fill = NA, color = 'black')) +
        scale_y_continuous(expand = c(0,0)) +
        labs(y = 'Expression', fill = legend.title) + color.values

      final.plot <- ggtree_plot_row + p +
        plot_layout(ncol = 2, widths = c(0.1, 0.9))
    }else
    {
      p <- ggplot(dat.2, aes(seurat.clusters, exp.value, fill = get(fill.var))) +
        geom_violin(scale = "width", adjust = 1, trim = TRUE) + p.add.points +
        facet_grid(rows =vars(Gene), scales = 'free_y' ) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = column.fontsize, angle = column.names.rotation, hjust = axis.text.hjust, vjust = axis.text.vjust),
              axis.title.y = element_text(size = row.title.fontsize), axis.text.y = element_text(size = row.fontsize),
              legend.title = element_text(size = legend.title.fontsize), legend.text = element_text(size = legend.fontsize),
              legend.position = legend.position,
              strip.background = element_blank(),
              strip.text.y = element_text(angle = 0,size = features.fontsize, hjust = 0),
              panel.background = element_rect(fill = NA, color = 'black')) +
        scale_y_continuous(expand = c(0,0)) +
        labs(y = 'Expression', fill = legend.title) + color.values

      if(fill.var == 'median.exp' | fill.var == 'mean.exp' | fill.var == 'Split')
      {
        legend.p <- plot_grid(get_legend(p + theme(legend.position="bottom",legend.title.align  = 0.5,
                                                   legend.title = element_text(size = legend.title.fontsize, hjust = 0.5),
                                                   legend.text = element_text(size = legend.fontsize))))
        p <- p + theme(legend.position = 'none')
      }

      feature.annotation$features <- factor(feature.annotation$features,levels = rev(levels(dat.2$Gene)))
      if(is.null(feature.annotation.cols))
      {
        set.seed(random.annotation.col.seed)
        feature.annotation.cols <- distinctColorPalette(length(unique(feature.annotation$annotation)))
      }

      labels <- ggplot(feature.annotation,aes(x = 1, y = features))+
        geom_tile(aes(fill = annotation)) +
        scale_fill_manual(values = feature.annotation.cols) +
        scale_y_discrete(expand = c(0,0)) +
        theme_nothing()

      legend <- plot_grid(get_legend(labels + theme(legend.position="bottom",
                                                    legend.title = element_text(size = legend.title.fontsize),
                                                    legend.text = element_text(size = legend.fontsize)) +
                                       guides(fill=guide_legend(nrow = 3))))

      if(fill.var == 'median.exp' | fill.var == 'mean.exp' | fill.var == 'Split')
      {
        final.plot <- ggtree_plot_row + p +  labels +
          plot_spacer() + legend.p + plot_spacer() +
          plot_spacer() + legend + plot_spacer()+
          plot_layout(ncol = 3, widths = c(0.1,0.7, 0.05,0.15), heights = c(0.8, 0.05,0.15))
      }else
      {
        final.plot <- ggtree_plot_row + p +  labels  +
          plot_spacer() + legend + plot_spacer()+
          plot_layout(ncol = 3, widths = c(0.1,0.85, 0.05), heights = c(0.8, 0.2))
      }
    }

  }else if(clust.row == FALSE & clust.column == TRUE)
  {
    if(is.null(split.by))
    {
      dat <- reshape2::dcast(dat.2,Gene ~ seurat.clusters,value.var = "mean.exp", fun.aggregate = mean)
    }else
    {
      dat <- reshape2::dcast(dat.2,Gene ~ seurat.clusters,value.var = "clus.mean.exp", fun.aggregate = mean)
    }
    rownames(dat) <- dat$Gene
    dat <- dat[,-1]

    v_clust <- hclust(dist(dat %>% as.matrix() %>% t())) # hclust with distance matrix
    ddgram_col <- as.dendrogram(v_clust)
    ggtree_plot_col <- ggtree(ddgram_col) + layout_dendrogram()

    dat.2$Gene <- factor(dat.2$Gene, levels = features)
    dat.2$seurat.clusters <- factor(dat.2$seurat.clusters, levels = v_clust$labels[v_clust$order])

    if(is.null(feature.annotation))
    {
      p <- ggplot(dat.2, aes(seurat.clusters, exp.value, fill = get(fill.var))) +
        geom_violin(scale = "width", adjust = 1, trim = TRUE) + p.add.points +
        facet_grid(rows =vars(Gene), scales = 'free_y' ) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = column.fontsize, angle = column.names.rotation, hjust = axis.text.hjust, vjust = axis.text.vjust),
              axis.title.y = element_text(size = row.title.fontsize), axis.text.y = element_text(size = row.fontsize),
              legend.title = element_text(size = legend.title.fontsize, hjust = 0.5), legend.text = element_text(size = legend.fontsize),
              legend.position = legend.position,legend.title.align  = 0.5,
              strip.background = element_blank(),
              strip.text.y = element_text(angle = 0,size = features.fontsize, hjust = 0),
              panel.background = element_rect(fill = NA, color = 'black')) +
        scale_y_continuous(expand = c(0,0)) +
        labs(y = 'Expression', fill = legend.title) + color.values

      final.plot <- ggtree_plot_col + p +
        plot_layout(ncol = 1, heights  = c(0.1, 0.9))
    }else
    {
      p <- ggplot(dat.2, aes(seurat.clusters, exp.value, fill = get(fill.var))) +
        geom_violin(scale = "width", adjust = 1, trim = TRUE) + p.add.points +
        facet_grid(rows =vars(Gene), scales = 'free_y' ) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = column.fontsize, angle = column.names.rotation, hjust = axis.text.hjust, vjust = axis.text.vjust),
              axis.title.y = element_text(size = row.title.fontsize), axis.text.y = element_text(size = row.fontsize),
              legend.title = element_text(size = legend.title.fontsize), legend.text = element_text(size = legend.fontsize),
              legend.position = legend.position,
              strip.background = element_blank(),
              strip.text.y = element_text(angle = 0,size = features.fontsize, hjust = 0),
              panel.background = element_rect(fill = NA, color = 'black')) +
        scale_y_continuous(expand = c(0,0)) +
        labs(y = 'Expression', fill = legend.title) + color.values

      if(fill.var == 'median.exp' | fill.var == 'mean.exp' | fill.var == 'Split')
      {
        legend.p <- plot_grid(get_legend(p + theme(legend.position="bottom",legend.title.align  = 0.5,
                                                   legend.title = element_text(size = legend.title.fontsize, hjust = 0.5),
                                                   legend.text = element_text(size = legend.fontsize))))
        p <- p + theme(legend.position = 'none')
      }

      feature.annotation$features <- factor(feature.annotation$features,levels = rev(levels(dat.2$Gene)))
      if(is.null(feature.annotation.cols))
      {
        set.seed(random.annotation.col.seed)
        feature.annotation.cols <- distinctColorPalette(length(unique(feature.annotation$annotation)))
      }

      labels <- ggplot(feature.annotation,aes(x = 1, y = features))+
        geom_tile(aes(fill = annotation)) +
        scale_fill_manual(values = feature.annotation.cols) +
        scale_y_discrete(expand = c(0,0)) +
        theme_nothing()

      legend <- plot_grid(get_legend(labels + theme(legend.position="bottom",
                                                    legend.title = element_text(size = legend.title.fontsize),
                                                    legend.text = element_text(size = legend.fontsize)) +
                                       guides(fill=guide_legend(nrow = 3))))

      if(fill.var == 'median.exp' | fill.var == 'mean.exp' | fill.var == 'Split')
      {
        final.plot <- ggtree_plot_col + plot_spacer() +
          p +  labels +
          legend.p + plot_spacer() +
          legend + plot_spacer()+
          plot_layout(ncol = 2, widths = c(0.95,0.05),heights = c(0.1,0.7,0.05,0.15))
      }else
      {
        final.plot <- ggtree_plot_col + plot_spacer() +
          p +  labels  +
          legend + plot_spacer()+
          plot_layout(ncol = 2, widths = c(0.95,0.05), heights = c(0.1, 0.75,0.15))
      }
    }
  }else if(clust.row == TRUE & clust.column == TRUE)
  {
    if(is.null(split.by))
    {
      dat <- reshape2::dcast(dat.2,Gene ~ seurat.clusters,value.var = "mean.exp", fun.aggregate = mean)
    }else
    {
      dat <- reshape2::dcast(dat.2,Gene ~ seurat.clusters,value.var = "clus.mean.exp", fun.aggregate = mean)
    }
    rownames(dat) <- dat$Gene
    dat <- dat[,-1]

    h_clust <- hclust(dist(dat %>% as.matrix())) # hclust with distance matrix
    ddgram_row <- as.dendrogram(h_clust)
    ggtree_plot_row <- ggtree(ddgram_row)

    v_clust <- hclust(dist(dat %>% as.matrix() %>% t())) # hclust with distance matrix
    ddgram_col <- as.dendrogram(v_clust)
    ggtree_plot_col <- ggtree(ddgram_col) + layout_dendrogram()

    dat.2$Gene <- factor(dat.2$Gene, levels = rev(h_clust$labels[h_clust$order]))
    dat.2$seurat.clusters <- factor(dat.2$seurat.clusters, levels = v_clust$labels[v_clust$order])

    if(is.null(feature.annotation))
    {
      p <- ggplot(dat.2, aes(seurat.clusters, exp.value, fill = get(fill.var))) +
        geom_violin(scale = "width", adjust = 1, trim = TRUE) +  p.add.points +
        facet_grid(rows =vars(Gene), scales = 'free_y' ) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = column.fontsize, angle = column.names.rotation, hjust = axis.text.hjust, vjust = axis.text.vjust),
              axis.title.y = element_text(size = row.title.fontsize), axis.text.y = element_text(size = row.fontsize),
              legend.title = element_text(size = legend.title.fontsize, hjust = 0.5), legend.text = element_text(size = legend.fontsize),
              legend.position = legend.position,legend.title.align  = 0.5,
              strip.background = element_blank(),
              strip.text.y = element_text(angle = 0,size = features.fontsize, hjust = 0),
              panel.background = element_rect(fill = NA, color = 'black')) +
        scale_y_continuous(expand = c(0,0)) +
        labs(y = 'Expression', fill = legend.title) + color.values

      final.plot <- plot_spacer() + ggtree_plot_col +
        ggtree_plot_row + p +
        plot_layout(ncol = 2, widths  = c(0.1, 0.9),heights = c(0.1,0.9))
    }else
    {
      p <- ggplot(dat.2, aes(seurat.clusters, exp.value, fill = get(fill.var))) +
        geom_violin(scale = "width", adjust = 1, trim = TRUE) + p.add.points +
        facet_grid(rows =vars(Gene), scales = 'free_y' ) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = column.fontsize, angle = column.names.rotation, hjust = axis.text.hjust, vjust = axis.text.vjust),
              axis.title.y = element_text(size = row.title.fontsize), axis.text.y = element_text(size = row.fontsize),
              legend.title = element_text(size = legend.title.fontsize), legend.text = element_text(size = legend.fontsize),
              legend.position = legend.position,
              strip.background = element_blank(),
              strip.text.y = element_text(angle = 0,size = features.fontsize, hjust = 0),
              panel.background = element_rect(fill = NA, color = 'black')) +
        scale_y_continuous(expand = c(0,0)) +
        labs(y = 'Expression', fill = legend.title) + color.values

      if(fill.var == 'median.exp' | fill.var == 'mean.exp' | fill.var == 'Split')
      {
        legend.p <- plot_grid(get_legend(p + theme(legend.position="bottom",legend.title.align  = 0.5,
                                                   legend.title = element_text(size = legend.title.fontsize, hjust = 0.5),
                                                   legend.text = element_text(size = legend.fontsize))))
        p <- p + theme(legend.position = 'none')
      }

      feature.annotation$features <- factor(feature.annotation$features,levels = rev(levels(dat.2$Gene)))
      if(is.null(feature.annotation.cols))
      {
        set.seed(random.annotation.col.seed)
        feature.annotation.cols <- distinctColorPalette(length(unique(feature.annotation$annotation)))
      }

      labels <- ggplot(feature.annotation,aes(x = 1, y = features))+
        geom_tile(aes(fill = annotation)) +
        scale_fill_manual(values = feature.annotation.cols) +
        scale_y_discrete(expand = c(0,0)) +
        theme_nothing()

      legend <- plot_grid(get_legend(labels + theme(legend.position="bottom",
                                                    legend.title = element_text(size = legend.title.fontsize),
                                                    legend.text = element_text(size = legend.fontsize))+
                                       guides(fill=guide_legend(nrow = 3))))

      if(fill.var == 'median.exp' | fill.var == 'mean.exp' | fill.var == 'Split')
      {
        final.plot <-  plot_spacer() + ggtree_plot_col + plot_spacer() +
          ggtree_plot_row +  p +  labels +
          plot_spacer() + legend.p + plot_spacer() +
          plot_spacer()+ legend + plot_spacer()+
          plot_layout(ncol = 3, widths = c(0.1,0.85,0.05),heights = c(0.1,0.7,0.05,0.15))
      }else
      {
        final.plot <- plot_spacer() + ggtree_plot_col + plot_spacer() +
          ggtree_plot_row +  p +  labels +
          plot_spacer()+ legend + plot_spacer()+
          plot_layout(ncol = 3, widths = c(0.1,0.85,0.05),heights = c(0.1,0.75,0.15))
      }
    }
  }

  return(final.plot)
}
