#' Expression Plot
#'
#'An internal function to plot expression of a gene groups by identity class or additional meta data. Use its wrapper function \code{\link[Ragas]{RunExpressionPlot}} for plotting.
#'
#'
#' @param object A Seurat object
#' @param feature Gene name
#' @param assay Name of the assay to use (default: "RNA")
#' @param group.by Name of metadata column to group cells by
#' @param split.by Name of the metadata column to split expression by (optional)
#' @param ident Cell identity variable (default: seurat_clusters)
#' @param res A list containing data.frame  with differential testing results. See \code{\link[muscat]{pbDS}}
#' @param plot.ds.p.value TRUE/FALSE if p-value from differential testing should be plotted (default: FALSE)
#' @param plot.ds.adjp.value TRUE/FALSE if adjusted p-value from differential testing should be plotted (default: FALSE)
#' @param colors Colors for plotting. Number of colors must be the equal to number of idents or the number of groups if split.by is used (optional)
#' @param random.col.seed Set seed to control colors
#' @param point.size Point size for box plot
#' @param axis.text.size,legend.text.size Text size for axis or legend (default: 12)
#' @param axis.title.size,legend.title.size Title size for axis or legend (default: 15)
#' @param axis.text.angle Rotation of axis text
#' @param axis.text.hjust Horizontal justification (in [0, 1]) for axis text
#' @param axis.text.vjust Vertical justification (in [0, 1]) for axis text
#' @param sig.size Asterisk size for significant tests
#' @param sig.color Asterisk color for significant tests
#' @param sig.bracket.color Bracket color for significant tests
#' @param sig.bracket.width Bracket width for significant tests
#'
#' @importFrom dplyr %>% summarise group_by
#' @importFrom reshape2 melt dcast
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom ggprism add_pvalue
#' @importFrom rstatix add_significance add_x_position add_y_position group_by
#' @importFrom ggplot2 ggplot aes geom_boxplot position_dodge geom_point position_jitterdodge scale_color_manual labs theme element_text
#'
#' @keywords internal



ExpressionPlot <- function(object, feature =NULL , assay = 'RNA',
                           group.by = NULL, split.by = NULL, ident = 'seurat_clusters',
                           res = NULL, plot.ds.p.value = FALSE, plot.ds.adjp.value = FALSE,
                           colors = NULL, random.col.seed = 1,
                           point.size = 1,
                           axis.text.size = 12,
                           axis.text.angle = 90, axis.text.hjust = 0, axis.text.vjust = 0,
                           axis.title.size = 15,
                           legend.text.size=12,legend.title.size=15,
                           sig.size = 8, sig.bracket.color = 'gray50', sig.color = 'gray50',
                           sig.bracket.width = 0.5)
{

  idx <- match(feature, rownames(object))
  if(sum(is.na(idx))>0)
  {
    print(paste(feature,"not found!",sep = ' '))
    stop()
  }else
  {
    dat <- GetAssayData(object, slot="data", assay = assay)[idx,]
  }

  if(is.null(group.by))
  {
    print("group.by must be assigned")
    stop()
  }

  clusters <- object@meta.data[,ident]
  groups <- object@meta.data[,group.by]

  if(is.null(split.by))
  {
    dat.1 <- data.frame(Clusters=clusters,Gene = feature, Group = groups,Exp= dat)
    dat.2 <- as.data.frame(dat.1 %>% group_by(Clusters,Group) %>%
                             dplyr::summarise(avg.exp = mean(Exp),.groups = 'drop'))

  }else
  {
    if(!is.factor(object@meta.data[,split.by]))
    {
      object@meta.data[,split.by] <- factor(object@meta.data[,split.by])
    }
    split.groups <- object@meta.data[,split.by]
    dat.1 <- data.frame(Clusters=clusters,Gene = feature, Group = groups,Split = split.groups,Exp= dat)
    dat.2 <- as.data.frame(dat.1 %>% group_by(Clusters, Group, Split) %>%
                             dplyr::summarise(avg.exp = mean(Exp),.groups = 'drop'))
  }

  if(is.null(split.by))
  {
    if(is.null(colors))
    {
      set.seed(random.col.seed)
      colors <- distinctColorPalette(length(levels(dat.2$Clusters)))
    }

    p <- ggplot(dat.2,aes(x = Clusters,y = avg.exp, color = Clusters)) +
      geom_boxplot(position = position_dodge(width = 0.8),outlier.shape = NA,aes(color = Clusters)) +
      geom_point(position = position_jitterdodge(seed = 2, dodge.width = 0.8),size=point.size,
                 aes(color = Clusters)) +
      scale_color_manual(values = colors) +
      labs(y = paste(feature,' - mean exp',sep='')) +
      theme(axis.text.y = element_text(size = axis.text.size),
            axis.text.x = element_text(size = axis.text.size, angle = axis.text.angle, hjust = axis.text.hjust, vjust = axis.text.vjust),
            axis.title = element_text(size = axis.title.size),
            legend.title = element_text(size = legend.title.size),
            legend.text = element_text(size = legend.text.size))

    final.plot <- p

  }else
  {
    if(is.null(colors))
    {
      set.seed(random.col.seed)
      colors <- distinctColorPalette(length(levels(dat.2$Split)))
    }

    p <- ggplot(dat.2,aes(x = Clusters,y = avg.exp, color = Split)) +
      geom_boxplot(position = position_dodge(width = 0.8),outlier.shape = NA,aes(color = Split)) +
      geom_point(position = position_jitterdodge(seed = 2, dodge.width = 0.8),size=point.size,
                 aes(color = Split)) +
      scale_color_manual(values = colors) +
      labs(y = paste(feature,' - mean exp',sep=''), color = 'Group') +
      theme(axis.text = element_text(size = axis.text.size),
            axis.text.x = element_text(size = axis.text.size, angle = axis.text.angle, hjust = axis.text.hjust, vjust = axis.text.vjust),
            axis.title = element_text(size = axis.title.size),
            legend.title = element_text(size = legend.title.size),
            legend.text = element_text(size = legend.text.size))

    if(plot.ds.p.value == TRUE & plot.ds.adjp.value == FALSE)
    {
      dat.p.value <- lapply(res, function(x){data.frame(p = do.call('rbind',lapply(x, function(y){y$p_val[match(feature,y$gene)]})))})
      dat.p.value <- lapply(dat.p.value, function(x){data.frame(Clusters = rownames(x), x)})

      myfun <- function(x,y){cbind(x,group1 =unlist(strsplit(y,split = '_'))[1],group2 =unlist(strsplit(y,split = '_'))[2])}
      dat.p.value <- do.call('rbind',mapply(myfun,dat.p.value,names(dat.p.value),SIMPLIFY = FALSE))
      dat.p.value$.y. <- 'avg.exp'
      rownames(dat.p.value) <- NULL

      fun.grpcheck <- function(x)
      {
        true.order <- levels(dat.2$Split)

        idx <- match(true.order,c(x['group1'],x['group2']));idx<- which(!is.na(idx))
        idx <- sort(idx)
        return(data.frame(group1 = true.order[idx[1]], group2 = true.order[idx[2]]))

      }
      dat.p.value[,c('group1','group2')] <- do.call('rbind',apply(dat.p.value,1, fun.grpcheck))

      dat.p.value <- dat.p.value %>% rstatix::add_significance(p.col = 'p', output.col = 'p.signif')

      dat.args <- list()
      dat.args$data <- dat.2
      attr(dat.p.value, "args") <- dat.args
      class(dat.p.value) <- c('rstatix_test',class(dat.p.value))

      dat.p.value <- add_x_position(dat.p.value, x = 'Clusters', group = 'Split', dodge = 0.8)
      dat.p.value <- dat.p.value[dat.p.value$p.signif!='ns',]

      dat.p.value <- dat.p.value[which(!is.na(dat.p.value$p)),]

      if(nrow(dat.p.value)!=0)
      {
        dat.p.value <- dat.2  %>%
                       group_by(Clusters) %>%
                       add_y_position(test = dat.p.value,formula = avg.exp ~ Split ,
                                      fun = 'max', step.increase = 0.05, scales = 'free')
        dat.p.value$Clusters <- factor(dat.p.value$Clusters, levels = levels(dat.2$Clusters))

        final.plot <- p + add_pvalue(dat.p.value,
                                     y.position = "y.position",xmin = 'xmin', xmax = 'xmax',
                                     label = 'p.signif', label.size = sig.size,
                                     bracket.color = sig.bracket.color,color = sig.color,
                                     bracket.size = sig.bracket.width, tip.length = 0.01)
      }else
      {
        final.plot <- p
      }
    }else if(plot.ds.p.value == FALSE & plot.ds.adjp.value == TRUE)
    {
      dat.p.value <- lapply(res, function(x){data.frame(p = do.call('rbind',lapply(x, function(y){y$p_adj.loc[match(feature,y$gene)]})))})
      dat.p.value <- lapply(dat.p.value, function(x){data.frame(Clusters = rownames(x), x)})

      myfun <- function(x,y){cbind(x,group1 =unlist(strsplit(y,split = '_'))[1],group2 =unlist(strsplit(y,split = '_'))[2])}
      dat.p.value <- do.call('rbind',mapply(myfun,dat.p.value,names(dat.p.value),SIMPLIFY = FALSE))
      dat.p.value$.y. <- 'avg.exp'
      rownames(dat.p.value) <- NULL

      fun.grpcheck <- function(x)
      {
        true.order <- levels(dat.2$Split)

        idx <- match(true.order,c(x['group1'],x['group2']));idx<- which(!is.na(idx))
        idx <- sort(idx)
        return(data.frame(group1 = true.order[idx[1]], group2 = true.order[idx[2]]))

      }
      dat.p.value[,c('group1','group2')] <- do.call('rbind',apply(dat.p.value,1, fun.grpcheck))


      dat.p.value <- dat.p.value %>% rstatix::add_significance(p.col = 'p', output.col = 'p.signif')

      dat.args <- list()
      dat.args$data <- dat.2
      attr(dat.p.value, "args") <- dat.args
      class(dat.p.value) <- c('rstatix_test',class(dat.p.value))

      dat.p.value <- add_x_position(dat.p.value, x = 'Clusters', group = 'Split', dodge = 0.8)
      dat.p.value <- dat.p.value[dat.p.value$p.signif!='ns',]

      dat.p.value <- dat.p.value[which(!is.na(dat.p.value$p)),]

      if(nrow(dat.p.value)!=0)
      {
        dat.p.value <- dat.2  %>%
          group_by(Clusters) %>%
          add_y_position(test = dat.p.value,formula = avg.exp ~ Split ,
                         fun = 'max', step.increase = 0.05, scales = 'free')
        dat.p.value$Clusters <- factor(dat.p.value$Clusters, levels = levels(dat.2$Clusters))

        final.plot <- p + add_pvalue(dat.p.value,
                                     y.position = "y.position",xmin = 'xmin', xmax = 'xmax',
                                     label = 'p.signif', label.size = sig.size,
                                     bracket.color = sig.bracket.color,color = sig.color,
                                     bracket.size = sig.bracket.width, tip.length = 0.01)
      }else
      {
        final.plot <- p
      }

    }else if(plot.ds.p.value == TRUE & plot.ds.adjp.value == TRUE)
    {
      print('both plot.ds.p.value & plot.ds.adjp.value cannot be set to TRUE')
      stop()
    }else
    {
      final.plot <- p
    }
  }

  return(final.plot)
}
