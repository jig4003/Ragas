#' Run Cell Proportion Analysis and Make Plots
#'
#' An internal function to run pooled or unpooled proportion analysis and make bar plot or box plot.
#' Use its wrapper function \code{\link[Ragas]{RunProportionPlot}} for plotting.
#'
#'
#' @param object A Seurat object
#' @param var.group Name of the metadata column to group cells by
#' @param method Method for cell proportion analysis: either "unpooled" or "pooled"
#' @param plot.all Whether to add proportion for all cells. Only applicable to pooled method (default: TRUE)
#' @param colors Colors for the plotting group
#' @param prop.by Choose one of the two methods to calculate pooled cell proportion: by "cluster" or "group"
#' @param split.by Name of the metadata column to split the pooled proportion plot by (optional)
#' @param parent.meta Meta data inherited from the parent object
#' @param unpool.plot.type Method for the unpooled plot: either boxplot or barplot
#' @param var.unpool Name of the secondary grouping variable from the metadata column, such as sample name
#' @param n.col Number of columns per plot (default: 5)
#' @param plot.p.sig Whether to highlight tests with significant raw p-values
#' @param plot.p.adj.sig Whether to highlight tests with significant adjusted p-values
#' @param title.text.size Text size for plot title
#' @param axis.title.size Axis title size
#' @param axis.text.size Axis text size
#' @param axis.text.angle Rotation of axis text
#' @param axis.text.hjust Horizontal justification (in [0, 1]) for axis text
#' @param axis.text.vjust Vertical justification (in [0, 1]) for axis text
#' @param legend.title.size Legend title size
#' @param legend.text.size Legend text size
#' @param point.size Point size for box plot
#' @param sig.size Asterisk size for significant tests
#' @param sig.color Asterisk color for significant tests
#' @param sig.bracket.color Bracket color for significant tests
#' @param sig.bracket.width Bracket width for significant tests
#' @param random.col.seed Set seed to control colors
#'
#' @importFrom scales hue_pal
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual theme element_blank element_text scale_y_continuous geom_hline facet_wrap
#' geom_boxplot position_dodge geom_point position_jitterdodge scale_color_manual expansion geom_errorbar labs
#' @importFrom reshape2 melt
#'
#' @importFrom dplyr %>% summarize group_by n
#'
#' @importFrom ggprism add_pvalue
#' @importFrom rstatix drop_na add_y_position
#' @importFrom randomcoloR distinctColorPalette
#'
#' @keywords internal
#'

ProportionPlot <- function(object, var.group, method = 'unpooled',plot.all = TRUE,
                           colors = NULL,
                           prop.by='cluster',  split.by = NULL,
                           parent.meta = parent.meta,
                           unpool.plot.type = 'boxplot',
                           var.unpool, n.col = 5,
                           plot.p.sig = TRUE, plot.p.adj.sig = FALSE,
                           title.text.size = 18, axis.title.size = 15, axis.text.size = 12,
                           axis.text.angle = 90, axis.text.hjust = 0, axis.text.vjust = 0,
                           legend.title.size = 15, legend.text.size = 12,
                           point.size = 1,
                           sig.size = 8, sig.bracket.color = 'gray50', sig.color = 'gray50',
                           sig.bracket.width = 0.5,
                           random.col.seed = 1)
{

  if(method == 'unpooled')
  {
    # Define modified t_test function
    t_test_mod <- function(t.dat,...)
    {
      out <- try(data.frame(t.dat %>%
                              rstatix::t_test(Proportion ~ Group,
                                              p.adjust.method = 'none', detailed = TRUE)),silent = TRUE)
      if(class(out)=='try-error')
      {
        if(!is.factor(t.dat$Group))
        {
          all.groups <- levels(factor(t.dat$Group))
        }else
        {
          all.groups <- levels(t.dat$Group)
        }

        final.out <- NULL
        for(i in 1:(length(all.groups)-1))
        {
          for(j in (i+1):length(all.groups))
          {
            t.dat.sub <- t.dat[which(t.dat$Group %in% c(all.groups[i],all.groups[j])),]
            t.dat.sub$Group <- droplevels(t.dat.sub$Group)

            out <- try(data.frame(t.dat.sub %>%
                                    rstatix::t_test(Proportion ~ Group,
                                                    p.adjust.method = 'none', detailed = TRUE)),silent = TRUE)
            if(class(out) == 'try-error')
            {
              e.out <- data.frame(t.dat.sub %>% group_by(Group) %>%  dplyr::summarise(mean = mean(Proportion),
                                                                                      n = length(Proportion),.groups = 'drop'))

              out <- data.frame(estimate = (e.out$mean[1] - e.out$mean[2]),
                                estimate1 = e.out$mean[1], estimate2 = e.out$mean[2],
                                .y.='Proportion',
                                group1 = e.out$Group[1],
                                group2 = e.out$Group[2],
                                n1 = e.out$n[1], n2 = e.out$n[2], statistic = NA, p = NA,
                                df = NA, conf.low = NA,conf.high = NA, method = NA, alternative = NA)
            }else
            {
              out <- data.frame(t.dat.sub %>%
                                  rstatix::t_test(Proportion ~ Group,
                                                  p.adjust.method = 'none', detailed = TRUE))
            }
            final.out <- rbind(final.out,out)
          }
        }

      }else
      {
        final.out <- data.frame(t.dat %>%
                                  rstatix::t_test(Proportion ~ Group,p.adjust.method = 'none', detailed = TRUE))
      }

      return(final.out)
    }

    # droplevels for var.unpool
    object@meta.data[,var.unpool] <- droplevels(object@meta.data[,var.unpool])

    dat_prop <- as.matrix(table(Idents(object),object@meta.data[,var.unpool]))

    if(is.null(parent.meta))
    {
      dat_prop <- prop.table(dat_prop,margin = 2)
    }else
    {
      if(!is.factor(object@meta.data[,var.unpool]))
      {
        object@meta.data[,var.unpool] <- factor(object@meta.data[,var.unpool])
      }
      if(!is.factor(parent.meta[,var.unpool]))
      {
        parent.meta[,var.unpool] <- factor(parent.meta[,var.unpool])
      }

      if(!identical(levels(parent.meta[,var.unpool]),levels(object@meta.data[,var.unpool])))
      {
        parent.meta[,var.unpool] <- factor(parent.meta[,var.unpool],levels = levels(object@meta.data[,var.unpool]))
      }

      dat_prop <- t(apply(dat_prop,1,function(x){x/as.numeric(table(parent.meta[,var.unpool]))}))
    }
    dat_prop_1 <- as.data.frame.matrix(t(dat_prop))

    if(!is.factor(object@meta.data[,var.group]))
    {
      object@meta.data[,var.group] <- factor(object@meta.data[,var.group])
    }

    idx <- match(rownames(dat_prop_1),object@meta.data[,var.unpool])
    dat_prop_1 <- cbind(dat_prop_1, Group =  object@meta.data[idx,var.group])
    dat_prop_1$Group <- factor(dat_prop_1$Group,levels = levels(object@meta.data[,var.group]))


    dat <- data.frame(Pool_Group = rownames(dat_prop_1),dat_prop_1,check.names = F)
    dat <- suppressMessages(melt(dat))
    colnames(dat) <- c('Pool_Group','Group','Cluster','Proportion')


    if(is.null(colors))
    {
      #colors <- hue_pal()(length(levels(factor(object@meta.data[,var.group]))))
      set.seed(random.col.seed)
      colors <- distinctColorPalette(length(levels(factor(object@meta.data[,var.group]))))
    }

    if(unpool.plot.type == 'boxplot')
    {
      # t-test
      if(!is.factor(dat$Cluster))
      {
        all.clusters <- levels(factor(dat$Cluster))
      }else
      {
        all.clusters <- levels(dat$Cluster)
      }

      test.out <- NULL
      for(i in 1:length(all.clusters))
      {
        clus.dat <- dat[which(dat$Cluster == all.clusters[i] ),]
        clus.dat$Group <- droplevels(clus.dat$Group)
        clus.dat$Cluster <- droplevels(clus.dat$Cluster)
        a <- as.data.frame(clus.dat %>%
                             t_test_mod(t.dat = clus.dat))
        test.out <- rbind(test.out,data.frame(Cluster = all.clusters[i],a))
      }

      test.out <- test.out %>%
        rstatix::adjust_pvalue(p.col = 'p', method = 'BH', output.col =  'p.adj') %>%
        rstatix::add_significance(p.col = 'p', output.col = 'p.signif') %>%
        rstatix::add_significance(p.col = 'p.adj', output.col = 'p.adj.signif')

      test.out <- dat  %>% rstatix::group_by(Cluster) %>%
        add_y_position(test = test.out, formula = Proportion ~ Group,
                       scales = 'free', fun = 'max',step.increase = 0.1)

      p <- ggplot(dat,aes(x = Group,y = Proportion)) +
        geom_boxplot(position = position_dodge(width = 0.8),outlier.shape = NA,aes(color = Group)) +
        facet_wrap(~Cluster,ncol = n.col,scales = 'free')+
        geom_point(position = position_jitterdodge(seed = 2, dodge.width = 0.8),size=point.size,
                   aes(color = Group)) +
        scale_color_manual(values = colors) +
        scale_y_continuous(expand = expansion(mult = c(.1, .15)))+
        theme(axis.text.y = element_text(size = axis.text.size),axis.text.x = element_blank(),
              axis.title = element_text(size = axis.title.size),
              legend.title = element_text(size = legend.title.size), legend.text = element_text(size = legend.text.size),
              strip.text = element_text(size = title.text.size),axis.ticks.x = element_blank())

      if(isTRUE(plot.p.sig) & isFALSE(plot.p.adj.sig))
      {
        test.out.1 <- test.out[which(test.out[,'p.signif'] != 'ns' &  test.out[,'p.signif'] != ''),]
        if(nrow(test.out.1)!=0)
        {
          test.out.1$Cluster <- factor(test.out.1$Cluster, levels = levels(dat$Cluster))
          final.plot <- p + add_pvalue(test.out.1,
                                       y.position = "y.position",
                                       label = 'p.signif', label.size = sig.size,
                                       step.increase = 0.1, step.group.by = 'Cluster',
                                       bracket.color = sig.bracket.color,color = sig.color,
                                       bracket.size = sig.bracket.width, bracket.nudge.y = 0.001)
        }else
        {
          final.plot <- p
        }
      }else if(isTRUE(plot.p.adj.sig))
      {
        test.out.1 <- test.out[which(test.out[,'p.adj.signif'] != 'ns' &  test.out[,'p.adj.signif'] != ''),]

        if(nrow(test.out.1)!=0)
        {
          test.out.1$Cluster <- factor(test.out.1$Cluster, levels = levels(dat$Cluster))
          final.plot <- p + add_pvalue(test.out.1,
                                       y.position = "y.position",
                                       label = 'p.adj.signif', label.size = sig.size,
                                       step.increase = 0.1, step.group.by = 'Cluster',
                                       bracket.color = sig.bracket.color,color = sig.color,
                                       bracket.size = sig.bracket.width, bracket.nudge.y = 0.001)
        }else
        {
          final.plot <- p
        }
      }else
      {
        final.plot <- p
      }
    }else if(unpool.plot.type == 'barplot')
    {
      dat_bar <- suppressMessages(as.data.frame(dat %>%
                                                  group_by(Cluster,Group) %>%
                                                  drop_na() %>%
                                                  summarize(mean = mean(Proportion),
                                                            se = (sd(Proportion))/sqrt(n()))))
      # t-test
      if(!is.factor(dat$Cluster))
      {
        all.clusters <- levels(factor(dat$Cluster))
      }else
      {
        all.clusters <- levels(dat$Cluster)
      }

      test.out <- NULL
      for(i in 1:length(all.clusters))
      {
        clus.dat <- dat[which(dat$Cluster == all.clusters[i] ),]
        clus.dat$Group <- droplevels(clus.dat$Group)
        clus.dat$Cluster <- droplevels(clus.dat$Cluster)
        a <- as.data.frame(clus.dat %>%
                             t_test_mod(t.dat = clus.dat))
        test.out <- rbind(test.out,data.frame(Cluster = all.clusters[i],a))
      }

      test.out <- test.out %>%
        rstatix::adjust_pvalue(p.col = 'p', method = 'BH', output.col =  'p.adj') %>%
        rstatix::add_significance(p.col = 'p', output.col = 'p.signif') %>%
        rstatix::add_significance(p.col = 'p.adj', output.col = 'p.adj.signif')

      test.out <- dat  %>% rstatix::group_by(Cluster) %>%
        add_y_position(test = test.out, formula = Proportion ~ Group,
                       scales = 'free', fun = 'max',step.increase = 0.1)

      p <- ggplot(dat_bar,aes(x = Group,y = mean)) +
        geom_bar(stat = 'identity', color = 'black', position = position_dodge(), aes(fill = Group)) +
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                      position=position_dodge(.9))  +
        facet_wrap(~Cluster,ncol = n.col,scales = 'free')+
        scale_fill_manual(values = colors) +  labs(y = 'Mean Proportion') +
        scale_y_continuous(expand = expansion(mult = c(0, .15)))+
        theme(axis.text.y = element_text(size = axis.text.size),axis.text.x = element_blank(),
              axis.title = element_text(size = axis.title.size),
              legend.title = element_text(size = legend.title.size), legend.text = element_text(size = legend.text.size),
              strip.text = element_text(size = title.text.size),axis.ticks.x = element_blank())

      if(isTRUE(plot.p.sig) & isFALSE(plot.p.adj.sig))
      {
        test.out.1 <- test.out[which(test.out[,'p.signif'] != 'ns' &  test.out[,'p.signif'] != ''),]
        if(nrow(test.out.1)!=0)
        {
          test.out.1$Cluster <- factor(test.out.1$Cluster, levels = levels(dat$Cluster))
          final.plot <- p + add_pvalue(test.out.1,
                                       y.position = "y.position",
                                       label = 'p.signif', label.size = sig.size,
                                       step.increase = 0.1, step.group.by = 'Cluster',
                                       bracket.color = sig.bracket.color,color = sig.color,
                                       bracket.size = sig.bracket.width, bracket.nudge.y = 0.001)


        }else
        {
          final.plot <- p
        }
      }else if(isTRUE(plot.p.adj.sig))
      {
        test.out.1 <- test.out[which(test.out[,'p.adj.signif'] != 'ns' &  test.out[,'p.adj.signif'] != ''),]

        if(nrow(test.out.1)!=0)
        {
          test.out.1$Cluster <- factor(test.out.1$Cluster, levels = levels(dat$Cluster))
          final.plot <- p + add_pvalue(test.out.1,
                                       y.position = "y.position",
                                       label = 'p.adj.signif', label.size = sig.size,
                                       step.increase = 0.1, step.group.by = 'Cluster',
                                       bracket.color = sig.bracket.color,color = sig.color,
                                       bracket.size = sig.bracket.width, bracket.nudge.y = 0.001)
        }else
        {
          final.plot <- p
        }
      }else
      {
        final.plot <- p
      }
    }


    final.stat.out <- test.out[,c('Cluster','group1','group2','estimate1','estimate2','estimate',
                                  'statistic','p','p.adj','p.signif','p.adj.signif')]

    final.stat.out <- list(data = dat, stats = final.stat.out)
    # return(list(plot = suppressMessages(final.plot), stat.out = final.stat.out))
    return(list(plot = final.plot, stat.out = final.stat.out))

  }else if(method == 'pooled')
  {
    if(prop.by=='cluster')
    {
      if(!is.null(split.by))
      {
        dat_prop <- table(Idents(object),object@meta.data[,var.group],object@meta.data[,split.by])

        split.by.total <- table(object@meta.data[,var.group],object@meta.data[,split.by])
        if(is.factor(unique(object@meta.data[,split.by])))
        {
          split.by.unique <- levels(unique(object@meta.data[,split.by]))
        }else
        {
          split.by.unique <- levels(as.factor(unique(object@meta.data[,split.by])))
        }

        dat_prop_1 <- NULL
        for(s in 1:length(split.by.unique))
        {
          if(plot.all == TRUE)
          {
            dat_prop_1[[s]] <- rbind(as.matrix(dat_prop[,,s]),All=colSums(dat_prop[,,s]))
          }else
          {
            dat_prop_1[[s]] <- rbind(as.matrix(dat_prop[,,s]))
          }

          dat_prop_1[[s]] <- prop.table(dat_prop_1[[s]],margin = 1)
        }
        names(dat_prop_1) <- split.by.unique
        dat_prop <- dat_prop_1

        dat <- reshape2::melt(dat_prop)
        colnames(dat) <- c("Cluster","Group","Proportion","Split_Group")
        if(!is.factor(dat$Group))
        {
          dat$Group <- factor(dat$Group)
        }

        if(plot.all == TRUE)
        {
          dat$Cluster <- factor(dat$Cluster,levels = c(levels(Idents(object)),'All'))
        }else
        {
          dat$Cluster <- factor(dat$Cluster,levels = levels(Idents(object)))
        }

        if(!is.factor(object@meta.data[,split.by]))
        {
          object@meta.data[,split.by] <- factor(object@meta.data[,split.by])
        }
        dat$Split_Group <- factor(dat$Split_Group, levels = levels(object@meta.data[,split.by]))

        if(is.null(colors))
        {
          #colors <- hue_pal()(length(levels(factor(object@meta.data[,var.group]))))
          set.seed(random.col.seed)
          colors <- distinctColorPalette(length(levels(factor(object@meta.data[,var.group]))))
        }


        p <- ggplot(dat, aes(fill=Group, y=Proportion, x=Cluster)) + facet_wrap(~Split_Group,nrow = length(split.by.unique)) +
          geom_bar(position="stack", stat="identity",color="grey50") +
          scale_fill_manual(values = colors) +
          theme(strip.text = element_text(size = title.text.size),axis.title.y = element_text(size = axis.title.size),axis.title.x = element_blank(),
                axis.text = element_text(size = axis.text.size),
                axis.text.x = element_text(angle = axis.text.angle, hjust = axis.text.hjust, vjust = axis.text.vjust),
                legend.title = element_text(size = legend.title.size), legend.text = element_text(size = legend.text.size))+
          scale_y_continuous(expand = c(0, 0))

        final.plot <- p
      }else
      {
        dat_prop <- as.matrix(table(Idents(object),object@meta.data[,var.group]))
        if(plot.all == TRUE)
        {
          dat_prop <- rbind(dat_prop,All=colSums(dat_prop))
        }else
        {
          dat_prop <- rbind(dat_prop)
        }

        dat_prop <- prop.table(dat_prop,margin = 1)
        dat <- reshape2::melt(dat_prop)
        colnames(dat) <- c("Cluster","Group","Proportion")
        if(!is.factor(dat$Group))
        {
          dat$Group <- factor(dat$Group)
        }

        if(plot.all == TRUE)
        {
          dat$Cluster <- factor(dat$Cluster,levels = c(levels(Idents(object)),'All'))
        }else
        {
          dat$Cluster <- factor(dat$Cluster,levels = levels(Idents(object)))
        }


        if(is.null(colors))
        {
          #colors <- hue_pal()(length(levels(factor(object@meta.data[,var.group]))))
          set.seed(random.col.seed)
          colors <- distinctColorPalette(length(levels(factor(object@meta.data[,var.group]))))
        }
        p <- ggplot(dat, aes(fill=Group, y=Proportion, x=Cluster)) +
          geom_bar(position="stack", stat="identity",color="grey50") +
          scale_fill_manual(values = colors) +
          theme(strip.text = element_text(size = title.text.size),axis.title.y = element_text(size = axis.title.size),axis.title.x = element_blank(),
                axis.text = element_text(size = axis.text.size),
                axis.text.x = element_text(angle = axis.text.angle, hjust = axis.text.hjust, vjust = axis.text.vjust),
                legend.title = element_text(size = legend.title.size), legend.text = element_text(size = legend.text.size))+
          scale_y_continuous(expand = c(0, 0))

        final.plot <- p
      }

    }else if(prop.by == 'group')
    {
      if(!is.null(split.by))
      {
        dat_prop <- table(Idents(object),object@meta.data[,var.group],object@meta.data[,split.by])

        if(is.null(parent.meta))
        {
          dat_prop <- prop.table(dat_prop,c(3,2))
        }else
        {
          if(!identical(levels(parent.meta[,var.group]),levels(object@meta.data[,var.group])))
          {
            parent.meta[,var.group] <- factor(parent.meta[,var.group],levels = levels(object@meta.data[,var.group]))
          }

          split.by.total <- table(parent.meta[,var.group],parent.meta[,split.by])
          if(is.factor(unique(object@meta.data[,split.by])))
          {
            split.by.unique <- levels(unique(object@meta.data[,split.by]))
          }else
          {
            split.by.unique <- levels(as.factor(unique(object@meta.data[,split.by])))
          }

          dat_prop_1 <- dat_prop
          for(s in 1:length(split.by.unique))
          {
            dat_prop_1[,,s] <- t(apply(dat_prop[,,split.by.unique[s]],1,function(x){x/ as.numeric(split.by.total[,split.by.unique[s]])}))
          }
          dat_prop <- dat_prop_1
        }

        dat <- melt(dat_prop)
        colnames(dat) <- c("Cluster","Group","Split_Group","Proportion")
        dat$Cluster <- factor(dat$Cluster,levels = levels(Idents(object)))

        dat <- dat[!is.nan(dat$Proportion),]

        if(is.null(colors))
        {
          #colors <- hue_pal()(length(levels(factor(object@meta.data[,var.group]))))
          set.seed(random.col.seed)
          colors <- distinctColorPalette(length(levels(factor(object@meta.data[,var.group]))))
        }
        p <- ggplot(dat, aes(fill=Group, y=Proportion, x=Cluster)) + facet_wrap(~Split_Group,nrow = length(unique(dat$Split_Group))) +
          geom_bar(position="dodge", stat="identity",color="grey50") +
          scale_fill_manual(values = colors) +
          theme(strip.text = element_text(size = title.text.size),axis.title.y = element_text(size = axis.title.size),axis.title.x = element_blank(),
                axis.text = element_text(size = axis.text.size),
                axis.text.x = element_text(angle = axis.text.angle, hjust = axis.text.hjust, vjust = axis.text.vjust),
                legend.title = element_text(size = legend.title.size), legend.text = element_text(size = legend.text.size))+
          scale_y_continuous(expand = c(0, 0),limits = c(0,(max(dat$Proportion)*1.1)))

        final.plot <- p

      }else
      {
        dat_prop <- as.matrix(table(Idents(object),object@meta.data[,var.group]))
        if(is.null(parent.meta))
        {
          dat_prop <- prop.table(dat_prop,margin = 2)
        }else
        {
          if(!identical(levels(parent.meta[,var.group]),levels(object@meta.data[,var.group])))
          {
            parent.meta[,var.group] <- factor(parent.meta[,var.group],levels = levels(object@meta.data[,var.group]))
          }

          dat_prop <- t(apply(dat_prop,1,function(x){x/as.numeric(table(parent.meta[,var.group]))}))
        }

        dat <- melt(dat_prop)
        colnames(dat) <- c("Cluster","Group","Proportion")
        dat$Cluster <- factor(dat$Cluster,levels = levels(Idents(object)))

        if(is.null(colors))
        {
          #colors <- hue_pal()(length(levels(factor(object@meta.data[,var.group]))))
          set.seed(random.col.seed)
          colors <- distinctColorPalette(length(levels(factor(object@meta.data[,var.group]))))
        }
        p <- ggplot(dat, aes(fill=Group, y=Proportion, x=Cluster)) +
          geom_bar(position="dodge", stat="identity",color="grey50") +
          scale_fill_manual(values = colors) +
          theme(strip.text = element_text(size = title.text.size),axis.title.y = element_text(size = axis.title.size),axis.title.x = element_blank(),
                axis.text = element_text(size = axis.text.size),
                axis.text.x = element_text(angle = axis.text.angle, hjust = axis.text.hjust, vjust = axis.text.vjust),
                legend.title = element_text(size = legend.title.size), legend.text = element_text(size = legend.text.size))+
          scale_y_continuous(expand = c(0, 0),limits = c(0,(max(dat$Proportion)*1.1)))

        final.plot <- p
      }
    }
    # return(suppressMessages(final.plot))
    #return(final.plot)
    final.stat.out <- list(data = dat, stats = NULL)
    return(list(plot = final.plot, stat.out = final.stat.out))

  }

}
