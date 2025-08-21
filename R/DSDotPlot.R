#' DotPlot for Differential State Analysis
#'
#' An internal function that visualizes differential state analysis results from \code{\link[Ragas]{RunPseudobulkAnalysis}}, which also
#' incorporates per gene expression frequencies calculated from \code{\link[Ragas]{CalculateExpFreqs}} as dot size.
#' Use its wrapper function \code{\link[Ragas]{RunDSDotPlot}} for plotting.
#'
#' @param object A Seurat object
#' @param res A list containing data.frame  with differential testing results. See \code{\link[muscat]{pbDS}}
#' @param exp.freq A list containing per gene expression frequencies. See \code{\link[Ragas]{CalculateExpFreqs}}
#' @param cpm.dat A data frame containing CPM
#' @param p.threshold P-Value cutoff
#' @param FC.threshold Fold-change cutoff
#' @param n.deg Number of genes to plot
#' @param p.adj.filter Whether to apply multiple testing correction
#' @param pct.threshold Cutoff for per gene expression frequencies
#' @param exp.filter Whether to perform additional expression filter based on median expression difference between groups (default: TRUE)
#' @param exp.filter.thr Threshold for difference between median expression from two groups being compared (default: 0)
#' @param group.col.name,sample.col.name Column names from the Seurat metadata used as group or sample variable by \code{\link[Ragas]{RunPseudobulkAnalysis}}
#' @param case.group,ref.group Names for the case and control/reference group, corresponding to group.1 and group.2 from \code{\link[Ragas]{RunPseudobulkAnalysis}}, respectively
#' @param direction Plot for either up or down-regulated genes, or both
#' @param bin.max Maximum exp proportion for binning pct expression data (default: 0.5)
#' @param bin.len Step size for binning (default: 0.1)
#' @param lfc.min,lfc.max Minimum/maximum logFC to plot
#' @param clust.row,clust.column Whether to cluster row and columns
#' @param heatmap.cols Three colors for heatmap (low, medium, high)
#' @param highlight.deg Whether to highlight significant genes
#' @param deg.col Color of the highlight
#' @param column.fontsize,row.fontsize Column and row font size
#' @param column.fontface,row.fontface Fontface for column/row labels either "plain", "bold", "italic" or "bold.italic" (default: "plain")
#' @param column.names.rotation Rotation for column names
#' @param axis.text.hjust,axis.text.vjust Horizontal/vertical justification (in [0, 1]) for axis text (default: 0.5)
#' @param legend.label.fontsize,legend.title.fontsize Font sizes for legend label and title
#'
#' @return A \code{\link[ggplot2]{ggplot}} object
#'
#' @references Hao and Hao et al. (2021). “Integrated analysis of multimodal single-cell data.” Cell. doi:10.1016/j.cell.2021.04.048
#'
#' @importFrom ggplot2 ggplot aes geom_point theme element_blank element_text element_rect element_line guides guide_legend scale_y_discrete scale_fill_gradient2 scale_color_manual
#' @importFrom scales squish
#' @importFrom dplyr  %>%
#' @importFrom ggtree ggtree layout_dendrogram
#' @importFrom patchwork plot_layout plot_spacer
#'
#' @keywords internal
#'

DSDotPlot <- function(object,
                      res,
                      exp.freq,
                      cpm.dat,
                      features = NULL,
                      p.threshold = 0.05,
                      FC.threshold = 1.5,
                      n.deg = Inf,
                      p.adj.filter = TRUE,
                      pct.threshold = 0.1,
                      direction = 'both',
                      exp.filter = TRUE,
                      exp.filter.thr = 0,
                      group.col.name = NULL,
                      sample.col.name = NULL,
                      case.group = NULL,
                      ref.group = NULL,
                      bin.max=0.5,
                      bin.len=0.1,
                      lfc.min=-2,
                      lfc.max=2,
                      clust.row=TRUE,
                      clust.column=FALSE,
                      heatmap.cols = NULL,
                      highlight.deg = TRUE,
                      deg.col = NULL,
                      column.fontsize=12,
                      row.fontsize=12,
                      column.fontface = 'plain',
                      row.fontface = 'plain',
                      column.names.rotation=90,
                      axis.text.hjust = 0.5,
                      axis.text.vjust = 0.5,
                      legend.label.fontsize=13,
                      legend.title.fontsize=15)
{

  if(is.null(features))
  {
    for(i in 1:length(res))
    {
      idx.pct <- match(res[[i]]$gene,rownames(exp.freq[[i]]))
      res[[i]]$pct.exp <- exp.freq[[i]][idx.pct,'AllCells']
    }


    if(p.adj.filter == TRUE & exp.filter == FALSE)
    {
      if(direction == 'up')
      {
        deg.genes <- lapply(res, function(x){x <- x[which(x$p_adj.loc < p.threshold & x$logFC > log2(FC.threshold) &
                                                            x$pct.exp > pct.threshold),];
        x <- x[order(abs(x$logFC),decreasing = TRUE),];
        if(nrow(x) > n.deg){x <- x[1:n.deg,]}else{x <- x}
        })
        deg.genes <- lapply(deg.genes,function(x){x$gene})
      }else if(direction == 'down')
      {
        deg.genes <- lapply(res, function(x){x <- x[which(x$p_adj.loc < p.threshold & x$logFC < -log2(FC.threshold) &
                                                            x$pct.exp > pct.threshold),];
        x <- x[order(abs(x$logFC),decreasing = TRUE),];
        if(nrow(x) > n.deg){x <- x[1:n.deg,]}else{x <- x}
        })
        deg.genes <- lapply(deg.genes,function(x){x$gene})
      }else if(direction == 'both')
      {
        deg.genes <- lapply(res, function(x){x <- x[which(x$p_adj.loc < p.threshold & abs(x$logFC) > log2(FC.threshold) &
                                                            x$pct.exp > pct.threshold),];
        x <- x[order(abs(x$logFC),decreasing = TRUE),];
        if(nrow(x) > n.deg){x <- x[1:n.deg,]}else{x <- x}
        })
        deg.genes <- lapply(deg.genes,function(x){x$gene})
      }
      all.deg.genes <- Reduce(union,deg.genes)
    }else if(p.adj.filter == FALSE & exp.filter == FALSE)
    {
      if(direction == 'up')
      {
        deg.genes <- lapply(res, function(x){x <- x[which(x$p_val < p.threshold & x$logFC > log2(FC.threshold) &
                                                            x$pct.exp > pct.threshold),];
        x <- x[order(abs(x$logFC),decreasing = TRUE),];
        if(nrow(x) > n.deg){x <- x[1:n.deg,]}else{x <- x}
        })
        deg.genes <- lapply(deg.genes,function(x){x$gene})
      }else if(direction == 'down')
      {
        deg.genes <- lapply(res, function(x){x <- x[which(x$p_val < p.threshold & x$logFC < -log2(FC.threshold) &
                                                            x$pct.exp > pct.threshold),];
        x <- x[order(abs(x$logFC),decreasing = TRUE),];
        if(nrow(x) > n.deg){x <- x[1:n.deg,]}else{x <- x}
        })
        deg.genes <- lapply(deg.genes,function(x){x$gene})
      }else if(direction == 'both')
      {
        deg.genes <- lapply(res, function(x){x <- x[which(x$p_val < p.threshold & abs(x$logFC) > log2(FC.threshold) &
                                                            x$pct.exp > pct.threshold),];
        x <- x[order(abs(x$logFC),decreasing = TRUE),];
        if(nrow(x) > n.deg){x <- x[1:n.deg,]}else{x <- x}
        })
        deg.genes <- lapply(deg.genes,function(x){x$gene})
      }
      all.deg.genes <- Reduce(union,deg.genes)
    }else if(p.adj.filter == TRUE & exp.filter == TRUE)
    {
      case.samples <- as.character(unique(object@meta.data[,sample.col.name][which(object@meta.data[,group.col.name]==case.group)]))
      ref.samples  <- as.character(unique(object@meta.data[,sample.col.name][which(object@meta.data[,group.col.name]==ref.group)]))

      cpm.dat$cluster_id <- factor(cpm.dat$cluster_id, levels = names(res))
      cpm.dat <- split(cpm.dat, f = cpm.dat$cluster_id)
      median.exp.dat <- lapply(cpm.dat, function(x){data.frame(gene = x[,'gene'],case.median = apply(x[,paste(case.samples,'.cpm',sep='')],1,
                                                                                                     FUN = function(y){median(y, na.rm = TRUE)}),
                                                               ref.median = apply(x[,paste(ref.samples,'.cpm',sep='')],1,
                                                                                  FUN = function(y){median(y, na.rm = TRUE)}))})
      median.exp.dat.1 <- lapply(median.exp.dat, function(x){data.frame(gene = x[,'gene'],
                                                                        abs.median.exp.diff = abs(x[,'case.median'] - x[,'ref.median']))})
      for(i in 1:length(res))
      {
        idx.median.exp <- match(res[[i]]$gene,median.exp.dat.1[[i]]$gene)
        res[[i]] <- cbind(res[[i]],median.exp.dat.1[[i]][idx.median.exp,])
      }

      if(direction == 'up')
      {
        deg.genes <- lapply(res, function(x){x <- x[which(x$p_adj.loc < p.threshold & x$logFC > log2(FC.threshold) &
                                                            x$pct.exp > pct.threshold &
                                                            x$abs.median.exp.diff > exp.filter.thr
        ),];
        x <- x[order(abs(x$logFC),decreasing = TRUE),];
        if(nrow(x) > n.deg){x <- x[1:n.deg,]}else{x <- x}
        })
        deg.genes <- lapply(deg.genes,function(x){x$gene})
      }else if(direction == 'down')
      {
        deg.genes <- lapply(res, function(x){x <- x[which(x$p_adj.loc < p.threshold & x$logFC < -log2(FC.threshold) &
                                                            x$pct.exp > pct.threshold &
                                                            x$abs.median.exp.diff > exp.filter.thr
        ),];
        x <- x[order(abs(x$logFC),decreasing = TRUE),];
        if(nrow(x) > n.deg){x <- x[1:n.deg,]}else{x <- x}
        })
        deg.genes <- lapply(deg.genes,function(x){x$gene})
      }else if(direction == 'both')
      {
        deg.genes <- lapply(res, function(x){x <- x[which(x$p_adj.loc < p.threshold & abs(x$logFC) > log2(FC.threshold) &
                                                            x$pct.exp > pct.threshold &
                                                            x$abs.median.exp.diff > exp.filter.thr
        ),];
        x <- x[order(abs(x$logFC),decreasing = TRUE),];
        if(nrow(x) > n.deg){x <- x[1:n.deg,]}else{x <- x}
        })
        deg.genes <- lapply(deg.genes,function(x){x$gene})
      }
      all.deg.genes <- Reduce(union,deg.genes)
    }else if(p.adj.filter == FALSE & exp.filter == TRUE)
    {
      case.samples <- as.character(unique(object@meta.data[,sample.col.name][which(object@meta.data[,group.col.name]==case.group)]))
      ref.samples  <- as.character(unique(object@meta.data[,sample.col.name][which(object@meta.data[,group.col.name]==ref.group)]))

      cpm.dat$cluster_id <- factor(cpm.dat$cluster_id, levels = names(res))
      cpm.dat <- split(cpm.dat, f = cpm.dat$cluster_id)
      median.exp.dat <- lapply(cpm.dat, function(x){data.frame(gene = x[,'gene'],case.median = apply(x[,paste(case.samples,'.cpm',sep='')],1,
                                                                                                     FUN = function(y){median(y, na.rm = TRUE)}),
                                                               ref.median = apply(x[,paste(ref.samples,'.cpm',sep='')],1,
                                                                                  FUN = function(y){median(y, na.rm = TRUE)}))})
      median.exp.dat.1 <- lapply(median.exp.dat, function(x){data.frame(gene = x[,'gene'],
                                                                        abs.median.exp.diff = abs(x[,'case.median'] - x[,'ref.median']))})
      for(i in 1:length(res))
      {
        idx.median.exp <- match(res[[i]]$gene,median.exp.dat.1[[i]]$gene)
        res[[i]] <- cbind(res[[i]],median.exp.dat.1[[i]][idx.median.exp,])
      }

      if(direction == 'up')
      {
        deg.genes <- lapply(res, function(x){x <- x[which(x$p_val < p.threshold & x$logFC > log2(FC.threshold) &
                                                            x$pct.exp > pct.threshold &
                                                            x$abs.median.exp.diff > exp.filter.thr
        ),];
        x <- x[order(abs(x$logFC),decreasing = TRUE),];
        if(nrow(x) > n.deg){x <- x[1:n.deg,]}else{x <- x}
        })
        deg.genes <- lapply(deg.genes,function(x){x$gene})
      }else if(direction == 'down')
      {
        deg.genes <- lapply(res, function(x){x <- x[which(x$p_val < p.threshold & x$logFC < -log2(FC.threshold) &
                                                            x$pct.exp > pct.threshold &
                                                            x$abs.median.exp.diff > exp.filter.thr),];
        x <- x[order(abs(x$logFC),decreasing = TRUE),];
        if(nrow(x) > n.deg){x <- x[1:n.deg,]}else{x <- x}
        })
        deg.genes <- lapply(deg.genes,function(x){x$gene})
      }else if(direction == 'both')
      {
        deg.genes <- lapply(res, function(x){x <- x[which(x$p_val < p.threshold & abs(x$logFC) > log2(FC.threshold) &
                                                            x$pct.exp > pct.threshold &
                                                            x$abs.median.exp.diff > exp.filter.thr),];
        x <- x[order(abs(x$logFC),decreasing = TRUE),];
        if(nrow(x) > n.deg){x <- x[1:n.deg,]}else{x <- x}
        })
        deg.genes <- lapply(deg.genes,function(x){x$gene})
      }
      all.deg.genes <- Reduce(union,deg.genes)
    }

  }else if(!is.null(features))
  {
    if(sum(duplicated(features))>=1)
    {
      print(paste('Features are not unique. Making them unique!'))
      features <- unique(features)
    }
    all.deg.genes <- features

    for(i in 1:length(res))
    {
      idx.pct <- match(res[[i]]$gene,rownames(exp.freq[[i]]))
      res[[i]]$pct.exp <- exp.freq[[i]][idx.pct,'AllCells']
    }

    res <- lapply(res, function(x){idx <- match(all.deg.genes,x$gene); idx <- idx[!is.na(idx)]; out <- x[idx,]; return(out)})

    if(p.adj.filter == TRUE & exp.filter == FALSE)
    {
      if(direction == 'up')
      {
        deg.genes <- lapply(res, function(x){x$gene[which(x$p_adj.loc < p.threshold & x$logFC > log2(FC.threshold) &
                                                            x$pct.exp > pct.threshold)]})
      }else if(direction == 'down')
      {
        deg.genes <- lapply(res, function(x){x$gene[which(x$p_adj.loc < p.threshold & x$logFC < -log2(FC.threshold) &
                                                            x$pct.exp > pct.threshold)]})
      }else if(direction == 'both')
      {
        deg.genes <- lapply(res, function(x){x$gene[which(x$p_adj.loc < p.threshold & abs(x$logFC) > log2(FC.threshold) &
                                                            x$pct.exp > pct.threshold)]})
      }
    }else if(p.adj.filter == FALSE & exp.filter == FALSE)
    {
      if(direction == 'up')
      {
        deg.genes <- lapply(res, function(x){x$gene[which(x$p_val < p.threshold & x$logFC > log2(FC.threshold) &
                                                            x$pct.exp > pct.threshold)]})
      }else if(direction == 'down')
      {
        deg.genes <- lapply(res, function(x){x$gene[which(x$p_val < p.threshold & x$logFC < -log2(FC.threshold) &
                                                            x$pct.exp > pct.threshold)]})
      }else if(direction == 'both')
      {
        deg.genes <- lapply(res, function(x){x$gene[which(x$p_val < p.threshold & abs(x$logFC) > log2(FC.threshold) &
                                                            x$pct.exp > pct.threshold)]})
      }
    }else if(p.adj.filter == TRUE & exp.filter == TRUE)
    {
      case.samples <- as.character(unique(object@meta.data[,sample.col.name][which(object@meta.data[,group.col.name]==case.group)]))
      ref.samples  <- as.character(unique(object@meta.data[,sample.col.name][which(object@meta.data[,group.col.name]==ref.group)]))

      cpm.dat$cluster_id <- factor(cpm.dat$cluster_id, levels = names(res))
      cpm.dat <- split(cpm.dat, f = cpm.dat$cluster_id)
      median.exp.dat <- lapply(cpm.dat, function(x){data.frame(gene = x[,'gene'],case.median = apply(x[,paste(case.samples,'.cpm',sep='')],1,
                                                                                                     FUN = function(y){median(y, na.rm = TRUE)}),
                                                               ref.median = apply(x[,paste(ref.samples,'.cpm',sep='')],1,
                                                                                  FUN = function(y){median(y, na.rm = TRUE)}))})
      median.exp.dat.1 <- lapply(median.exp.dat, function(x){data.frame(gene = x[,'gene'],
                                                                        abs.median.exp.diff = abs(x[,'case.median'] - x[,'ref.median']))})
      for(i in 1:length(res))
      {
        idx.median.exp <- match(res[[i]]$gene,median.exp.dat.1[[i]]$gene)
        res[[i]] <- cbind(res[[i]],median.exp.dat.1[[i]][idx.median.exp,])
      }

      if(direction == 'up')
      {
        deg.genes <- lapply(res, function(x){x$gene[which(x$p_adj.loc < p.threshold & x$logFC > log2(FC.threshold) &
                                                            x$pct.exp > pct.threshold &
                                                            x$abs.median.exp.diff > exp.filter.thr)]})
      }else if(direction == 'down')
      {
        deg.genes <- lapply(res, function(x){x$gene[which(x$p_adj.loc < p.threshold & x$logFC < -log2(FC.threshold) &
                                                            x$pct.exp > pct.threshold &
                                                            x$abs.median.exp.diff > exp.filter.thr)]})
      }else if(direction == 'both')
      {
        deg.genes <- lapply(res, function(x){x$gene[which(x$p_adj.loc < p.threshold & abs(x$logFC) > log2(FC.threshold) &
                                                            x$pct.exp > pct.threshold &
                                                            x$abs.median.exp.diff > exp.filter.thr)]})
      }
    }else if(p.adj.filter == FALSE & exp.filter == TRUE)
    {
      case.samples <- as.character(unique(object@meta.data[,sample.col.name][which(object@meta.data[,group.col.name]==case.group)]))
      ref.samples  <- as.character(unique(object@meta.data[,sample.col.name][which(object@meta.data[,group.col.name]==ref.group)]))

      cpm.dat$cluster_id <- factor(cpm.dat$cluster_id, levels = names(res))
      cpm.dat <- split(cpm.dat, f = cpm.dat$cluster_id)
      median.exp.dat <- lapply(cpm.dat, function(x){data.frame(gene = x[,'gene'],case.median = apply(x[,paste(case.samples,'.cpm',sep='')],1,
                                                                                                     FUN = function(y){median(y, na.rm = TRUE)}),
                                                               ref.median = apply(x[,paste(ref.samples,'.cpm',sep='')],1,
                                                                                  FUN = function(y){median(y, na.rm = TRUE)}))})
      median.exp.dat.1 <- lapply(median.exp.dat, function(x){data.frame(gene = x[,'gene'],
                                                                        abs.median.exp.diff = abs(x[,'case.median'] - x[,'ref.median']))})
      for(i in 1:length(res))
      {
        idx.median.exp <- match(res[[i]]$gene,median.exp.dat.1[[i]]$gene)
        res[[i]] <- cbind(res[[i]],median.exp.dat.1[[i]][idx.median.exp,])
      }

      if(direction == 'up')
      {
        deg.genes <- lapply(res, function(x){x$gene[which(x$p_val < p.threshold & x$logFC > log2(FC.threshold) &
                                                            x$pct.exp > pct.threshold &
                                                            x$abs.median.exp.diff > exp.filter.thr)]})
      }else if(direction == 'down')
      {
        deg.genes <- lapply(res, function(x){x$gene[which(x$p_val < p.threshold & x$logFC < -log2(FC.threshold) &
                                                            x$pct.exp > pct.threshold &
                                                            x$abs.median.exp.diff > exp.filter.thr)]})
      }else if(direction == 'both')
      {
        deg.genes <- lapply(res, function(x){x$gene[which(x$p_val < p.threshold & abs(x$logFC) > log2(FC.threshold) &
                                                            x$pct.exp > pct.threshold &
                                                            x$abs.median.exp.diff > exp.filter.thr)]})
      }
    }
  }

  if(length(all.deg.genes) == 0)
  {
    print('No genes passed cut-off criteria')
    stop()
  }
  idx.res <- lapply(res,function(x){out <- match(all.deg.genes,x$gene); out <- out[!is.na(out)]; return(out)})

  cluster.ids <- names(idx.res)
  dat <- NULL
  for(i in 1:length(cluster.ids))
  {
    tmp <- res[[cluster.ids[i]]][idx.res[[cluster.ids[i]]],c('gene','logFC','cluster_id','pct.exp')]
    tmp$deg <- 'no'
    idx.deg <- match(deg.genes[[cluster.ids[i]]],tmp$gene)
    if(length(idx.deg)!=0)
    {
      tmp$deg[idx.deg] <- 'yes'
    }
    dat[[i]] <- tmp
  }
  names(dat) <- cluster.ids
  dat.1 <- do.call('rbind',dat)
  if(sum(cluster.ids %in% c('0','1')) > 0)
  {
    dat.1$cluster_id <- factor(dat.1$cluster_id,levels = sort(as.numeric(unique(dat.1$cluster_id)),decreasing = F))
  }

  bin.seq <- seq(0,bin.max,by=bin.len)

  tmp <- cut(dat.1$pct.exp,breaks = bin.seq,include.lowest = TRUE,labels = paste(bin.seq[1:(length(bin.seq)-1)]*100,'%',sep = ''))
  tmp <- as.character(tmp)
  tmp[which(is.na(tmp))] <- paste('>',bin.seq[length(bin.seq)]*100,'%',sep = '')

  dat.1$pct.exp.bin <- factor(tmp,levels = c(paste(bin.seq[1:(length(bin.seq)-1)]*100,'%',sep = ''),
                                             paste('>',bin.seq[length(bin.seq)]*100,'%',sep = '')))

  dat.1$deg.highlight <- ifelse(dat.1$deg=='yes',dat.1$pct.exp+10,0)

  dat.2 <- reshape2::dcast(dat.1,gene ~ cluster_id,value.var = "logFC")
  rownames(dat.2) <- dat.2$gene
  dat.2 <- dat.2[,-1]

  dat.2[is.na(dat.2)] <- 0

  if(is.null(heatmap.cols))
  {
    heatmap.cols <- c('blue','white','red')
  }

  if(is.null(deg.col))
  {
    deg.col <- c('white','gray66')
  }else
  {
    deg.col <- c('white',deg.col[1])
  }

  if(clust.row == TRUE & clust.column == TRUE)
  {
    h_clust <- hclust(dist(dat.2 %>% as.matrix())) # hclust with distance matrix
    ddgram_row <- as.dendrogram(h_clust)
    ggtree_plot_row <- ggtree(ddgram_row)

    v_clust <- hclust(dist(dat.2 %>% as.matrix() %>% t())) # hclust with distance matrix
    ddgram_col <- as.dendrogram(v_clust)
    ggtree_plot_col <- ggtree(ddgram_col) + layout_dendrogram()

    dat.1$gene <- factor(dat.1$gene,levels = h_clust$labels[h_clust$order])
    dat.1$cluster_id <- factor(dat.1$cluster_id,levels = v_clust$labels[v_clust$order])

    if(highlight.deg == TRUE)
    {
      p <- ggplot(dat.1, mapping = aes(x = cluster_id, y = gene)) +
        geom_point(aes( size = pct.exp.bin,fill=logFC),shape=21,color='black',alpha=1) +
        geom_point(aes(color = deg),shape=22,size = 8)+
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
              axis.line = element_blank(),
              axis.text.x = element_text(size = column.fontsize, angle = column.names.rotation,
                                         hjust = axis.text.hjust, vjust = axis.text.vjust,
                                         face = column.fontface),
              axis.text.y = element_text(size = row.fontsize, face = row.fontface),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              panel.background = element_blank(),
              panel.grid.major = element_line(colour = 'gray70',size = 0.2,linetype = 'dashed'),
              panel.grid.minor = element_line(colour = 'gray70',size = 0.2,linetype = 'dashed'),
              legend.title = element_text(size = legend.title.fontsize), legend.text = element_text(size = legend.label.fontsize)) +
        guides(size = guide_legend(title = "Percent Expressed",nrow=2),color = guide_legend("Differentially Expressed")) +
        scale_y_discrete(position = "left") +
        scale_fill_gradient2(low = heatmap.cols[1],mid = heatmap.cols[2],high = heatmap.cols[3],limits=c(lfc.min, lfc.max), oob=squish) +
        scale_color_manual(values = deg.col)
    }else if(highlight.deg == FALSE)
    {
      p <- ggplot(dat.1, mapping = aes(x = cluster_id, y = gene)) +
        geom_point(aes( size = pct.exp.bin,fill=logFC),shape=21,color='black',alpha=1) +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
              axis.line = element_blank(),
              axis.text.x = element_text(size = column.fontsize, angle = column.names.rotation,
                                         hjust = axis.text.hjust, vjust = axis.text.vjust,
                                         face = column.fontface),
              axis.text.y = element_text(size = row.fontsize, face = row.fontface),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              panel.background = element_blank(),
              panel.grid.major = element_line(colour = 'gray70',size = 0.2,linetype = 'dashed'),
              panel.grid.minor = element_line(colour = 'gray70',size = 0.2,linetype = 'dashed'),
              legend.title = element_text(size = legend.title.fontsize), legend.text = element_text(size = legend.label.fontsize)) +
        guides(size = guide_legend(title = "Percent Expressed",nrow=2)) +
        scale_y_discrete(position = "left") +
        scale_fill_gradient2(low = heatmap.cols[1],mid = heatmap.cols[2],high = heatmap.cols[3],limits=c(lfc.min, lfc.max), oob=squish) +
        scale_color_manual(values = deg.col)
    }


    final.plot <- plot_spacer() + ggtree_plot_col + ggtree_plot_row+ p +
      plot_layout(nrow = 2,widths  = c(4, 10),heights = c(1,10))

  }else if (clust.row == TRUE & clust.column == FALSE)
  {
    h_clust <- hclust(dist(dat.2 %>% as.matrix())) # hclust with distance matrix
    ddgram_row <- as.dendrogram(h_clust)
    ggtree_plot_row <- ggtree(ddgram_row)

    dat.1$gene <- factor(dat.1$gene,levels = h_clust$labels[h_clust$order])

    if(highlight.deg == TRUE)
    {
      p <- ggplot(dat.1, mapping = aes(x = cluster_id, y = gene)) +
        geom_point(aes( size = pct.exp.bin,fill=logFC),shape=21,color='black',alpha=1) +
        geom_point(aes(color = deg),shape=22,size = 8)+
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
              axis.line = element_blank(),
              axis.text.x = element_text(size = column.fontsize, angle = column.names.rotation,
                                         hjust = axis.text.hjust, vjust = axis.text.vjust,
                                         face = column.fontface),
              axis.text.y = element_text(size = row.fontsize, face = row.fontface),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              panel.background = element_blank(),
              panel.grid.major = element_line(colour = 'gray70',size = 0.2,linetype = 'dashed'),
              panel.grid.minor = element_line(colour = 'gray70',size = 0.2,linetype = 'dashed'),
              legend.title = element_text(size = legend.title.fontsize), legend.text = element_text(size = legend.label.fontsize)) +
        guides(size = guide_legend(title = "Percent Expressed", nrow = 2),color = guide_legend("Differentially Expressed")) +
        scale_y_discrete(position = "left") +
        scale_fill_gradient2(low = heatmap.cols[1],mid = heatmap.cols[2],high = heatmap.cols[3],limits=c(lfc.min, lfc.max), oob=squish) +
        scale_color_manual(values = deg.col)
    }else if(highlight.deg == FALSE)
    {
      p <- ggplot(dat.1, mapping = aes(x = cluster_id, y = gene)) +
        geom_point(aes( size = pct.exp.bin,fill=logFC),shape=21,color='black',alpha=1) +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
              axis.line = element_blank(),
              axis.text.x = element_text(size = column.fontsize, angle = column.names.rotation,
                                         hjust = axis.text.hjust, vjust = axis.text.vjust,
                                         face = column.fontface),
              axis.text.y = element_text(size = row.fontsize, face = row.fontface),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              panel.background = element_blank(),
              panel.grid.major = element_line(colour = 'gray70',size = 0.2,linetype = 'dashed'),
              panel.grid.minor = element_line(colour = 'gray70',size = 0.2,linetype = 'dashed'),
              legend.title = element_text(size = legend.title.fontsize), legend.text = element_text(size = legend.label.fontsize)) +
        guides(size = guide_legend(title = "Percent Expressed", nrow = 2)) +
        scale_y_discrete(position = "left") +
        scale_fill_gradient2(low = heatmap.cols[1],mid = heatmap.cols[2],high = heatmap.cols[3],limits=c(lfc.min, lfc.max), oob=squish) +
        scale_color_manual(values = deg.col)
    }
    final.plot <- ggtree_plot_row+ p +
      plot_layout(nrow = 1,widths  = c(4, 10))
  }else if (clust.row == FALSE & clust.column == TRUE)
  {
    v_clust <- hclust(dist(dat.2 %>% as.matrix() %>% t())) # hclust with distance matrix
    ddgram_col <- as.dendrogram(v_clust)
    ggtree_plot_col <- ggtree(ddgram_col) + layout_dendrogram()

    dat.1$cluster_id <- factor(dat.1$cluster_id,levels = v_clust$labels[v_clust$order])

    if(highlight.deg == TRUE)
    {
      p <- ggplot(dat.1, mapping = aes(x = cluster_id, y = gene)) +
        geom_point(aes( size = pct.exp.bin,fill=logFC),shape=21,color='black',alpha=1) +
        geom_point(aes(color = deg),shape=22,size = 8)+
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
              axis.line = element_blank(),
              axis.text.x = element_text(size = column.fontsize, angle = column.names.rotation,
                                         hjust = axis.text.hjust, vjust = axis.text.vjust,
                                         face = column.fontface),
              axis.text.y = element_text(size = row.fontsize, face = row.fontface),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              panel.background = element_blank(),
              panel.grid.major = element_line(colour = 'gray70',size = 0.2,linetype = 'dashed'),
              panel.grid.minor = element_line(colour = 'gray70',size = 0.2,linetype = 'dashed'),
              legend.title = element_text(size = legend.title.fontsize), legend.text = element_text(size = legend.label.fontsize)) +
        guides(size = guide_legend(title = "Percent Expressed",nrow=2),color = guide_legend("Differentially Expressed")) +
        scale_y_discrete(position = "left") +
        scale_fill_gradient2(low = heatmap.cols[1],mid = heatmap.cols[2],high = heatmap.cols[3],limits=c(lfc.min, lfc.max), oob=squish) +
        scale_color_manual(values = deg.col)
    }else if(highlight.deg == FALSE)
    {
      p <- ggplot(dat.1, mapping = aes(x = cluster_id, y = gene)) +
        geom_point(aes( size = pct.exp.bin,fill=logFC),shape=21,color='black',alpha=1) +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
              axis.line = element_blank(),
              axis.text.x = element_text(size = column.fontsize, angle = column.names.rotation,
                                         hjust = axis.text.hjust, vjust = axis.text.vjust,
                                         face = column.fontface),
              axis.text.y = element_text(size = row.fontsize, face = row.fontface),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              panel.background = element_blank(),
              panel.grid.major = element_line(colour = 'gray70',size = 0.2,linetype = 'dashed'),
              panel.grid.minor = element_line(colour = 'gray70',size = 0.2,linetype = 'dashed'),
              legend.title = element_text(size = legend.title.fontsize), legend.text = element_text(size = legend.label.fontsize)) +
        guides(size = guide_legend(title = "Percent Expressed",nrow=2)) +
        scale_y_discrete(position = "left") +
        scale_fill_gradient2(low = heatmap.cols[1],mid = heatmap.cols[2],high = heatmap.cols[3],limits=c(lfc.min, lfc.max), oob=squish) +
        scale_color_manual(values = deg.col)
    }
    final.plot <-  ggtree_plot_col + p +
      plot_layout(nrow = 2, heights = c(1,10))
  }else if (clust.row == FALSE & clust.column == FALSE)
  {
    if(highlight.deg == TRUE)
    {
      p <- ggplot(dat.1, mapping = aes(x = cluster_id, y = gene)) +
        geom_point(aes( size = pct.exp.bin,fill=logFC),shape=21,color='black',alpha=1) +
        geom_point(aes(color = deg),shape=22,size = 8)+
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
              axis.line = element_blank(),
              axis.text.x = element_text(size = column.fontsize, angle = column.names.rotation,
                                         hjust = axis.text.hjust, vjust = axis.text.vjust,
                                         face = column.fontface),
              axis.text.y = element_text(size = row.fontsize, face = row.fontface),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              panel.background = element_blank(),
              panel.grid.major = element_line(colour = 'gray70',size = 0.2,linetype = 'dashed'),
              panel.grid.minor = element_line(colour = 'gray70',size = 0.2,linetype = 'dashed'),
              legend.title = element_text(size = legend.title.fontsize), legend.text = element_text(size = legend.label.fontsize)) +
        guides(size = guide_legend(title = "Percent Expressed",nrow=2),color = guide_legend("Differentially Expressed")) +
        scale_y_discrete(position = "left") +
        scale_fill_gradient2(low = heatmap.cols[1],mid = heatmap.cols[2],high = heatmap.cols[3],limits=c(lfc.min, lfc.max), oob=squish) +
        scale_color_manual(values = deg.col)
    }else if(highlight.deg == FALSE)
    {
      p <- ggplot(dat.1, mapping = aes(x = cluster_id, y = gene)) +
        geom_point(aes( size = pct.exp.bin,fill=logFC),shape=21,color='black',alpha=1) +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
              axis.line = element_blank(),
              axis.text.x = element_text(size = column.fontsize, angle = column.names.rotation,
                                         hjust = axis.text.hjust, vjust = axis.text.vjust,
                                         face = column.fontface),
              axis.text.y = element_text(size = row.fontsize, face = row.fontface),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              panel.background = element_blank(),
              panel.grid.major = element_line(colour = 'gray70',size = 0.2,linetype = 'dashed'),
              panel.grid.minor = element_line(colour = 'gray70',size = 0.2,linetype = 'dashed'),
              legend.title = element_text(size = legend.title.fontsize), legend.text = element_text(size = legend.label.fontsize)) +
        guides(size = guide_legend(title = "Percent Expressed",nrow=2)) +
        scale_y_discrete(position = "left") +
        scale_fill_gradient2(low = heatmap.cols[1],mid = heatmap.cols[2],high = heatmap.cols[3],limits=c(lfc.min, lfc.max), oob=squish) +
        scale_color_manual(values = deg.col)
    }
    final.plot <-  p
  }
  # return(suppressWarnings(print(final.plot)))
  return(final.plot)

}
