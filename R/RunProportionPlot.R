#' Run ProportionPlot
#'
#' A wrapper function for \code{\link[Ragas]{ProportionPlot}}
#'
#' This function supports two methods for single-cell proportion analysis, namely "pooled" and "unpooled". The pooled method
#' will pool all cells across all samples and calculate cell proportion per cluster per group (as determined by the argument \code{group.by}).
#' Only bar plot is applicable for the pooled method.
#'
#' For the pooled method, proportions can be calculated by two ways, either by cluster or by group, given by the \code{pooled.prop.by} argument.
#' If \code{pooled.prop.by} is set to "cluster", proportions of cells per cluster per group will be "normalized" by the total number of cells in each cluster.
#' In this case, for a given cluster, the proportions of cells from all groups will sum up to 1. This will lead to stacked bar plot for each cluster.
#' If \code{pooled.prop.by} is set to "group", proportion of cells for per cluster per group will be "normalized" by the total number of cells in each group, combining all cell clusters.
#' This will lead to grouped bar plots.
#'
#' The "unpooled" method is suitable for multi-sample, multi-group study design, where an additional grouping argument \code{unpool.by} will decide how to split cells
#' into biological replicates. A typical choice for \code{unpool.by} is the metadata column that contains the sample/replicate IDs/names.
#' Two plotting methods are supported under the unpooled mode: box plot and bar plot. For box plot, T test of cell proportion differences will be performed between all groups
#' using the \code{\link[rstatix]{t_test}} function and significant proportional differences will be highlighted between corresponding groups; for bar plot,
#' standard errors will be plotted.
#'
#' If \code{use.parent.as.ref} is set to TRUE, cell numbers from the parent object will be use as reference (i.e., denominator) for cell proportion calculation.
#' This will give users the flexibility to investigate cell proportions for any cluster relative to a given choice of "total population". For example, when
#' comparing cell proportions of memory CD4 T cells, if \code{use.parent.as.ref} is FALSE, then the derived cell proportion per cluster means percentage of
#' cells within the total memory CD4 cells; if \code{use.parent.as.ref} is TRUE and a PBMC object is assigned as \code{parent.object} when running
#' \code{\link[Ragas]{CreatePostIntegrationObject}}, then the calculated percentage means percentage of cell in total PBMC. The \code{use.parent.as.ref}
#' argument only affect unpooled methods and pooled method with \code{pooled.prop.by} set to "group".
#'
#' @param object A Seurat or \code{\link[Ragas]{Pi}} object
#' @param group.by Name of the metadata column to group cells by
#' @param ident Cell identity variable (default: seurat_clusters)
#' @param method Method for cell proportion analysis: either "unpooled" or "pooled" (default: "unpooled")
#' @param use.parent.as.ref Whether to use cell counts from the parent object to calculate cell proportions (default: FALSE)
#' @param parent.meta.data.key A unique identifier to retrieve parent metadata
#' @param return.value Type of the return value (default: "stats")
#' @param group.colors Colors for the plotting group
#' @param pooled.prop.by Choose one of the two methods to calculate pooled cell proportion: by "cluster" or "group" (default: "cluster")
#' @param pooled.split.by Name of the metadata column to split the pooled proportion plot by (optional)
#' @param pooled.add.all Whether to add proportion for all cells (default: TRUE)
#' @param unpool.by Name of the secondary grouping variable from the metadata column, such as sample name
#' @param unpool.plot.type Method for the unpooled plot: either boxplot or barplot (default: "boxplot")
#' @param unpool.ncol Number of columns per plot (default: 5)
#' @param unpool.plot.p.sig Whether to highlight tests with significant raw p-values (default: TRUE)
#' @param unpool.plot.p.adj.sig Whether to highlight tests with significant adjusted p-values (default: FALSE)
#' @param unpool.pt.size Point size for box plot (default: 1)
#' @param unpool.sig.size,unpool.sig.bracket.color,unpool.sig.color,unpool.sig.bracket.width Additional aesthetic parameters that control highlighting of significant tests
#' @param title.text.size,axis.title.size,axis.text.size,legend.title.size,legend.text.size Text size for plot title, axis, or legend (default: 12)
#' @param axis.text.angle Rotation for axis text (default: 90)
#' @param axis.text.hjust,axis.text.vjust Horizontal/vertical justification (in [0, 1]) for axis text (default: 1)
#' @param random.col.seed Set seed to control colors
#'
#' @return When \code{return.value} is "ggplot", a ggplot object will be returned; when \code{return.value} is "stats" and \code{object} is a Seurat object,
#' a data summary list is returned. If \code{method} is "pooled", the data summary list only contains cell proportions used for the pooled plot; If \code{method} is "unpooled",
#' the summary list encloses an additional data.frame summarizing test statistics for cell proportions,
#' which contains the following 11 columns for all pairwise cell proportion comparisons:
#' \itemize{
#' \item{Cluster} : Name of cell cluster/identity
#' \item{group1} : Name of group 1
#' \item{group2} : Name of group 2
#' \item{estimate1} : Estimated proportion of group1
#' \item{estimate2} : Estimated proportion of group2
#' \item{estimate} : Estimated proportion difference, i.e., estimate1 - estimate2
#' \item{statistic} : Test statistic
#' \item{p} : P-value
#' \item{p.adj} : Adjusted P-value (method: "BH")
#' \item{p.signif} : A Boolean variable indicating whether the a test is significant based on P-value
#' \item{p.adj.signif} : A Boolean variable indicating whether the a test is significant based on adjusted P-value
#' }
#' When \code{return.value} is "stats" and \code{object} is a Pi object, an updated Pi object will be returned: a \code{\link[Ragas]{PiCellPropData}}
#' object that encloses the summary list will be created and added to the \code{cell.prop} field.
#'
#' @references Hao and Hao et al. (2021). “Integrated analysis of multimodal single-cell data.” Cell. doi:10.1016/j.cell.2021.04.048
#'
#' @examples
#' \dontrun{
#' ## Create a Pi object for b cells and assign parent object to the PBMC
#' my.pi <- CreatePostIntegrationObject(object = csle.bcell.small,
#'                                      parent.object = csle.pbmc.small)
#'
#' ## pooled, proportion by cluster
#' my.pi <- RunProportionPlot(my.pi,
#'                            ident = "cluster.annotation",
#'                            group.by = "Groups",
#'                            method = "pooled")
#'
#' ## pooled and split, proportion by group
#' my.pi <- RunProportionPlot(my.pi,
#'                            ident = "cluster.annotation",
#'                            group.by = "Groups",
#'                            method = "pooled",
#'                            pooled.prop.by = "group",
#'                            pooled.split.by = "Ethnicity")
#'
#' ## unpooled, boxplot, use parent (pbmc) to normalize cell counts
#' my.pi <- RunProportionPlot(my.pi,
#'                            ident = "cluster.annotation",
#'                            group.by = "Groups",
#'                            method = "unpooled",
#'                            unpool.by = "Names",
#'                            use.parent.as.ref = TRUE)
#'
#' ## unpooled, barplot
#' my.pi <- RunProportionPlot(my.pi,
#'                            ident = "cluster.annotation",
#'                            group.by = "Groups",
#'                            method = "unpooled",
#'                            unpool.by = "Names",
#'                            unpool.plot.type = "barplot")
#'
#' ## Input is a Seurat object, output is a ggplot
#' RunProportionPlot(csle.bcell.small,
#'                   ident = "cluster.annotation",
#'                   group.by = "Groups",
#'                   method = "unpooled",
#'                   unpool.by = "Names",
#'                   return.value = "ggplot")
#' }
#'
#' @export
#'
RunProportionPlot <- function(object,
                              group.by,
                              ident = "seurat_clusters",
                              method = c("unpooled", "pooled"),
                              group.colors = NULL,
                              use.parent.as.ref = FALSE,
                              parent.meta.data.key = NULL,
                              return.value = c("stats", "ggplot"),

                              ## parameters for pooled method
                              pooled.prop.by = c("cluster", "group"),
                              pooled.split.by = NULL,
                              pooled.add.all = TRUE,

                              ## parameters for unpooled method
                              unpool.by = NULL,
                              unpool.plot.type = c("boxplot", "barplot"),
                              unpool.ncol = 5,
                              unpool.plot.p.sig = TRUE,
                              unpool.plot.p.adj.sig = FALSE,
                              unpool.pt.szie = 1,
                              unpool.sig.size = 6,
                              unpool.sig.bracket.color = "black",
                              unpool.sig.color = "black",
                              unpool.sig.bracket.width = 0.5,

                              ## change plot aesthetics
                              title.text.size = 12,
                              axis.title.size = 12,
                              axis.text.size = 12,
                              legend.title.size = 12,
                              legend.text.size = 12,
                              axis.text.angle = 90,
                              axis.text.hjust = 1,
                              axis.text.vjust = 1,
                              random.col.seed = 1
                              ){
  if(!(is(object, "Pi") || is(object, "Seurat"))){stop("Invalid input argument \"object\". Either Seurat or Pi object is required.")}
  if(is(object, "Pi")){seurat.obj <- object$seurat.obj}
  if(is(object, "Seurat")){seurat.obj <- object}

  ## check input args
  if(!(group.by %in% names(seurat.obj[[]]))){stop("argument \"group.by\" not found in the meta data")}

  if(!(ident %in% names(seurat.obj[[]]))){
    stop("argument \"ident\" not found in the meta data")
  }else{
    Idents(seurat.obj) <- ident
    }



  if(!is.null(pooled.split.by)){
    if(!(pooled.split.by %in% names(seurat.obj[[]]))){stop("argument \"pooled.split.by\" not found in the meta data")}
  }

  if(use.parent.as.ref){
    if(is(object, "Seurat")){stop("use.parent.as.ref=TRUE is only supported for Pi input.")}
    if(length(object$parent.meta.data) == 0){stop("use.parent.as.ref is set to TRUE but parent.meta.data field of the Pi object is empty. Did you set the \"parent\" argument properly when running CreatePostIntegrationObject?")}
    if(is.null(parent.meta.data.key)){
      if(length(object$parent.meta.data) == 1){
        ## unique parent metadata
        parent.meta.data.key <- names(object$parent.meta.data)
      }else{
        stop("Multiple parent metadata exist but \"parent.meta.data.key\" is not set")
      }
    }
    if(!parent.meta.data.key %in% names(object$parent.meta.data)){
      stop("Invalid key \"",parent.meta.data.key, "\" for the parent.meta.data object. The following key(s) are available:\n    ",
           paste(names(object$parent.meta.data), collapse = "\n    "),"\nChoose an availble key.")
    }else{
      parent.meta <- object$parent.meta.data[[parent.meta.data.key]][["data"]]
    }
  }else{
    parent.meta <- NULL
  }

  method <- match.arg(method)
  return.value <- match.arg(return.value)
  pooled.prop.by <- match.arg(pooled.prop.by)

  if(!is.null(group.colors)){
    if(length(group.colors) < length(levels(factor(seurat.obj@meta.data[,group.by])))){
      stop("The length of user-provided \"group.colors\" (",length(group.colors),") is smaller than the number of levels of metadata \"", group.by, "\" (", length(levels(factor(seurat.obj@meta.data[,group.by]))),") in the Seurat object.")
    }else if(length(group.colors) > length(levels(factor(seurat.obj@meta.data[,group.by])))){
      warning("The length of user-provided \"group.colors\" (",length(group.colors),") is greater than the number of levels of metadata \"", group.by, "\" (", length(levels(factor(seurat.obj@meta.data[,group.by]))),") in the Seurat object.\n",
              "Taking the first ", length(levels(factor(seurat.obj@meta.data[,group.by]))), " color(s) for the plot.")

    }
  }

  if(method == "unpooled" & is.null(unpool.by)){stop("argument \"unpool.by\" cannot be NULL when \"method\" is \"unpooled\"")}
  if(method == "unpooled"){
    if(!(unpool.by %in% names(seurat.obj[[]]))){stop("argument \"unpool.by\" not found in the meta data")}
    unpool.plot.type <- match.arg(unpool.plot.type)
  }

  stopifnot(is.numeric(axis.text.hjust),is.numeric(axis.text.vjust),
            axis.text.hjust >=0 && axis.text.hjust <=1, axis.text.vjust >=0 && axis.text.vjust <=1)

  if(method =="pooled"){
    out <- ProportionPlot(object = seurat.obj,
                   var.group = group.by,
                   method = method,
                   plot.all = pooled.add.all,
                   colors = group.colors,
                   prop.by= pooled.prop.by,
                   split.by = pooled.split.by,

                   parent.meta = parent.meta,
                   title.text.size = title.text.size,
                   axis.title.size = axis.title.size,
                   axis.text.size = axis.text.size,
                   legend.title.size = legend.title.size,
                   legend.text.size = legend.text.size,
                   axis.text.angle = axis.text.angle,
                   axis.text.hjust = axis.text.hjust,
                   axis.text.vjust = axis.text.vjust
                   )
  }else{
    out <- ProportionPlot(seurat.obj,
                   var.group = group.by,
                   method = method,
                   colors = group.colors ,
                   parent.meta = parent.meta,

                   var.unpool = unpool.by,
                   n.col = unpool.ncol,
                   unpool.plot.type = unpool.plot.type,
                   plot.p.sig = unpool.plot.p.sig,
                   plot.p.adj.sig = unpool.plot.p.adj.sig,
                   point.size = unpool.pt.szie,
                   sig.size = unpool.sig.size,
                   sig.bracket.color = unpool.sig.bracket.color,
                   sig.color = unpool.sig.color,
                   sig.bracket.width = unpool.sig.bracket.width,

                   title.text.size = title.text.size,
                   axis.title.size = axis.title.size,
                   axis.text.size = axis.text.size,
                   legend.title.size = legend.title.size,
                   legend.text.size = legend.text.size,
                   axis.text.angle = axis.text.angle,
                   axis.text.hjust = axis.text.hjust,
                   axis.text.vjust = axis.text.vjust)
  }

  if(is(object, "Seurat")){
    if(return.value == "stats"){
      plot(out$plot)
      return(out$stat.out)
    }else{
      return(out$plot)
    }

  }else{
    if(return.value == "stats"){
      plot(out$plot)
      cell.prop <- out$stat.out

      if(!use.parent.as.ref){
        parent.info <- NULL
      }else{
        parent.info <- object[[6]][[parent.meta.data.key]][["param"]]$parent
      }

      my.method <- paste(method, ifelse(method == 'pooled',
                                        ifelse(is.null(pooled.split.by), '', paste('|split_by_',pooled.split.by, sep = '')),
                                               paste('|unpool_by_',unpool.by, sep = '')),
                                        sep = '')
      cp <- PiCellPropData(data = cell.prop,
                           ident = ident,
                           method = my.method,
                           group = group.by,
                           parent = parent.info
      )

      CheckPiData(cp, seurat.obj = object$seurat.obj)

      object <- AddPiData(object, cp)
      return(object)
    }else{
      return(out$plot)
    }

  }


}
