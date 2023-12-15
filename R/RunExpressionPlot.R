#' Run Expression Plot
#'
#' A wrapper function for \code{\link[Ragas]{ExpressionPlot}}
#'
#' Plot expression levels of gene per identity class. The average expression level is calculated by the \code{group.by} argument. If \code{split.by} is set, plot will be
#' further separated by data defined by the \code{split.by} argument.
#'
#' To provide user-defined colors, a vector of colors is required. Number of colors must be equal to number of idents or
#' number of groups if \code{split.by} is used.
#'
#' @param object A Seurat object or a \code{\link[Ragas]{Pi}} object
#' @param feature Gene name
#' @param assay Name of the assay to use (default: "RNA")
#' @param group.by Name of metadata column to group cells by
#' @param split.by Name of the metadata column to split expression by (optional)
#' @param ident Cell identity variable (default: seurat_clusters)
#' @param ds.key A named list with unique identifier to retrieve data from the \code{ds} object. List names must be of format group1_group2 where group1 and group2 are groups in the split.by argument. See  \code{\link[Ragas]{RunPseudobulkAnalysis}} for more details to create ds object.
#' @param plot.ds.p.value TRUE/FALSE if p-value from differential testing should be plotted (default: FALSE)
#' @param plot.ds.adjp.value TRUE/FALSE if adjusted p-value from differential testing should be plotted (default: FALSE)
#' @param colors Colors for plotting. Number of colors must be the equal to number of idents or the number of groups if split.by is used (optional)
#' @param random.col.seed Set seed to control colors
#' @param point.size Point size for box plot
#' @param axis.text.size,legend.text.size Text size for axis or legend (default: 12)
#' @param axis.title.size,legend.title.size Title size for axis or legend (default: 15)
#' @param sig.size Asterisk size for significant tests
#' @param sig.color Asterisk color for significant tests
#' @param sig.bracket.color Bracket color for significant tests
#' @param sig.bracket.width Bracket width for significant tests
#'
#' @return A \code{\link[ggplot2]{ggplot}} object will be returned
#'
#'
#' @examples
#' \dontrun{
#' # Seurat input and user defined colors
#' library(scales)
#' cols <-  hue_pal()(length(levels(csle.pbmc.small$Groups)))
#' RunExpressionPlot(object = csle.pbmc.small,
#'                   feature = "ISG15",
#'                   ident = "cluster.annotation",
#'                   group.by = "Names",
#'                   split.by = "Groups",
#'                   colors = cols)
#'
#' # Pi input without ds p-value
#' RunExpressionPlot(object = csle.pbmc.pi,
#'                   feature = "ISG15",
#'                   ident = "subcluster_idents",
#'                   group.by = "Names",
#'                   split.by = "Groups")
#'
#' # Pi input with ds p-value
#' RunExpressionPlot(object = csle.pbmc.pi,
#'                   feature = "ISG15",
#'                   ident = "subcluster_idents",
#'                   group.by = "Names",
#'                   ds.key = list(cSLE_cHD = "DS|subcluster_idents|edgeR|group=Groups;sample=Names;gp1=cSLE;gp2=cHD;contrast=cSLE-cHD"),
#'                   plot.ds.p.value = TRUE,
#'                   plot.ds.adjp.value = FALSE,
#'                   split.by = "Groups")
#' }
#' @export

RunExpressionPlot <- function(object,
                              feature = NULL,
                              assay = 'RNA',
                              group.by = NULL,
                              split.by = NULL,
                              ident = 'seurat_clusters',
                              ds.key = NULL,
                              plot.ds.p.value = FALSE, plot.ds.adjp.value = FALSE,
                              colors = NULL,
                              random.col.seed = 1,
                              point.size = 1,
                              axis.text.size = 12,
                              axis.text.angle = 90,
                              axis.text.hjust = 0,
                              axis.text.vjust = 0,
                              axis.title.size = 15,
                              legend.text.size=12,
                              legend.title.size=15,
                              sig.size = 8,
                              sig.bracket.color = 'gray50',
                              sig.color = 'gray50',
                              sig.bracket.width = 0.5){
  if(!(is(object, "Pi") || is(object, "Seurat"))){stop("Invalid input argument \"object\". Either Seurat or Pi object is required.")}
  if(is(object, "Pi")){seurat.obj <- object$seurat.obj}
  if(is(object, "Seurat")){seurat.obj <- object}

  if(is.null(assay)){
    assay <- DefaultAssay(seurat.obj)
  }else{
    if(assay != DefaultAssay(seurat.obj)){
      DefaultAssay(seurat.obj) <- assay
    }
  }

  # check feature
  if(is.null(feature)){
    stop("Feature must be provided.")
  }

  if(ident != 'seurat_clusters')
  {
    if(!(ident %in% names(seurat.obj[[]]))){
      stop("argument \"ident\" not found in the meta data")
    }
  }

  if(is.null(group.by))
  {
    stop(" argument \"group.by\" is required.")
  }else
  {
    if(!(group.by %in% names(seurat.obj[[]]))){
      stop("argument \"group.by\" not found in the meta data")
    }
  }

  if(!is.null(split.by))
  {
    if(!(split.by %in% names(seurat.obj[[]]))){
      stop("argument \"split.by\" not found in the meta data")
    }
  }

  # plotting ds p-value
  if(!is.null(split.by) & is.null(ds.key) & (plot.ds.p.value == TRUE | plot.ds.adjp.value == TRUE))
  {
    stop("argument \"ds.key\" should be provided to plot p-value")
  }else if(!is.null(ds.key) & (plot.ds.p.value == FALSE & plot.ds.adjp.value == FALSE))
  {
    stop("argument \"plot.ds.p.value\" or \"plot.ds.adjp.value\ should be provided to plot p-value")
  }else if(!is.null(split.by) & !is.null(ds.key) & (plot.ds.p.value == TRUE | plot.ds.adjp.value == TRUE))
  {
    # check if object is a Pi object
    if(!(is(object, "Pi")))
    {
      stop("argument \"object\" should be a Pi object")
    }else if(is(object,"Pi") & is.null(names(object$ds)))
    {
      stop("The \"ds\" field of the Pi object is empty. Run RunPseudobulkAnalysis on the Pi object first.")
    }else if(is(object,"Pi") & !is.null(names(object$ds)))
    {
      key.check <- unlist(lapply(ds.key, function(x){x %in% names(object$ds)}))
      if(sum(key.check)!=length(ds.key))
      {
        key.idx <- which(key.check == FALSE)
        stop(paste("Invalid key:",paste(unlist(ds.key[key.idx]), sep = "\n"),"for the ds object"))
      }else
      {
        # check if  groups in ds.key(s) match with groups in split.by
        key.grps <- lapply(ds.key, function(x){c(object$ds[[x]]$param$gp1,object$ds[[x]]$param$gp2)})
        if(!is.factor(object$seurat.obj@meta.data[,split.by]))
        {
          object$seurat.obj@meta.data[,split.by] <- factor(object$seurat.obj@meta.data[,split.by])
        }
        key.grps <- lapply(key.grps, function(x){x %in% levels(object$seurat.obj@meta.data[,split.by])})
        key.grps <- unlist(lapply(key.grps, function(x){sum(x)}))
        key.idx <- which(key.grps !=2)
        if(length(key.idx)!=0)
        {
          stop("The groups in argument \"",ds.key[key.idx],"\" does not match groups in argument \"split.by\".")
        }

        # check if  names of ds.key(s) match with groups in split.by
        ds.key.names <- unique(unlist(lapply(names(ds.key), function(x){unlist(strsplit(x, split = '_'))})))
        idx.key.names <- match(ds.key.names, levels(object$seurat.obj@meta.data[,split.by]))
        if(sum(is.na(idx.key.names))>=1)
        {
          stop("The names of argument \"ds.key\" should match the groups in argument \"split.by\".")
        }

        # check if ident match with the ident in ds ds.key
        ds.idents <- unique(unlist(lapply(ds.key, function(x){c(object$ds[[x]]$ident)})))
        idx.ds.idents <- match(ds.idents, ident)
        if(sum(is.na(idx.ds.idents))>=1 )
        {
          stop("The ident in \"ds.key\" should match the argument \"ident\".")
        }else
        {
          res <- lapply(ds.key, function(x){object$ds[[x]][["data"]]$res$table[[1]]})
        }
      }
    }
  }else if(is.null(split.by) & !is.null(ds.key) & (plot.ds.p.value == TRUE | plot.ds.adjp.value == TRUE))
  {
    stop("argument \"split.by\" should be provided to plot p-value")
  }


  p <- ExpressionPlot(object = seurat.obj,
                      feature = feature,
                      assay = assay,
                      group.by = group.by,
                      split.by = split.by,
                      ident = ident,
                      res = res,
                      plot.ds.p.value = plot.ds.p.value, plot.ds.adjp.value = plot.ds.adjp.value,
                      colors = colors,
                      random.col.seed = random.col.seed,
                      point.size = point.size,
                      axis.text.size = axis.text.size,
                      axis.text.angle = axis.text.angle,
                      axis.text.hjust = axis.text.hjust,
                      axis.text.vjust = axis.text.vjust,
                      axis.title.size = axis.title.size,
                      legend.text.size=legend.text.size,
                      legend.title.size=legend.title.size)

  return(p)
}
