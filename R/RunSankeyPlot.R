#' Run Sankey Plot
#'
#' A wrapper function for \code{\link[Ragas]{SankeyPlot}}
#'
#' This function takes 2 metadata column from seurat object and creates a sankey diagram
#'
#' @param object  A \code{\link[Ragas]{Pi}} object or Seurat object
#' @param meta.1 Name of metadata column to map from
#' @param meta.2 Name of metadata column to map to
#' @param node.colors A character vector with colors for the node. The length must be the same as the unique number of items in meta.1 and meta.2
#' @param node.outline.color Color for node outline (default: black)
#' @param node.text.size Text size for node (default: 5)
#' @param node.text.color Color for node text (default: black)
#' @param node.text.hjust Horizontal justification (in [0, 1]) for node text
#' @param node.text.vjust Vertical justification (in [0, 1]) for node text
#' @param axis.text.size Axis text size
#' @param random.col.seed Set seed to control colors
#'
#' @return A \code{\link[ggplot2]{ggplot}} object
#' @export
#'
#' @examples
#' \dontrun{
#' ## Sankey plot
#' RunSankeyPlot(csle.pbmc.pi,
#'               meta.1 = "cluster.annotation", meta.2 = "subcluster_idents")
#' }
#'
RunSankeyPlot <- function(object, meta.1, meta.2,
                          node.colors = NULL, node.outline.color = 'black',
                          node.text.size = 5, node.text.color = 'black',
                          node.text.hjust = 0,node.text.vjust = 0.5,
                          axis.text.size = 20, random.col.seed = 1)
{
  if(!(is(object, "Pi") || is(object, "Seurat"))){stop("Invalid input argument \"object\". Either Seurat or Pi object is required.")}
  if(is(object, "Pi")){seurat.obj <- object$seurat.obj}
  if(is(object, "Seurat")){seurat.obj <- object}

  if(!(meta.1 %in% names(seurat.obj[[]])) & !(meta.2 %in% names(seurat.obj[[]]))){
    stop("argument \"meta.1\" and \"meta.2\" not found in the meta data")
  }else if(!(meta.1 %in% names(seurat.obj[[]])) & (meta.2 %in% names(seurat.obj[[]])))
  {
    stop("argument \"meta.1\" not found in the meta data")
  }else if((meta.1 %in% names(seurat.obj[[]])) & !(meta.2 %in% names(seurat.obj[[]])))
  {
    stop("argument \"meta.2\" not found in the meta data")
  }

  SankeyPlot(object = seurat.obj, meta.1, meta.2,
             node.colors = node.colors, node.outline.color = node.outline.color,
              node.text.size = node.text.size, node.text.color = node.text.color,
              node.text.hjust = node.text.hjust,node.text.vjust = node.text.vjust,
              axis.text.size = axis.text.size, random.col.seed = random.col.seed)

}
