#' Make Sankey Plot
#'
#' A function to make Sankey plot between metadata columns for Seurat or Pi object.
#'
#' @param object A \code{\link[Ragas]{Pi}} object or Seurat object
#' @param meta.1 Name of metadata column to map from
#' @param meta.2 Name of metadata column to map to
#' @param node.colors A character vector with colors for the node. The length must be the same as the unique number of items in meta.1 and meta.2
#' @param node.outline.color Color for node outline (default: black)
#' @param node.text.size Text size for node (default: 5)
#' @param node.text.color Color for node text (default: black)
#' @param node.text.hjust Horizontal justification for node text (default: -0.05)
#' @param node.text.vjust Vertical justification for node text (default: 0.5)
#' @param axis.text.size Axis text size (default: 20)
#' @param random.col.seed Set seed to control colors (default: 1)
#'
#' @importFrom ggsankey make_long geom_sankey geom_sankey_text theme_sankey
#' @importFrom ggplot2 ggplot aes scale_fill_manual theme element_blank element_text position_nudge
#' @importFrom randomcoloR distinctColorPalette
#' @return A \code{\link[ggplot2]{ggplot}} object
#' @examples
#' \dontrun{
#' SankeyPlot(object = csle.pbmc.pi, meta.1 = 'cluster.annotation', meta.2 = 'subcluster_idents')
#' }
#' @export
#'

SankeyPlot <- function(object, meta.1, meta.2,
                       node.colors = NULL, node.outline.color = 'black',
                       node.text.size = 5, node.text.color = 'black',
                       node.text.hjust = (-0.05),node.text.vjust = 0.5,
                       axis.text.size = 20, random.col.seed = 1)
{
  if(!(is(object, "Pi") || is(object, "Seurat"))){stop("Invalid input argument \"object\". Either Seurat or Pi object is required.")}
  if(is(object, "Pi")){object <- object$seurat.obj}
  dat <- object@meta.data[,c(meta.1,meta.2)]
  dat <- suppressWarnings(dat %>% make_long(c(meta.1,meta.2)))

  if(is.null(node.colors))
  {
    set.seed(random.col.seed)
    node.colors <- distinctColorPalette(length(unique(dat$node)))
  }else if(!is.null(node.colors))
  {
    if(length(node.colors) != length(unique(dat$node)))
    {
      stop("The length of \"node.colors\" must be the same as the unique number of items in \"meta.1\" and \"meta.2\" .")
    }
  }

  final.plot <- ggplot(dat, aes(x = x,
                       next_x = next_x,
                       node = node,
                       next_node = next_node,
                       fill = factor(node),
                       label = node)) +
                geom_sankey(flow.alpha = 0.8, node.color = node.outline.color, type = 'sankey') +
                geom_sankey_text(size = node.text.size, color = node.text.color,
                                 hjust = node.text.hjust,vjust = node.text.vjust,
                                 position = position_nudge(x = 0.05)) +
                scale_fill_manual(values = node.colors)  +
                theme_sankey(base_size = 11) +
                theme(legend.position = "none", axis.title.x  = element_blank(),
                      axis.text.x = element_text(size = axis.text.size))

  return(final.plot)
}
