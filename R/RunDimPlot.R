#' Run DimPlot
#'
#' Run DimPlot for Seurat or post-integration (Pi) object
#'
#' Wrapper functiuon for \code{\link[Seurat]{DimPlot}} that is enhanced to: (1) uses \code{\link[randomcoloR]{distinctColorPalette}} to generate random colors;
#' (2) allows adjustment of overlap among identity labels. By default, column \code{seurat_clusters} from the metadata
#' will be by assigned as \code{group.by}.
#'
#' @param object A Seurat or \code{\link[Ragas]{Pi}} object
#' @param group.by Name of the metadata column to group (color) cells by (default: seurat_clusters)
#' @param reduction Which dimension reduction to use (default: umap)
#' @param label Whether to show the label of clusters (default: TRUE)
#' @param repel Whether to repel labels (default: TRUE)
#' @param pt.size Plot point size (default: 0.1)
#' @param label.size Label size (default: 3)
#' @param random.col.seed Set seed to control colors
#' @param cols User-defined colors. When cols is set, random.col.seed will be ignored
#' @param ggrepel.max.overlaps Adjust overlap for ggrepel (default: Inf)
#' @param ... Extra parameters passed to DimPlot
#'
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom ggplot2 scale_color_manual
#'
#' @references Hao and Hao et al. (2021). “Integrated analysis of multimodal single-cell data.” Cell. doi:10.1016/j.cell.2021.04.048
#' @references H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
#' @references Ron Ammar (2019). randomcoloR: Generate Attractive Random Colors. R package version 1.1.0.1.https://CRAN.R-project.org/package=randomcoloR
#' @return A \code{\link[ggplot2]{ggplot}} object
#' @examples
#' \dontrun{
#' ## default, use "seurat_clusters" to group cells
#' RunDimPlot(object = csle.pbmc.pi)
#'
#' ## use subcluster identities and re-projected UMAP
#' RunDimPlot(object = csle.pbmc.pi, group.by = "subcluster_idents", reduction = "rp", label.size = 5, pt.size = 1)
#' }
#' @export

RunDimPlot <- function(object,
                       group.by = "seurat_clusters",
                       reduction = "umap",
                       label = TRUE,
                       repel = TRUE,
                       pt.size = .1,
                       label.size = 3,
                       random.col.seed = 42,
                       cols = NULL,
                       ggrepel.max.overlaps = Inf,
                       ...){
  if(!(is(object, "Pi") || is(object, "Seurat"))){stop("Invalid input argument \"object\". Either Seurat or Pi object is required.")}
  if(is(object, "Pi")){seurat.obj <- object$seurat.obj}
  if(is(object, "Seurat")){seurat.obj <- object}

  ## manipulate grouping variable
  if(length(group.by) > 1){stop("Invalid \"group.by\" argument: more than 1 metadata column provided. Only supports 1 metadata column.")}
  if(group.by == "seurat_clusters" & !group.by %in% names(seurat.obj[[]])){stop("Cannot find \"seurat_clusters\" (default value for argument \"group.by\") in the metadata. Consider changing value for \"group.by\".")}
  if(!group.by %in% names(seurat.obj[[]])){stop("Invalid argument \"group.by\": ", group.by)}
    group.dat <- seurat.obj[[group.by, drop = TRUE]]
  if(!is.factor(group.dat)){group.dat <- as.factor(group.dat)}
  group.dat <- droplevels(x = group.dat)


  if(is.null(cols)){
    set.seed(random.col.seed)
    cols <- distinctColorPalette(length(levels(group.dat)))
    names(cols) <- levels(group.dat)
  }


  p <- DimPlot(object = seurat.obj, reduction = reduction, group.by = group.by, pt.size = pt.size,repel = repel, ...) +
    scale_color_manual(values = cols)
  if(label){
    options(ggrepel.max.overlaps = ggrepel.max.overlaps)
    p <- LabelClusters(plot = p, id = group.by, size = label.size)
  }

  return(p)

}
