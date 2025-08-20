#' Run Matrix Plot
#'
#' A wrapper function for \code{\link[Ragas]{MatrixPlot}}
#'
#' This function takes pre-calculated results from the \code{markers} field of the Pi object and plots heatmap of
#' top genes for each identity group (e.g., cluster). Running \code{\link[Ragas]{RunFindAllMarkers}} is always required prior to run this function,
#' which will store new  \code{\link[Seurat]{FindAllMarkers}} results in the \code{markers} object with a unique identifier.
#' This unique identifier needs to be passed to the \code{markers.key} argument so the function can pull relevant differential
#' expression data for plotting.
#'
#' A new assay with assay name determined by argument \code{alt.assay.name} will be created, followed by
#' RC normalization and scaling. When scaling, no variables are regressed out so the scaled expression will be more consistent with that from
#' \code{\link[Seurat]{FeaturePlot}}. The new assay only contains \code{top.n} genes from each identity group (e.g., clusters, in most cases) to
#' reduce computational cost.
#'
#'
#' @param object A \code{\link[Ragas]{Pi}} object
#' @param markers.key A unique identifier for \code{\link[Ragas]{RunFindAllMarkers}} to pull results from the \code{markers} field
#' @param alt.assay.name Name of the alternative assay to store RC normalized data
#' @param top.n Number of top markers to plot (per identity group)
#' @param up.genes Whether to plot the up-regulated genes (default: TRUE)
#' @param heatmap.cols Colors to plot the heatmap (optional). A character vector with colors for low, medium and high expression.
#' @param min.exp The minimum expression value to plot (default: -1.5)
#' @param max.exp The maximum expression value to plot (default: 1.5)
#' @param column.fontsize Font size for column texts (default: 6)
#' @param row.fontsize Font size for row texts (default: 8)
#' @param column.fontface Fontface for column labels, either "plain", "bold", "italic", "oblique" or "bold.italic" (default: "plain")
#' @param row.fontface Fontface for column labels, either "plain", "bold", "italic", "oblique" or "bold.italic" (default: "plain")
#' @param column.anno.cols Named character vector with colors for column annotations (default: NULL). Names should match Idents of object.
#' @param column.anno.name.fontsize Font size for column annotations, such as clusters (default: 10)
#' @param column.anno.name.rot Rotation of column annotations (default: 0)
#' @param column.anno.name.just Alignment of column annotations (default: "left")
#' @param legend.label.fontsize Font size for legend label (default: 10)
#' @param legend.title.fontsize Font size for legend title (default: 10)
#' @param heatmap.width Width of the heatmap (default: 20)
#' @param heatmap.height Height of the heatmap (default: 8)
#'
#' @return A \code{\link[ComplexHeatmap]{Heatmap-class}} object
#'
#' @references Hao and Hao et al. (2021). “Integrated analysis of multimodal single-cell data.” Cell. doi:10.1016/j.cell.2021.04.048
#' @references Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics.
#' @export
#'
#' @examples
#' \dontrun{
#' ## Run matrix plot for pre-calculated markers in the example data
#' RunMatrixPlot(csle.pbmc.pi,
#'               markers.key = "Markers|subcluster_idents|AllMarkers|test.use=wilcox",
#'               column.anno.name.rot = 45)
#' }
#'
RunMatrixPlot <- function(object,
                          markers.key,
                          alt.assay.name = "RNAalt",
                          top.n = 5,
                          up.genes = TRUE,
                          heatmap.cols = NULL,
                          min.exp = -1.5,
                          max.exp = 1.5,
                          column.fontsize = 6,
                          row.fontsize = 8,
                          column.fontface = "plain",
                          row.fontface = "plain",
                          column.anno.cols = NULL,
                          column.anno.name.fontsize = 10,
                          column.anno.name.rot = 0,
                          column.anno.name.just="left",
                          legend.label.fontsize = 10,
                          legend.title.fontsize = 10,
                          heatmap.width = 20,
                          heatmap.height = 8
                          ){
  if(is(object, "Pi")){
    object <- unclass(object)
    if(length(object$markers) == 0){
      stop("The \"markers\" field of the Pi object is empty. Run RunFindAllMarkers on the Pi object first.")
    }else{
      if(!markers.key %in% names(object$markers)){
        stop("Invalid key \"",markers.key, "\" for the markers object. The following key(s) are available:\n    ", paste(names(object$markers), collapse = "\n    "),"\nEither choose an availble key or rerun RunFindAllMarkers.")
      }else{
        my.markers <- object$markers[[markers.key]][["data"]]
        my.marker.type <- unlist(strsplit(markers.key, split = "|", fixed = TRUE))[3]
        my.marker.idents <- unlist(strsplit(markers.key, split = "|", fixed = TRUE))[2]

        if(!identical(Idents(object$seurat.obj), object$seurat.obj[[my.marker.idents, drop = TRUE]])){
          Idents(object$seurat.obj) <- my.marker.idents
          message("Set active identity to ", my.marker.idents)
        }
      }
    }

    stopifnot(top.n > 0, is.numeric(top.n))
    ## only run re-normalization on top markers
    markers.split0 <- split(my.markers, f = my.markers$cluster)
    if(up.genes){
      markers.split <- lapply(markers.split0, function(x) x %>%
                                dplyr::filter(avg_log2FC > 0) %>%
                                dplyr::arrange(desc(avg_log2FC))
      )
    }else{
      markers.split <- lapply(markers.split0, function(x) x %>%
                                dplyr::filter(avg_log2FC < 0) %>%
                                dplyr::arrange(avg_log2FC)
      )
    }
    n.list <- unlist(lapply(markers.split, function(x) dim(x)[1]))

    n.deg.list <- pmin(n.list, top.n)

    top.markers.list <- lapply(seq_along(markers.split), function(i) markers.split[[i]][1:n.deg.list[i], "gene"])

    top.markers.list[which(is.na(top.markers.list))] <- NULL
    top.markers <- unique(unlist(top.markers.list))

    ## Re-run normalization using RC method
    dat <- GetAssayData(object$seurat.obj,assay = 'RNA',slot = 'counts')
    RNA.alt <- CreateAssayObject(count = dat)
    object$seurat.obj[[alt.assay.name]] <- RNA.alt
    DefaultAssay(object$seurat.obj) <- alt.assay.name
    object$seurat.obj <- NormalizeData(object$seurat.obj,normalization.method = 'RC')
    object$seurat.obj <- ScaleData(object$seurat.obj, features = top.markers)

    all.marker.types <- c('all', 'conserved')
    names(all.marker.types) <- c('AllMarkers','AllConservedMarkers')

    p <- MatrixPlot(object = object$seurat.obj,
                     markers.list = my.markers,
                     n.genes = top.n,
                     type = all.marker.types[my.marker.type],
                     up.genes = up.genes,
                     down.genes = !up.genes,
                     heatmap.cols = heatmap.cols,
                     min.exp = min.exp,
                     max.exp = max.exp,
                     column.fontsize = column.fontsize,
                     row.fontsize = row.fontsize,
                     column.fontface = column.fontface,
                     row.fontface = row.fontface,
                     column.anno.cols = column.anno.cols,
                     column.anno.name.fontsize = column.anno.name.fontsize,
                     column.anno.name.rot = column.anno.name.rot,
                     column.anno.name.just = column.anno.name.just,
                     legend.label.fontsize = legend.label.fontsize,
                     legend.title.fontsize = legend.title.fontsize,
                     heatmap.width = heatmap.width,
                     heatmap.height = heatmap.height)


    return(p)

  }else{
    stop('Invalid input object. Post Integration (Pi) object is required.')
  }
}
