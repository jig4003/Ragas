#' Run Annotated DotPlot for Seurat or post-integration (Pi) object
#'
#' A wrapper function for \code{\link[Ragas]{AnnotatedDotPlot}}
#'
#'
#' @param object A Seurat or \code{\link[Ragas]{Pi}} object
#' @param group.by A factor that contains identity information to group the cells by (default: seurat_clusters)
#' @param assay Name of the assay to use (default: the active assay as indicated by  \code{\link[Seurat]{DefaultAssay}})
#' @param annotations User defined annotations in one of the following two formats:
#' \itemize{
#' \item{list} : a named list object. Each list element is a character vector of a group of features; the name of the list element will be extracted to represent the feature group
#' \item{data.frame} : a data.frame with two columns: "features" and "annotation"
#' }
#' @param annotation.cols Colors for feature annotations. Needs to be a character vector (preferably a named vector with its names matching the "annotations" argument. See example below) of the same length as annotation feature groups; will be automatically assigned with random colors if set to NULL
#' @param features a list of features/genes to plot without annotations. Will be ignored if argument annotations is set
#' @param cols Colors to plot. Same as the "cols" argument of \code{\link[Seurat]{DotPlot}} (default: c("lightgrey", "blue"))
#' @param split.by Factor to split the groups by. Same as the "split.by" argument of \code{\link[Seurat]{DotPlot}}
#' @param clust.row Whether to cluster the rows/identities (default: TRUE)
#' @param clust.column Whether to cluster the columns/features (default: FALSE)
#' @param column.fontsize Size of column text (default: 12)
#' @param row.fontsize Size of row text (default: 12)
#' @param column.fontface,row.fontface Fontface for column/row labels either "plain", "bold", "italic" or "bold.italic" (default: "plain")
#' @param legend.label.fontsize Size of the legend labels (default: 13)
#' @param legend.title.fontsize Size of the legend title (default: 15)
#' @param random.col.seed Set seed to control colors
#' @param ... Extra parameters passed to DotPlot
#'
#' @return A \code{\link[ggplot2]{ggplot}} object
#'
#' @references Hao and Hao et al. (2021). “Integrated analysis of multimodal single-cell data.” Cell. doi:10.1016/j.cell.2021.04.048
#'
#' @examples
#' \dontrun{
#' my.list <- list(Monocytes = c("CD14","FCGR3A"),
#'                 B = "MS4A1",
#'                 T = c("CCR7", "CD8A"),
#'                 NK = "NKG7")
#' RunAnnotatedDotPlot(object = csle.pbmc.pi,
#'                     annotation.cols = c( 'Monocytes'= 'tomato', 'B' = 'seagreen', 'T' = 'steelblue', 'NK' = 'purple'),
#'                     annotations = my.list)
#' }
#'
#' @export

RunAnnotatedDotPlot <- function(object,
                                assay = NULL,
                                group.by = "seurat_clusters",
                                annotations = NULL,
                                annotation.cols = NULL,
                                features = NULL,
                                cols = c("lightgrey", "blue"),
                                split.by = NULL,
                                clust.row = TRUE,
                                clust.column = FALSE,
                                column.fontsize = 12,
                                row.fontsize = 12,
                                column.fontface = "plain",
                                row.fontface = "plain",
                                legend.label.fontsize=13,
                                legend.title.fontsize=15,
                                random.col.seed = 42,
                                ...
                                ){
  if(!(is(object, "Pi") || is(object, "Seurat"))){stop("Invalid input argument \"object\". Either Seurat or Pi object is required.")}
  if(is(object, "Pi")){seurat.obj <- object$seurat.obj}
  if(is(object, "Seurat")){seurat.obj <- object}

  if(is.null(assay)){
    assay <- DefaultAssay(seurat.obj)
  }

  ## check input args
  if(!(group.by %in% names(seurat.obj[[]]))){stop("argument \"group.by\" not found in the meta data")}

  if(is.null(annotations) && is.null(features)){stop("Arguments \"annotations\" and \"features\" are both NULL. One of them has to be set.")}

  if (!is.null(annotations)){
    if(is.list(annotations) & !is.data.frame(annotations)){
      if(is.null(names(annotations))){stop("Invalid annotations in list format: a named list is required.")}
      column.annotation <- stack(annotations)
      colnames(column.annotation) <- c("features", "annotation")
    }else if(is.data.frame(annotations)){
      if(!identical(colnames(annotations), c("features", "annotation"))){stop("Invalid column names for annotation: must be have two columns named \"features\" and \"annotation\".")}
      column.annotation <- annotations
    }else{stop("Argument \"annotation\" must be a named list or data.frame.")}
    column.annotation$features <- as.character(column.annotation$features)
    column.annotation$annotation <- as.character(column.annotation$annotation)


    if(length(which(nchar(column.annotation$annotation) == 0)) > 0){print(column.annotation); stop("Name of annotation cannot be empty character.")}
    # print(column.annotation)

    feature_exist <- column.annotation$features %in% rownames(seurat.obj)
    if(length(feature_exist) == length(which(feature_exist == FALSE))){
      stop("None of the features provided in the \"annotation\" argument exist in the Seurat object.")
    }else{
      if(length(which(feature_exist == FALSE)) > 0){
        warning("Warning: the following feature(s) ", paste(column.annotation$features[!feature_exist], collapse = ", "), " do(es) not exist in the Seurat object.")
      }
    }

    ## assign annotation color
    if(is.null(annotation.cols)){
      set.seed(random.col.seed)
      annotation.cols <- distinctColorPalette(length(unique(column.annotation$annotation)))
    }else{
      if(length(annotation.cols) != length(unique(column.annotation$annotation))){
        stop("The length of user provided \"annotation.cols\" does not match the number of feature-annotation groups")
      }
      if(!is.null(names(annotation.cols))){
        if(!all(names(annotation.cols) %in% unique(column.annotation$annotation))){
          stop("Names of \"annotation.cols\" do not match \"annotations\"!")
        }

      }
    }

    features = column.annotation$features
  }else{
    names(features) <- NULL
    feature_exist <- features %in% rownames(seurat.obj)
    if(length(feature_exist) == length(which(feature_exist == FALSE))){
      stop("None of the features provided in the \"features\" argument exist in the Seurat object.")
    }else{
      if(length(which(feature_exist == FALSE)) > 0){
        warning("Warning: the following feature(s) ", paste(features[!feature_exist], collapse = ", "), " do(es) not exist in the Seurat object.")
      }
    }
    column.annotation <- NULL
  }


  p <- AnnotatedDotPlot(object = seurat.obj,
                        features = features,
                        assay = assay,
                        split.by = split.by,
                        group.by = group.by,
                        clust.row = clust.row,
                        clust.column = clust.column,
                        cols = cols,
                        column.annotation = column.annotation,
                        column.annotation.cols = annotation.cols,
                        column.fontsize = column.fontsize,
                        row.fontsize = row.fontsize,
                        column.fontface = column.fontface,
                        row.fontface = row.fontface,
                        legend.label.fontsize = legend.label.fontsize,
                        legend.title.fontsize = legend.title.fontsize,
                        ...)
  return(p)
}
