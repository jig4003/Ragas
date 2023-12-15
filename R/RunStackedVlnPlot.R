#' Run Stacked VlnPLot for Seurat or post-integration (Pi) object
#'
#' A wrapper function for the \code{\link[Ragas]{StackedVlnPlot}}
#'
#' The Stacked VlnPlot provides several new functionalities compared to \code{\link[Seurat]{VlnPlot}}:
#' (1) stacked output when multiple features are provided;(2) allows feature annotation; (3) allows clustering
#' of rows and/or columns; (4) provides multiple way to color volin plots, see argument \code{color.by} for more details.
#'
#' @param object A Seurat or \code{\link[Ragas]{Pi}} object
#' @param assay Name of the assay to use (default: "RNA")
#' @param ident Name of the metadata column to be used as cell identity (default: seurat_clusters)
#' @param features User defined feature list in one of the following three formats:
#' \itemize{
#' \item{list} : a named list object. Each list element is a character vector of a group of features; the name of the list element will be extracted to represent the feature group
#' \item{data.frame}: a data.frame with two columns: "features" and "annotation"
#' \item{vector} : a character vector of features
#' }
#' @param feature.annotation.cols A character vector of colors for feature annotation. Its length must equal the number of feature groups provided by the "features" argument.
#' If the vector is named, its names should match the "features" argument (default: NULL)
#' @param split.by Factor to split the columns by (default: NULL)
#' @param color.by One of the five ways to color for the violin plot: "features", "clusters", "median.exp" (median expression), "mean.exp" (mean expression), "split.var" (the variable that splits the violin plot).
#' When \code{split.by} is set (not to NULL), \code{color.by} must set to "split.var"; and vice versa. (default: "features")
#' @param clust.row Whether to cluster the rows (default: FALSE)
#' @param clust.column Whether to cluster the columns (default: TRUE)
#' @param add.points Whether to add data points to the plot (default: FALSE)
#' @param points.size Point size (default: 0.1)
#' @param column.fontsize Column font size (default: 12)
#' @param row.fontsize Row font size (default: 8)
#' @param row.title.fontsize Row title font size (default: 15)
#' @param legend.fontsize Legend font size (default: 12)
#' @param legend.title.fontsize Legend title font size (default: 15)
#' @param features.fontsize Feature font size (default: 12)
#' @param column.names.rotation Rotation for column names (default: 0)
#' @param axis.text.hjust,axis.text.vjust Horizontal/vertical justification (in [0, 1]) for axis text (default: 1)
#' @param fill.colors Colors to fill the violin plot when \code{color.by} is set to "median.exp" or "mean.exp"
#' @param random.annotation.col.seed Seed to generate random colors for feature annotations
#'
#' @return A \code{\link[ggplot2]{ggplot}} object
#' @export
#' @references Hao and Hao et al. (2021). “Integrated analysis of multimodal single-cell data.” Cell. doi:10.1016/j.cell.2021.04.048
#'
#' @examples
#' \dontrun{
#' my.list <- list(Monocytes = c("CD14","FCGR3A"),
#'                 B = "MS4A1",
#'                 T = c("CCR7", "CD8A"),
#'                 NK = "NKG7")
#'
#' ## Seurat input; split.by is NULL; color by mean expression
#' RunStackedVlnPlot(csle.pbmc.small,
#'                   ident = "cluster.annotation",
#'                   features = my.list,
#'                   feature.annotation.cols = c( 'Monocytes'= 'tomato', 'B' = 'seagreen', 'T' = 'steelblue', 'NK' = 'purple'),
#'                   color.by = "mean.exp",
#'                   column.names.rotation = 90)
#' ## Pi input; split.by is not NULL
#' RunStackedVlnPlot(csle.pbmc.pi,
#'                   ident = "cluster.annotation",
#'                   features = my.list,
#'                   split.by = "Groups",
#'                   color.by = "split.var",
#'                   column.names.rotation = 90)
#' }
#'
RunStackedVlnPlot <- function(object,
                              assay = "RNA",
                              ident = "seurat_clusters",
                              features = NULL,
                              feature.annotation.cols = NULL,
                              split.by = NULL,
                              color.by = "features",
                              clust.row = FALSE,
                              clust.column = TRUE,
                              add.points = FALSE,
                              points.size = 0.1,
                              column.fontsize = 12,
                              row.fontsize = 8,
                              row.title.fontsize = 15,
                              legend.fontsize = 12,
                              legend.title.fontsize = 15,
                              features.fontsize = 12,
                              column.names.rotation = 0,
                              axis.text.hjust = 1,
                              axis.text.vjust = 1,
                              fill.colors = NULL,
                              random.annotation.col.seed = 1){
  if(!(is(object, "Pi") || is(object, "Seurat"))){stop("Invalid input argument \"object\". Either Seurat or Pi object is required.")}
  if(is(object, "Pi")){seurat.obj <- object$seurat.obj}
  if(is(object, "Seurat")){seurat.obj <- object}

  if(!(ident %in% names(seurat.obj[[]]))){
    stop("argument \"ident\" not found in the meta data")
  }else{
    Idents(seurat.obj) <- ident
  }

  #color.by <- match.arg(color.by)

  if(!is.list(features)){## without annotations
    stopifnot(is.character(features))
    all.features <- features
    feature.annotation = NULL
  }else{ ## with annotations
    if(is.list(features) & !is.data.frame(features)){ ## annotation is named list
      if(is.null(names(features))){stop("Invalid \"features\" in list format: a named list is required.")}
      feature.annotation <- stack(features)
      colnames(feature.annotation) <- c("features", "annotation")
    }else{
      if(!identical(colnames(features), c("features", "annotation"))){stop("Invalid column names for \"features\": the data.frame must be have two columns named \"features\" and \"annotation\".")}
      feature.annotation <- features
    }

    feature.annotation$features <- as.character(feature.annotation$features)
    feature.annotation$annotation <- as.character(feature.annotation$annotation)
    all.features <- feature.annotation[, "features"]

    if(length(which(nchar(feature.annotation$annotation) == 0)) > 0){print(feature.annotation); stop("Name of annotation cannot be empty character.")}
  }


  ## Check whether features in seurat.obj
  feature.exist <- all.features %in% rownames(seurat.obj)
  if(length(feature.exist) == length(which(feature.exist == FALSE))){
    stop("None of the features provided in the \"features\" argument exist in the Seurat object.")
  }else{
    if(length(which(feature.exist == FALSE)) > 0){ ## at least one missing feature
      warning("Warning: the following feature(s) ", paste(all.features[!feature.exist], collapse = ", "), " do(es) not exist in the Seurat object.")
      all.valid.features <- all.features[feature.exist]

      if(!is.null(feature.annotation)){
        feature.annotation <- feature.annotation[feature.exist,,drop = FALSE]
      }
    }else{## no missing feature
      all.valid.features <- all.features
    }
  }
  # if(!is.null(feature.annotation)){
  #   print(feature.annotation)
  # }

  if(!is.null(feature.annotation.cols)){
    if(length(feature.annotation.cols) != length(unique(feature.annotation$annotation))){
      stop("The length of user provided \"feature.annotation.cols\" does not match the number of feature-annotation groups")
    }
    if(!is.null(names(feature.annotation.cols))){
      if(!all(names(feature.annotation.cols) %in% unique(feature.annotation$annotation))){
        stop("Names of \"feature.annotation.cols\" do not match \"feature.annotation\"!")
        }
    }
  }

  if(!is.null(split.by))
  {
    if(is.null(color.by))
    {
      print("Warning: setting color.by = 'split.var'")
      color.by <- 'split.var'
    }else if (color.by %in% c("features", "clusters", "median.exp", "mean.exp"))
    {
      stop("Error: split.by is compatible only with color.by = split.var")
    }
  }else if(is.null(split.by) & color.by == 'split.var')
  {
    stop("Error: color.by = split.var is compatible only with split.by")
  }



  StackedVlnPlot(object = seurat.obj,
                 features = all.valid.features,
                 feature.annotation = feature.annotation,
                 feature.annotation.cols = feature.annotation.cols,
                 assay = assay,
                 split.by = split.by,
                 color.by = color.by,
                 clust.row = clust.row,
                 clust.column = clust.column,
                 add.points = add.points,
                 points.size = points.size,
                 column.fontsize = column.fontsize,
                 row.fontsize = row.fontsize,
                 row.title.fontsize = row.title.fontsize,
                 legend.fontsize = legend.fontsize,
                 legend.title.fontsize = legend.title.fontsize,
                 features.fontsize = features.fontsize,
                 column.names.rotation = column.names.rotation,
                 colors = fill.colors,
                 axis.text.hjust = axis.text.hjust,
                 axis.text.vjust = axis.text.vjust,
                 random.annotation.col.seed = random.annotation.col.seed)

}
