#' Run Summarized Heatmap for Seurat or post-integration (Pi) object
#'
#' A wrapper function for \code{\link[Ragas]{SummarizedHeatmap}}
#'
#' Plot heatmap summarizing expression levels per identity class. If \code{split.by} is set, the columns of the heatmap will be
#' further separated by data defined by the \code{split.by} argument. Additional metadata columns are also allowed to be included for column
#' annotation, provided that their levels are consistent with those of the \code{split.by} argument in a one/multiple (data to split) to one
#' (additional metadata) relationship. For example, if metadata column "patient" is used to split the columns of the heatmap, additional metadata like "age",
#' "gender", or "disease group" are allowed since one patient typically has unique values for age, gender, and the disease group; on the other hand,  if metadata
#' "disease group" is used to split the columns, one can no longer add "age", "gender", or "patient" as additional metadata since one level of the disease group
#' (i.e., healthy) can correspond to multiple patients with difference values for age and gender.
#'
#' To provide user-defined colors for column annotations, a named list of color vectors is required, with each of the list element corresponds to annotation data
#' passed by arguments "ident", "split.by", or "additional.metadata". The list element name for argument "ident" should always be "Cluster", and use the corresponding
#' metadata column names for annotations provided by "split.by" or "additional.metadata" (see example).
#'
#' @param object A Seurat or \code{\link[Ragas]{Pi}} object
#' @param features User defined feature list in one of the following two formats:
#' \itemize{
#' \item{list} : a named list object. Each list element is a character vector of a group of features; the name of the list element will be extracted to represent the feature group
#' \item{vector} : a character vector of features; no row annotation is applicable so leave \code{row.annotation.cols} to NULL
#' }
#' @param assay Name of the assay to use (default: the active assay as indicated by  \code{\link[Seurat]{DefaultAssay}})
#' @param ident Name of the metadata column to be used as cell identity (default: seurat_clusters)
#' @param split.by Factor to split the columns by
#' @param additional.metadata Names of additional metadata used for column annotation
#' @param heatmap.cols Colors for filling the heatmap (default: c("cadetblue","black","yellow"))
#' @param clust.row,clust.column Whether to cluster the rows/columns (default: TRUE)
#' @param row.annotation.cols A named list containing colors for column annotations; If set to NULL, \code{\link[randomcoloR]{distinctColorPalette}} will be used to generate random colors
#' @param column.annotation.cols A character vector containing colors for row annotations; If set to NULL, \code{\link[randomcoloR]{distinctColorPalette}} will be used to generate random colors
#' @param show.row.names,show.column.names Whether to show row/column names (default: TRUE)
#' @param min.exp,max.exp Minimum/maximum scaled expression value to plot (default: -2/2)
#' @param column.names.rotation Rotation for column names (default: 90)
#' @param row.fontsize,column.fontsize Row/column font size (default: 10/8)
#' @param row.fontface,column.fontface Fontface for row/column labels, either "plain", "bold", "italic", "oblique" or "bold.italic" (default: "plain")
#' @param legend.label.fontsize,legend.title.fontsize,annotation.name.fontsize Font sizes for legend label/title and annotation name (default: 10)
#' @param heatmap.width,heatmap.height Width and height of heatmap (default: 12/10)
#' @param random.col.seed Set seed to control row and column annotation colors (default: 1)
#' @param ... Extra parameters passed to \code{\link[ComplexHeatmap]{Heatmap}}
#'
#' @return A \code{\link[ComplexHeatmap]{Heatmap-class}} object
#'
#' @references Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics.
#'
#' @examples
#' \dontrun{
#' # Example 1: Seurat input and named list for features
#' features <- list("T cell" = c("CD8A", "CD8B"),
#'                  "ISG" = c("ISG15","ISG20","IFI44L"),
#'                  "Housekeeping" = c("ACTB"))
#'
#' RunSummarizedHeatmap(object = csle.pbmc.small,
#'                      ident = "cluster.annotation",
#'                      features = features,
#'                      split.by = "Groups")
#'
#' # Example 2: a more complex example with more metadata and user-defined column and row annotation color
#' ## make a small example data set
#' my.bcell <- subset(csle.bcell.small, subset = Names %in% c("cHD1", "CHD2", "cHD3","cSLE1", "cSLE2", "cSLE3"))
#' my.bcell$Names <- droplevels(my.bcell$Names)
#'
#' ## user-defined column annotation color
#' library(randomcoloR)
#' cols.cluster <- distinctColorPalette(length(levels(my.bcell$cluster.annotation))); names(cols.cluster) <- levels(my.bcell$cluster.annotation)
#' cols.name <- distinctColorPalette(length(levels(my.bcell$Names))); names(cols.name) <- levels(my.bcell$Names)
#' cols.group <- distinctColorPalette(length(levels(my.bcell$Groups))); names(cols.group) <- levels(my.bcell$Groups)
#' cols.batch <- distinctColorPalette(length(levels(factor(my.bcell$Batch)))); names(cols.batch) <- levels(factor(my.bcell$Batch))
#' cols.column.anno <- list("Cluster" = cols.cluster,
#'                          "Names" = cols.name,
#'                          "Groups" = cols.group,
#'                          "Batch" = cols.batch)
#' ## user-defined row/feature annotation color
#' features <- list("B cell" = c("TCL1A","ZEB2","CD27","ITGAX"),
#'                  "ISG" = c("ISG15","ISG20","IFI44L"),
#'                  "Housekeeping" = c("ACTB"))
#' cols.features <- distinctColorPalette(length(features))
#'
#' RunSummarizedHeatmap(object = my.bcell,
#'                      ident = "cluster.annotation",
#'                      features = features,
#'                      split.by = "Names",
#'                      additional.metadata = c("Groups","Batch"),
#'                      column.annotation.cols = cols.column.anno,
#'                      row.annotation.cols = cols.features )
#'
#' # Example 3: Pi input and vector input for features
#' features <- c("CD8A","CD8B","ISG15","ISG20","IFI44L","ACTB")
#' RunSummarizedHeatmap(object = csle.pbmc.pi,
#'                      ident = "subcluster_idents",
#'                      features = features,
#'                      split.by = "Groups")
#' }
#' @export
#'
# require(dplyr)
# require(Seurat)
# require(randomcoloR)
# require(ComplexHeatmap)
# require(circlize)
# require(gplots)
RunSummarizedHeatmap <- function(object,
                                 features,
                                 assay = NULL,
                                 ident = "seurat_clusters",
                                 split.by = NULL,
                                 additional.metadata = NULL,
                                 heatmap.cols = c("cadetblue","black","yellow"),
                                 clust.row = TRUE,
                                 clust.column = TRUE,
                                 row.annotation.cols = NULL,
                                 column.annotation.cols = NULL,
                                 show.row.names = TRUE,
                                 show.column.names = TRUE,
                                 min.exp=-2,
                                 max.exp=2,
                                 column.names.rotation=90,
                                 row.fontsize=10,
                                 column.fontsize=8,
                                 row.fontface = "plain",
                                 column.fontface = "plain",
                                 legend.label.fontsize=10,
                                 legend.title.fontsize=10,
                                 annotation.name.fontsize=10,
                                 heatmap.width=12,
                                 heatmap.height=10,
                                 random.col.seed = 1,
                                 ...){
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


  row.annotation <- NULL
  if(is.list(features)){ ##
    if(is.null(names(features))){stop("Invalid features in list format: a named list is required.")}
    features.df <- stack(features)
    colnames(features.df) <- c("features", "Group")
    row.annotation <- features.df[, "Group", drop = FALSE]

    if(!is.null(row.annotation.cols)){
      if(length(row.annotation.cols) != length(unique(row.annotation$Group))){
        stop("The length of user provided \"row.annotation.cols\" does not match the number of feature-annotation groups.")
      }else{
        names(row.annotation.cols) <- levels(row.annotation$Group)
        row.annotation.cols <- list(Group = row.annotation.cols)
      }
    }
    all.features <- features.df[,"features"]
  }else{
    all.features <- features
    if(!is.null(row.annotation.cols)){
      warning("User provides \"row.annotation.cols\", yet no annotation is availble from \"features\". Consider changing the \"features\" argument to a named list")
    }
  }

  ## evaluate features
  feature.exist <- all.features %in% rownames(seurat.obj)
  if(length(feature.exist) == length(which(feature.exist == FALSE))){
    stop("None of the features provided in the \"features\" argument exist in the Seurat object.")
  }else{
    if(length(which(feature.exist == FALSE)) > 0){
      warning("Warning: the following feature(s) ", paste(all.features[!feature.exist], collapse = ", "), " not found in the Seurat object.")
      ##if(length(which(feature.exist == TRUE)) == 1){stop("Need at least two valid features to make the heatmap.")}
      all.valid.features <- all.features[feature.exist]
      if(!is.null(row.annotation)){row.annotation <- row.annotation[feature.exist,,drop = FALSE]}
    }else{
      all.valid.features <- all.features
    }
  }

  if(!(ident %in% names(seurat.obj[[]]))){
    stop("argument \"ident\" not found in the meta data")
  }else{
    Idents(seurat.obj) <- ident
  }

  if(!is.null(split.by)){
    data.split.by <- seurat.obj[[split.by, drop = TRUE]]
    if(!is.factor(data.split.by)){data.split.by <- as.factor(data.split.by)}
    column.annotation <- data.frame(Cluster = rep(levels(Idents(seurat.obj)),
                                                  each = length(levels(data.split.by))),
                                    Split = rep(levels(data.split.by),
                                                length(levels(Idents(seurat.obj)))))
    colnames(column.annotation) <- c("Cluster", split.by)

    if(!is.null(additional.metadata)){
      if(!(all(additional.metadata %in% names(seurat.obj[[]])))){stop("Not all names in \"additional.metadata\" are found in the Seurat object.")}

      ## check consistency between split.by and additional metadata
      for(i in 1:length(additional.metadata)){
        .CheckMetadataSummarizedHeatmap(seurat.obj, split.by, additional.metadata[i])
      }

      for(i in 1:length(additional.metadata)){
        my.metadata1 <- seurat.obj[[split.by, drop = TRUE]]
        my.metadata2 <- seurat.obj[[additional.metadata[i], drop = TRUE]]
        m1 <- match(column.annotation[,2], my.metadata1)
        column.annotation[[2+i]] <- my.metadata2[m1]
      }
      colnames(column.annotation) <- c("Cluster", split.by, additional.metadata)

    }

  }else{
    if(!is.null(additional.metadata)){message("Ignore \"additional.metadata\" becasue \"split.by\" is NULL.")}
    column.annotation <- data.frame(Cluster = levels(Idents(seurat.obj)))
  }

  # print(column.annotation)
  if(!is.null(column.annotation.cols)){
    if(!all(sort(colnames(column.annotation))==sort(names(column.annotation.cols)))){
      stop("Names of \"column.annotation.cols\" do not match column annotation. \nRun ? RunSummarizedHeatmap to see an example of how to provide user-defined colors for column annotation.")
    }
  }

  p <- SummarizedHeatmap(object = seurat.obj,
                    features = all.valid.features,
                    assay = assay,
                    heatmap.cols = heatmap.cols,
                    split.by = split.by,
                    clust.row = clust.row,
                    clust.column = clust.column,
                    row.annotation = row.annotation,
                    row.annotation.cols = row.annotation.cols,
                    column.annotation = column.annotation,
                    column.annotation.cols = column.annotation.cols,
                    show.column.names = show.column.names,
                    show.row.names = show.row.names,
                    min.exp = min.exp,
                    max.exp = max.exp,
                    column.names.rotation = column.names.rotation,
                    column.fontsize = column.fontsize,
                    row.fontsize = row.fontsize,
                    column.fontface = column.fontface,
                    row.fontface = row.fontface,
                    legend.label.fontsize = legend.label.fontsize,
                    legend.title.fontsize = legend.title.fontsize,
                    annotation.name.fontsize = annotation.name.fontsize,
                    heatmap.width = heatmap.width,
                    heatmap.height = heatmap.height,
                    random.col.seed = random.col.seed,
                    ...

  )
  return(p)
}
