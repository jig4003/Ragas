## S3 classes Pi, PiDataList and PiData, and their methods
## constructors

#' The Post Integration (Pi) Class.
#'
#' An object of class Pi is a list of six fields that covers a variety of analyses for single cell RNA-seq.
#' These six fields are: (1) "seurat.obj", a \code{\link[Seurat]{Seurat}} object that contains processed data (i.e., data that are filtered, normalized,
#' batch-integrated, clustered and dimension-reduced); (2) "exp.freq", a list that contains per-gene expression frequency, which is the output of
#' \code{\link[Ragas]{CalculateExpFreqs}}; (3) "markers", a data.frame that contains results from differential gene expression analysis, such as markers
#' between different clusters, by running \code{\link[Ragas]{RunFindAllMarkers}}; (4) "ds", differential state analysis results for multi-sample, multi-group study
#' design by running \code{\link[Ragas]{RunPseudobulkAnalysis}}; (5) "cell.prop", differential cell proportion analysis results from \code{\link[Ragas]{RunProportionPlot}};
#' (6) "Parent.meta.data", a data.frame containing metadata of from the parent Seurat object.
#'
#' @param object A Seurat object
#' @param exp.freq A \code{\link[Ragas]{PiDataList}} containing calculated per gene expression frequency
#' @param markers A \code{\link[Ragas]{PiDataList}} containing expression markers
#' @param ds A \code{\link[Ragas]{PiDataList}} containing differential state analysis results by pseudobulk method
#' @param cell.prop A \code{\link[Ragas]{PiDataList}} containing differential analysis results on cell proportions
#' @param parent.meta.data A \code{\link[Ragas]{PiDataList}} containing Parent metadata
#' @seealso \code{\link[Ragas]{CreatePostIntegrationObject}} for a formal way to initialize a Pi object.
#'
#' Other related classes: \code{\link[Ragas]{PiData}}, \code{\link[Ragas]{PiExpFreqData}}, \code{\link[Ragas]{PiMarkerData}},
#' \code{\link[Ragas]{PiDSData}}, \code{\link[Ragas]{PiCellPropData}}, \code{\link[Ragas]{PiParentMetaData}}, and
#' \code{\link[Ragas]{PiDataList}}.
#' @export
#' @examples
#' ## Create a minimum Pi object from a Seurat object
#' Pi(object = csle.pbmc.small)
#'

Pi <- function(object,
               exp.freq = NULL,
               markers = NULL,
               ds = NULL,
               cell.prop = NULL,
               parent.meta.data = NULL){
  stopifnot(inherits(object, "Seurat"))
  pi.obj <- list(seurat.obj = object,
                 exp.freq = list(),
                 markers = list(),
                 ds = list(),
                 cell.prop = list(),
                 parent.meta.data = list()
  )

  if(!is.null(exp.freq)){
    pi.obj[["exp.freq"]] <- exp.freq
  }

  if(!is.null(markers)){
    pi.obj[["markers"]] <- markers
  }

  if(!is.null(ds)){
    pi.obj[["ds"]] <- ds
  }

  if(!is.null(cell.prop)){
    pi.obj[["cell.prop"]] <- cell.prop
  }

  if(!is.null(parent.meta.data)){
    pi.obj[["parent.meta.data"]] <- parent.meta.data
  }

  class(pi.obj) <- "Pi"

  return(pi.obj)

}

#' PiDataList Class
#'
#' PiDataList is a simple list of \code{\link[Ragas]{PiData}} objects
#' @param ... \code{\link[Ragas]{PiData}} objects or additional parameters to \code{PiDataList} methods
#'
#' @export
#'

PiDataList <- function(...){
  x <- list(...)
  if(length(x) == 0){stop("Empty PiDataList")}
  if(!all(unlist(lapply(x, function(u) is.PiData(u))))){stop("Input for PiDataList must be an object of class PiData")}
  if(length(x) > 1){
    type <- x[[1]][["type"]]
    if(!all(unlist(lapply(x, function(u) u[["type"]] == type)))){stop("Cannot create PiDataList from different types of PiData")}
  }
  names(x) <- unlist(lapply(x, function(u) GetPiDataName(u)))
  class(x) <- "PiDataList"
  x
}

#' PiData Class
#'
#' PiData class is the most fundamental class for Ragas analysis, which is the parent class for \code{\link[Ragas]{PiExpFreqData}}, \code{\link[Ragas]{PiMarkerData}},
#' \code{\link[Ragas]{PiDSData}}, \code{\link[Ragas]{PiCellPropData}}, and \code{\link[Ragas]{PiParentMetaData}}. When multiple analyses of the same type have been performed,
#' A \code{\link[Ragas]{PiDataList}} object will be created, which assembles several \code{PiData} objects of the same analysis/data type into one list.
#' A \code{\link[Ragas]{Pi}} object consists of different types of \code{PiDataList}.
#'
#' @param data Analysis data. Depending on the type of the analysis, \code{data} can be a list or data.frame
#' @param type A character string indicating the type of the analysis ("un-specified", "exp.freq", "markers", "ds", "cell.prop", or "parent.meta.data")
#' @param ident Name for the identity data
#' @param method A character string describing the analysis method
#' @param param Additional descriptions of the analysis parameters
#'
#' @seealso Other related classes: \code{\link[Ragas]{PiExpFreqData}}, \code{\link[Ragas]{PiMarkerData}},
#' \code{\link[Ragas]{PiDSData}}, \code{\link[Ragas]{PiCellPropData}}, \code{\link[Ragas]{PiParentMetaData}}, and
#' \code{\link[Ragas]{PiDataList}}.
#'
#' @export
#'

PiData <- function(data,
                   type = c("un-specified", "exp.freq", "markers", "ds", "cell.prop", "parent.meta.data"),
                   ident,
                   method = NULL,
                   param = NULL){
  type <- match.arg(type)
  output <- list(data = data,
                 type = type,
                 ident = ident,
                 method = method,
                 param = param)
  class(output) <- "PiData"
  return(output)
}

## classes that inherit PiData
### PiExpFreqData

#' PiExpFreqData Class
#'
#' The PiExpFreqData inherts \code{\link[Ragas]{PiData}} class, which is used to store per gene expression frequency data calculated by \code{\link[Ragas]{CalculateExpFreqs}}.
#'
#' @param data A list containing per-gene expression frequencies for each identity class
#' @param group Factor name to group cells by
#' @param cutoff Log-expression cutoff (default: 0)
#' @inheritParams PiData
#' @export
#' @examples
#' \dontrun{
#' ## Calculate gene expression frequency and construct a PiExpFreqData object
#' res <- CalculateExpFreqs(csle.bcell.small,
#'                          ident = "seurat_clusters",
#'                          group.by = NULL,
#'                          cutoff = 0)
#' ef <- PiExpFreqData(data = res,
#'                     ident = "seurat_clusters",
#'                     group = NULL,
#'                     cutoff = 0)
#'
#' ## Print object
#' ef
#'
#' ## Check object validity
#' CheckPiData(ef, seurat.obj = csle.bcell.small)
#'
#' ## Add PiExpFreqData object to an existing Pi object
#' my.pi <- CreatePostIntegrationObject(csle.bcell.small) ## a minimum Pi object
#' my.pi <- AddPiData(my.pi, ef)
#' }
#'

PiExpFreqData <- function(data,
                          ident,
                          group = NULL,
                          cutoff = 0){
  if(is.null(group)){
    output <- list(data = data,
                   type = "exp.freq",
                   ident = ident,
                   method = NULL,
                   param = list("cutoff" = cutoff))
  }else{
    output <- list(data = data,
                   type = "exp.freq",
                   ident = ident,
                   method = NULL,
                   param = list("group" = group,
                                "cutoff" = cutoff))
  }

  class(output) <- c("PiExpFreqData", "PiData")
  return(output)
}


### PiMarkerData
#' PiMarkerData Class
#'
#' The PiMarkerData inherts \code{\link[Ragas]{PiData}} class, which stores marker data from \code{\link[Ragas]{RunFindAllMarkers}}.
#' @param data A data.frame containing output from \code{\link[Ragas]{RunFindAllMarkers}}
#' @param test.use Which test to use. See details in \code{\link[Seurat]{FindAllMarkers}}
#' @param latent.vars  Latent variables to test. See details in \code{\link[Seurat]{FindAllMarkers}}
#' @inheritParams PiData
#' @export
#' @examples
#' \dontrun{
#' ## Run marker analysis and construct a PiMarkerData object
#' res <- RunFindAllMarkers(csle.bcell.small,
#'                          ident = "seurat_clusters",
#'                          test.use = "wilcox",
#'                          latent.vars = NULL)
#' mk <- PiMarkerData(data = res,
#'                    ident = "seurat_clusters",
#'                    method = "AllMarkers",
#'                    test.use = "wilcox",
#'                    latent.vars = NULL)
#' ## Print object
#' mk
#'
#' ## Check object validity
#' CheckPiData(mk, seurat.obj = csle.bcell.small)
#'
#' ## Add PiMarkerData object to an existing Pi object
#' my.pi <- CreatePostIntegrationObject(csle.bcell.small) ## a minimum Pi object
#' my.pi <- AddPiData(my.pi, mk)
#' }
#'

PiMarkerData <- function(data,
                         ident,
                         method = c("AllMarkers"),
                         test.use,
                         latent.vars = NULL){
  method <- match.arg(method)
  if(is.null(latent.vars)){
    output <- list(data = data,
                   type = "markers",
                   ident = ident,
                   method = method,
                   param = list("test.use" = test.use))
  }else{
    output <- list(data = data,
                   type = "markers",
                   ident = ident,
                   method = method,
                   param = list("test.use" = test.use,
                                "latent.vars" = paste(latent.vars, collapse = "+")))
  }

  class(output) <- c("PiMarkerData", "PiData")
  return(output)
}


### PiDSData
#' PiDSData Class
#'
#' The PiDSData inherts \code{\link[Ragas]{PiData}} class, which stores differential state analysis results from \code{\link[Ragas]{RunPseudobulkAnalysis}}.
#' @param data Output from \code{\link[Ragas]{RunPseudobulkAnalysis}}
#' @param group Name of the group from Seurat metadata for differential state analysis
#' @param sample Name of the metadata column to identify samples
#' @param group.1 Group name for case
#' @param group.2 Group name for control
#' @param contrast A character string summarizing the comparison
#' @param block Blocking variable, such as Batch
#' @inheritParams PiData
#' @export
#' @examples
#' \dontrun{
#' ## Run pseudobulk analysis and construct a PiDSData object
#' gp1 <- "cSLE"; gp2 <- "cHD";
#' res <- RunPseudobulkAnalysis(csle.bcell.small,
#'                              ident.var = "seurat_clusters",
#'                              group.var = "Groups",
#'                              sample.var = "Names",
#'                              group.1 = gp1,
#'                              group.2 = gp2,
#'                              pbDS.method = "edgeR")
#' ds <- PiDSData(data = res,
#'                ident = "seurat_clusters",
#'                method = "edgeR",
#'                group = "Groups",
#'                sample = "Names",
#'                group.1 = gp1,
#'                group.2 = gp2,
#'                contrast = paste(gp1, "-",  gp2, sep = ""),
#'                block = NULL)
#'
#' ## Print object
#' ds
#'
#' ## Check object validity
#' CheckPiData(ds, seurat.obj = csle.bcell.small)
#'
#' ## Add PiDSData object to an existing Pi object
#' my.pi <- CreatePostIntegrationObject(csle.bcell.small) ## a minimum Pi object
#' my.pi <- AddPiData(my.pi, ds)
#' }
#'
PiDSData <- function(data,
                     ident,
                     group,
                     sample,
                     block = NULL,
                     method,
                     group.1,
                     group.2,
                     contrast){
  if(is.null(block)){
    output <- list(data = data,
                   type = "ds",
                   ident = ident,
                   method = method,
                   param = list("group" = group,
                                "sample" = sample,
                                "gp1" = group.1,
                                "gp2" = group.2,
                                "contrast" = contrast))
  }else{
    output <- list(data = data,
                   type = "ds",
                   ident = ident,
                   method = method,
                   param = list("group" = group,
                                "sample" = sample,
                                "gp1" = group.1,
                                "gp2" = group.2,
                                "contrast" = contrast,
                                "block" = block))
  }

  class(output) <- c("PiDSData", "PiData")
  return(output)
}

### PiCellPropData
#' PiCellPropData Class
#'
#' The PiCellPropData inherts \code{\link[Ragas]{PiData}} class, which stores differential cell proportion analysis results from \code{\link[Ragas]{RunProportionPlot}}.
#' @param data A data.frame containing test results from \code{\link[Ragas]{RunProportionPlot}}
#' @param group Name of the metadata column to group cells by
#' @param parent Name of the parent object
#' @inheritParams PiData
#' @export
#' @examples
#' \dontrun{
#' ## Run cell proportion analysis on a Seurat object and construct a PiCellPropData object
#' ## In this case, use.parent.as.ref = TRUE is not supported.
#' ## Run ?RunProportionPlot to see how to use parent metadata for proportion analysis
#' res <- RunProportionPlot(csle.bcell.small,
#'                          ident = "seurat_clusters",
#'                          group.by = "Groups",
#'                          method = "unpooled",
#'                          unpool.by = "Names")
#' cp <- PiCellPropData(data = res,
#'                      ident = "seurat_clusters",
#'                      method = "unpooled",
#'                      group = "Groups",
#'                      parent = NULL)
#' ## Print object
#' cp
#'
#' ## Check object validity
#' CheckPiData(cp, seurat.obj = csle.bcell.small)
#'
#' ## Add PiCellPropData object to an existing Pi object
#' my.pi <- CreatePostIntegrationObject(csle.bcell.small) ## a minimum Pi object
#' my.pi <- AddPiData(my.pi, cp)
#' }
#'

PiCellPropData <- function(data,
                           ident,
                           method,
                           group,
                           by,
                           parent = NULL){
  if(is.null(parent)){
    output <- list(data = data,
                   type = "cell.prop",
                   ident = ident,
                   method = method,
                   param = list("group" = group))
  }else{
    output <- list(data = data,
                   type = "cell.prop",
                   ident = ident,
                   method = method,
                   param = list("group" = group,
                                "parent" = parent))
  }


  class(output) <- c("PiCellPropData", "PiData")
  return(output)
}


### PiParentMetaData
#' PiParentMetaData Class
#'
#' The PiParentMetaData inherts \code{\link[Ragas]{PiData}} class, which stores metadata of the parent Seurat object.
#' @param data A data.frame for metadata
#' @param parent.name Name of the parent object
#' @export
#' @examples
#' \dontrun{
#' ## Construct a PiParentMetaData object from the parent Seurat object csle.pbmc.small
#' pm <- PiParentMetaData(data = csle.pbmc.small[[]],
#'                        parent.name = "csle.pbmc.small")
#'
#' ## Print object
#' pm
#'
#' ## Check object validity
#' ## Cells in Seurat object csle.bcell.small is a subset of cells from its PBMC parent
#' CheckPiData(pm, csle.bcell.small)
#'
#' ## Add PiParentMetaData object to an existing Pi object
#' my.pi <- CreatePostIntegrationObject(csle.bcell.small) ## a minimum Pi object
#' my.pi <- AddPiData(my.pi, pm)
#' }
#'
PiParentMetaData <- function(data,
                             parent.name){
  output <- list(data = data,
                 type = "parent.meta.data",
                 ident = NULL,
                 method = NULL,
                 param = list("parent" = parent.name))
  class(output) <- c("PiParentMetaData", "PiData")
  return(output)
}



## check and validate class
is.Pi <- function(x){
  ## TBA
  return("TBA")
}


is.PiDataList <- function(x,
                          type = NULL){
  data <- unclass(x)

  if(inherits(x, "PiDataList") &&
     all(unlist(lapply(data, function(x) is.PiData(x))))){
    if(is.null(type)){
      return(TRUE)
    }else{
      if(all(unlist(lapply(data, function(x) inherits(x, type))))){
        return(TRUE)
      }else{
        return(FALSE)
      }
    }

  }else{
    return(FALSE)
  }
}

is.PiData <- function(x){
  if(inherits(x, "PiData") &&
     identical(sort(names(x)), sort(c("data", "type", "ident", "method", "param")))
  ){
    return(TRUE)
  }else{
    return(FALSE)
  }
}


is.PiExpFreqData <- function(x){
  if(inherits(x, "PiData") &&
     inherits(x, "PiExpFreqData") &&
     identical(sort(names(x)), sort(c("data", "type", "ident", "method", "param")))
  ){
    return(TRUE)
  }else{
    return(FALSE)
  }
}


is.PiMarkerData <- function(x){
  if(inherits(x, "PiData") &&
     inherits(x, "PiMarkerData") &&
     identical(sort(names(x)), sort(c("data", "type", "ident", "method", "param")))
  ){
    return(TRUE)
  }else{
    return(FALSE)
  }
}


is.PiDSData <- function(x){
  if(inherits(x, "PiData") &&
     inherits(x, "PiDSData") &&
     identical(sort(names(x)), sort(c("data", "type", "ident", "method", "param")))
  ){
    return(TRUE)
  }else{
    return(FALSE)
  }
}


is.PiCellPropData <- function(x){
  if(inherits(x, "PiData") &&
     inherits(x, "PiCellPropData") &&
     identical(sort(names(x)), sort(c("data", "type", "ident", "method", "param")))
  ){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

is.PiParentMetaData <- function(x){
  if(inherits(x, "PiData") &&
     inherits(x, "PiParentMetaData") &&
     identical(sort(names(x)), sort(c("data", "type", "ident", "method", "param")))
  ){
    return(TRUE)
  }else{
    return(FALSE)
  }
}


#' Check Validity for PiData
#'
#' Check the validity of PiData/PiDataList/Pi object
#'
#' @param x A PiData object
#' @param ... Additional arguments
#'
#' @export
#' @seealso Usage of \code{CheckPiData} for different classes: \code{\link[Ragas]{PiExpFreqData}}, \code{\link[Ragas]{PiMarkerData}},
#' \code{\link[Ragas]{PiDSData}}, \code{\link[Ragas]{PiCellPropData}}, \code{\link[Ragas]{PiParentMetaData}},
#' \code{\link[Ragas]{PiDataList}}, \code{\link[Ragas]{Pi}}.
#'
CheckPiData <- function(x, ...){
  UseMethod("CheckPiData", x)
}

#' @rdname CheckPiData
#' @param seurat.obj A Seurat object
#' @param verbose Whether to print messages (default: FALSE)
#' @method CheckPiData default
#' @export
#'
CheckPiData.default <- function(x,
                                seurat.obj,
                                verbose = FALSE){
  stopifnot(inherits(seurat.obj, "Seurat"))
  if(verbose && length(x) == 0){
    cat("(skipping check for empty list)")
  }
  if(length(x) > 0){stop("Unexpected list size")}

  return(invisible(x = NULL))
}



#' @rdname PiExpFreqData
#' @param x A PiExpFreqData object
#' @param seurat.obj A Seurat object
#' @method CheckPiData PiExpFreqData
#' @export
#'
CheckPiData.PiExpFreqData <- function(x,
                                      seurat.obj){
  stopifnot(x$type == "exp.freq", ## whether type description matches
            x$ident %in% names(seurat.obj[[]]), ## whether ident is in seurat.obj metadata column
            inherits(seurat.obj, "Seurat"))

  if(!is.factor(seurat.obj[[x$ident]])){
    seurat.obj[[x$ident]] <- factor(seurat.obj[[x$ident, drop = TRUE]])
  }

  ## whether names of list match identity levels
  stopifnot(identical(sort(names(x$data)), sort(levels(seurat.obj[[x$ident, drop = TRUE]]))))

  ## check per identity matrix size, column names, row names
  stopifnot(all(unlist(lapply(x$data, function(x) dim(x)[1] == dim(seurat.obj)[1]))),
            all(unlist(lapply(x$data, function(x) length(colnames(x)) > 0))),
            all(unlist(lapply(x$data, function(x) colnames(x)[1] == "AllCells"))),
            all(unlist(lapply(x$data, function(x) identical(rownames(x), rownames(seurat.obj)))))
  )


  return(invisible(x = NULL))
}

#' @rdname PiMarkerData
#' @param x A PiMarkerData object
#' @param seurat.obj A Seurat object
#' @method CheckPiData PiMarkerData
#' @export
#'
CheckPiData.PiMarkerData <- function(x,
                                     seurat.obj){
  stopifnot(x$type == "markers", ## whether type description matches
            x$ident %in% names(seurat.obj[[]]), ## whether ident is in seurat.obj metadata column
            inherits(seurat.obj, "Seurat"))

  if(!is.factor(seurat.obj[[x$ident]])){
    seurat.obj[[x$ident]] <- factor(seurat.obj[[x$ident, drop = TRUE]])
  }

  data <- x$data

  if(x$method == "AllMarkers"){
    if(!is.data.frame(data)){stop("Data must be a data.frame")}
    if(!(all(colnames(data) %in% c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene")))){
      stop("Invalid column names")
    }
    if(!(all(data$cluster %in% levels(seurat.obj[[x$ident, drop = TRUE]])))){
      stop("Invalid cluster name")
    }

    if(!(all(data$gene %in% rownames(seurat.obj)))){
      stop("Invalid gene name")
    }

  }else{
    stop("Invalid analysis method")
  }

  return(invisible(x = NULL))
}


#' @rdname PiDSData
#' @param x A PiDSData object
#' @param seurat.obj A Seurat object
#' @method CheckPiData PiDSData
#' @export
#'
CheckPiData.PiDSData <- function(x,
                                 seurat.obj){
  stopifnot(x$type == "ds", ## whether type description matches
            x$ident %in% names(seurat.obj[[]]), ## whether ident is in seurat.obj metadata column
            inherits(seurat.obj, "Seurat"))
  if(!is.factor(seurat.obj[[x$ident]])){
    seurat.obj[[x$ident]] <- factor(seurat.obj[[x$ident, drop = TRUE]])
  }
  data <- x$data

  # check res object
  if((!all(names(data$res$table[[1]]) %in% levels(seurat.obj[[x$ident, drop = TRUE]]))) ||
     (!all(unique(data$tbl$cluster_id) %in% levels(seurat.obj[[x$ident, drop = TRUE]])))){
    cat(sort(names(data$res$table[[1]])),"\n")
    cat(sort(unique(data$tbl$cluster_id)),"\n")
    cat(sort(levels(seurat.obj[[x$ident, drop = TRUE]])), "\n")
    stop("Cluster identity mismatch")
  }

  all.genes <- Reduce(union, lapply(data$res$table[[1]], function(u) u$gene))
  if(!(all(all.genes %in% rownames(seurat.obj))) ||
     !( all(data$tbl$gene %in% rownames(seurat.obj)))){
    stop("Genes don't match")
    }

  return(invisible(x = NULL))
}

#' @rdname PiCellPropData
#' @param x A PiCellPropData object
#' @param seurat.obj A Seurat object
#' @method CheckPiData PiCellPropData
#' @export
#'
CheckPiData.PiCellPropData <- function(x,
                                       seurat.obj){
  stopifnot(x$type == "cell.prop", ## whether type description matches
            x$ident %in% names(seurat.obj[[]]), ## whether ident is in seurat.obj metadata column
            inherits(seurat.obj, "Seurat"))
  if(!is.factor(seurat.obj[[x$ident]])){
    seurat.obj[[x$ident]] <- factor(seurat.obj[[x$ident, drop = TRUE]])
  }
  freqs <- x$data$data
  stats <- x$data$stats

  if(length(grep('unpooled', x$method)) == 1){
    # check freqs (TBA)
    # check stats
    my.columns <- c("Cluster", "group1", "group2", "estimate1", "estimate2","estimate", "statistic", "p", "p.adj", "p.signif", "p.adj.signif")
    if(!all(colnames(stats) == my.columns)){
      stop("T-test result for cell proportions must have the following columns: ", paste(my.columns, collapse = ", "))
    }
    my.params <- x$param
    my.group <- my.params[["group"]]
    if(!(my.group %in% names(seurat.obj[[]]))){
      stop("Cannot find grouping variable ", my.group, " in the seurat object")
    }
    if(!(all(stats$group1 %in% unique(seurat.obj[[my.group, drop = TRUE]]))) |
       !(all(stats$group2 %in% unique(seurat.obj[[my.group, drop = TRUE]]))) ){
      stop("Inconsistent group levels")
    }

  }

  return(invisible(x = NULL))
}

#' @rdname PiParentMetaData
#' @param x A PiParentMetaData object
#' @param seurat.obj A Seurat object
#' @method CheckPiData PiParentMetaData
#' @export
#'
CheckPiData.PiParentMetaData <- function(x,
                                         seurat.obj){
  stopifnot(x$type == "parent.meta.data", ## whether type description matches
            inherits(seurat.obj, "Seurat"))
  data <- x$data
  if(!(all(colnames(seurat.obj) %in% rownames(data)))){
    stop("Not all cells in the Seurat object are from the parent oject: ",x$param[["parent"]])
  }

  return(invisible(x = NULL))
}


#' @rdname PiDataList
#' @param x A PiDataList object
#' @param seurat.obj A Seurat object
#' @method CheckPiData PiDataList
#' @export
#'
CheckPiData.PiDataList <- function(x,
                                   seurat.obj,
                                   ...){
  stopifnot(inherits(seurat.obj, "Seurat"))
  x <- unclass(x)
  lapply(x, function(u) CheckPiData(u, seurat.obj))

  return(invisible(x = NULL))
}

#' @rdname Pi
#' @param x A Pi object
#' @param seurat.obj A Seurat object
#' @param verbose Whether to print messages (default: FALSE)
#' @method CheckPiData Pi
#'
#' @examples
#' \dontrun{
#' ## Check the validity of Pi object
#' CheckPiData(csle.pbmc.pi, verbose = TRUE)}
#'
#' @export
#'
CheckPiData.Pi <- function(x,
                           verbose = FALSE,
                           ...){
  x <- unclass(x)
  seurat.obj <- x$seurat.obj
  pidatalist <- x[2:6]
  lapply(1:length(pidatalist), function(i) {
    if(verbose){cat("Checking",names(pidatalist)[i], "... ")}
    CheckPiData(pidatalist[[i]], seurat.obj, verbose = verbose)
    if(verbose){cat("\n")}
    }
    )
  if(verbose){cat("Check for PiData successful!\n")}

  return(invisible(x = NULL))
}

## print

#' @rdname Pi
#' @param x A Pi object
#' @param ... Additional arguments
#' @method print Pi
#' @export
#'
#' @importFrom crayon bold red yellow green blue magenta cyan
#' @examples
#' ## print the Pi object
#' csle.pbmc.pi
#'

print.Pi <- function(x, ...){
  cat("An object of class", class(x),"\n")
  x <- unclass(x)
  cat(length(x), " fields in the object: ", paste(names(x), collapse = ", "),".\n", sep = "")
  x.is.not.empty.list <- unlist(lapply(x, function(x) length(x) != 0))
  idx.is.not.empty.list <- which(x.is.not.empty.list)
  if(length(idx.is.not.empty.list) > 1){
    cat("The following fields have been processed:\n")
  }else{
    cat("The following field has been processed:\n")
  }
  for(i in 1:length(idx.is.not.empty.list)){
    idx <- idx.is.not.empty.list[i]
    obj <- x[[idx]]
    if(names(x)[idx] == "seurat.obj"){
      cat("\t", bold(red("seurat.obj")), ": A Seurat object of ", dim(obj)[1], " features and ", dim(obj)[2],
          " cells.\n\t\t", length(Seurat::Assays(obj)), ifelse(length(Seurat::Assays(obj)) == 1, " assay: ", " assays: "), paste(Seurat::Assays(obj), collapse = ", "), ", and ",
          length(Reductions(obj)), ifelse(length(Reductions(obj)) == 1, " reduction: ", " reductions: "), paste(Reductions(obj), collapse = ", "),"\n", sep = "")
    }
    if(names(x)[idx] == "exp.freq"){
      cat("\t", bold(yellow("exp.freq")), ": A list of numeric matrices containing per gene expression frequencies\n\t\t",
          length(obj), " analysis ", ifelse(length(obj) == 1, "run", "runs"), ": \n\t\t  ", paste(names(obj), collapse = "\n\t\t  "), "\n", sep = "")
    }

    if(names(x)[idx] == "markers"){
      cat("\t", bold(green("markers")), ": A list of data frames containing marker results.\n\t\t",
          length(obj), " analysis ", ifelse(length(obj) == 1, "run", "runs"), ": \n\t\t  ", paste(names(obj), collapse = "\n\t\t  "), "\n", sep = "")
    }

    if(names(x)[idx] == "ds"){
      cat("\t", bold(blue("ds")), ": A list of lists and data frames containing pseudobulk analysis results.\n\t\t",
          length(obj), " analysis ", ifelse(length(obj) == 1, "run", "runs"), ": \n\t\t  ", paste(names(obj), collapse = "\n\t\t  "), "\n", sep = "")
    }

    if(names(x)[idx] == "cell.prop"){
      cat("\t", bold(magenta("cell.prop")), ": A list of data frames containing cell proportion analysis results\n\t\t",
          length(obj), " analysis ", ifelse(length(obj) == 1, "run", "runs"), ": \n\t\t  ", paste(names(obj), collapse = "\n\t\t  "), "\n", sep = "")
    }

    if(names(x)[idx] == "parent.meta.data"){
      cat("\t", bold(cyan("parent.meta.data")), ": A list of data frames containing parent metadata\n\t\t",
          length(obj), " metadata: \n\t\t  ", paste(names(obj), collapse = "\n\t\t  "), "\n", sep = "")
    }
  }
  if(length(x$parent.meta.data) == 0){
    cat("Metadata from the parent object provided?", red("No"), "\n")
  }else{
    cat("Metadata from the parent object provided?", green("Yes"), "\n")
  }
  if(!("subcluster_idents" %in% names(x$seurat.obj[[]]))){
    cat("Subclusters integrated?", red("No"))
  }else{
    cat("Subclusters integrated?", green("Yes"))
  }

}


#' @rdname PiDataList
#' @param x A PiDataList object
#' @method print PiDataList
#' @export
print.PiDataList <- function(x, ...){
  x <- unclass(x)
  type <- x[[1]]$type
  data.type <- c("analysis", "metadata")
  cat(length(x), type, "object(s):\n")
  for(i in 1:length(x)){
    cat(" +", ifelse(type == "parent.meta.data", data.type[2], data.type[1]), i,"\n")
    writeLines(paste("  ", capture.output(print(x[[i]])), sep=""))
  }
}


#' @rdname PiData
#' @param x A PiData object
#' @param ... Additional arguments
#' @method print PiData
#' @export
print.PiData <- function(x, ...){
  cat("type:", x$type, "\n")
  if(!is.null(x$ident)){
    cat("identity:", x$ident, "\n")
  }

  if(!is.null(x$method)){
    cat("method:", x$method, "\n")
  }
  if(!is.null(x$param)){
    param.str <- ParamListToStr(x$param)
    cat("param:", param.str, "\n")
  }
}

#' @rdname PiExpFreqData
#' @param x A PiExpFreqData object
#' @param ... Additional arguments
#' @method print PiExpFreqData
#' @export
print.PiExpFreqData <- function(x, ...){
  cat("A PiExpFreqData object\n")
  cat(length(x$data), "lists,",  dim(x$data[[1]])[1], "features,", dim(x$data[[1]])[2], "column(s)\n")
  NextMethod()
}

#' @rdname PiMarkerData
#' @param x A PiMarkerData object
#' @param ... Additional arguments
#' @method print PiMarkerData
#' @export
print.PiMarkerData <- function(x, ...){
  cat("A PiMarkerData object\n")
  if(x$method == "AllMarkers"){
    cat(dim(x$data)[1], "markers for",  length(unique(x$data$cluster)), "clusters\n")
  }
  NextMethod()
}

#' @rdname PiDSData
#' @param x A PiDSData object
#' @param ... Additional arguments
#' @method print PiDSData
#' @export
print.PiDSData <- function(x, ...){
  cat("A PiDSData object\n")
  cat(length(x$data$res$table[[1]]), " lists for contrast \"", names(x$data$res$table), "\"\n", sep = "")
  NextMethod()
}

#' @rdname PiCellPropData
#' @param x A PiCellPropData object
#' @param ... Additional arguments
#' @method print PiCellPropData
#' @export
print.PiCellPropData <- function(x, ...){
  cat("A PiCellPropData object\n")
  cat("Comparisons of cell proportions in", length(unique(x$data$stats$Cluster)),"clusters between", paste(unique(c(x$data$group1, x$data$group2)), collapse = ", "), "\n")
  NextMethod()
}

#' @rdname PiParentMetaData
#' @param x A PiParentMetaData object
#' @param ... Additional arguments
#' @method print PiParentMetaData
#' @export
print.PiParentMetaData <- function(x, ...){
  cat("A PiParentMetaData object\n")
  cat(dim(x$data)[2], "metadata columns for", dim(x$data)[1], "cells\n")
  NextMethod()
}


## extract and replacement

#' @rdname Pi
#' @param x A Pi object
#' @param i Index for list element to extract or replace
#' @method [[ Pi
#' @export
#' @examples
#' ## Get a sub-object
#' csle.pbmc.pi[["exp.freq"]]
#' csle.pbmc.pi[["ds"]]
#'
`[[.Pi` <- function(x,
                    i){
  if(length(i) > 1){stop("Index i must be a length-one numeric or character")}
  if(is.character(i)){
    if(!(i %in% c("seurat.obj", "exp.freq", "markers", "ds", "cell.prop", "parent.meta.data"))){
      stop("Invalid Pi object field name. Choose one of the following: seurat.obj, exp.freq, markers, ds, cell.prop, parent.meta.data.")
    }
  }
  # if(is.numeric(i)){
  #   if(i < 1 | i > 6){stop("Invlaid index i. The index i should be integer between 1 and 6.")}
  # }
  NextMethod()
}

#' @rdname PiDataList
#' @param  i Index for the list element to be extracted
#' @method [[ PiDataList
#' @export
`[[.PiDataList` <- function(x,
                            i){
  if(length(i) > 1){stop("Index i must be a length-one numeric or character")}
  if(is.character(i) && !(i %in% names(unclass(x)))){
    stop("Invalid character index.")
  }
  NextMethod()
}

#' @rdname PiData
#' @param  i Index for the list element to be extracted
#' @method [[ PiData
#' @export
`[[.PiData` <- function(x, i){
  if(length(i) > 1){stop("Index i must be a length-one numeric or character")}
  if(is.character(i)){
    if(!(i %in% c("data", "type", "ident", "method", "param" ))){
      stop("Invalid Pi object field name. Choose one of the following: data, type, ident, method, param.")
    }
  }
  # if(is.numeric(i)){
  #   if(i < 1 | i > 5){stop("Invlaid index i. The index i should be integer between 1 and 5.")}
  # }
  NextMethod()

}

#' @rdname PiExpFreqData
#' @param  i Index for the list element to be extracted
#' @method [[ PiExpFreqData
#' @export
`[[.PiExpFreqData` <- function(x, i){
  NextMethod()
}

#' @rdname PiMarkerData
#' @param  i Index for the list element to be extracted
#' @method [[ PiMarkerData
#' @export
`[[.PiMarkerData` <- function(x, i){
  NextMethod()
}

#' @rdname PiDSData
#' @param  i Index for the list element to be extracted
#' @method [[ PiDSData
#' @export
`[[.PiDSData` <- function(x, i){
  NextMethod()
}


#' @rdname PiCellPropData
#' @param  i Index for the list element to be extracted
#' @method [[ PiCellPropData
#' @export
`[[.PiCellPropData` <- function(x, i){
  NextMethod()
}

#' @rdname PiParentMetaData
#' @param  i Index for the list element to be extracted
#' @method [[ PiParentMetaData
#' @export
`[[.PiParentMetaData` <- function(x, i){
  NextMethod()
}



#' @rdname Pi
#' @param x A Pi object
#' @param i index for list element to extract or replace
#' @param value An object of class \code{PiData}, \code{PiDataList}, or \code{Seurat} to replace its current value
#' @method [[<- Pi
#' @export
`[[<-.Pi` <- function(x, i, value){
  if(length(i) > 1){stop("Index i must be a length-one numeric or character")}
  valid.types <- c("seurat.obj", "exp.freq", "markers", "ds", "cell.prop", "parent.meta.data")
  valid.classes <- c("Seurat", "PiExpFreqData", "PiMarkerData", "PiDSData", "PiCellPropData", "PiParentMetaData")
  names(valid.classes) <- c("seurat.obj", "exp.freq", "markers", "ds", "cell.prop", "parent.meta.data")
  if(is.character(i)){
    if(!(i %in% valid.types)){
      stop("Invalid Pi object field name. Choose one of the following: seurat.obj, exp.freq, markers, ds, cell.prop, parent.meta.data.")
    }
  }
  i <- ifelse(is.character(i), i, valid.types[i])

  cl <- class(x)
  x <- unclass(x)
  if(i == "seurat.obj"){
    if(inherits(value, "Seurat")){
      x[[i]] <- value
    }else{
      stop("Overwriting a seurat.obj field with an object of class ", class(value)[1])
    }
  }else{
    if(is.PiDataList(value, type = valid.classes[i])){## if input is a PiDataList, it will overwirte the current PiDataList
      x[[i]] <- value
    }else if(is.PiData(value)){
      if(!inherits(value, valid.classes[i])){
        stop("Incosistent data provided: trying to modify field ", i, " but the provdied object is of class ", class(value)[1])
      }else{
        x[[i]] <- PiDataList(value)
      }
    }else{
      stop("Invalid input. Argument \"value\" must be one of the following classes: Seurat, PiDataList or PiData ")
    }
  }


  class(x) <- cl
  return(x)

}


`[[<-.PiDataList` <- function(x, i, value){
  stop(MsgStopListReplacement())

}


`[[<-.PiData` <- function(x, i, value){
  stop(MsgStopListReplacement())
}


`[<-.Pi` <- function(x, i, value){
  stop(MsgStopListReplacement())
}


`[<-.PiDataList` <- function(x, i, value){
  stop(MsgStopListReplacement())
}


`[<-.PiData` <- function(x, i, value){
  stop(MsgStopListReplacement())
}


MsgStopListReplacement <- function(){
  msg1 <- "Direct replacement/modification of Pi/PiDataList/PiData list elements is prohibited. \n"
  msg2 <- "  Run one of the following functions on the Pi object to add new analysis results:\n"
  msg3 <- "\tCalculateExpFreqs: calculate per gene expression frequency\n"
  msg4 <- "\tRunFindAllMarkers: find markers for each cluster\n"
  msg5 <- "\tRunPseudobulkAnalysis: differential state analysis\n"
  msg6 <- "\tRunProportionPlot: calculate cell proportions and test differences\n"
  return(paste(msg1, msg2, msg3, msg4, msg5, msg6, sep = ""))
}



## Add PiData
#' Add PiData
#'
#' Add PiData object to an existing PiData/PiDataList/Pi object. This function is intended to combine multiple analysis results of the same type.
#' @param x A PiData/PiDataList/Pi or list object
#' @param ... Additional arguments
#'
#' @return A PiDataList or Pi object.
#' @export
#' @examples
#' \dontrun{
#' ## Create two PiData objects for markers analysis
#' res1 <- RunFindAllMarkers(csle.bcell.small,
#'                          ident = "seurat_clusters",
#'                          test.use = "wilcox",
#'                          latent.vars = NULL)
#' mk1 <- PiMarkerData(data = res1,
#'                    ident = "seurat_clusters",
#'                    method = "AllMarkers",
#'                    test.use = "wilcox",
#'                    latent.vars = NULL)
#' res2 <- RunFindAllMarkers(csle.bcell.small,
#'                          ident = "main.cluster.anno",
#'                          test.use = "wilcox",
#'                          latent.vars = NULL)
#' mk2 <- PiMarkerData(data = res2,
#'                    ident = "main.cluster.anno",
#'                    method = "AllMarkers",
#'                    test.use = "wilcox",
#'                    latent.vars = NULL)
#'
#' ## Add a PiData object to another PiData object of the same type returns a PiDataList of two elements
#' mk.list <- AddPiData(mk1, mk2)
#'
#' ## Add a PiData object to an empty list returns a PiDataList
#' mk.list <- AddPiData(list(), mk1)
#'
#' ## Add a PiData object to a PiDataList. If the same PiData object already exists in the list, it will be overwritten
#' mk.list1 <- AddPiData(mk.list, mk1)
#'
#' ## Add a PiData object to a Pi object
#' my.pi <- CreatePostIntegrationObject(csle.bcell.small) ## a minimum Pi object
#' my.pi <- AddPiData(my.pi, mk1)
#' }
#'
AddPiData <- function(x, ...){
  UseMethod("AddPiData", x)
}

#' @rdname AddPiData
#' @param value A \code{PiData} object
#' @method AddPiData default
#' @export
#'
AddPiData.default <- function(x, value){
  if(!inherits(value, "PiData")){
    stop("Argument value must be an object of class PiData")
  }
  x <- list(value)
  names(x) <- GetPiDataName(value)
  class(x) <- "PiDataList"
  return(x)
}

#' @rdname PiData
#' @param value A \code{PiData} object
#' @method AddPiData PiData
#' @export
AddPiData.PiData <- function(x, value){
  if(!inherits(value, "PiData")){
    stop("Argument value must be an object of class PiData")
  }
  return(PiDataList(x, value))
}


#' @rdname PiDataList
#' @param value A \code{PiData} object
#' @method AddPiData PiDataList
#' @export
AddPiData.PiDataList <- function(x, value){
  if(!inherits(value, "PiData")){
    stop("Argument value must be an object of class PiData")
  }

  type <- value[["type"]]
  if(!all(unlist(lapply(x, function(u) u[["type"]] == type)))){stop("Type of existing elements and type of new PiData must be the same")}

  cl <- class(x)
  x <- unclass(x)
  x.names <- names(x)
  if(GetPiDataName(value) %in% x.names){
    message("PiData ",GetPiDataName(value), " already exisits. Overwriting...")
    idx <- which(x.names == GetPiDataName(value))
    x[[idx]] <- value
  }else{
    x[[length(x) + 1]] <- value
    names(x) <- c(x.names, GetPiDataName(value))
  }

  class(x) <- cl
  return(x)
}

#' @rdname Pi
#' @method AddPiData Pi
#' @export
AddPiData.Pi <- function(x, value){
  if(!is.PiData(value)){
    stop("Argument value must be an object of class PiData")
  }
  # valid.types <- c("exp.freq", "markers", "ds", "cell.prop", "parent.meta.data")
  # names(valid.types) <- c("PiExpFreqData", "PiMarkerData", "PiDSData", "PiCellPropData", "PiParentMetaData")

  cl <- class(x)
  x <- unclass(x)
  obj <- x[[value[["type"]]]]
  x[[value[["type"]]]] <- AddPiData(obj, value)

  class(x) <- cl
  return(x)
}

## other utils for Pi/PiDataList/PiData
ParamListToStr <- function(x){
  if(length(x) == 0){
    warning("Unexpected empty param list")
    return("")
  }else{
      my.param.str <- ""
      for(i in 1:length(x)){
        if(i == 1){
          my.param.str <- paste(my.param.str, names(x)[i], "=", x[i], sep = "" )
        }else{
          my.param.str <- paste(my.param.str, ";", names(x)[i], "=", x[i], sep = "" )
        }
      }
      return(my.param.str)
    }
}

GetPiDataName <- function(x){
  stopifnot(inherits(x, "PiData"))
  valid.class.short.names <- c("ExpFreq", "Markers", "DS", "CellProp", "ParentMetaData")
  names(valid.class.short.names) <- c("PiExpFreqData", "PiMarkerData", "PiDSData", "PiCellPropData", "PiParentMetaData")
  value <- valid.class.short.names[class(x)[1]]

  value <- ifelse(is.null(x$ident), value, paste(value, x$ident, sep = "|"))
  value <- ifelse(is.null(x$method), value, paste(value, x$method, sep = "|"))
  value <- ifelse(is.null(x$param), value, paste(value, ParamListToStr(x$param), sep = "|"))
  return(value)
}
