#' Run Pseudobulk Analysis
#'
#' Run differential state analysis with the \code{muscat} package on
#' pseudobulk counts for multi-sample, multi-group study design
#'
#' @param object A Seurat or \code{\link[Ragas]{Pi}} object
#' @param ident.var Cell identity variable
#' @param group.var Group variable, used to test differential state
#' @param sample.var Sample variable
#' @param blocking.var Additional blocking variable used for batch correction (default: NULL)
#' @param group.1 A character string for the case group, must be consistent with the levels of \code{group.var}
#' @param group.2 A character string for the ctrl group, must be consistent with the levels of \code{group.var}
#' @param pbDS.method Name of the method for pseudobulk analysis (default: edgeR)
#' @param pbDS.verbose Verbosity of the \code{\link[muscat]{pbDS}} function
#' @param pbDS.min.cell A numeric number that specifies the minimum number of cells required in a given cluster-sample.
#' See \code{\link[muscat]{pbDS}} for details (default 1)
#' @param pbDS.filter Whether to filter on genes and/or samples. See \code{\link[muscat]{pbDS}} for details (default: "none")
#' @param resDS.bind,resDS.frq,resDS.cpm,resDS.digits,resDS.sep  parameters for \code{\link[muscat]{resDS}}
#' @param calcExprFreqs.assay,calcExprFreqs.th parameters for \code{\link[muscat]{calcExprFreqs}}
#' @param auto.make.names Whether to automatically check and make syntactically valid names (default: TRUE)
#' @param verbose Verbosity (default: TRUE)
#' @param ... Extra parameters passed to \code{\link[muscat]{pbDS}}
#'
#' @references Hao and Hao et al. (2021). “Integrated analysis of multimodal single-cell data.” Cell. doi:10.1016/j.cell.2021.04.048
#' @references Helena L. Crowell, Pierre-Luc Germain, Charlotte Soneson, Anthony Sonrel and Mark D. Robinson (2020). muscat: Multi-sample
#' multi-group scRNA-seq data analysis tools. R package version 1.4.0. https://github.com/HelenaLC/muscat
#' @references Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth, G.K. (2015). limma powers differential expression analyses
#' for RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47.
#'
#' @return When \code{object} is a Seurat object, a list containing the output of \code{\link[muscat]{pbDS}} and  \code{\link[muscat]{resDS}} will be returned;
#' If \code{object} is a Pi object (recommended), an updated Pi object will be returned:  a \code{\link[Ragas]{PiDSData}} object that encloses pseudobulk analysis results
#' will be created and added to the \code{ds} field.
#'
#' @importFrom muscat prepSCE aggregateData pbDS resDS
#' @importFrom SingleCellExperiment colData
#' @importFrom limma makeContrasts
#' @examples
#' \dontrun{
#' my.pi <- CreatePostIntegrationObject(object = csle.pbmc.small) ## a minimum Pi object
#' my.pi <- RunPseudobulkAnalysis(object = my.pi,
#'                                ident.var = "cluster.annotation",
#'                                group.var = "Groups",
#'                                sample.var = "Names",
#'                                group.1 = "cSLE",
#'                                group.2 = "cHD")
#' }
#' @export
#require(muscat)
#require(limma)
RunPseudobulkAnalysis <- function(object,
                                  ident.var,
                                  group.var,
                                  sample.var,
                                  blocking.var = NULL,
                                  group.1,
                                  group.2,
                                  pbDS.method = c("edgeR", "DESeq2", "limma-trend", "limma-voom"),
                                  pbDS.verbose = TRUE,
                                  pbDS.min.cell = 1,
                                  pbDS.filter = "none",
                                  resDS.bind = "row",
                                  resDS.frq = TRUE,
                                  resDS.cpm = TRUE,
                                  resDS.digits = Inf,
                                  resDS.sep = "__",
                                  calcExprFreqs.assay = "counts",
                                  calcExprFreqs.th = 0,
                                  auto.make.names = TRUE,
                                  verbose = TRUE,
                                  ...
){
  if(is(object, "Pi") || is(object, "Seurat")){
    if(is(object, "Pi")){seurat.obj <- object$seurat.obj}
    if(is(object, "Seurat")){seurat.obj <- object}

    pbDS.method <- match.arg(pbDS.method)
    ## check input args
    if(!(ident.var %in% names(seurat.obj[[]]))){stop("ident.var not found in the meta data")}
    if(!(group.var %in% names(seurat.obj[[]]))){stop("group.var not found in the meta data")}
    if(!(sample.var %in% names(seurat.obj[[]]))){stop("sample.var not found in the meta data")}
    if(!is.null(blocking.var)){
      if(!(blocking.var %in% names(seurat.obj[[]]))){stop("blocking.var not found in the meta data")}

      blocking.var.dat <- seurat.obj[[blocking.var, drop = TRUE]]
      # message("Blocking var: ",paste(unique(blocking.var.dat), collapse = ","))
      if(length(unique(blocking.var.dat)) < 2){
        stop("blocking.var should have at least 2 levels, otherwsie leave it NULL")
      }
    }

    if(!(group.1 %in% unique(seurat.obj[[group.var, drop = TRUE]]))){stop("Invalid group.1")}
    if(!(group.2 %in% unique(seurat.obj[[group.var, drop = TRUE]]))){stop("Invalid group.2")}

    ## check whether levels of group_id and block_id are valid names
    if(verbose){message("Checking validity of names...")}
    if(is.null(blocking.var)){
      seurat.obj <- .PrepareSeuratObjectPseudobulkAnalysis(object = seurat.obj,
                                             var.list = c(ident.var, group.var, sample.var),
                                             check.names = c(FALSE, TRUE, FALSE),
                                             auto.make.names = auto.make.names,
                                             verbose = verbose)
    }else{
      seurat.obj <- .PrepareSeuratObjectPseudobulkAnalysis(object = seurat.obj,
                                             var.list = c(ident.var, group.var, sample.var, blocking.var),
                                             check.names = c(FALSE, TRUE, FALSE, TRUE),
                                             auto.make.names = auto.make.names,
                                             verbose = verbose)
    }

    if(auto.make.names){
      if(!(group.1 %in% unique(seurat.obj[[group.var, drop = TRUE]]) & group.2 %in% unique(seurat.obj[[group.var, drop = TRUE]]))){
        group.1 <- make.names(group.1)
        group.2 <- make.names(group.2)
        message("Changed arguments group.1 and group.2 to ", group.1, " and ", group.2, " accordingly.")
      }

    }

    my.contrast <- paste("group_id", group.1, "-", "group_id", group.2, sep = "")
    if(verbose){
      message("Contrast created: ", my.contrast)
    }
    # if(is.null(blocking.var)){
    #   run.identifier <- paste(group.1, "-",  group.2, sep = "")
    # }else{
    #   run.identifier <- paste(group.1, "-",  group.2, ";", blocking.var, sep = "")
    # }


    ## create single cell experiment assay
    if(verbose){message("Creating single cell experiment object...")}
    sce <- as.SingleCellExperiment(seurat.obj)
    sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
    sce <- prepSCE(sce,
                   kid = ident.var, # subpopulation assignments
                   gid = group.var,
                   sid = sample.var,
                   drop = FALSE)
    ## assign block_id if needed
    if(!is.null(blocking.var)){
      sce$block_id <- sce[[blocking.var]]
      }

    ## aggregate counts across cells by cluster
    if(verbose){message("Aggregating count data...")}
    pb <- aggregateData(sce,
                        assay = "counts", fun = "sum",
                        by = c("cluster_id", "sample_id"))
    cd <- as.data.frame(colData(pb))
    if(!is.null(blocking.var)){
      design <- model.matrix(~ 0 + group_id + block_id, data = cd)
    }else{
      design <- model.matrix(~ 0 + group_id, data = cd)
    }


    contrast <- makeContrasts(contrasts = my.contrast, levels = design)
    ## print(contrast)
    .CheckDesign(sce,design)

    if(verbose){message("Running differential state analysis...")}
    res <- pbDS(pb,
                method = pbDS.method,
                design = design,
                contrast = contrast,
                verbose = pbDS.verbose,
                min_cells = pbDS.min.cell,
                filter = pbDS.filter,
                ...)

    ## formatting output
    if(verbose){message("\nFormatting output...")}
    tbl <- resDS(sce,
                 res,
                 bind = resDS.bind,
                 frq = resDS.frq,
                 cpm = resDS.cpm,
                 digits = resDS.digits,
                 sep = resDS.sep,
                 assay = calcExprFreqs.assay,
                 th = calcExprFreqs.th
                 )

    if(verbose){message("Pseudobulk-analysis completed!")}
    if(is(object, "Pi")){
      ds <- PiDSData(data = list(res = res,tbl = tbl),
                     ident = ident.var,
                     method = pbDS.method,
                     group = group.var,
                     sample = sample.var,
                     group.1 = group.1,
                     group.2 = group.2,
                     contrast = paste(group.1, "-",  group.2, sep = ""),
                     block = blocking.var)

      CheckPiData(ds, seurat.obj = object$seurat.obj)

      object <- AddPiData(object, ds)
      return(object)
    }else{
      return(list(res = res,
                  tbl = tbl))
    }

  }else{
    stop("Invalid input object. Either Seurat or PI object is required.")
    }
}
