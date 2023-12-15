#' Run DSDotPlot
#'
#' A wrapper function for \code{\link[Ragas]{DSDotPlot}}
#'
#' The input Pi object must have analyzed results in the \code{ds} and \code{exp.freq} fields, which can be accomplished by running \code{\link[Ragas]{RunPseudobulkAnalysis}}
#' and \code{\link[Ragas]{CalculateExpFreqs}}, respectively.
#'
#' Gene filtering will be performed based on significance (i.e., \code{p.filter} and \code{FC.filter}) and expression frequency (\code{exp.freq.filter}).
#' The purpose of significance filtering is twofold: (1) only features that passed the filtering criteria in at least one of the identity groups (e.g., cell clusters,
#' in most cases) will be included in the plot; (2) significant features that passed the filtering criteria will be highlighted on their corresponding clusters.
#' The \code{FC.filter} argument expects un-logged fold-change and a \code{FC.filter} of 2 can mean a two-fold expression change for both up or down-regulation.
#' If \code{top.n} is set to a number other than Inf, only the top genes will be plotted.
#'
#' If the argument \code{features} is set, all valid features (i.e., exist in the Seurat object) will be plotted in the figure regardless of their significance. However, only features
#' passed the filtering criteria will be highlighted.
#'
#' @param object A \code{\link[Ragas]{Pi}} object
#' @param features A list of features/genes to plot
#' @param exp.freq.key A unique identifier to retrieve data from the \code{exp.freq} object. See  \code{\link[Ragas]{CalculateExpFreqs}} for more details
#' @param ds.key A unique identifier to retrieve data from the \code{ds} object. See  \code{\link[Ragas]{RunPseudobulkAnalysis}} for more details
#' @param p.filter P-Value cutoff for gene filtering (default: 0.05)
#' @param FC.filter Fold change cutoff for gene filtering (default: 2)
#' @param exp.freq.filter Gene expression frequency filtering cutoff (default: 0.1)
#' @param filter.by.median.diff Whether to perform additional filtering based on median CPM differences between groups. This will remove differential genes due to outliers.
#' @param median.diff.filter Median difference filtering cutoff (default: 0)
#' @param to.adjust Whether to perform adjustment for multiple testing (default: FALSE)
#' @param top.n Number of top genes to plot (default: Inf, i.e., all genes passing filtering criteria will be plotted)
#' @param direction Whether to plot "up" or "down" regulated genes, or "both" (default: "both")
#' @param to.highlight Whether to highlight feature/identity pairs that passed the filtering criteria (default: TRUE)
#' @param clust.row Whether to cluster rows (default: TRUE)
#' @param clust.column Whether to cluster columns (default: FALSE)
#' @param verbose Verbosity (default: TRUE)
#' @param ... Extra parameters passed to \code{\link[Ragas]{DSDotPlot}} to change plot font and color, etc.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object
#'
#' @references Hao and Hao et al. (2021). “Integrated analysis of multimodal single-cell data.” Cell. doi:10.1016/j.cell.2021.04.048
#'
#' @examples
#' \dontrun{
#' RunDSDotPlot(object = csle.bcell.pi,
#'              exp.freq.key = "ExpFreq|cluster.annotation|cutoff=0",
#'              ds.key = "DS|cluster.annotation|edgeR|group=Groups;sample=Names;gp1=cSLE;gp2=cHD;contrast=cSLE-cHD",
#'              p.filter = 0.05)
#'
#' }
#' @export
#'
RunDSDotPlot <- function(object,
                         features = NULL,
                         exp.freq.key,
                         ds.key,
                         p.filter = 0.05,
                         FC.filter = 2,
                         exp.freq.filter = 0.1,
                         filter.by.median.diff = TRUE,
                         median.diff.filter = 0,
                         to.adjust = FALSE,
                         top.n = Inf,
                         direction = "both",
                         to.highlight = TRUE,
                         clust.row = TRUE,
                         clust.column = FALSE,
                         verbose = TRUE,
                         ...
                         ){
  if(!is(object, "Pi")){
    if(is(object, "Seurat")){
      message("RunDEDotPlot requires Pi object as input, but Seurat object is provided. Run \"CreatePostIntegrationObject\" first.")
    }else{
      stop("Invalid input argument: object. Post-integration object is required.")
      }
  }

  object <- unclass(object)

  ## Check and load exp.freq
  if(length(object$exp.freq) == 0){
    stop("The \"exp.freq\" field of the Pi object is empty. Run CalculateExpFreqs on the Pi object first.")
  }else{
    if(!exp.freq.key %in% names(object$exp.freq)){
      stop("Invalid key \"",exp.freq.key, "\" for the exp.freq object. The following key(s) are available:\n    ", paste(names(object$exp.freq), collapse = "\n    "),"\nEither choose an availble key or rerun CalculateExpFreqs.")
    }else{
      exp.freq <- object$exp.freq[[exp.freq.key]][["data"]]
    }
  }

  ## Check and load pseudobulk analysis results
  if(length(object$ds) == 0){
    stop("The \"ds\" field of the Pi object is empty. Run RunPseudobulkAnalysis on the Pi object first.")
  }else{
    if(!ds.key %in% names(object$ds)){
      stop("Invalid key \"",ds.key ,"\" for the ds object. The following key(s) are available:\n    ", paste(names(object$ds), collapse = "\n    "),"\nEither choose an availble key or rerun RunPseudobulkAnalysis")
    }else{
      res <- object$ds[[ds.key]][["data"]]$res$table[[1]]
    }
  }
  if(!(length(names(res)) == length(names(exp.freq)))){
    cat("Idents from the exp.freq object: ",names(exp.freq),"\n")
    cat("Idents from the ds object: ",names(res),"\n")
    stop("Identity names of the exp.freq and ds objects do not match")
  }else{
    if(!all(names(res) == names(exp.freq))){
      cat("Idents from the exp.freq object: ",names(exp.freq),"\n")
      cat("Idents from the ds object: ",names(res),"\n")
      stop("Identity names of the exp.freq and ds objects do not match")
    }
  }

  ## Check features
  if(!is.null(features)){
    if(length(features) != length(unique(features))){
      features <- unique(features)
      if(verbose){message("The features provided are not unique. Taking unique.")}
    }
    features.valid <- intersect(features, rownames(object$seurat.obj))
    features.invalid <- setdiff(features, features.valid)
    if(length(features.invalid) > 0){
      message("The following feature(s): ", paste(features.invalid, collapse = ", "), " not found in the Pi object")
    }
    if(length(features.valid) == 0){
      stop("None of the provided features are found in the Pi object")
    }
  }

  ## Check direction
  if(!direction %in% c("up", "down", "both")){stop("The \"direction argument\" must be one of the following: up, down, or both!")}

  ## Additional filtering when applicable
  if(filter.by.median.diff){
    cpm.dat <- object$ds[[ds.key]][["data"]]$tbl
    group.col.name <- object$ds[[ds.key]]$param$group
    sample.col.name <- object$ds[[ds.key]]$param$sample
    case.group <- object$ds[[ds.key]]$param$gp1
    ref.group <- object$ds[[ds.key]]$param$gp2
  }else{
    cpm.dat <- NULL
  }

  p <- DSDotPlot(object = object$seurat.obj,
                 res = res,
                 exp.freq = exp.freq,
                 features = features,
                 p.threshold = p.filter,
                 FC.threshold = FC.filter,
                 pct.threshold = exp.freq.filter,
                 n.deg = top.n,
                 p.adj.filter = to.adjust,
                 highlight.deg = to.highlight,
                 direction = direction,
                 clust.row = clust.row,
                 clust.column = clust.column,
                 cpm.dat = cpm.dat,
                 exp.filter = filter.by.median.diff,
                 exp.filter.thr = median.diff.filter,
                 group.col.name = group.col.name,
                 sample.col.name = sample.col.name,
                 case.group = case.group,
                 ref.group = ref.group,
                 ...)
  return(suppressWarnings(print(p)))
  # return(p)
}
