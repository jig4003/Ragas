#' Calculate Expression Frequencies
#'
#' Calculate single-cell expression frequency for each gene per identity class.
#'
#' The function will return a list object, either directly as output or indirectly as part of the updated Pi object, with each list element
#' being a data.frame to store the expression frequencies for all genes per \code{ident}. If \code{group.by} is NULL, the data.frame will have only one column called "AllCells" to store gene expression frequencies
#' across all the cells for each identity; if \code{group.by} is set, additional columns storing gene expression frequencies in cells groups defined by \code{group.by} will be calculated and stored.
#'
#' @param object A Seurat or \code{\link[Ragas]{Pi}} object
#' @param ident Name of the metadata column as cell identity (default: seurat_clusters)
#' @param group.by Additional factor to group the cells by
#' @param cutoff Log-expression cutoff (default: 0)
#' @param verbose Verbosity (default: TRUE)
#'
#' @return When \code{object} is a Seurat object, a list containing expression frequencies per \code{ident} will be returned;
#' If \code{object} is a Pi object, an updated Pi object will be returned: a \code{\link[Ragas]{PiExpFreqData}} object that encloses the list
#' containing expression frequencies per \code{ident} will be created and added to the \code{exp.freq} field.
#'
#' @importFrom SingleCellExperiment counts
#' @importFrom Matrix rowSums
#' @examples
#' \dontrun{
#' my.pi <- CreatePostIntegrationObject(object = csle.pbmc.small) ## a minimum Pi object
#' my.pi <- CalculateExpFreqs(my.pi)
#' }
#' @export

CalculateExpFreqs <- function(object,
                              ident = "seurat_clusters",
                              group.by = NULL,
                              cutoff = 0,
                              verbose = TRUE
                               ){
  if(is(object, "Pi") || is(object, "Seurat")){
    if(is(object, "Pi")){seurat.obj <- object$seurat.obj}
    if(is(object, "Seurat")){seurat.obj <- object}

    ## check input args
    if(!(ident %in% names(seurat.obj[[]]))){
      stop("The ident argument ", ident," not found in the meta data")
    }else{
        ident.dat <- seurat.obj[[ident, drop = TRUE]]
        if(!is.factor(ident.dat)){ident.dat <- as.factor(ident.dat)}
        ident.dat <- droplevels(x = ident.dat)
      }
    if(!is.null(group.by)){
      if(!(group.by %in% names(seurat.obj[[]]))){
        stop("The group.by argument ", group.by, " not found in the meta data")
      }else{
          group.dat <- seurat.obj[[group.by, drop = TRUE]]
          if(!is.factor(group.dat)){group.dat <- as.factor(group.dat)}
          if(length(levels(group.dat)) < 2){stop("The group.by argument should at least have 2 levels. Otherwise set is to NULL")}
          group.dat <- droplevels(x = group.dat)
        }
    }

    ## retrieve counts and calculate cell frequencies
    sce <- as.SingleCellExperiment(seurat.obj)
    dat <- counts(sce)

    ident.levels <- levels(ident.dat)
    group.levels <- NULL
    if(!is.null(group.by)){
      group.levels <- levels(group.dat)
      n.group <- length(group.levels) + 1
    }else{
        n.group <- 1
    }

    exp.freq <- list()
    for(i in 1:length(ident.levels))
    {
      if(verbose){message("Calculating expression frequencies for ",ident, ": ", ident.levels[i], appendLF = TRUE)}
      exp.freq.i <- matrix(0, dim(dat)[1], n.group)
      colnames(exp.freq.i) <- c('AllCells',group.levels)
      rownames(exp.freq.i) <- rownames(dat)

      dat.sub <- dat[, which(ident.dat == ident.levels[i]), drop = FALSE]
      exp.freq.i[,'AllCells'] <- rowSums(dat.sub > cutoff) / dim(dat.sub)[2]

      if(!is.null(group.by)){
        group.dat.sub <- group.dat[which(ident.dat == ident.levels[i])]
        for(j in 1:length(group.levels)){
          idx <- which(group.dat.sub == group.levels[j])
          if(length(idx) > 0){
            dat.sub1 <- dat.sub[, idx, drop = FALSE]
            exp.freq.i[,group.levels[j]] <- rowSums(dat.sub1 > cutoff) / dim(dat.sub1)[2]
          }else{
            exp.freq.i[,group.levels[j]] <- 0
            }
          }
        }

      exp.freq[[i]] <- exp.freq.i
    }
    names(exp.freq) <- ident.levels

    if(verbose){message("Done!")}
    if(is(object, "Pi")){
      ef <- PiExpFreqData(data = exp.freq,
                          ident = ident,
                          group = group.by,
                          cutoff = cutoff)
      CheckPiData(ef, seurat.obj = object$seurat.obj)
      object <- AddPiData(object, ef)
      return(object)
    }else{
      return(exp.freq)
    }

  }else{
    stop("Invalid input object. Either Seurat or Pi object is required.")
  }
}
