#' Run FindAllMarkers from Seurat
#'
#' This is a wrapper function to run \code{\link[Seurat]{FindAllMarkers}} function on Seurat or Pi object
#'
#' @param object A Seurat or \code{\link[Ragas]{Pi}} object
#' @param ident Name of the metadata column used as cell identity (default: seurat_clusters)
#' @param test.use Which test to use. See details in \code{\link[Seurat]{FindAllMarkers}} (default: wilcox)
#' @param latent.vars Latent variables to test. See details in \code{\link[Seurat]{FindAllMarkers}}
#' @param future.strategy Parallelization strategies. See \code{\link[future]{plan}} for more details (default: "multisession)
#' @param future.workers Number of parallelizations (default: 1)
#' @param future.globals.maxSize.MB Maximum allowed total size (in Megabyte) of global variables. See \code{\link[future]{future.options}} (default: 1000)
#' @param ... Additional parameters passed to \code{\link[Seurat]{FindAllMarkers}}
#'
#' @return If \code{object} is a Seurat object, a matrix containing markers and their associated statistics will be returned
#' (see \code{\link[Seurat]{FindMarkers}} for more details); If \code{object} is a Pi object, a \code{\link[Ragas]{PiMarkerData}} object containing
#' the marker result matrix will be created and attached to the \code{markers} field.
#'
#' @references Hao and Hao et al. (2021). “Integrated analysis of multimodal single-cell data.” Cell. doi:10.1016/j.cell.2021.04.048
#'
#' @importFrom future plan
#' @export
#'
#' @examples
#' \dontrun{
#' my.pi <- CreatePostIntegrationObject(object = csle.pbmc.small) ## a minimum Pi object
#' my.pi <- RunFindAllMarkers(object = my.pi)
#' }
#'
RunFindAllMarkers <- function(object,
                              ident = "seurat_clusters",
                              test.use = "wilcox",
                              latent.vars = NULL,
                              future.strategy = "multisession",
                              future.workers = 1,
                              future.globals.maxSize.MB = 1000,
                              ...){
  if(is(object, "Pi") || is(object, "Seurat")){
    if(is(object, "Pi")){seurat.obj <- object$seurat.obj}
    if(is(object, "Seurat")){seurat.obj <- object}

    run.type <- "AllMarkers"

    if(!ident %in% names(seurat.obj[[]])){
      stop("\"ident\" not found in the meta data")
    }else{
      Idents(seurat.obj) <- ident
      }

    if(future.workers > 1){
      message("Running multiple workers...")
      options(future.globals.maxSize = future.globals.maxSize.MB * 10^6)
      plan(future.strategy, workers = future.workers)
      plan()
    }
    markers <- FindAllMarkers(seurat.obj, test.use = test.use, latent.vars = latent.vars, ...)
    ## need to implement DeFindAllMarkers

    plan("default")

    if(is(object, "Pi")){
      mk <- PiMarkerData(data = markers,
                         ident = ident,
                         method = run.type,
                         test.use = test.use,
                         latent.vars = latent.vars)
      CheckPiData(mk, seurat.obj = object$seurat.obj)

      object <- AddPiData(object, mk)
      return(object)
    }
    if(is(object, "Seurat")){return(markers)}
  }else{
    stop('Invalid input object. Either Seurat or Pi object is required.')
  }
}
