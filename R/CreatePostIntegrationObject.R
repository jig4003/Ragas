#' Create a Post Integration Object
#'
#' Initiates a \code{\link[Ragas]{Pi}} object that collects information after single-cell integration at either main level or sub-cluster level.
#'
#' This function loads Seurat object(s) after basic analytic steps (i.e., normalization, dimension reduction, single-cell integration)
#' and prepares for downstream analysis. When sub-cluster analysis has been performed, the function provides an interface to link
#' parent and child Seurat objects so seamless analysis can be orchestrated.
#'
#' @param object A Seurat object
#' @param child.object.list A list of Seurat/Pi objects that are the children of \code{object} from sub-cluster analysis
#' @param keep.child.object.name Whether to keep the child object name as the prefix for subcluster identities. See Example 2 below
#' @param parent.object A Seurat object that is the parent of \code{object} or a data frame of the metadata from the parent Seurat object (i.e., seurat_object[[]]). Use data frame as input can significantly reduce use of memory.
#' @param parent.key An optional character to describe the parent object. If not provided, the name of the parent object will be used. If parent.object is a data frame, assign parent.key becomes mandatory.
#' @param rp.main.cluster.anno Additional cluster annotation from the metadata column, such as user-annotated cluster identities.
#' Only applies to re-projection of subclusters
#' @param rp.main.cluster.umap.config A list containing UMAP run parameters for main clusters. See \code{\link[Ragas]{ConfigureReprojection}} for more details. When this argument is set to NULL,
#' defaults of \code{\link[Ragas]{ConfigureReprojection}} will be used when re-running UMAP
#' @param rp.main.cluster.to.preserve Name of the metadata column containing main cluster identities that needs to be preserved during re-projection (default: seurat_clusters)
#' @param rp.subcluster.colname Name of the metadata column of the child object contains subcluster identities (default: seurat_clusters). Two formats are accepted:
#' \itemize{
#' \item{a single character} : name of the (common) metadata column containing cell identity/cluster information for one or more child objects
#' \item{a character vector} : must be a named vector with its elements correspond to the child objects. Its names must match those of the \code{child.object.list}
#' }
#' @param rp.subcluster.umap.config A list containing UMAP run parameters for subclusters per child object. Same as the \code{rp.main.cluster.umap.config}, when set to NULL,
#' defaults of \code{\link[Ragas]{ConfigureReprojection}} will be used when re-running UMAP. Names of the list must match those of the \code{child.object.list}
#' @param rp.weight A numeric number between 0 and 1 indicating the weight to reserve nearest neighbor structures from the parent object.
#' Only applies to re-projection of subclusters (default: 1)
#' @param rp.reduction.name Name to store the re-projected dimension reduction object (default: rp)
#' @param rp.reduction.key Key for the re-projected dimension reduction object (default: RPUMAP_)
#' @param rp.seed Random seed for re-projected UMAP
#' @param verbose Verbosity (default: TRUE)
#'
#' @references Hao and Hao et al. (2021). “Integrated analysis of multimodal single-cell data.” Cell. doi:10.1016/j.cell.2021.04.048
#' @return A \code{\link[Ragas]{Pi}} object includes the following components:
#' \itemize{
#'  \item{"seurat.obj}: a Seurat object
#'  \item{"exp.freq}: cell expression frequency
#'  \item{"markers}: results of differential expression analysis to identify markers for each cell identity
#'  \item{"ds}: results of differential state analysis by pseudobulk method
#'  \item{"cell.prop}: results of differential analysis on cell proportions
#'  \item{"parent.meta.data}: metadata from the parent object, if applicable
#' }
#' @examples
#' \dontrun{
#' ## Example 1: create a Pi object without integrating subclusters
#' my.pi <- CreatePostIntegrationObject(object = csle.pbmc.small)
#'
#' ## Example 2: add parent object
#' my.pi <- CreatePostIntegrationObject(object = csle.bcell.small, parent.object = csle.pbmc.small)
#'
#' ## a more efficient alternative
#' metadata <- csle.pbmc.small[[]] ## get metadata from Seurat object
#' my.pi <- CreatePostIntegrationObject(object = csle.bcell.small, parent.object = metadata, parent.key = 'pbmc') ## a key name is required
#'
#' ## Example 3: create a Pi object and integrate subclusters with default UMAP run parameters
#' subclusters <- list(T = "csle.tcell.small")
#' my.pi <- CreatePostIntegrationObject(object = csle.pbmc.small,
#'                                      keep.child.object.name = FALSE,
#'                                      child.object.list = subclusters,
#'                                      rp.main.cluster.anno = "cluster.annotation",
#'                                      rp.subcluster.colname = "cluster.annotation")
#' ### A new UMAP called "rp" will be created; by default, cell cluster information is stored in a Seurat metadata column called "subcluster_idents"
#' RunDimPlot(my.pi, reduction = "rp", group.by = "subcluster_idents",label = TRUE, raster = FALSE)
#'
#' ## Example 4: create a Pi object with user-defined parameters for re-projection (multi-level reprojection)
#' ### level 1: integrate T cell and its subclusters
#' subclusters <- list("CD4mem" = "csle.cd4.mem.small", "Treg" = "csle.treg.small")
#' tcell.pi <- CreatePostIntegrationObject(object = csle.tcell.small, rp.weight = 0,
#'                                         child.object.list = subclusters,
#'                                         rp.main.cluster.anno = "cluster.annotation",
#'                                         rp.subcluster.colname = "cluster.annotation")
#'
#' ### level 2: integrate PBMC with T and B subclusters
#' subclusters <- list("B cell" = "csle.bcell.small", "T cell" = "tcell.pi")
#' sc.config <- ConfigureReprojection(type = "sc", "B cell") ## default
#' sc.config <- ConfigureReprojection(type = "sc", "T cell", umap.name = "rp", append.to = sc.config) ## use UMAP named "rp"
#' sc.colnames <- c("B cell" = "cluster.annotation", "T cell" = "subcluster_idents")
#' pbmc.pi <- CreatePostIntegrationObject(object = csle.pbmc.small,
#'                                        child.object.list = subclusters,
#'                                        rp.subcluster.umap.config = sc.config,
#'                                        rp.main.cluster.anno = "cluster.annotation",
#'                                        rp.subcluster.colname = sc.colnames)
#' RunDimPlot(pbmc.pi, reduction = "rp", group.by = "subcluster_idents",label = TRUE, raster = FALSE)
#' }
#'
#' @export
#'
CreatePostIntegrationObject <- function(object,
                                        child.object.list = NULL,
                                        keep.child.object.name = TRUE,
                                        parent.object = NULL,
                                        parent.key = NULL,
                                        rp.main.cluster.anno = NULL,
                                        rp.main.cluster.umap.config = NULL,
                                        rp.main.cluster.to.preserve = "seurat_clusters",
                                        rp.subcluster.colname = NULL,
                                        rp.subcluster.umap.config = NULL,
                                        rp.weight = 1,
                                        rp.reduction.name = "rp",
                                        rp.reduction.key = "RPUMAP_",
                                        rp.seed = 42L,
                                        verbose = TRUE
                                        ){
  if(!is(object, "Seurat")){stop("Invalid input argument: object. Seurat object is required.")}

  if(!is.null(rp.subcluster.umap.config)){
    if(length(unique(unlist(lapply(rp.subcluster.umap.config, function(x) x$umap.n.neighbors)))) > 1){stop("UMAP analysis for all subclusters should have the same umap.n.neighbors!")}
    if(!is.null(rp.main.cluster.umap.config)){
      if(rp.main.cluster.umap.config[[1]]$umap.n.neighbors != unique(unlist(lapply(rp.subcluster.umap.config, function(x) x$umap.n.neighbors))) ){
        stop("Main cluster and subcluster analysis should have the same umap.n.neighbors!\t")
      }
    }
  }
  if((!is.null(child.object.list)) && (!is.null(rp.subcluster.umap.config))){
    if(!identical(sort(names(child.object.list)), sort(names(rp.subcluster.umap.config)))){
      stop("Names of child.object.list do not match rp.subcluster.umap.config!")
    }
  }
  if(!is.null(child.object.list)){
    if(is.null(rp.subcluster.colname)){
      rp.subcluster.colname <- rep("seurat_clusters", length(child.object.list))
      names(rp.subcluster.colname) <- names(child.object.list)
    }else{
      if(!is.character(rp.subcluster.colname)){
        stop("rp.subcluster.colname is not character!")
      }else{
        if(length(rp.subcluster.colname) == 1){
          rp.subcluster.colname <- rep(rp.subcluster.colname, length(child.object.list))
          names(rp.subcluster.colname) <- names(child.object.list)
        }else{
          if(is.null(names(rp.subcluster.colname))){
            stop("rp.subcluster.colname must be a named character vector when its length  is greater than 1!")
          }
          if(!identical(sort(names(child.object.list)), sort(names(rp.subcluster.colname)))){
            stop("Names of child.object.list do not match rp.subcluster.colname!")
          }
        }
      }
    }
  }


  if(!is.null(child.object.list)){
    seurat.obj <- RunSubclusterReprojection(object,
                                            subclusters.list = child.object.list,
                                            keep.child.object.name = keep.child.object.name,
                                            weight = rp.weight,
                                            main.cluster.anno = rp.main.cluster.anno,
                                            main.cluster.umap.config = rp.main.cluster.umap.config,
                                            main.cluster.to.preserve = rp.main.cluster.to.preserve,
                                            subcluster.colname = rp.subcluster.colname,
                                            subcluster.umap.config = rp.subcluster.umap.config,
                                            rp.reduction.name = rp.reduction.name,
                                            rp.reduction.key = rp.reduction.key,
                                            rp.seed = rp.seed,
                                            verbose = verbose)
  }else{
    seurat.obj <- object
  }


  exp.freq <- NULL
  markers <- NULL
  ds <- NULL
  cell.prop <-NULL
  if(!is.null(parent.object)){
    if(!inherits(parent.object, "Seurat") & !(is.data.frame(parent.object))){stop("Parent object must be a Seurat object or a data frame containing metadata from Seurat object!")}
    if(inherits(parent.object, "Seurat")){
      parent.object.name = deparse(substitute(parent.object))
      parent.key <- ifelse(is.null(parent.key), parent.object.name, parent.key)
      parent.meta.data <- PiParentMetaData(data = parent.object[[]],
                                           parent.name = parent.key)
      CheckPiData(parent.meta.data, seurat.obj = object)
    }
    if(is.data.frame(parent.object)){
      if(is.null(parent.key)){stop('Argument \'parent.key\' cannot be NULL while \'parent.object\' is a data frame!')}
      parent.meta.data <- PiParentMetaData(data = parent.object,
                                           parent.name = parent.key)
      CheckPiData(parent.meta.data, seurat.obj = object)
    }

  }else{
    parent.meta.data <- NULL
  }

  pi.obj <- Pi(object = seurat.obj)

  if(!is.null(parent.meta.data)){
    pi.obj <- AddPiData(pi.obj, parent.meta.data)
  }


  if(verbose == TRUE){message('Post-integration object created')}
  return(pi.obj)
}
