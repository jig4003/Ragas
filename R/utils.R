#' Run UMAP
#'
#' Run UMAP and return nearest neighbor information
#'
#' This function re-uses code from the Seurat interface to run UMAP from uwot package
#' and return nearest neighbor data for downstream analysis
#' @importFrom uwot umap
#' @importFrom future nbrOfWorkers
#' @noRd

.RunUMAP <- function(object = NULL,
                     dims = 1:20,
                     reduction = 'harmony',
                     features = NULL,
                     graph = NULL,
                     assay = DefaultAssay(object = object),
                     nn.name = NULL,
                     slot = 'data',
                     umap.method = 'uwot',
                     reduction.model = NULL,
                     return.model = FALSE,
                     n.neighbors = 30L,
                     n.components = 2L,
                     metric = 'cosine',
                     n.epochs = NULL,
                     learning.rate = 1,
                     min.dist = 0.3,
                     spread = 1,
                     set.op.mix.ratio = 1,
                     local.connectivity = 1L,
                     repulsion.strength = 1,
                     negative.sample.rate = 5L,
                     a = NULL,
                     b = NULL,
                     uwot.sgd = FALSE,
                     seed.use = 42L,
                     metric.kwds = NULL,
                     angular.rp.forest = FALSE,
                     verbose = FALSE,
                     nn_method = NULL,
                     y = NULL,
                     target_weight = 0.5){

  if(!is.null(object)){
    data.use <- Embeddings(object[[reduction]])[, dims]
    assay <- DefaultAssay(object = object[[reduction]])

    set.seed(seed.use)
    ## 3. uwot umap defaults
    my.umap <- umap(
      X = data.use,
      n_threads = nbrOfWorkers(),
      n_neighbors = as.integer(x = n.neighbors),
      n_components = as.integer(x = n.components),
      metric = metric,
      n_epochs = n.epochs,
      learning_rate = learning.rate,
      min_dist = min.dist,
      spread = spread,
      set_op_mix_ratio = set.op.mix.ratio,
      local_connectivity = local.connectivity,
      repulsion_strength = repulsion.strength,
      negative_sample_rate = negative.sample.rate,
      a = a,
      b = b,
      fast_sgd = uwot.sgd,
      verbose = verbose,
      ret_model = return.model,
      ret_nn = TRUE,
      y = y,
      target_weight = target_weight
    )
  }else{
    set.seed(seed.use)
    my.umap <- umap(
      X = NULL,
      n_threads = nbrOfWorkers(),
      n_neighbors = as.integer(x = n.neighbors),
      n_components = as.integer(x = n.components),
      metric = metric,
      n_epochs = n.epochs,
      learning_rate = learning.rate,
      min_dist = min.dist,
      spread = spread,
      set_op_mix_ratio = set.op.mix.ratio,
      local_connectivity = local.connectivity,
      repulsion_strength = repulsion.strength,
      negative_sample_rate = negative.sample.rate,
      a = a,
      b = b,
      fast_sgd = uwot.sgd,
      verbose = verbose,
      ret_model = return.model,
      ret_nn = TRUE,
      nn_method = nn_method,
      y = y,
      target_weight = target_weight
    )
  }


  return(my.umap)
}

#' Run Single cell Re-projection
#'
#' Combine child Seurat objects from sub-cluster analysis and their re-project cell embedding on parent Seurat object
#' @noRd
#'
.RunReprojection <- function(mc,
                             sc.list,
                             mc.umap.config = NULL,
                             sc.umap.config = NULL,
                             w = 1,
                             mc.cluster.to.preserve = "seurat_clusters",
                             reduction.name = "rp",
                             reduction.key = "RPUMAP_",
                             seed = 42L,
                             verbose = TRUE){

  if(is.null(mc.umap.config)){
    mc.umap.name <- "umap"
    mc.umap.input.reduction.name = "harmony"
    mc.umap.dims = deparse(c(1:20))
    mc.umap.n.neighbors <- 30
    mc.umap.metric <- "cosine"
  }else{
    mc.umap.name <- mc.umap.config[[1]][["umap.name"]]
    mc.umap.input.reduction.name = mc.umap.config[[1]][["umap.input.reduction.name"]]
    mc.umap.dims = mc.umap.config[[1]][["umap.dims"]]
    mc.umap.n.neighbors <- mc.umap.config[[1]][["umap.n.neighbors"]]
    mc.umap.metric <-  mc.umap.config[[1]][["umap.metric"]]

  }
  # cat(mc.umap.name, mc.umap.input.reduction.name, mc.umap.dims, mc.umap.n.neighbors, mc.umap.metric,'\n')

  if(is.null(Misc(mc[[mc.umap.name]], slot = mc.umap.metric))){
    if(verbose){message("Re-running UMAP for the main object...")}
    mc.umap <- .RunUMAP(mc,
                        reduction = mc.umap.input.reduction.name,
                        dims = eval(parse(text = mc.umap.dims)),
                        n.neighbors = mc.umap.n.neighbors,
                        metric = mc.umap.metric)
    mc.nn <- mc.umap$nn[[mc.umap.metric]]
  }else{
    if(verbose){message("KNN for main cluster reduction ", mc.umap.name, " recycled.")}
    mc.nn <- Misc(mc[[mc.umap.name]], slot = mc.umap.metric)
  }


  sc.nn.list <- list()
  for(i in 1:length(sc.list)){
    sc <- sc.list[[i]]
    if(is.null(sc.umap.config)){
      umap.name <- "umap"
      umap.input.reduction.name <- "harmony"
      umap.dims = deparse(c(1:20))
      umap.n.neighbors <- 30
      umap.metric <- "cosine"
    }else{
      umap.name <- sc.umap.config[[names(sc.list)[i]]][["umap.name"]]
      umap.input.reduction.name = sc.umap.config[[names(sc.list)[i]]][["umap.input.reduction.name"]]
      umap.dims = sc.umap.config[[names(sc.list)[i]]][["umap.dims"]]
      umap.n.neighbors <- sc.umap.config[[names(sc.list)[i]]][["umap.n.neighbors"]]
      umap.metric <-  sc.umap.config[[names(sc.list)[i]]][["umap.metric"]]
      # cat(umap.name, umap.input.reduction.name, umap.n.neighbors, umap.metric,'\n')
    }

    if(is.null(Misc(sc[[umap.name]], slot = umap.metric))){
      if(verbose){message("Child object for ", names(sc.list)[i], " does not contain KNN. Re-running UMAP...")}
      new.umap <- .RunUMAP(sc,
                           reduction = umap.input.reduction.name,
                           dims = eval(parse(text = umap.dims)),
                           n.neighbors = umap.n.neighbors,
                           metric = umap.metric
                           )
      sc.nn.list[[names(sc.list)[i]]] <- new.umap$nn[[umap.metric]]
    }else{
      if(verbose){message("KNN for child object for ", names(sc.list)[i], " recycled.")}
      sc.nn.list[[names(sc.list)[i]]] <- Misc(sc[[umap.name]], slot = umap.metric)
    }
  }


  ## cell info
  mc.cells <- colnames(mc)
  mc.idx <- c(1:length(mc.cells))

  sc.cells.list <- lapply(sc.list, function(x) colnames(x))
  sc.cells <- Reduce(base::union, sc.cells.list)
  is.sc <- mc.cells %in% sc.cells

  ## Make intra-cluster nn for each sub-cluster
  mc.nn.new.intra <- mc.nn
  for (i in 1:length(sc.nn.list)){
    nn.dist.new <- sc.nn.list[[i]]$dist
    nn.idx.new <- sc.nn.list[[i]]$idx
    cells.new <- sc.cells.list[[i]]

    idx <- match(cells.new, mc.cells)
    nn.dist.old <- mc.nn$dist[idx,]

    nn.dist.new.adj <- nn.dist.new / median(nn.dist.new) * median(nn.dist.old)

    mc.nn.new.intra$dist[idx, ] <- nn.dist.new.adj

    mc.nn.new.intra$idx[idx, ] <- matrix(match(cells.new[nn.idx.new], mc.cells), nrow = length(cells.new))

  }

  # mc.clusters <- mc$seurat_clusters
  mc.clusters <- mc[[mc.cluster.to.preserve, drop = TRUE]]
  K <- dim(mc.nn$idx)[2]

  ## new code with weight
  mc.nn.new.w <- list()
  mc.nn.new.w$idx <- matrix(0, nrow = dim(mc.nn$idx)[1], ncol = K)
  mc.nn.new.w$dist <- matrix(0, nrow = dim(mc.nn$idx)[1], ncol = K)

  for(i in 1:dim(mc.nn.new.w$idx)[1]){
    if(!is.sc[i]){
      ## use nn from main cluster if no subclustering is perfomed for a cell
      mc.nn.new.w$idx[i,] <- mc.nn$idx[i,]
      mc.nn.new.w$dist[i,] <- mc.nn$dist[i,]
    }else{
      is.reserved <- (mc.clusters[mc.nn$idx[i,]] != mc.clusters[i])
      k.res <- floor(length(which(is.reserved)) * w)
      if(k.res == (K - 1)){
        d <- mc.nn$dist[i,which(is.reserved)]
        idx <- mc.nn$idx[i,which(is.reserved)]
        mc.nn.new.w$idx[i,] <- c(i,idx)
        mc.nn.new.w$dist[i,] <- c(0,d)
      }else if(k.res > 0){
        # d <- c(mc.nn.new.intra$dist[i,1:(K - k.res)],mc.nn$dist[i,which(is.reserved)[1:k.res]])
        # idx <- c(mc.nn.new.intra$idx[i,1:(K - k.res)],mc.nn$idx[i,which(is.reserved)[1:k.res]])
        mc.nn.new.intra.d <- mc.nn.new.intra$dist[i,]
        mc.nn.new.intra.idx <- mc.nn.new.intra$idx[i,]
        intra.idx.diff <- base::setdiff(mc.nn.new.intra.idx, mc.nn$idx[i,which(is.reserved)[1:k.res]]) ## if a neighbor presents in both mc and sc, keep the mc
        m1 <- base::match(intra.idx.diff, mc.nn.new.intra.idx)
        intra.d.diff <- mc.nn.new.intra.d[m1]

        d <- c(intra.d.diff[1:(K - k.res)],mc.nn$dist[i,which(is.reserved)[1:k.res]])
        idx <- c(intra.idx.diff[1:(K - k.res)],mc.nn$idx[i,which(is.reserved)[1:k.res]])

        my.order <- order(d, decreasing = F)
        d1 <- d[my.order]
        idx1 <- idx[my.order]
        mc.nn.new.w$idx[i,] <- idx1
        mc.nn.new.w$dist[i,] <- d1
      }else if(k.res == 0){
        ## if no nn from mc to be reserved, simply use nn from sc
        mc.nn.new.w$idx[i,] <- mc.nn.new.intra$idx[i,]
        mc.nn.new.w$dist[i,] <- mc.nn.new.intra$dist[i,]
      }

    }
  }

  mc.umap.new <- .RunUMAP(nn_method = mc.nn.new.w, seed.use = seed)
  # plot(mc.umap.new$embedding[,1],mc.umap.new$embedding[,2], pch = 19, cex = .2)

  ## Create dimension reduction object named "RP"
  rp <- mc.umap.new$embedding
  colnames(rp) <- paste0("RP_", 1:2)
  rownames(rp) <- mc.cells

  mc[[reduction.name]] <- CreateDimReducObject(embeddings = rp, key = reduction.key, assay = DefaultAssay(mc))
  Misc(mc[[reduction.name]], slot = mc.umap.metric) <- mc.nn.new.w

  # print(table(apply(mc@reductions$rp@misc$cosine$idx, 1, function(x) length(unique(x))), is.sc))

  return(mc)

}


#' Prepare Seurat Object for Pseudobulk Analysis
#'
#' Convert columns to factors, drop unused levels, perform name check and make valid names
#' @noRd
#'
.PrepareSeuratObjectPseudobulkAnalysis <- function(object,
                                                   var.list,
                                                   check.names,
                                                   auto.make.names = FALSE,
                                                   verbose){
  if(is(object, "Seurat")){
    for(i in 1:length(var.list)){
      if(!is.factor(object[[var.list[i], drop = TRUE]])){
        if(verbose == TRUE){message("Converting ", var.list[i], " to factor")}
        object[[var.list[i]]] <- as.factor(object[[var.list[i], drop = TRUE]])
      }else{
        object[[var.list[i]]] <- droplevels(object[[var.list[i], drop = TRUE]])
      }
      if(check.names[i]){
        my.dat <- object[[var.list[i], drop = TRUE]]
        my.levels <- levels(my.dat)

        if(!all(make.names(my.levels) == my.levels)){
          my.levels.cat <- paste(my.levels, collapse = ',')
          if(auto.make.names == TRUE){
            message("The levels: ",my.levels.cat," in ", var.list[i], " are not synatically valid. ")
            object[[var.list[i]]] <- factor(make.names(my.dat), levels = make.names(my.levels))
            my.new.levels <- levels(object[[var.list[i], drop = TRUE]])
            my.new.levels.cat <- paste(my.new.levels, collapse = ',')
            message("Levels of ", var.list[i], " are automatically converted to ", my.new.levels.cat)

          }else{
            stop("The levels: ", my.levels.cat," in ", var.list[i], " are not synatically valid. Either manually change the names or set auto.make.names to TRUE")

          }
        }
      }

    }

    return(object)

  }else{
    stop("Invalid input object. Seurat object is required.")
  }
}

#' Check design
#'
#' @noRd

.CheckDesign <- function(sce,
                         design,
                         n = 1
){
  tbl.sample.cluster <- table(sce$sample_id, sce$cluster_id)
  ## print(tbl.sample.cluster)
  ## check overall design
  # print(design)
  # message(nrow(design), " ", ncol(design), " ", qr(design)$rank)
  if(qr(design)$rank == nrow(design)){stop("Invalid design: rank of the design matrix equals the number of samples.")}
  if(qr(design)$rank < ncol(design)){stop("Invalid design: design matrix not full column rank.")}

  ## check final design per identity class
  kids <- colnames(tbl.sample.cluster)
  for(i in 1:length(kids)){
    idx <- which(tbl.sample.cluster[,i] >= n)
    design.i <- design[idx, , drop = FALSE]
    ## print(design.i)
    # message(i, ':', dim(design.i)[2], ' ' ,qr(design.i)$rank, ' ',  dim(design.i)[1])
    if(qr(design.i)$rank == nrow(design.i)){stop(paste("Invalid design for cluster ", kids[i], ": rank of the design matrix equals the number of samples.", sep = ""))}
    if(qr(design.i)$rank < ncol(design.i)){stop(paste("Invalid design for cluster ", kids[i], ": design matrix not full column rank.", sep = ""))}
  }

  return(invisible(x = NULL))
}

#' Check consistency between two metadata columns for SummarizedHeatmap
#'
#' @noRd

.CheckMetadataSummarizedHeatmap <- function(object,
                                            var1,
                                            var2){
  mdat1 <- object[[var1, drop = TRUE]]
  mdat2 <- object[[var2, drop = TRUE]]
  ct <- table(mdat1, mdat2)
  if(any(rowSums(ct>0) >1))
  {
    stop("The argument used to split the columns (\"",var1,"\") is not compatible with addtional metadata (\"",var2,"\") for column annotation.\n? RunSummarizedHeatmap for more details.")
  }

  return(invisible(x = NULL))
}

#' Create UMAP configuration for re-projection
#'
#' Allow users to provide custom parameters to retrieve KNN data or re-run UMAP for re-projection
#'
#' @param type Type of UMAP configure, either "main" for main cluster or "sc" for subcluster
#' @param sc.name Name of the child object for subcluster analysis, which should be consistent with the \code{child.object.list}
#'  argument of the \code{\link[Ragas]{CreatePostIntegrationObject}} function
#' @param umap.name Name of the UMAP object (default: umap)
#' @param umap.input.reduction.name Name of the reduction as input to re-run UMAP (default: harmony)
#' @param umap.dims Input dimensions to run UMAP with. Same as the \code{dims} argument of the \code{\link[Seurat]{RunUMAP}} function (default: 1:20).
#' @param umap.n.neighbors Number of neighbors for the KNN graph. Same as the \code{n.neighbors} argument of the \code{\link[Seurat]{RunUMAP}} function (default: 30)
#' @param umap.metric Distance metric. Same as the \code{metric} argument of the \code{\link[Seurat]{RunUMAP}} function (default: "cosine")
#' @param append.to Existing subcluster configuration to append to
#'
#' @return A list containing configurations for UAMP re-projection
#' @export
#'
#' @examples
#' \dontrun{
#' my.config <- ConfigureReprojection(type = "sc", sc.name = "Bcell")
#' my.config <- ConfigureReprojection(type = "sc", sc.name = "Tcell", append.to = my.config)
#' }
#'
ConfigureReprojection <- function(type = c("main","sc"),
                                  sc.name = NULL,
                                  umap.name = "umap",
                                  umap.input.reduction.name = "harmony",
                                  umap.dims = 1:20,
                                  umap.n.neighbors = 30,
                                  umap.metric = "cosine",
                                  append.to = NULL
){
  type <- match.arg(type)
  if(is.null(sc.name) && type == "sc"){stop("sc.name cannot be NULL for subcluster analysis!")}
  if(!is.null(append.to) && type == "main"){
    warning("Argument append.to not applicable when type = main! Ignored.")
    append.to <- NULL}
  if(is.null(append.to)){my.config <- list()}else{my.config = append.to}
  if(type == "main"){
    my.config[["main.cluster"]] <- list("umap.name" = umap.name,
                                "umap.input.reduction.name" = umap.input.reduction.name,
                                "umap.dims" = deparse(umap.dims),
                                "umap.n.neighbors" = umap.n.neighbors,
                                "umap.metric" = umap.metric)
  }else{
    my.config[[sc.name]] <- list("umap.name" = umap.name,
                                 "umap.input.reduction.name" = umap.input.reduction.name,
                                 "umap.dims" = deparse(umap.dims),
                                 "umap.n.neighbors" = umap.n.neighbors,
                                 "umap.metric" = umap.metric)
  }
  all.n.neighbors <- unlist(lapply(my.config, function(x) x$umap.n.neighbors ))
  if(length(unique(all.n.neighbors)) > 1){
    stop("Multiple n.neighbors provided (", paste(all.n.neighbors, collapse = ","),"), but current algorithm only supports integration of nearest-neighbor graphs with the same n.neighbors!")
  }
  return(my.config)
}
