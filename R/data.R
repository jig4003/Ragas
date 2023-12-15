#' A small example of PBMC data from child SLE
#'
#' This object contains 3000 down-sampled PBMC cells from Nehar-Belaid et al., 2020
#'
#' @format A Seurat object
#' @source \url{https://www.nature.com/articles/s41590-020-0743-0}
#'
"csle.pbmc.small"


#' A small example of B cell data from child SLE
#'
#' This object contains down-sampled cells from Nehar-Belaid et al., 2020
#'
#' @format A Seurat object
#' @source \url{https://www.nature.com/articles/s41590-020-0743-0}
#'
"csle.bcell.small"

#' A small example of T cell data from child SLE
#'
#' This object contains down-sampled cells from Nehar-Belaid et al., 2020
#'
#' @format A Seurat object
#' @source \url{https://www.nature.com/articles/s41590-020-0743-0}
#'
"csle.tcell.small"

#' A small example of Treg data from child SLE
#'
#' This object contains down-sampled cells from Nehar-Belaid et al., 2020
#'
#' @format A Seurat object
#' @source \url{https://www.nature.com/articles/s41590-020-0743-0}
#'
"csle.treg.small"

#' A small example of CD4+ memory T cell data from child SLE
#'
#' This object contains down-sampled cells from Nehar-Belaid et al., 2020
#'
#' @format A Seurat object
#' @source \url{https://www.nature.com/articles/s41590-020-0743-0}
#'
"csle.cd4.mem.small"

#' An example Pi object built from csle.pbmc.small
#'
#' This data contained 3000 down-sampled PBMC cells from Nehar-Belaid et al., 2020.
#' Run the example code to re-produce the Pi object from Seurat object \code{\link[Ragas]{csle.pbmc.small}}:
#'
#' @examples
#' \dontrun{
#' ## A list of subclusters to combine
#' subclusters <- list("B cell" = "csle.bcell.small", "T cell" = "csle.tcell.pi")
#'
#' ## UMAP configuration
#' sc.config <- ConfigureReprojection(type = "sc", "B cell")
#' sc.config <- ConfigureReprojection(type = "sc", "T cell", umap.name = "rp", append.to = sc.config)
#'
#' ## select columns for cell identities
#' sc.colnames <- c("B cell" = "cluster.annotation", "T cell" = "subcluster_idents")
#'
#' ## integrate subclusters
#' library(dplyr)
#' csle.pbmc.pi <- csle.pbmc.small %>%
#' CreatePostIntegrationObject(child.object.list = subclusters,
#'                             rp.subcluster.umap.config = sc.config,
#'                             rp.main.cluster.anno = "cluster.annotation",
#'                             rp.subcluster.colname = sc.colnames) %>%
#' CalculateExpFreqs(ident = "subcluster_idents") %>%
#' RunFindAllMarkers(ident = "subcluster_idents") %>%
#' RunPseudobulkAnalysis(ident.var = "subcluster_idents",
#'                       group.var = "Groups",
#'                       sample.var = "Names",
#'                       group.1 = "cSLE",
#'                       group.2 = "cHD") %>%
#' RunProportionPlot(group.by = "Groups",
#'                   ident = "subcluster_idents",
#'                   method = "unpooled",
#'                   unpool.by = "Names"
#'                   )
#'
#'}
#' @format A \code{\link[Ragas]{Pi}} object
#' @source \url{https://www.nature.com/articles/s41590-020-0743-0}
#'
"csle.pbmc.pi"

#' An example Pi object built from csle.bcell.small
#'
#' This data contained down-sampled B cells from Nehar-Belaid et al., 2020.
#' Run the example code to re-produce the Pi object from Seurat object \code{\link[Ragas]{csle.bcell.small}}:
#'
#' @examples
#' \dontrun{
#' csle.bcell.pi <- csle.bcell.small %>%
#' CreatePostIntegrationObject(parent.object = csle.pbmc.small) %>%
#' CalculateExpFreqs(ident = "cluster.annotation") %>%
#' RunFindAllMarkers(ident = "cluster.annotation") %>%
#' RunPseudobulkAnalysis(ident.var = "cluster.annotation",
#'                       group.var = "Groups",
#'                       sample.var = "Names",
#'                       group.1 = "cSLE",
#'                       group.2 = "cHD") %>%
#' RunProportionPlot(group.by = "Groups",
#'                   ident = "cluster.annotation",
#'                   method = "unpooled",
#'                   unpool.by = "Names",
#'                   use.parent.as.ref = TRUE
#'                   )
#'}
#' @format A \code{\link[Ragas]{Pi}} object
#' @source \url{https://www.nature.com/articles/s41590-020-0743-0}
#'
"csle.bcell.pi"


#' An example Pi object built from csle.tcell.small
#'
#' This data contained down-sampled T cells from Nehar-Belaid et al., 2020.
#' Run the example code to re-produce the Pi object from Seurat object \code{\link[Ragas]{csle.tcell.small}}:
#'
#' @examples
#' \dontrun{
#' ## A list of T cell subclusters to combine
#' subclusters <- list("CD4mem" = "csle.cd4.mem.small", "Treg" = "csle.treg.small")
#'
#' csle.tcell.pi <- csle.tcell.small %>%
#' CreatePostIntegrationObject(child.object.list = subclusters,
#'                             rp.main.cluster.anno = "cluster.annotation",
#'                             rp.subcluster.colname = "cluster.annotation",
#'                             parent.object = csle.pbmc.small) %>%
#' CalculateExpFreqs(ident = "subcluster_idents") %>%
#' RunFindAllMarkers(ident = "subcluster_idents") %>%
#' RunPseudobulkAnalysis(ident.var = "subcluster_idents",
#'                       group.var = "Groups",
#'                       sample.var = "Names",
#'                       group.1 = "cSLE",
#'                       group.2 = "cHD") %>%
#' RunProportionPlot(ident = "subcluster_idents",
#'                   group.by = "Groups",
#'                   method = "unpooled",
#'                   unpool.by = "Names",
#'                   use.parent.as.ref = TRUE
#'                   )
#'
#'}
#' @format A \code{\link[Ragas]{Pi}} object
#' @source \url{https://www.nature.com/articles/s41590-020-0743-0}
#'
"csle.tcell.pi"

#' An example Pi object built from csle.treg.small
#'
#' This data contained down-sampled Treg cells from Nehar-Belaid et al., 2020.
#' Run the example code to re-produce the Pi object from Seurat object \code{\link[Ragas]{csle.treg.small}}:
#'
#' @examples
#' \dontrun{
#' csle.treg.pi <- csle.treg.small %>%
#' CreatePostIntegrationObject(parent.object = csle.tcell.small) %>% ## use total T cell as reference for later proportion analysis
#' CalculateExpFreqs(ident = "cluster.annotation") %>%
#' RunFindAllMarkers(ident = "cluster.annotation") %>%
#' RunPseudobulkAnalysis(ident.var = "cluster.annotation",
#'                       group.var = "Groups",
#'                       sample.var = "Names",
#'                       group.1 = "cSLE",
#'                       group.2 = "cHD") %>%
#' RunProportionPlot(group.by = "Groups",
#'                   method = "unpooled",
#'                   unpool.by = "Names",
#'                   use.parent.as.ref = TRUE
#'                   )
#'
#'}
#' @format A \code{\link[Ragas]{Pi}} object
#' @source \url{https://www.nature.com/articles/s41590-020-0743-0}
#'
"csle.treg.pi"

#' An example Pi object built from csle.cd4.mem.small
#'
#' This data contained down-sampled CD4+ memory T cells from Nehar-Belaid et al., 2020.
#' Run the example code to re-produce the Pi object from Seurat object \code{\link[Ragas]{csle.cd4.mem.small}}:
#'
#' @examples
#' \dontrun{
#' csle.cd4.mem.pi <- csle.cd4.mem.small %>%
#' CreatePostIntegrationObject(parent.object = csle.tcell.small) %>% ## use total T cell as reference for later proportion analysis
#' CalculateExpFreqs(ident = "cluster.annotation") %>%
#' RunFindAllMarkers(ident = "cluster.annotation") %>%
#' RunPseudobulkAnalysis(ident.var = "cluster.annotation",
#'                       group.var = "Groups",
#'                       sample.var = "Names",
#'                       group.1 = "cSLE",
#'                       group.2 = "cHD") %>%
#' RunProportionPlot(group.by = "Groups",
#'                   method = "unpooled",
#'                   unpool.by = "Names",
#'                   use.parent.as.ref = TRUE
#'                   )
#'
#'}
#' @format A \code{\link[Ragas]{Pi}} object
#' @source \url{https://www.nature.com/articles/s41590-020-0743-0}
#'
"csle.cd4.mem.pi"
