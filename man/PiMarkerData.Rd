% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Pi.R
\name{PiMarkerData}
\alias{PiMarkerData}
\alias{CheckPiData.PiMarkerData}
\alias{print.PiMarkerData}
\alias{[[.PiMarkerData}
\title{PiMarkerData Class}
\usage{
PiMarkerData(
  data,
  ident,
  method = c("AllMarkers"),
  test.use,
  latent.vars = NULL
)

\method{CheckPiData}{PiMarkerData}(x, seurat.obj)

\method{print}{PiMarkerData}(x, ...)

\method{[[}{PiMarkerData}(x, i)
}
\arguments{
\item{data}{A data.frame containing output from \code{\link[Ragas]{RunFindAllMarkers}}}

\item{ident}{Name for the identity data}

\item{method}{A character string describing the analysis method}

\item{test.use}{Which test to use. See details in \code{\link[Seurat]{FindAllMarkers}}}

\item{latent.vars}{Latent variables to test. See details in \code{\link[Seurat]{FindAllMarkers}}}

\item{x}{A PiMarkerData object}

\item{seurat.obj}{A Seurat object}

\item{...}{Additional arguments}

\item{i}{Index for the list element to be extracted}
}
\description{
The PiMarkerData inherts \code{\link[Ragas]{PiData}} class, which stores marker data from \code{\link[Ragas]{RunFindAllMarkers}}.
}
\examples{
\dontrun{
## Run marker analysis and construct a PiMarkerData object
res <- RunFindAllMarkers(csle.bcell.small,
                         ident = "seurat_clusters",
                         test.use = "wilcox",
                         latent.vars = NULL)
mk <- PiMarkerData(data = res,
                   ident = "seurat_clusters",
                   method = "AllMarkers",
                   test.use = "wilcox",
                   latent.vars = NULL)
## Print object
mk

## Check object validity
CheckPiData(mk, seurat.obj = csle.bcell.small)

## Add PiMarkerData object to an existing Pi object
my.pi <- CreatePostIntegrationObject(csle.bcell.small) ## a minimum Pi object
my.pi <- AddPiData(my.pi, mk)
}

}
