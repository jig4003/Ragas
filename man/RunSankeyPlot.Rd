% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RunSankeyPlot.R
\name{RunSankeyPlot}
\alias{RunSankeyPlot}
\title{Run Sankey Plot}
\usage{
RunSankeyPlot(
  object,
  meta.1,
  meta.2,
  node.colors = NULL,
  node.outline.color = "black",
  node.text.size = 5,
  node.text.color = "black",
  node.text.hjust = 0,
  node.text.vjust = 0.5,
  axis.text.size = 20,
  random.col.seed = 1
)
}
\arguments{
\item{object}{A \code{\link[Ragas]{Pi}} object or Seurat object}

\item{meta.1}{Name of metadata column to map from}

\item{meta.2}{Name of metadata column to map to}

\item{node.colors}{A character vector with colors for the node. The length must be the same as the unique number of items in meta.1 and meta.2}

\item{node.outline.color}{Color for node outline (default: black)}

\item{node.text.size}{Text size for node (default: 5)}

\item{node.text.color}{Color for node text (default: black)}

\item{node.text.hjust}{Horizontal justification (in [0, 1]) for node text}

\item{node.text.vjust}{Vertical justification (in [0, 1]) for node text}

\item{axis.text.size}{Axis text size}

\item{random.col.seed}{Set seed to control colors}
}
\value{
A \code{\link[ggplot2]{ggplot}} object
}
\description{
A wrapper function for \code{\link[Ragas]{SankeyPlot}}
}
\details{
This function takes 2 metadata column from seurat object and creates a sankey diagram
}
\examples{
\dontrun{
## Sankey plot
RunSankeyPlot(csle.pbmc.pi,
              meta.1 = "cluster.annotation", meta.2 = "subcluster_idents")
}

}
