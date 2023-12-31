% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RunStackedVlnPlot.R
\name{RunStackedVlnPlot}
\alias{RunStackedVlnPlot}
\title{Run Stacked VlnPLot for Seurat or post-integration (Pi) object}
\usage{
RunStackedVlnPlot(
  object,
  assay = "RNA",
  ident = "seurat_clusters",
  features = NULL,
  feature.annotation.cols = NULL,
  split.by = NULL,
  color.by = "features",
  clust.row = FALSE,
  clust.column = TRUE,
  add.points = FALSE,
  points.size = 0.1,
  column.fontsize = 12,
  row.fontsize = 8,
  row.title.fontsize = 15,
  legend.fontsize = 12,
  legend.title.fontsize = 15,
  features.fontsize = 12,
  column.names.rotation = 0,
  axis.text.hjust = 1,
  axis.text.vjust = 1,
  fill.colors = NULL,
  random.annotation.col.seed = 1
)
}
\arguments{
\item{object}{A Seurat or \code{\link[Ragas]{Pi}} object}

\item{assay}{Name of the assay to use (default: "RNA")}

\item{ident}{Name of the metadata column to be used as cell identity (default: seurat_clusters)}

\item{features}{User defined feature list in one of the following three formats:
\itemize{
\item{list} : a named list object. Each list element is a character vector of a group of features; the name of the list element will be extracted to represent the feature group
\item{data.frame}: a data.frame with two columns: "features" and "annotation"
\item{vector} : a character vector of features
}}

\item{feature.annotation.cols}{A character vector of colors for feature annotation. Its length must equal the number of feature groups provided by the "features" argument.
If the vector is named, its names should match the "features" argument (default: NULL)}

\item{split.by}{Factor to split the columns by (default: NULL)}

\item{color.by}{One of the five ways to color for the violin plot: "features", "clusters", "median.exp" (median expression), "mean.exp" (mean expression), "split.var" (the variable that splits the violin plot).
When \code{split.by} is set (not to NULL), \code{color.by} must set to "split.var"; and vice versa. (default: "features")}

\item{clust.row}{Whether to cluster the rows (default: FALSE)}

\item{clust.column}{Whether to cluster the columns (default: TRUE)}

\item{add.points}{Whether to add data points to the plot (default: FALSE)}

\item{points.size}{Point size (default: 0.1)}

\item{column.fontsize}{Column font size (default: 12)}

\item{row.fontsize}{Row font size (default: 8)}

\item{row.title.fontsize}{Row title font size (default: 15)}

\item{legend.fontsize}{Legend font size (default: 12)}

\item{legend.title.fontsize}{Legend title font size (default: 15)}

\item{features.fontsize}{Feature font size (default: 12)}

\item{column.names.rotation}{Rotation for column names (default: 0)}

\item{axis.text.hjust, axis.text.vjust}{Horizontal/vertical justification (in [0, 1]) for axis text (default: 1)}

\item{fill.colors}{Colors to fill the violin plot when \code{color.by} is set to "median.exp" or "mean.exp"}

\item{random.annotation.col.seed}{Seed to generate random colors for feature annotations}
}
\value{
A \code{\link[ggplot2]{ggplot}} object
}
\description{
A wrapper function for the \code{\link[Ragas]{StackedVlnPlot}}
}
\details{
The Stacked VlnPlot provides several new functionalities compared to \code{\link[Seurat]{VlnPlot}}:
(1) stacked output when multiple features are provided;(2) allows feature annotation; (3) allows clustering
of rows and/or columns; (4) provides multiple way to color volin plots, see argument \code{color.by} for more details.
}
\examples{
\dontrun{
my.list <- list(Monocytes = c("CD14","FCGR3A"),
                B = "MS4A1",
                T = c("CCR7", "CD8A"),
                NK = "NKG7")

## Seurat input; split.by is NULL; color by mean expression
RunStackedVlnPlot(csle.pbmc.small,
                  ident = "cluster.annotation",
                  features = my.list,
                  feature.annotation.cols = c( 'Monocytes'= 'tomato', 'B' = 'seagreen', 'T' = 'steelblue', 'NK' = 'purple'),
                  color.by = "mean.exp",
                  column.names.rotation = 90)
## Pi input; split.by is not NULL
RunStackedVlnPlot(csle.pbmc.pi,
                  ident = "cluster.annotation",
                  features = my.list,
                  split.by = "Groups",
                  color.by = "split.var",
                  column.names.rotation = 90)
}

}
\references{
Hao and Hao et al. (2021). “Integrated analysis of multimodal single-cell data.” Cell. doi:10.1016/j.cell.2021.04.048
}
