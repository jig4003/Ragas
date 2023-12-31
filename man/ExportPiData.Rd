% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExportPiData.R
\name{ExportPiData}
\alias{ExportPiData}
\title{Export Pi Data}
\usage{
ExportPiData(
  object,
  field = NULL,
  key = NULL,
  file.prefix = NULL,
  ds.export.cpm = FALSE,
  ds.export.freq = FALSE,
  export.path = NULL
)
}
\arguments{
\item{object}{A \code{\link[Ragas]{Pi}} object}

\item{field}{Name of pi object component to be exported}

\item{key}{A unique identifier to retrieve data from the \code{field} parameter}

\item{ds.export.cpm}{Option to export CPM data. Default: FALSE}

\item{ds.export.freq}{Option to export frequency data. Default: FALSE}

\item{export.path}{path to the directory where data should be exported. Default to current working directory}
}
\description{
Function to export data from \code{\link[Ragas]{Pi}} object
}
\details{
This function takes in a Pi object, field and key parameters to export data
}
\examples{
\dontrun{
ExportPiData(object = csle.pbmc.pi,
             field = 'exp.freq',
             key = 'ExpFreq|subcluster_idents|cutoff=0'
           )

}
}
