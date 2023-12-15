#' Export Pi Data
#'
#' Function to export data from \code{\link[Ragas]{Pi}} object
#'
#' This function takes in a Pi object, field and key parameters to export data
#'
#' @param object  A \code{\link[Ragas]{Pi}} object
#' @param field Name of pi object component to be exported
#' @param key A unique identifier to retrieve data from the \code{field} parameter
#' @param ds.export.cpm Option to export CPM data. Default: FALSE
#' @param ds.export.freq Option to export frequency data. Default: FALSE
#' @param export.path path to the directory where data should be exported. Default to current working directory
#' @examples
#' \dontrun{
#' ExportPiData(object = csle.pbmc.pi,
#'              field = 'exp.freq',
#'              key = 'ExpFreq|subcluster_idents|cutoff=0'
#'            )
#'
#' }
#' @export
#'
ExportPiData <- function(object, field = NULL, key = NULL, file.prefix = NULL, ds.export.cpm = FALSE, ds.export.freq = FALSE, export.path = NULL)
{
  if(!(is(object, "Pi"))){stop("Invalid input argument \"object\". Pi object is required.")}
  if(is.null(file.prefix)){stop('The argument file.prefix must be provided. It is recommended to use informative prefix such as \'pbmc.ds\'.')}
  if(is.null(export.path))
  {
    export.path <- getwd()
  }
  if(!(dir.exists(export.path)))
  {
    dir.create(export.path)
  }

  if(is.null(field))
  {
    stop("field cannot be NULL. The following field(s) are available:\n    ", paste(c("exp.freq","markers","ds"), collapse = "\n    "))
  }else if(field == 'exp.freq')
  {
    if(length(object$exp.freq) == 0)
    {
      stop("The \"exp.freq\" field of the Pi object is empty. Run CalculateExpFreqs on the Pi object first.")
    }else
    {
      if(is.null(key))
      {
        stop("key cannot be NULL. The following key(s) are available:\n    ", paste(names(object$exp.freq), collapse = "\n    "),"\nEither choose an availble key or rerun CalculateExpFreqs.")
      }else if(!key %in% names(object$exp.freq)){
        stop("Invalid key \"",key, "\" for the exp.freq object. The following key(s) are available:\n    ", paste(names(object$exp.freq), collapse = "\n    "),"\nEither choose an availble key or rerun CalculateExpFreqs.")
      }else
      {
        export.data <- object$exp.freq[[key]][["data"]]
        names(export.data) <- gsub('/','_',names(export.data))

        dirName <- file.path(export.path,'exp.freq')
        dir.create(path = dirName)
        my.write.fun <- function(x,y)
        {
          fileName <- paste(file.prefix, '_', y,'.csv',sep = '')
          write.csv(x, file = file.path(dirName,fileName),row.names = TRUE)
        }
        invisible(mapply(my.write.fun, export.data, names(export.data) ))
      }
    }
  }else if(field == 'markers')
  {
    if(length(object$markers) == 0)
    {
      stop("The \"markers\" field of the Pi object is empty. Run RunFindAllMarkers on the Pi object first.")
    }else
    {
      if(is.null(key))
      {
        stop("key cannot be NULL. The following key(s) are available:\n    ", paste(names(object$markers), collapse = "\n    "),"\nEither choose an availble key or rerun RunFindAllMarkers.")
      }else if(!key %in% names(object$markers)){
        stop("Invalid key \"",key, "\" for the markers object. The following key(s) are available:\n    ", paste(names(object$markers), collapse = "\n    "),"\nEither choose an availble key or rerun RunFindAllMarkers.")
      }else
      {
        export.data <- object$markers[[key]][["data"]]
        fileName <- paste(file.prefix, '.csv',sep = '')

        dirName <- file.path(export.path,'markers')
        dir.create(path = dirName)

        write.csv(export.data, file = file.path(dirName,fileName),row.names = TRUE)
      }
    }
  }else if(field == 'ds')
  {
    if(length(object$ds) == 0)
    {
      stop("The \"ds\" field of the Pi object is empty. Run RunPseudobulkAnalysis on the Pi object first.")
    }else
    {
      if(is.null(key))
      {
        stop("key cannot be NULL. The following key(s) are available:\n    ", paste(names(object$ds), collapse = "\n    "),"\nEither choose an availble key or rerun RunPseudobulkAnalysis")
      }else if(!key %in% names(object$ds)){
        stop("Invalid key \"",key, "\" for the ds object. The following key(s) are available:\n    ", paste(names(object$ds), collapse = "\n    "),"\nEither choose an availble key or rerun RunPseudobulkAnalysis")
      }else
      {
        export.data <- object$ds[[key]][["data"]]$tbl
        fileName <- paste(file.prefix, '.csv',sep = '')

        dirName <- file.path(export.path,'ds')
        dir.create(path = dirName)

        if(isTRUE(ds.export.cpm) & isTRUE(ds.export.freq))
        {
          write.csv(export.data, file = file.path(dirName,fileName),row.names = FALSE)
        }else if(isTRUE(ds.export.cpm) & isFALSE(ds.export.freq))
        {
          idx <- grep('frq',colnames(export.data))
          export.data <- export.data[,-idx]
          write.csv(export.data, file = file.path(dirName,fileName),row.names = FALSE)
        }else if(isFALSE(ds.export.cpm) & isTRUE(ds.export.freq))
        {
          idx <- grep('cpm',colnames(export.data))
          export.data <- export.data[,-idx]
          write.csv(export.data, file = file.path(dirName,fileName),row.names = FALSE)
        }else if(isFALSE(ds.export.cpm) & isFALSE(ds.export.freq))
        {
          idx <- c(grep('cpm',colnames(export.data)), grep('frq',colnames(export.data)))
          export.data <- export.data[,-idx]
          write.csv(export.data, file = file.path(dirName,fileName),row.names = FALSE)
        }
      }
    }
  }else if(field == 'cell.prop')
  {
    if(length(object$cell.prop) == 0)
    {
      stop("The \"cell.prop\" field of the Pi object is empty. Run RunProportionPlot on the Pi object first.")
    }else
    {
      if(is.null(key))
      {
        stop("key cannot be NULL. The following key(s) are available:\n    ", paste(names(object$cell.prop), collapse = "\n    "),"\nEither choose an availble key or rerun RunProportionPlot.")
      }else if(!key %in% names(object$cell.prop)){
        stop("Invalid key \"",key, "\" for the cell.prop object. The following key(s) are available:\n    ", paste(names(object$cell.prop), collapse = "\n    "),"\nEither choose an availble key or rerun RunProportionPlot")
      }else
      {
        # export data
        export.data <- object$cell.prop[[key]]$data$data

        fileName <- paste(file.prefix, '.data.csv',sep = '')

        dirName <- file.path(export.path,'cell.prop')
        dir.create(path = dirName)

        write.csv(export.data, file = file.path(dirName,fileName),row.names = FALSE)

        # export stats (if applicable)
        export.stats <- object$cell.prop[[key]]$data$stats
        if(!is.null(export.stats))
        {
          fileName <- paste(file.prefix, '.stats.csv',sep = '')
          write.csv(export.stats, file = file.path(dirName,fileName),row.names = FALSE)
        }
      }
    }
  }else
  {
    stop("Invalid field \"",field, "\". The following field(s) are available:\n    ", paste(c("exp.freq","markers","ds"), collapse = "\n    "))
  }

}
