\name{ee_plot}
\alias{ee_plot}
\title{Plot EarthExplorer scene list}
\usage{
ee_plot(x, start_date, end_date, min_clear = 0.7, exclude = list(),
  normalize = FALSE)
}
\arguments{
  \item{x}{a \code{data.frame} with a list of Landsat
  scenes as output by \code{\link{ee_read}}}

  \item{start_date}{starting date as a \code{Date} object}

  \item{end_date}{end date as a \code{Date} object}

  \item{min_clear}{the minimum percent clear to plot
  (calculated as 1 - percent cloud cover). Images with less
  than \code{min_clear} fraction of the image area free of
  clouds will be ignored.}

  \item{exclude}{a list of sensors to exclude (for example,
  set \code{exclude=c('LE7', 'LT4')} to exclude Landsat 7
  ETM+ and Landsat 4 TM images.}

  \item{normalize}{if \code{TRUE}, plot as a normalized
  line plot}
}
\value{
used for side effect of producing a plot
}
\description{
This function can produce two different types of plots from
a USGS EarthExplorer Landsat CDR Surface Reflectance scene
list.
}

