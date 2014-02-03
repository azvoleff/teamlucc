\name{plot_ee_line}
\alias{plot_ee_line}
\title{Plot EarthExplorer metadata file}
\usage{
plot_ee_line(x, start_year, end_year, min_clear = 0.7, exclude = list())
}
\arguments{
  \item{x}{a \code{data.frame} with a list of Landsat
  scenes as output from the save metadata function on
  http://earthexplorer.usgs.gov}

  \item{start_year}{an integer indicating the first year to
  plot}

  \item{end_year}{an integer indicating the last year to
  plot}

  \item{min_clear}{the minimum percent clear to plot
  (calculated as 1 - percent cloud cover). Images with less
  than \code{min_clear} fraction of the image area free of
  clouds will be ignored.}

  \item{exclude}{a list of sensors to exclude (for example,
  set \code{exclude=c('LE7', 'LT4')} to exclude Landsat 7
  ETM+ and Landsat 4 TM images.}
}
\value{
used for side effect of producing a plot
}
\description{
Plot EarthExplorer metadata file
}

