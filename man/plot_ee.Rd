\name{plot_ee}
\alias{plot_ee}
\title{Plot EarthExplorer metadata file}
\usage{
plot_ee(x, start_year, end_year, min_clear = 0.7)
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
}
\value{
used for side effect of producing a plot
}
\description{
Plot EarthExplorer metadata file
}

