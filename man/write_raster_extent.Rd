\name{write_raster_extent}
\alias{write_raster_extent}
\title{Makes a polygon shapefile of a \code{Raster*} extent.}
\usage{
  write_raster_extent(x, out_file,
    fields = c("t1_class", "t2_class"))
}
\arguments{
  \item{x}{A \code{Raster*} object}

  \item{out_file}{Filename for the output shapefile}

  \item{fields}{list of fields to include in the output
  shapefile}
}
\description{
  Makes a polygon shapefile of a \code{Raster*} extent.
}

