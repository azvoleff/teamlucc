\name{class_statistics}
\alias{class_statistics}
\title{Exports statistics on pixels within each of a set of land cover classes}
\usage{
class_statistics(x, y, classcol)
}
\arguments{
  \item{x}{A \code{RasterLayer} from which class statistics
  will be calculated.}

  \item{y}{A \code{SpatialPolygonsDataFrame} with cover
  class polygons}

  \item{classcol}{the name of the column containing the
  response variable (for example the land cover type of
  each pixel)}
}
\value{
A data.frame of class statistics.
}
\description{
Exports statistics on pixels within each of a set of land
cover classes
}
\examples{
\dontrun{
# Don't run due to raster bug that will be fixed in next version of raster
#package (as of 1/8/2014). See: http://bit.ly/1fe5x0m
class_statistics(L5TSR_1986, L5TSR_1986_2001_training, "class_1986")
}
}

