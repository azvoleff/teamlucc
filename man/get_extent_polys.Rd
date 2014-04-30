\name{get_extent_polys}
\alias{get_extent_polys}
\title{Generate a SpatialPolygonDataFrame of raster extents}
\usage{
get_extent_polys(x)
}
\arguments{
  \item{x}{a \code{Raster*} object, or \code{list}
  of\code{Raster*} objects}
}
\value{
\code{SpatialPolygonDataFrame} with the extent of each
raster object as a polygon, with a "filename" attribute
giving the filename for the raster object from with each
extent is derived.
}
\description{
Also includes the filename associated with each raster
object. Useful for providing the \code{dem_extents}
argument to the \code{\link{auto_setup_dem}} function.
}

