\name{simplify_polygon}
\alias{simplify_polygon}
\title{Simplify a polygon to contain less than a certain number of vertices}
\usage{
simplify_polygon(poly_obj, max_vertices, maxit = 100, multiplier = 1.05,
  initial_tolerance = "dynamic")
}
\arguments{
  \item{poly_obj}{a polygon as an sp object}

  \item{max_vertices}{the maximum number of vertices to
  allow in the simplified polygon}

  \item{maxit}{the maximum number of iterations to use to
  simplify the polygon}

  \item{multiplier}{a number used to increase the tolerance
  for each subsequent iteration, if the number of vertices
  in the simplified polygon is not less than
  /code{max_vertices}}

  \item{initial_tolerance}{initial value for tolerance used
  to remove vertices from polygon.  If set to the default
  option, "dynamic", the code will automatically set the
  \code{initial_tolerance} to .01 * the length of the
  diagonal of the bounding box of the polygon.
  \code{initial_tolerance} can also be set to an arbitrary
  value, in the same units as the polygon object.}
}
\value{
polygon with less than \code{max_vertices} vertices
}
\description{
Useful for simplifying area of interest (AOI) polygons for
use with online data portals (like USGS EarthExplorer) that
limit the number of vertices allowed in uploaded AOI
shapefiles.
}
\examples{
# TODO: add an example
}

