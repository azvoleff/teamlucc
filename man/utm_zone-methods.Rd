\docType{methods}
\name{utm_zone}
\alias{utm_zone}
\alias{utm_zone,Spatial,ANY-method}
\alias{utm_zone,Spatial,logical-method}
\alias{utm_zone,numeric,numeric,logical-method}
\alias{utm_zone,numeric,numeric-method}
\title{Given a spatial object, calculate the UTM zone of the centroid}
\usage{
utm_zone(x, y, ...)

\S4method{utm_zone}{numeric,numeric}(x, y, proj4string = FALSE)

\S4method{utm_zone}{Spatial,ANY}(x, proj4string = FALSE)
}
\arguments{
  \item{proj4string}{if FALSE (default) return the UTM zone
  as a string (for example "34S" for UTM Zone 34 South). If
  TRUE, return a proj4string using the EPSG code as an
  initialization string.}

  \item{x}{a longitude (with western hemisphere longitudes
  negative)}

  \item{y}{a latitude (with southern hemisphere latitudes
  negative)}

  \item{x}{a Spatial object (UTM zone will be given for the
  object's centroid)}
}
\description{
For a line or polygon, the UTM zone of the centroid is
given, after reprojecting the object into WGS-84.
}
\details{
Based on the code on gis.stackexchange.com at
http://bit.ly/17SdcuN.
}
\examples{
utm_zone(45, 10)
utm_zone(45, -10)
utm_zone(45, 10, proj4string=TRUE)
}

