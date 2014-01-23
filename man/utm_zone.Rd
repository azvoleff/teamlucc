\name{utm_zone}
\alias{utm_zone}
\title{Given a spatial object, calculate the UTM zone of the centroid}
\usage{
utm_zone(long, lat, proj4string = FALSE)
}
\arguments{
  \item{long}{a longitude (with western hemisphere
  longitudes negative)}

  \item{lat}{a latitude (with southern hemisphere latitudes
  negative)}

  \item{proj4string}{if FALSE (default) return the UTM zone
  as a string (for example "34S" for UTM Zone 34 South). If
  TRUE, return a proj4string using the EPSG code as an
  initialization string.}
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

