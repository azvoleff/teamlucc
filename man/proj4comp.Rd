\name{proj4comp}
\alias{proj4comp}
\title{Performs a rough comparison of two proj4strings to see if they match}
\usage{
proj4comp(x, y)
}
\arguments{
  \item{x}{}

  \item{y}{}
}
\value{
TRUE or FALSE depending on if the projections match
}
\description{
Compares the proj, ellps, zone (if applicable), units, and
datum tags in two proj4strings to determine if two
projections match. Requires proj and ellps to be present.
If present, zone, units, and datum must match in both
strings.
}

