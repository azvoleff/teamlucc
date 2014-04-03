\docType{methods}
\name{adj_areas}
\alias{adj_areas}
\alias{adj_areas,accuracy,missing-method}
\alias{adj_areas,numeric,matrix-method}
\alias{adj_areas,numeric,missing-method}
\alias{adj_areas,numeric,table-method}
\title{Calculated adjusted class areas for an image classification}
\usage{
adj_areas(x, y)

\S4method{adj_areas}{numeric,table}(x, y)

\S4method{adj_areas}{numeric,matrix}(x, y)

\S4method{adj_areas}{accuracy,missing}(x)
}
\arguments{
  \item{x}{an \code{accuracy} object or a list of
  populations as a \code{numeric}}

  \item{y}{missing, or a contingency table}
}
\description{
Calculates the adjusted areas of each class in an image
after taking account of omission and commission errors. For
unbiased adjustments, error rates should be calculated
using a population sample matrix (see
\code{\link{accuracy}}.
}
\details{
Standard errors for the adjusted areas are calculated as in
Olofsson et al. (2013).
}
\references{
Olofsson, P., G. M. Foody, S. V. Stehman, and C. E.
Woodcock. 2013. Making better use of accuracy data in land
change studies: Estimating accuracy and area and
quantifying uncertainty using stratified estimation. Remote
Sensing of Environment 129:122-131.
}

