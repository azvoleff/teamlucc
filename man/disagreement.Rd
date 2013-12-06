\name{disagreement}
\alias{disagreement}
\title{Calculate quantity and allocation disagreement for two maps}
\usage{
  disagreement(ct, pop = NULL)
}
\arguments{
  \item{an}{observed sample matrix as calculated by
  \code{accuracy}}

  \item{pop}{(optional) either 1) a list of length equal to
  the number of classes in ct giving the total number of
  pixels in the entire population for each class, or 2) the
  predicted cover map as a \code{RasterLayer}, from which
  the population frequencies will be derived.}
}
\value{
  list with two entries: Q (quantity disagreement) and A
  (allocation disagreement)
}
\description{
  This function calculate quantity disagreement \code{Q}
  and allocation disagreement \code{A} for use in accuracy
  assessment and map comparison. These quantities are
  calculated based on Pontius and Millones (2011).
}
\references{
  Pontius, R. G., and M. Millones. 2011. Death to Kappa:
  birth of quantity disagreement and allocation
  disagreement for accuracy assessment. International
  Journal of Remote Sensing 32:4407-4429.
}

