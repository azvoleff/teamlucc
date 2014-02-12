\name{team_QA_stats}
\alias{team_QA_stats}
\title{Calculate statistics on imagery within an AOI}
\usage{
team_QA_stats(image_dirs, aoi)
}
\arguments{
  \item{image_dirs}{list of paths to a set of Landsat CDR
  image files in ENVI format as output by the
  \code{unstack_ledapscdr} function.}

  \item{aoi}{an area of interest (AOI) to crop from each
  image}
}
\value{
a \code{data.frame}
}
\description{
Calculate statistics on imagery within an AOI
}

