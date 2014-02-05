\name{team_fill_gaps}
\alias{team_fill_gaps}
\title{Perform SLC-off gap fill of Landsat 7 ETM+ image}
\usage{
team_fill_gaps(sample_size = 20, size_wind = 12, class_num = 4,
  DN_min = 0, DN_max = 1, patch_long = 500)
}
\arguments{
  \item{sample_size}{the sample size of sample pixels}

  \item{size_wind}{the maximum window size}

  \item{class_num}{the estimated number of classes}

  \item{DN_min}{the minimum DN value of the image}

  \item{DN_max}{the maximum DN value of the image}

  \item{patch_long}{the size of block, to process whole ETM
  scene, set to 1000}
}
\description{
Calls GNSPI.pro IDL script by Zhu Xiaolin to fill gaps in
SLC-off Landsat 7 ETM+ image.
}
\references{
Zhu, X., Liu, D., Chen, J., 2012. A new geostatistical
approach for filling gaps in Landsat ETM+ SLC-off images.
Remote Sensing of Environment 124, 49-60.
}

