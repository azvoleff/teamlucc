\name{cloud_remove}
\alias{cloud_remove}
\title{R/C++ implementation of Xiaolin Zhu's CLOUD_REMOVE IDL script}
\usage{
cloud_remove(cloudy, clear, cloud_mask, out_name, num_class = 1,
  min_pixel = 20, cloud_nbh = 1, DN_min = 0, DN_max = 255, n_cpus = 1,
  ...)
}
\arguments{
  \item{cloudy}{the cloudy image (base image) as a
  \code{Raster*}}

  \item{clear}{the clear image as a \code{Raster*} to use
  for filling \code{img_cloudy}}

  \item{cloud_mask}{the cloud mask as a \code{RasterLayer},
  with each cloud patch assigned a unique integer code.
  Areas that are clear in both \code{cloudy_rast} and
  \code{clear_rast} should be coded 0, while areas that are
  clouded in \code{clear_rast} should be coded -1.}

  \item{out_name}{filename for cloud filled image}

  \item{num_class}{set the estimated number of classes in
  image}

  \item{min_pixel}{the sample size of similar pixels}

  \item{cloud_nbh}{the range of cloud neighborhood (in
  pixels)}

  \item{DN_min}{the minimum valid DN value}

  \item{DN_max}{the maximum valid DN value}

  \item{n_cpus}{the number of CPUs to use for processes
  that can run in parallel}

  \item{...}{additional arguments to pass to rasterEngine}
}
\description{
R/C++ implementation of Xiaolin Zhu's CLOUD_REMOVE IDL
script
}
\references{
Zhu, X., Gao, F., Liu, D., Chen, J., 2012. A modified
neighborhood similar pixel interpolator approach for
removing thick clouds in Landsat images. Geoscience and
Remote Sensing Letters, IEEE 9, 521--525.
}

