\name{cloud_remove}
\alias{cloud_remove}
\title{Remove clouds using Xiaolin Zhu's NSPI algorithm}
\usage{
cloud_remove(cloudy, clear, cloud_mask, out_name, use_IDL = TRUE,
  fast = TRUE, num_class = 1, min_pixel = 20, max_pixel = 1000,
  cloud_nbh = 1, DN_min = 0, DN_max = 255,
  idl = "C:/Program Files/Exelis/IDL83/bin/bin.x86_64/idl.exe")
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

  \item{use_IDL}{whether to use Xiaolin Zhu's original IDL
  scripts (\code{use_IDL=TRUE}) or the C++ implementation
  native to \code{teamlucc} (\code{use_IDL=FALSE})}

  \item{fast}{if \code{TRUE}, use the CLOUD_REMOVE_FAST.pro
  script. If \code{FALSE}, use the CLOUD_REMOVE.pro
  script.}

  \item{num_class}{set the estimated number of classes in
  image}

  \item{min_pixel}{the sample size of similar pixels}

  \item{max_pixel}{the maximum sample size to search for
  similar pixels}

  \item{cloud_nbh}{the range of cloud neighborhood (in
  pixels)}

  \item{DN_min}{the minimum valid DN value}

  \item{DN_max}{the maximum valid DN value}

  \item{idl}{path to the IDL binary on your machine (on
  Windows, the path to idl.exe)}

  \item{...}{additional arguments to pass to
  \code{rasterEngine}}
}
\value{
\code{Raster*} with cloud-filled image
}
\description{
This code uses the Neighborhood Similar Pixel Interpolator
(NSPI) algorithm by Xiaolin Zhu to fill heavy clouds in a
Landsat image.
}
\details{
This code can use either a Xiaolin Zhu's original IDL code,
or an R/C++ implementation native to \code{teamlucc}. The
\code{use_IDL} parameter (defaults to TRUE) decides whether
to use Xiaolin's code \code{use_IDL=TRUE) or the or the R
and C++ implementation native to \code{teamlucc}
(\code{use_IDL=FALSE}).

The results from running the two alternative versions of
the cloud removal algorithm should be identical. However,
there is one important difference to note: the R/C++
implementation of cloud removal does not currently
implement the \code{fast=TRUE} option.
}
\references{
Zhu, X., Gao, F., Liu, D., Chen, J., 2012. A modified
neighborhood similar pixel interpolator approach for
removing thick clouds in Landsat images. Geoscience and
Remote Sensing Letters, IEEE 9, 521--525.
}

