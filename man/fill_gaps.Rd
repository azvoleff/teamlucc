\name{fill_gaps}
\alias{fill_gaps}
\title{Perform SLC-off gap fill of Landsat 7 ETM+ image}
\usage{
fill_gaps(slc_off, fill, timeseries = c(), out_base = NULL,
  use_IDL = TRUE, sample_size = 20, size_wind = 12, class_num = 4,
  DN_min = 0, DN_max = 1, patch_long = 1000,
  idl = "C:/Program Files/Exelis/IDL83/bin/bin.x86_64/idl.exe",
  verbose = FALSE)
}
\arguments{
  \item{slc_off}{the SLC-off Landsat 7 file to gap fill, as
  a \code{Raster*}}

  \item{fill}{the first TM image to use to fill in the
  gaps, as a \code{Raster*}}

  \item{timeseries}{a timeseries of TM images as
  \code{Raster*} objects to use as additional inputs to the
  gap fill algorithm (optional)}

  \item{out_base}{path and base filename for the output
  file. The script will save the output files by appending
  "_GNSPI.envi" and "_GNSPI_uncertainty.envi" to this base
  filename.}

  \item{use_IDL}{whether to use Xiaolin Zhu's original IDL
  scripts (\code{use_IDL=TRUE}) or the C++ implementation
  native to \code{teamlucc} (\code{use_IDL=FALSE})}

  \item{sample_size}{the sample size of sample pixels}

  \item{size_wind}{the maximum window size}

  \item{class_num}{the estimated number of classes}

  \item{DN_min}{the minimum DN value of the image}

  \item{DN_max}{the maximum DN value of the image}

  \item{patch_long}{the size of block, to process whole ETM
  scene, set to 1000}

  \item{idl}{path to the IDL binary}

  \item{verbose}{whether to print detailed status messages}
}
\value{
a list of two rasters: 1) "filled", the gap filled image,
and 2) "uncertainty", the uncertainty image.
}
\description{
Calls GNSPI.pro IDL script by Xiaolin Zhu to fill gaps in
SLC-off Landsat 7 ETM+ image. The script requires
\code{fill} to be a TM (Landsat 5) image. \code{slc_off}
must be a Landsat 7 SLC-off image.
}
\details{
If supplied, \code{timeseries} should be a list of TM
images.

Pixels in gaps, background, and/or clouds in
\code{slc_off}, \code{input_image}, and the images in
\code{timeseries} should be coded as 0.
}
\references{
Zhu, X., Liu, D., Chen, J., 2012. A new geostatistical
approach for filling gaps in Landsat ETM+ SLC-off images.
Remote Sensing of Environment 124, 49--60.
}

