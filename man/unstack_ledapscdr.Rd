\name{unstack_ledapscdr}
\alias{unstack_ledapscdr}
\title{Convert LEDAPS CDR Landsat format images from HDF4 to ENVI format}
\usage{
unstack_ledapscdr(x, out, overwrite = FALSE)
}
\arguments{
  \item{x}{input HDF4 file}

  \item{out}{output file}

  \item{overwrite}{whether to overwrite existing files}
}
\description{
This function converts a Landsat surface reflectance (SR)
image from the Landsat Climate Data Record (CDR) archive
into a series of single band images in ENVI format.
}

