\name{unstack_ledapscdr}
\alias{unstack_ledapscdr}
\title{Convert LEDAPS CDR Landsat format images from HDF4 to ENVI format}
\usage{
unstack_ledapscdr(x, output_folder = NULL, overwrite = FALSE,
  rmhdf = FALSE)
}
\arguments{
  \item{x}{input HDF4 file}

  \item{out_folder}{output folder (if \code{NULL}, defaults
  to input folder)}

  \item{overwrite}{whether to overwrite existing files}

  \item{rmhdf}{whether to remove hdf files after unstacking
  them}
}
\description{
This function converts a Landsat surface reflectance (SR)
image from the Landsat Climate Data Record (CDR) archive
into a series of single band images in ENVI format.
}
\details{
Note that this function requires a local GDAL installation.
See http://trac.osgeo.org/gdal/wiki/DownloadingGdalBinaries
or http://trac.osgeo.org/osgeo4w/ to download the
appropriate installer for your operating system.
}

