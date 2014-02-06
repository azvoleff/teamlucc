\name{fill_gaps}
\alias{fill_gaps}
\title{Perform SLC-off gap fill of Landsat 7 ETM+ image}
\usage{
fill_gaps(slc_off_file, input_file, timeseries_files,
  script_path = "GNSPI.sav", out_base = NULL, sample_size = 20,
  size_wind = 12, class_num = 4, DN_min = 0, DN_max = 1,
  patch_long = 1000,
  idl = "C:/Program Files/Exelis/IDL83/bin/bin.x86_64/idl.exe")
}
\arguments{
  \item{slc_off_file}{the SLC-off Landsat 7 file to gap
  fill}

  \item{input_file}{the first file to use to fill in the
  gaps}

  \item{timeseries_files}{a timeseries of files to use as
  additional inputs to the gap fill algorithm}

  \item{script_path}{the path to the GNSPI.sav file}

  \item{out_base}{path and base filename for the output
  file. The script will save the output files by appending
  "_GNSPI.envi" and "_GNSPI_uncertainty.envi" to this base
  filename.}

  \item{sample_size}{the sample size of sample pixels}

  \item{size_wind}{the maximum window size}

  \item{class_num}{the estimated number of classes}

  \item{DN_min}{the minimum DN value of the image}

  \item{DN_max}{the maximum DN value of the image}

  \item{patch_long}{the size of block, to process whole ETM
  scene, set to 1000}

  \item{idl}{path to the IDL binary}
}
\description{
Calls GNSPI.pro IDL script by Zhu Xiaolin to fill gaps in
SLC-off Landsat 7 ETM+ image. The GNSPI script should be
provided as a .sav file.
}
\references{
Zhu, X., Liu, D., Chen, J., 2012. A new geostatistical
approach for filling gaps in Landsat ETM+ SLC-off images.
Remote Sensing of Environment 124, 49-60.
}

