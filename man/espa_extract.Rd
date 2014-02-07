\name{espa_extract}
\alias{espa_extract}
\title{Extract a set of Landsat tarballs into a folder tree organized by image date}
\usage{
espa_extract(in_folder, out_folder)
}
\arguments{
  \item{in_folder}{Path to a folder of .tar.gz Landsat
  surface reflectance images}

  \item{out_folder}{output folder}
}
\value{
used for side effect of unzipping Landsat tarballs
}
\description{
Each image .tar.gz file will be extracted into a subfolder
within \code{output_folder}. The subfolders will be named
according to the year and Julian date of capture, and the
sensor type (LT4, LT5 or LE7 for Landsat 4 TM, Landsat 5
TM, and Landsat 7 ETM+ respectively). For example,
"LT50150531986037-SC20130816144215.tar.gz" would be
extracted into a subfolder named "1986_037_LT5', for 1986,
Julian day 37, and Landsat 5 TM.
}

