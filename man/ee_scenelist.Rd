\name{ee_scenelist}
\alias{ee_scenelist}
\title{Save scenelist from EarthExplorer metadata for upload to ESPA}
\usage{
ee_scenelist(x, start_date, out_file, num_months = 12, min_clear = 0.7,
  exclude = list())
}
\arguments{
  \item{x}{a \code{data.frame} with a list of Landsat
  scenes as output from the save metadata function on
  http://earthexplorer.usgs.gov}

  \item{start_date}{starting date as a \code{Date} object}

  \item{out_file}{filename for output text file for later
  upload to ESPA}

  \item{num_months}{number of months to include}

  \item{min_clear}{the minimum percent clear to plot
  (calculated as 1 - percent cloud cover). Images with less
  than \code{min_clear} fraction of the image area free of
  clouds will be ignored.}

  \item{exclude}{a list of sensors to exclude (for example,
  set \code{exclude=c('LE7', 'LT4')} to exclude Landsat 7
  ETM+ and Landsat 4 TM images.}
}
\value{
used for side effect of producing ESPA scene list
}
\description{
Save scenelist from EarthExplorer metadata for upload to
ESPA
}

