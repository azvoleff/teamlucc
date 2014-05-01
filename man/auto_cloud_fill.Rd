\name{auto_cloud_fill}
\alias{auto_cloud_fill}
\title{Automated removal of clouds using Xiaolin Zhu's NSPI algorithm}
\usage{
auto_cloud_fill(data_dir, wrspath, wrsrow, start_date, end_date,
  base_date = NULL, tc = TRUE, threshold = 1, max_iter = 5,
  n_cpus = 1, notify = print, verbose = TRUE, ...)
}
\arguments{
  \item{data_dir}{folder where input images are located,
  with filenames as output by the
  \code{\link{auto_preprocess_landsat}} function. This
  folder will be searched recursively for images (taking
  the below path/row, date, and topographic correction
  options into account).}

  \item{wrspath}{World Reference System (WRS) path}

  \item{wrsrow}{World Reference System (WRS) row}

  \item{start_date}{start date of period from which images
  will be chosen to fill cloudy areas in the base image (as
  \code{Date} object)}

  \item{end_date}{end date of period from which images will
  be chosen to fill cloudy areas in the the base image (as
  \code{Date} object)}

  \item{base_date}{ideal date for base image (base image
  will be chosen as the image among the available images
  that is closest to this date). If NULL, then the base
  image will be the image with the lowest cloud cover.}

  \item{tc}{if \code{TRUE}, use topographically corrected
  imagery as output by \code{auto_preprocess_landsat}. IF
  \code{FALSE} use bands 1-5 and 7 surface reflectance as
  output by \code{unstack_ledaps} or
  \code{auto_preprocess_landsat} (if
  \code{auto_preprocess_landsat} was also run with
  tc=FALSE).}

  \item{threshold}{maximum percent cloud cover allowable in
  base image (cloud fill iterate until percent cloud cover
  in base image is below this value) or until more than
  \code{max_iter} iterations have been run}

  \item{max_iter}{maximum number of times to run cloud fill
  script}

  \item{n_cpus}{the number of CPUs to use for processes
  that can run in parallel}

  \item{notify}{notifier to use (defaults to \code{print}
  function).  See the \code{notifyR} package for one way of
  sending notifications from R.  The \code{notify} function
  should accept a string as the only argument.}

  \item{verbose}{whether to print detailed status messages}

  \item{...}{additional arguments passed to
  \code{\link{cloud_remove}}, such as \code{DN_min},
  \code{DN_max}, \code{use_IDL}, \code{fast},
  \code{verbose}, etc. See \code{\link{cloud_remove}}.}
}
\value{
\code{Raster*} object with cloud filled image.
}
\description{
Uses the NSPI algorithm from Zhu et al. See
\code{\link{cloud_remove}} for details. In hilly areas,
cloud fill should be done after topographic correction.
}
\details{
The \code{auto_cloud_fill} function allows an analyst to
automatically construct a cloud-filled image after
specifying: \code{data_dir} (a folder of Landsat images),
\code{wrspath} and \code{wrsrow} (the WRS-2 path/row to
use), and \code{start_date} and \code{end_date} (a start
and end date limiting the images to use in the algorithm).
The analyst can also optionally specify a \code{base_date},
and the \code{auto_cloud_fill} function will automatically
pick the image closest to that date to use as the base
image.

As the \code{auto_cloud_fill} function automatically
chooses images for inclusion in the cloud fill process, it
relies on having images stored on disk in a particular way,
and currently only supports cloud fill for Landsat CDR
surface reflectance images. To ensure that images are
correctly stored on your hard disk, use the
\code{\link{auto_preprocess_landsat}} function to extract
the original Landsat CDR hdf files from the USGS archive.
The \code{auto_preprocess_landsat} function will ensure
that images are extracted and renamed properly so that they
can be used with the \code{auto_cloud_fill} script.
}
\references{
Zhu, X., Gao, F., Liu, D., Chen, J., 2012. A modified
neighborhood similar pixel interpolator approach for
removing thick clouds in Landsat images.  Geoscience and
Remote Sensing Letters, IEEE 9, 521--525.
doi:10.1109/LGRS.2011.2173290
}

