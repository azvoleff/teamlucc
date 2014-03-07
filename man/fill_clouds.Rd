\name{fill_clouds}
\alias{fill_clouds}
\title{Perform heavy cloud filling}
\usage{
fill_clouds(img_cloudy, img_clear, img_cloud_mask, out_name = NULL,
  fast = TRUE, num_class = 1, min_pixel = 20, extent1 = 1, DN_min = 0,
  DN_max = 255, patch_long = 1000,
  idl = "C:/Program Files/Exelis/IDL83/bin/bin.x86_64/idl.exe")
}
\arguments{
  \item{img_cloudy}{the cloudy image (base image)}

  \item{img_clear}{the clear image to use for filling
  \code{img_cloudy}}

  \item{img_cloud_mask}{cloud mask, with cloud patches
  assigned unique integer codes, and all background coded
  as zero. The \code{ConnCompLabel} function in the
  \code{SDMTools} package can automate this process.}

  \item{out_name}{name for output image, or NULL. If null
  an output name will be automatically assigned based on
  \code{img_clear}}

  \item{fast}{if \code{TRUE}, use the CLOUD_REMOVE_FAST.pro
  script. If \code{FALSE}, use the CLOUD_REMOVE.pro
  script.}

  \item{num_class}{set the estimated number of classes in
  image}

  \item{min_pixel}{the sample size of similar pixels}

  \item{extent1}{set the range of cloud neighborhood}

  \item{DN_min}{minimum of DN}

  \item{DN_max}{maximum of DN}

  \item{patch_long}{block size when process large image}

  \item{idl}{path to the IDL binary}
}
\description{
Calls either of the CLOUD_REMOVE.pro or
CLOUD_REMOVE_FAST.pro IDL scripts by Xiaolin Zhu to fill
heavy clouds in a Landsat image.
}
\references{
Zhu, X., Gao, F., Liu, D., Chen, J., 2012. A modified
neighborhood similar pixel interpolator approach for
removing thick clouds in Landsat images. Geoscience and
Remote Sensing Letters, IEEE 9, 521-525.
}

