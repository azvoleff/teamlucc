\name{fill_clouds}
\alias{fill_clouds}
\title{Perform heavy cloud filling}
\usage{
fill_clouds(img_cloudy, img_clear, img_cloud_mask, out_name = NULL,
  script_path = "CLOUD_REMOVE.sav", num_class = 1, extent1 = 1,
  DN_min = 0, DN_max = 255, patch_long = 1000,
  idl = "C:/Program Files/Exelis/IDL83/bin/bin.x86_64/idl.exe")
}
\arguments{
  \item{img_cloudy}{the clear image (base image)}

  \item{img_clear}{the cloudy image}

  \item{img_cloud_mask}{cloud mask, with clouded areas set
  to}

  \item{out_name}{name for output image, or NULL. If null
  an output name will be automatically assigned based on
  \code{img_clear}}

  \item{script_path}{the path to the CLOUD_REMOVE.sav file}

  \item{num_class}{set the estimated number of classes in
  image}

  \item{extent1}{set the range of cloud neighborhood}

  \item{DN_min}{minimum of DN}

  \item{DN_max}{maximum of DN}

  \item{patch_long}{block size when process large image}

  \item{idl}{path to the IDL binary}
}
\description{
Calls CLOUD_REMOVE.pro IDL script by Zhu Xiaolin to fill
heavy clouds in Landsat image. The CLOUD_REMOVE script
should be provided as a .sav file.
}
\references{
Zhu, X., Gao, F., Liu, D., Chen, J., 2012. A modified
neighborhood similar pixel interpolator approach for
removing thick clouds in Landsat images. Geoscience and
Remote Sensing Letters, IEEE 9, 521-525.
}

