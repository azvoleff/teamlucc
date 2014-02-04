\name{espa_download}
\alias{espa_download}
\title{Download a completed ESPA order}
\usage{
espa_download(email, order_ID, output_folder)
}
\arguments{
  \item{email}{address used to place the order}

  \item{order_ID}{the ESPA order ID}

  \item{output_folder}{the folder to save output data in}
}
\value{
used for the side effect of downloading Landsat scenes
}
\description{
Function to download a set of Landsat images from ESPA
given a valid order ID and the email that placed the ESPA
order.
}

