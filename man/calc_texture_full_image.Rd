\name{calc_texture_full_image}
\alias{calc_texture_full_image}
\title{Calculates a glcm texture for use in the glcm.R script}
\usage{
calc_texture_full_image(rast, n_grey, window_dims, shift, statistics)
}
\arguments{
  \item{rast}{a matrix containing the pixels to be used in
  the texture calculation}

  \item{n_grey}{number of grey levels to use in texture
  calculation}

  \item{window_dims}{2 element list with row and column
  dimensions of the texture window}

  \item{shift}{a length 2 vector with the number of cells
  to shift when computing co-ocurrency matrices}

  \item{statistics}{a list of strings naming the texture
  statistics to calculate}
}
\value{
a list of length equal to the length of the
\code{statistics} input parameter, containing the selected
textures measures
}
\description{
This function is called by the \code{\link{glcm}} function.
It is not intended to be used directly.
}

