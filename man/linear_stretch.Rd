\name{linear_stretch}
\alias{linear_stretch}
\title{Apply a linear stretch to an image}
\usage{
linear_stretch(x, pct = 2, max_val = 1)
}
\arguments{
  \item{x}{image to stretch}

  \item{pct}{percent stretch}

  \item{max_val}{maximum value of final output (image will
  be rescaled to range from 0 - \code{max_val})}
}
\value{
image with stretch applied
}
\description{
Applies a linear stretch to an image (default linear 2%
stretch), and returns the image with each band individually
stretched and rescaled to range between zero and
\code{max_val} (default of \code{max_val} is 1).
}

