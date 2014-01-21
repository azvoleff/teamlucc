\name{glcm}
\alias{glcm}
\title{Image texture measures from grey-level co-occurrence matrices (GLCM)}
\usage{
glcm(x, n_grey = 32, window = c(3, 3), shift = c(1, 1), statistics =
c("mean", "variance", "homogeneity", "contrast", "dissimilarity", "entropy",
"second_moment", "correlation"),
min_x=NULL, max_x=NULL, na_opt="any", na_val=NA, scale_factor=1,
asinteger=FALSE)
}
\arguments{
  \item{x}{a /code{RasterLayer} or matrix}

  \item{n_grey}{number of grey levels to use in texture
  calculation}

  \item{window}{the window size to consider for texture
  calculation as a two element integer vector}

  \item{shift}{a 2 element integer vector giving the shift
  (Q in Gonzalez and Woods, 2008).}

  \item{statistics}{A list of GLCM texture measures to
  calculate.  Can include any (one or more) of the
  following: 'mean', 'mean_ENVI', 'variance',
  'variance_ENVI', 'homogeneity', 'contrast',
  'dissimilarity', 'entropy', 'second_moment', and/or
  'correlation'.}

  \item{min_x}{minimum value of input \code{RasterLayer}
  (optional, \code{glcm} will calculate if not supplied).
  Useful when running \code{glcm} over blocks of a raster.}

  \item{max_x}{maximum value of input \code{RasterLayer}
  (optional, \code{glcm} will calculate if not supplied).
  Useful when running \code{glcm} over blocks of a raster.}

  \item{na_opt}{How to handle NA values in \code{x}. Can be
  set to "ignore", "any" or "center". If set to "any", all
  textures statistics for a given pixel will be set to NA
  if there are any NA values in the \code{window} around
  that pixel. If set to "center" this will only occur if
  the center value is an NA. If set to "ignore", NA values
  in \code{window} will be ignored.}

  \item{na_val}{the value to use to fill NA values on edges
  of \code{x} where textures cannot be calculated due to
  the window falling outside of the image, and as necessary
  depending on the chosen \code{na_opt}.}

  \item{scale_factor}{factor by which to multiply results.
  Useful if rounding results to integers (see
  \code{asinteger} argument).}

  \item{asinteger}{whether to round results to nearest
  integer. Can be used to save space by saving results as,
  for example, an 'INT2S' \code{raster}.}
}
\value{
A /code{RasterLayer} with the calculated requested GLCM
texture measures.
}
\description{
Image texture measures from grey-level co-occurrence
matrices (GLCM)
}
\examples{
\dontrun{
textures <- glcm(raster(L5TSR_1986, layer=1))
plot(textures)
}
}
\references{
Lu, D., and M. Batistella. 2005. Exploring TM image texture
and its relationships with biomass estimation in
Rond\^{o}nia, Brazilian Amazon. Acta Amazonica 35:249-257.

Gonzalez, R. C. 2008. Digital image processing. 3rd ed.
Prentice Hall, Upper Saddle River, N.J, pages 830-836.

Haralick, R. M., K. Shanmugam, and I. Dinstein. 1973.
Textural features for image classification. IEEE
Transactions on Systems, Man and Cybernetics SMC-3:610-621.

Pratt, W. K. 2007. Digital image processing: PIKS
Scientific inside. 4th ed. Wiley-Interscience, Hoboken, N.J
pages 540-541, 563-566.
}

