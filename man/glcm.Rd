\name{glcm}
\alias{glcm}
\title{Image texture measures from grey-level co-occurrence matrices (GLCM)}
\usage{
  glcm(layer, n_grey = 32, window = c(3, 3),
    shift = c(1, 1),
    statistics = c("mean", "variance", "homogeneity", "contrast", "dissimilarity", "entropy", "second_moment", "correlation"),
    ...)
}
\arguments{
  \item{layer}{a /code{RasterLayer}}

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

  \item{...}{additional parameters to pass to rasterEngine}
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
L5TSR_1986 <- stack(system.file('extdata/L5TSR_1986.dat', package='teamr'))
textures <- glcm(raster(L5TSR_1986, layer=1))
plot(textures)
}
\references{
  Lu, D., and M. Batistella. 2005. Exploring TM image
  texture and its relationships with biomass estimation in
  Rond\^{o}nia, Brazilian Amazon. Acta Amazonica
  35:249-257.

  Gonzalez, R. C. 2008. Digital image processing. 3rd ed.
  Prentice Hall, Upper Saddle River, N.J, pages 830-836.

  Haralick, R. M., K. Shanmugam, and I. Dinstein. 1973.
  Textural features for image classification. IEEE
  Transactions on Systems, Man and Cybernetics
  SMC-3:610-621.

  Pratt, W. K. 2007. Digital image processing: PIKS
  Scientific inside. 4th ed. Wiley-Interscience, Hoboken,
  N.J pages 540-541, 563-566.
}

