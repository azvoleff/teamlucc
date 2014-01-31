\name{accuracy}
\alias{accuracy}
\title{Calculate statistics summarizing classification accuracy}
\usage{
accuracy(model, test_data = NULL, pop = NULL, reclass_mat = NULL)
}
\arguments{
  \item{model}{a classification model with a \code{predict}
  method}

  \item{test_data}{a \code{link{Training_data}} object, or
  NULL. If test_data is NULL, \code{accuracy} will try to
  use the training data included in the \code{model} object
  (this will only work if the model object is of class
  \code{train} from the \code{caret} package. If test_data
  is a \code{Training_data} object, \code{accuracy} will
  use the included \code{training_flag} indicator to
  separate testing and training data.}

  \item{pop}{(optional) used to convert from sample matrix
  to population matrix as in Pontius and Millones 2011.Can
  be: 1) NULL, in which case the sample frequencies will be
  used as estimates of the population frequencies, 2) a
  list of length equal to the number of classes in the map
  giving the total number of pixels in the entire
  population for each class, or 3) the predicted cover map
  from \code{model} as a \code{RasterLayer}, from which the
  population frequencies will be tabulated.}

  \item{reclass_mat}{a reclassification matrix to be used
  in the case of a model fit by \code{classify_image} with
  the \code{do_split} option selected}
}
\value{
\code{\link{accuracy-class}} instance
}
\description{
Calculates a contingency table and various statistics for
use in image classification accuracy assessment and map
comparison. Contingency table includes user's, producer's,
and overall accuracies for an image classification, and
quantity disagreement \code{Q} and allocation disagreement
\code{A}. Q and A are calculated based on Pontius and
Millones (2011). 95% confidence intervals for the user's,
producer's, and overall accuracies are calculated as in
Olofsson et al. (2013).
}
\details{
To avoid bias due to the use of a sample contingency table,
the contingency table can be converted to a population
contingency table, if the variable 'pop' is provided. For
an accuracy assessment based on a random sample, 'pop' does
not need to be provided.
}
\examples{
accuracy(classified_LT5SR_1986$model)
}
\references{
Pontius, R. G., and M. Millones. 2011. Death to Kappa:
birth of quantity disagreement and allocation disagreement
for accuracy assessment. International Journal of Remote
Sensing 32:4407-4429. Olofsson, P., G. M. Foody, S. V.
Stehman, and C. E. Woodcock.  2013. Making better use of
accuracy data in land change studies: Estimating accuracy
and area and quantifying uncertainty using stratified
estimation.  Remote Sensing of Environment 129:122-131.
}

