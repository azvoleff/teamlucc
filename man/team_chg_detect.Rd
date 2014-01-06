\name{team_chg_detect}
\alias{team_chg_detect}
\title{Perform change detection for two Landsat CDR surface reflectance images}
\usage{
team_chg_detect(t1_predictions, t2_predictions)
}
\arguments{
  \item{t1_predictions}{output from
  \code{team_classify_image} for time 1 image}

  \item{t2_predictions}{output from
  \code{team_classify_image} for time 2 image}
}
\description{
First the images should be classified using the
\code{team_classify} function.
}

