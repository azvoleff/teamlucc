\name{tpa2df}
\alias{tpa2df}
\title{Function to convert TIMESAT .tpa binary format file to an R dataframe.}
\usage{
tpa2df(x)
}
\arguments{
  \item{x}{A string giving the location of a .tpa file
  output by TIMESAT}
}
\value{
A data.frame containing 14 columns: row, col, season,
start, end, length, base_value, peak_time, peak_value, amp,
left_deriv, right_deriv, large_integ, and small_integ
}
\description{
Function to convert TIMESAT .tpa binary format file to an R
dataframe.
}
\examples{
# TODO: Need to add examples here, and need to include a sample TIMESAT tpa
# file in the package data.
}

