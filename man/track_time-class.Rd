\docType{class}
\name{Track_time-class}
\alias{Track_time-class}
\title{A class for tracking running time of individual sections of an R script}
\description{
A class for tracking running time of individual sections of
an R script
}
\section{Slots}{

\describe{
\item{\code{timers}}{a \code{data.frame} tracking timer names and start times}

\item{\code{notify}}{function to use for outputting timers (defaults to
\code{\link{print}}}
}}

