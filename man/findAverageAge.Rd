% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{findAverageAge}
\alias{findAverageAge}
\title{Find average age of first infection/hospitalisation, second+
infection/hospitalisation for all, term and preterm}
\usage{
findAverageAge(modOut, ageYears, startT, endT, risk = 1)
}
\arguments{
\item{modOut}{model output returned from aggregate_output function}

\item{ageYears}{lower bound of age groups in years}

\item{startT}{start t}

\item{endT}{end t}

\item{risk}{if 1, then risk model, 0 for base}
}
\value{
function
}
\description{
Find average age of first infection/hospitalisation, second+
infection/hospitalisation for all, term and preterm
}
