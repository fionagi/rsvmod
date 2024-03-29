% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{create_likelihood}
\alias{create_likelihood}
\title{Set up likelihood function in a closure environment for use with lazymcmc}
\usage{
create_likelihood(parTab, data, PRIOR_FUNC, ...)
}
\arguments{
\item{parTab}{parameter table setup for use in lazymcmc fitting functions}

\item{data}{observational incidence data}

\item{PRIOR_FUNC}{optional prior function}
}
\value{
function
}
\description{
Set up likelihood function in a closure environment for use with lazymcmc
}
