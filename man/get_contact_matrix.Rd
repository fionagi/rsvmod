% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{get_contact_matrix}
\alias{get_contact_matrix}
\title{Get contact matrix
Using R package conmat, population for greater Perth and assumed
age groupings for model}
\usage{
get_contact_matrix(year, sym = 0, save = 1, fileName = NA)
}
\arguments{
\item{year}{year of ABS population data data}

\item{sym}{return a symmetric matrix}

\item{save}{1 to save matrix, 0 to return}

\item{fileName}{if save, give file name}
}
\value{
matrix
}
\description{
Get contact matrix
Using R package conmat, population for greater Perth and assumed
age groupings for model
}
