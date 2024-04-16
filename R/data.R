# data.R - This file includes documentation of the data used.

#' ABS population of Greater Perth by 5-year age group - 2011
#'
#' Data is formatted to be easily used with the R package conmat
#' to find the applicable contact matrix
#'
#' @format A data table containing 18 observations of 4 variables.
#' \describe{
#'   \item{lga}{Local Government Area classification}
#'   \item{lower.age.limit}{}
#'   \item{year}{census data year}
#'   \item{population}{population in age group}
#' }
#' @source \url{https://abs.gov.au/census/find-census-data/quickstats/2021/5GPER}
"data.greaterPerth.2011"

#' Monthly contact matrix
#'
#' Monthly number of contacts between five-year age groups within the population.
#' This contact matrix was generated using the R package conmat based on the
#' POLYMOD study results for the UK, normalised to 2011 metropolitan Perth
#' population demographics. The matrix was made symmetric by averaging the
#' non-diagonal values of the matrix. The <5 year age group is divided into
#' monthly age groups in the model, with contacts being uniformly distributed
#' into these resulting age groups.
#'
#' @format A 75 x 75 data matrix
#' @source \url{https://rdrr.io/github/njtierney/conmat/}
"data.contact.monthly"
