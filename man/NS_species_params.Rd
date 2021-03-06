% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_therMizer.R
\docType{data}
\name{NS_species_params}
\alias{NS_species_params}
\title{Example parameter set based on the North Sea}
\format{A data frame with 12 rows and 7 columns. Each row is a species.
\describe{
\item{species}{Name of the species}
\item{w_inf}{The von Bertalanffy Winfinity parameter}
\item{w_mat}{Size at maturity}
\item{beta}{Size preference ratio}
\item{sigma}{Width of the size-preference}
\item{r_max}{Maxium recruitment}
\item{k_vb}{The von Bertalanffy k parameter}
}}
\source{
{Blanchard et al.}
}
\usage{
NS_species_params
}
\description{
This data set is based on species in the North Sea (Blanchard et al.).
It is a data.frame that contains all the necessary information to be used by the
\code{\link{therMizerParams}} constructor.
As there is no gear column, each species is asssumed to be fished by a seperate gear.
}
\keyword{datasets}
