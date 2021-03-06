% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/project_methods_therMizer.R
\docType{methods}
\name{getEReproAndGrowth}
\alias{getEReproAndGrowth}
\alias{getEReproAndGrowth,therMizerParams,matrix,numeric,numeric,matrix-method}
\title{getEReproAndGrowth method for the size based model}
\usage{
getEReproAndGrowth(object, n, n_pp, ocean_temp, feeding_level)


  \S4method{getEReproAndGrowth}{therMizerParams,matrix,numeric,numeric,matrix}(object,
  n, n_pp, ocean_temp, feeding_level)
}
\arguments{
\item{object}{A \code{therMizerParams} object.}

\item{n}{A matrix of species abundance (species x size).}

\item{n_pp}{A vector of the background abundance by size.}

\item{ocean_temp}{A numeric vector of the ocean temperature by realm or a 
single numeric ocean_temp value which is used for all realms.}

\item{feeding_level}{The current feeding level (optional). A matrix of size
no. species x no. size bins. If not supplied, is calculated internally
using the \code{getFeedingLevel()} method.}
}
\value{
A two dimensional array (species x size)
}
\description{
Calculates the energy available by species and size for reproduction and
growth after metabolism and movement have been accounted for. Used by the
\code{project_therMizer} method for performing simulations.

\code{getEReproAndGrowth} method with \code{feeding_level} argument.
}
\examples{
\dontrun{
data(NS_species_params_gears)
data(inter)
params <- therMizerParams(NS_species_params_gears, inter)
# Project with constant fishing effort for all gears for 20 time steps
sim <- project_therMizer(params, t_max = 20, effort = 0.5)
# Get the energy at a particular time step
getEReproAndGrowth(params,sim@n[21,,],sim@n_pp[21,])
}
}
\seealso{
\code{\link{project_therMizer}} and \code{\link{getFeedingLevel}}.
}
