#' therMizer: mizer with the effects of temperature on metabolism and encounter rate
#' Mizer: Multi-species size-based modelling in R
#'
#' The mizer package implements multi-species size-based modelling in R. It has 
#' been designed for modelling marine ecosystems.
#'
#' therMizer expands on the mizer package.
#'
#' Using \pkg{therMizer} is relatively simple. There are three main stages: 
#' \enumerate{
#'
#' \item Setting the model parameters. This is done by creating an object of 
#' class \code{therMizerParams}. This includes model parameters such as the life 
#' history parameters of each species, and the range of the size spectrum.
#'
#' \item Running a simulation. This is done by calling the \code{project_therMizer()} 
#' method on the model parameters. This produces an object of \code{therMizerSim} 
#' which contains the results of the simulation.
#'
#' \item Exploring results. After a simulation has been run, the results can be 
#' explored using a range of plots and summaries. 
#' }
#'
#' See the accompanying vignette for full details of the principles behind mizer
#' and how the package can be used to peform size-based modelling.
#'
#' See Woodworth-Jefcoats et al. in prep for details on therMizer.
#'
#' @import plyr ggplot2 grid methods
#' @importFrom reshape2 melt
#'
#' @docType package
#' @name therMizer
#' @aliases therMizer, therMizer-package
NULL