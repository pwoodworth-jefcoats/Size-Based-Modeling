# Class specification and constructors for the simulation class

# Copyright 2012 Finlay Scott and Julia Blanchard. 
# Distributed under the GPL 2 or later
# Maintainer of mizer package: Finlay Scott, CEFAS. finlay.scott@cefas.co.uk

# Modified for therMizer by Phoebe.Woodworth-Jefcoats@noaa.gov
# Added ability to include:
# ocean temperature
# plankton densities input at each time step to set background resource

# Hackiness to get past the 'no visible binding ... ' warning when running check
utils::globalVariables(c("lm"))

# Validity check
valid_therMizerSim <- function(object){
	errors <- character()
	validObject(object@params)
	# array dimensions
	if(length(dim(object@n)) != 3){
		msg <- "n slot must have three dimensions"
		errors <- c(errors, msg)
	}
	if(length(dim(object@effort)) != 2){
		msg <- "effort slot must have two dimensions"
		errors <- c(errors, msg)
	}
	if(length(dim(object@ocean_temp)) != 2){ 
		msg <- "ocean_temp slot must have two dimension" 
		errors <- c(errors, msg) 
	} 
	if(length(dim(object@n_pp)) != 2){
		msg <- "n_pp slot must have two dimensions"
		errors <- c(errors, msg)
	}
	if(length(dim(object@n_plank)) != 2){ 
		msg <- "n_plank slot must have two dimensions" 
		errors <- c(errors, msg) 
	} 
	# Check time dimension is good - size, dim name, and names
	if(!all(c(dim(object@n)[1]-1,dim(object@n_pp)[1]-1,dim(object@ocean_temp)[1],dim(object@n_plank[1])) == dim(object@effort)[1])){ 
		msg <- "First dimension of effort, ocean_temp, n_plank, n and n_pp slots must be the same length (n and n_pp are + 1 longer)" 
		errors <- c(errors, msg)
	}
	if(!all(c(names(dimnames(object@n))[1], names(dimnames(object@n_pp))[1], names(dimnames(object@effort))[1], names(dimnames(object@ocean_temp))[1], names(dimnames(object@n_plank))[1]) == "time")){ 
		msg <- "First dimension of effort, ocean_temp, n, n_plank, and n_pp slots must be called 'time'" 
		errors <- c(errors, msg)
	}
	if(!all(c(names(dimnames(object@n))[1], names(dimnames(object@n_pp))[1], names(dimnames(object@ocean_temp))[1], names(dimnames(object@n_plank))[1]) == 
		names(dimnames(object@effort))[1])){
		msg <- "First dimension of effort, ocean_temp, n, and n_pp slots must have the same names" 
		errors <- c(errors, msg)
	}
	# species dimension of n
	if(dim(object@n)[2] != dim(object@params@psi)[1]){
		msg <- "Second dimension of n slot must have same length as the species names in the params slot"
		errors <- c(errors, msg)
	}
	if(names(dimnames(object@n))[2] != "sp"){
		msg <- "Second dimension of n slot must be called 'sp'"
		errors <- c(errors, msg)
	}
	if(!all(names(dimnames(object@n))[2] == names(dimnames(object@params@psi))[1])){
		msg <- "Second dimension of n slot must have same species names as in the params slot"
		errors <- c(errors, msg)
	}
	# w dimension of n
	if(dim(object@n)[3] != length(object@params@w)){
		msg <- "Third dimension of n slot must have same length as w in the params slot"
		errors <- c(errors, msg)
	}
	if(names(dimnames(object@n))[3] != "w"){
		msg <- "Third dimension of n slot must be called 'w'"
		errors <- c(errors, msg)
	}
	if(!all(names(dimnames(object@n))[3] == names(dimnames(object@params@psi))[2])){
		msg <- "Third dimension of n slot must have same size names as in the params slot"
		errors <- c(errors, msg)
	}
	# w dimension of n_pp
	if(dim(object@n_pp)[2] != length(object@params@w_full)){
		msg <- "Second dimension of n_pp slot must have same length as w_full in the params slot"
		errors <- c(errors, msg)
	}
	if(names(dimnames(object@n_pp))[2] != "w"){
		msg <- "Second dimension of n_pp slot must be called 'w'"
		errors <- c(errors, msg)
	}
	if(!all(dimnames(object@n_pp)$w == names(object@params@rr_pp))){
		msg <- "Second dimension of n_pp slot must have same size names as rr_pp in the params slot"
		errors <- c(errors, msg)
	}
	# w dimension of n_plank 
	if(dim(object@n_plank)[2] != (length(object@params@w_full))){ 
		msg <- "Second dimension of n_plank slot must have same length as w_fullin the params slot" 
		errors <- c(errors, msg) 
	} 
	if(names(dimnames(object@n_plank))[2] != "w"){ 
		msg <- "Second dimension of n_plank slot must be called 'w'" 
		errors <- c(errors, msg) 
	}
	if(!all(dimnames(object@n_plank)$w == names(object@params@rr_pp))){ 
		msg <- "Second dimension of n_plank slot must have same size names as plankton rr_pp in the params slot" 
		errors <- c(errors, msg) 
	} 
	# gear dimension of effort
	if(dim(object@effort)[2] != dim(object@params@catchability)[1]){
		msg <- "Second dimension of effort slot must have same number of gears as in the params slot"
		errors <- c(errors, msg)
	}
	if(names(dimnames(object@effort))[2] != "gear"){
		msg <- "Second dimension of effort slot must be called 'gear'"
		errors <- c(errors, msg)
	}
	if(!all(names(dimnames(object@effort))[2] == names(dimnames(object@params@catchability)[1]))){
		msg <- "Second dimension of effort slot must have same gear names as in the params slot"
		errors <- c(errors, msg)
	}
	# realm dimension of ocean_temp 
	if(dim(object@ocean_temp)[2] != dim(object@params@exposure)[1]){ 
		msg <- "Second dimension of ocean_temp slot must have the same number of realms as in the params slot" 
		errors <- c(errors, msg) 
	} 
	if(names(dimnames(object@ocean_temp))[2] != "realm"){ 
		msg <- "Second dimension of ocean_temp slot must be called 'realm'" 
		errors <- c(errors, msg) 
	} 
	if(!all(names(dimnames(object@ocean_temp))[2] == names(dimnames(object@params@exposure)[1]))){ 
		msg <- "Second dimension of ocean_temp slot must have same realms names as in the params slot" 
		errors <- c(errors, msg) 
	} 
	if (length(errors) == 0) TRUE else errors
}

# Soundtrack: Yob - Quantum Mystic

#' therMizerSim
#' 
#' A class that holds the results of projecting a \linkS4class{therMizerParams}
#' object through time.
#' 
#' \code{therMizerSim} objects are created by using the \code{\link{project_therMizer}} method
#' on an object of type \code{therMizerParams}.
#' 
#' There are several plotting methods available to explore the contents of a
#' \code{therMizerSim} object. See the package vignette for more details.
#' 
#' @slot params An object of type \linkS4class{therMizerParams}.
#' @slot n Array that stores the projected community population abundances by
#'   time, species and size
#' @slot effort Array that stores the fishing effort through time by time and
#'   gear
#' @slot ocean_temp Array that stores the ocean temperature through time by
#'   time and realm 
#' @slot n_pp Array that stores the projected background population by time and
#'   size
#' 
#' @seealso \code{\link{project_therMizer}} \code{\link{therMizerParams}}
#' @export
setClass(
	"therMizerSim",
	representation(
		params = "therMizerParams",
		n = "array",
		effort = "array",
		ocean_temp = "array", 
		n_pp = "array",
		n_plank = "array" 
	),
	prototype = prototype(
		params = new("therMizerParams"),
		n = array(
			NA,dim = c(1,1,1), dimnames = list(time = NULL, sp = NULL, w = NULL)
		),
		effort = array(
			NA,dim = c(1,1), dimnames = list(time = NULL, gear = NULL)
		),
		ocean_temp = array( 
			NA,dim = c(1,1), dimnames = list(time = NULL, realm = NULL) 
		), 
		n_pp = array(
			NA,dim = c(1,1), dimnames = list(time = NULL, w = NULL)
		),
		n_plank = array( 
			NA,dim = c(1,1), dimnames = list(time = NULL, w = NULL) 
		) 
	),
	validity = valid_therMizerSim
)

setValidity("therMizerSim", valid_therMizerSim)
remove(valid_therMizerSim)


# Constructors

#' Constructor for the \code{therMizerSim} class
#' 
#' A constructor for the \code{therMizerSim} class. This is used by the
#' \code{project_therMizer} method to create \code{therMizerSim} objects of the right
#' dimensions. It is not necessary for users to use this constructor.
#' 
#' @param object a \linkS4class{therMizerParams} object
#' @param t_dimnames Numeric vector that is used for the time dimensions of the
#'   slots. Default = NA.
#' @param t_max The maximum time step of the simulation. Only used if t_dimnames
#'   = NA. Default value = 100.
#' @param t_save How often should the results of the simulation be stored. Only
#'   used if t_dimnames = NA. Default value = 1.
#' @param ... Other arguments (currently not used).
#' 
#' @return An object of type \linkS4class{therMizerSim}
#' @seealso \code{\link{project_therMizer}} \linkS4class{therMizerParams}
#'   \linkS4class{therMizerSim}
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- therMizerParams(NS_species_params_gears, inter)
#' sim <- project_therMizer(params)
#' }
setGeneric('therMizerSim', function(object, ...)
	standardGeneric('therMizerSim'))

#' therMizerSim constructor taking only a \code{therMizerParams} object.
#' @rdname therMizerSim
setMethod('therMizerSim', signature(object='therMizerParams'),
	function(object, t_dimnames = NA, t_max = 100, t_save=1){
		# If the dimnames for the time dimension not passed in, calculate them
		# from t_max and t_save
		if (any(is.na(t_dimnames))){
			if((t_max %% t_save) != 0)
				stop("t_max must be divisible by t_save with no remainder")
			t_dimnames <- seq(from = t_save, to = t_max, by = t_save)
		}
		if (is.character(t_dimnames)){
			stop("The t_dimnames argument must be numeric.")
		}
		no_sp <- nrow(object@species_params)
		species_names <- dimnames(object@psi)$sp
		no_w <- length(object@w)
		w_names <- dimnames(object@psi)$w
		t_dimnames_n <-
			c(t_dimnames[1] - (t_dimnames[2] - t_dimnames[1]),t_dimnames) 
			# N is 1 bigger because it holds the initial population
		t_dim_n <- length(t_dimnames_n) 
		t_dim_effort <- length(t_dimnames)
		t_dim_ocean_temp <- length(t_dimnames) 
		t_dim_n_plank <- length(t_dimnames) 
		array_n <- array(NA, dim = c(t_dim_n, no_sp, no_w), 
						 dimnames = list(time = t_dimnames_n, 
						 				 sp = species_names, w = w_names))
		
		no_gears <- dim(object@selectivity)[1]
		gear_names <- dimnames(object@selectivity)$gear
		array_effort <- array(NA, dim = c(t_dim_effort, no_gears), 
							dimnames = list(time = t_dimnames, 
											gear = gear_names))

		no_realms <- dim(object@ontogenetic_migration)[1] 
		realm_names <- dimnames(object@ontogenetic_migration)$realm 
		array_ocean_temp <- array(NA, dim = c(t_dim_ocean_temp, no_realms), 
								dimnames = list(time = t_dimnames, 
												realm = realm_names)) 
		
		no_w_full <- length(object@w_full)
		w_full_names <- names(object@rr_pp)
		array_n_pp <- array(NA, dim = c(t_dim_n, no_w_full), 
							dimnames = list(time=t_dimnames_n, 
											w = w_full_names))
											
		no_n_plank <- length(object@w_full) 
		w_plank_names <- names(object@rr_pp) 
		array_n_plank <- array(NA, dim = c(t_dim_n_plank, no_n_plank), 
								dimnames = list(time=t_dimnames, 
												w = w_plank_names)) 

		sim <- new('therMizerSim',
				n = array_n, 
				effort = array_effort,
				ocean_temp = array_ocean_temp,
				n_pp = array_n_pp,
				n_plank = array_n_plank,
				params = object)
		return(sim)
		}
)