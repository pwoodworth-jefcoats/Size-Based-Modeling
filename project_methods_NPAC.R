# Methods used for project_NPACing for the size based modelling package

# Copyright 2012 Finlay Scott and Julia Blanchard. Distributed under the GPL 2 or later
# Maintainer: Finlay Scott, CEFAS
# Modified to include temperature effect by:
# Phoebe.Woodworth-Jefcoats@noaa.gov
# September 2017

# Calculate the amount of food exposed to each predator by predator size

#' getPhiPrey method for the size based model
#' 
#' Calculates the amount \eqn{E_{a,i}(w)} of food exposed to each predator by
#' predator size. This method is used by the \code{\link{project_NPAC}} method for
#' performing simulations.
#' @param object An \linkS4class{MizerParams_NPAC} object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the background abundance by size
#' @param ... Other arguments (currently unused)
#' 
#' @return A two dimensional array (predator species x predator size)
#' @seealso \code{\link{project_NPAC}}
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears_NPAC)
#' data(inter)
#' params <- MizerParams_NPAC(NS_species_params_gears_NPAC, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project_NPAC(params, t_max = 20, effort = 0.5)
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' getPhiPrey(params,n,n_pp)
#' }
setGeneric('getPhiPrey', function(object, n, n_pp,...)
	standardGeneric('getPhiPrey'))

#' @rdname getPhiPrey
setMethod('getPhiPrey', signature(object='MizerParams_NPAC', n = 'matrix', n_pp='numeric'),
	function(object, n, n_pp, ...){
# 		 cat("In getPhiPrey\n")
		# Check n dims
		if(dim(n)[1] != dim(object@interaction)[1])
			stop("n does not have the right number of species (first dimension)")
		if(dim(n)[2] != length(object@w))
			stop("n does not have the right number of size groups (second dimension)")
		if(length(n_pp) != length(object@w_full))
			stop("n_pp does not have the right number of size groups")
		# n_eff_prey is the total prey abundance by size exposed to each predator (prey
		# not broken into species - here we are just working out how much a predator
		# eats - not which species are being eaten - that is in the mortality calculation
	n_eff_prey <- sweep(object@interaction %*% n, 2, object@w * object@dw, "*", check.margin=FALSE) 
		# Quick reference to just the fish part of the size spectrum
		idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
		# pred_kernel is predator x predator size x prey size
		# So multiply 3rd dimension of pred_kernel by the prey abundance
		# Then sum over 3rd dimension to get total eaten by each predator by predator size
	# This line is a bottle neck
		phi_prey_species <- rowSums(sweep(object@pred_kernel[,,idx_sp,drop=FALSE],c(1,3),n_eff_prey,"*", check.margin=FALSE),dims=2)
		# Eating the background
	# This line is a bottle neck
		phi_prey_background <- rowSums(sweep(object@pred_kernel,3,object@dw_full*object@w_full*n_pp,"*", check.margin=FALSE),dims=2)
		return(phi_prey_species+phi_prey_background)
})

# Temperature effect by realm (PW)
# Gets temperature effect in each realm (PW)

#' getTempEffectRealm method for the size based model (PW)
#' 
#' Calcuates the temperature effect by realm, species, and size at each time step
#' in the \code{ocean_temp} argument. Used by the \code{project_NPAC} method to 
#' perform simulations (PW)
#'
#' @param object A \code{MizerParams_NPAC} object or a \code{MizerSim_NPAC}
#'   object.
#' @param ocean_temp The ocean temperature in each realm. Only needed if the object 
#'   argument is of class \code{MizerParams_NPAC}. See notes below.
#' @param time_range Subset the returned temperature effects by time. The time
#'   range is either a vector of values, a vector of min and max time, or a 
#'   single value. Default is the whole time range. Only used if the \code{object}
#'   argument is of type \code{MizerSim_NPAC}. 
#' @param ... Other arguments (currently unused)/ (PW)
#' 
#' @return An array. If the ocean_temp argument has a time dimension, or a 
#'   \code{MizerSim_NPAC} is passed in, the output array has four dimensions
#'   (time x realm x species x size). If the ocean_temp argument does not have a 
#'   time dimension (i.e., it is a vector or a single numeric), the output 
#'   array has three dimensions (realm x species x size).
#' @note Here: Temperature effect results from exposure x ontogenetic migration
#'   x ocean temperature. (PW)
#'
#' The \code{ocean_temp} argument is only used of a \code{MizerParams_NPAC} 
#'   object is passed in. The \code{ocean_temp} argument can be a two dimensional
#'   array (time x realm), a vector of length equal to the number of realms 
#'   (each realm has a different temperature that is constant in time), or a
#'   single numeric value (each realm has the same temperature that is constant
#'   in time). The order of realms in the \code{ocean_temp} argument must be the
#'   same as in the \code{MizerParams_NPAC} object. (PW)
#'
#' If the object argument is of class \code{MizerSim_NPAC} then the 
#'   ocean_temp slot of the \code{MizerSim_NPAC} object is used and the 
#'   \code{ocean_temp} argument is not used.
setGeneric('getTempEffectRealm', function(object, ocean_temp, ...) # PW
	standardGeneric('getTempEffectRealm')) # PW

#' \code{getTempEffectRealm} method for \code{MizerParams_NPAC} object with constant temperature. (PW)
#' @rdname getTempEffectRealm
# Ocean temperature is a single value or a numeric vector (PW)
# Ocean temperature has no time dimension (PW)
setMethod('getTempEffectRealm', signature(object='MizerParams_NPAC', ocean_temp='numeric'), # PW
	function(object, ocean_temp, ...){ # PW
		no_realm <- dim(object@exposure)[1] # PW
		# If a single value, just repeat it for all realms (PW)
		if(length(ocean_temp) == 1) # PW
			ocean_temp <- rep(ocean_temp, no_realm) # PW
		if(length(ocean_temp) != no_realm) # PW
			stop("Ocean_temp must be a single value or a vector as long as the number of realms\n") # PW
		out <- object@ontogenetic_migration # PW
		out[] <- (exp(25.22 - (0.63/((8.62*10^-5)*(273+ocean_temp))))) * c(object@exposure) * c(object@ontogenetic_migration) # PW
		return(out) # PW
	} # PW
) # PW

#' \code{getTempEffectRealm} method for \code{MizerParams_NPAC} object with time-changing temp.
#' @rdname getTempEffectRealm
# Always returns a 4D array: time x realm x species x size (PW)
setMethod('getTempEffectRealm',signature(object='MizerParams_NPAC', ocean_temp = 'matrix'), # PW
	function(object, ocean_temp, ...){ # PW
		no_realm <- dim(object@exposure)[1] # PW
		if (dim(ocean_temp)[2] != no_realm) # PW
			stop("Ocean_temp array must have a single value or a vector as long as the number of realms for each time step\n") # PW
		# Make the output array- note that we put time as last dimention and then aperm before returning (PW)
		# This is because of the order of the values when we al the other getTempEffectRealm method (PW)
		# Fill it up by calling the other method and passing in each line of the ocean_temp matrix (PW)
		out <- array(NA, dim=c(dim(object@ontogenetic_migration), dim(ocean_temp)[1]), dimnames= c(dimnames(object@ontogenetic_migration), list(time = dimnames(ocean_temp)[[1]]))) # PW
		out[] <- apply(ocean_temp, 1, function(x) getTempEffectRealm(object, x)) # PW
		out <- aperm(out, c(4,1,2,3)) # PW
	} # PW
) # PW

# Returns the temperature: time * realm * species * size (PW)
#' \code{getTempEffectRealm} method for \code{MizerSim_PW} object.
#' @rdname getTempEffectRealm
setMethod('getTempEffectRealm', signature(object='MizerSim_NPAC', ocean_temp='missing'), # PW
	function(object,ocean_temp, time_range=dimnames(object@ocean_temp)$time, ...){ # PW
		time_elements <- get_time_elements(object,time_range, slot_name="ocean_temp") # PW
		temp_effect_realm <- getTempEffectRealm(object@params, object@ocean_temp, ...) # PW
		return(temp_effect_realm[time_elements,,,,drop=FALSE]) # PW
}) # PW

# Total temperature effect across all realms (PW)
# species x size and maybe also by time if ocean_temp is time based (PW)

#' Get the total temperature effect across all realms by time, species, and size.
#' 
#' Calculates the temperature effect across all realms by species and size at 
#' each time step in teh \code{ocean_temp} argument.
#' The total temperature effect is just the sum of the temperature effects in
#' each realm (PW)
#'
#' @param object A \code{MizerParams_NPAC} object or a \code{MizerSim_NPAC} object
#' @param ocean_temp The ocean temperature in each realm. Only needed if the object
#'   argument is of class \code{MizerParams_NPAC}. See notes below.
#' @param time_range Subset the returned temperature effects by time. The time
#'   range is either a vector of values, a vector of min and max time, or a
#'   single value. Default is the whole time range. Only used if the
#'   \code{object} argument is of type \code{MizerSim_NPAC}.
#' @param drop Only used when object is of type \code{MizerSim_NPAC}. Should
#'   dimensions of length 1 be dropped, e.g., if your community has only one
#'   species it might make presentation of results easier. Default is TRUE.
#' @param ... Other arguments passed to \code{getTempEffect} method.
#'
#' @return An array. If the effort argument has a time dimension, or object is
#'   of class \code{MizerSim_PW}, the output array has three dimensions (time x 
#'   species x size). If the ocean_temp argument does not have a time dimension, the
#'   output array has two dimensions (species x size).
#' @note Here: temperature effect is determined from exposure x ontogenetic_migration
#'   x ocean_temp. (PW)
#' 
#' The \code{ocean_temp} argument is only used if a \code{MizerParams_NPAC} object is
#'   passed in. The \code{ocean_temp} argument can be a two dimension array (time x 
#'   realm), a vector of length equal to the number of realms (each realm has a 
#'   different temperature that is constant in time), or a single numeric value 
#'   (each realm has the same temperature that is constant in time). The order of 
#'   realms in the \code{ocean_temp} argument must be the same as in the 
#'   \code{MizerParams_NPAC} object. (PW)
#'
#' If the object argument is of class \code{MizerSim_NPAC} then the ocean_temp slot
#'   of the \code{MizerSim_NPAC} object is used and the \code{ocean_temp} argument 
#'   is not used. (PW)
setGeneric('getTempEffect', function(object, ocean_temp, ...) # PW
	standardGeneric('getTempEffect')) # PW

#' \code{getTempEffect} method for \code{MizerParams_PW} object with constant temperature.
#' @rdname getTempEffect
# Called from project_NPAC -> getFeedingLevel & getEReproAndGrowth
setMethod('getTempEffect', signature(object='MizerParams_NPAC', ocean_temp='numeric'), # PW
	function(object, ocean_temp, ...){ # PW
		tempEffectRealm <- getTempEffectRealm(object, ocean_temp, ...) # PW
		tempEffect <-colSums(tempEffectRealm) # PW
		return(tempEffect) # PW
}) # PW

#' \code{getTempEffect} method for \code{MizerParams_NPAC} object with time changing effort.
#' @rdname getTempEffect
setMethod('getTempEffect', signature(object='MizerParams_NPAC', ocean_temp='matrix'), # PW
	function(object, ocean_temp, ...){ # PW
		tempEffectRealm <- getTempEffectRealm(object, ocean_temp, ...) # PW
		tempEffect <- apply(tempEffectRealm, c(1,3,4), sum) # PW
		return(tempEffect) # PW
}) # PW

#' \code{getTempEffect} method for \code{MizerSim_NPAC} object.
#' @rdname getTempEffect
setMethod('getTempEffect', signature(object='MizerSim_NPAC', ocean_temp='missing'), # PW
	function(object, ocean_temp, time_range=dimnames(object@ocean_temp)$time, drop=TRUE, ...){ # PW
		time_elements <- get_time_elements(object,time_range, slot_name="ocean_temp") # PW
		tempEffect <- getTempEffect(object@params, object@ocean_temp, ...) # PW
		return(tempEffect[time_elements,,,drop=drop]) # PW
	} # PW
) # PW


# Feeding level
# The amount of food consumed by a predator, by each predator size

#' getFeedingLevel method for the size based model
#' 
#' Calculates the amount of food \eqn{f_i(w)} consumed by a predator by predator
#' size based on food availability, search volume and maximum intake. This
#' method is used by the \code{\link{project_NPAC}} method for performing
#' simulations.
#' @param object A \code{MizerParams_NPAC} or \code{MizerSim_NPAC} object
#' @param n A matrix of species abundance (species x size). Only used if 
#'   \code{object} argument is of type \code{MizerParams_NPAC}.
#' @param n_pp A vector of the background abundance by size. Only used if 
#'   \code{object} argument is of type \code{MizerParams_NPAC}.
#' @param ocean_temp A numeric vectory of the ocean temperature by realm or a 
#'   single numeric temperature value which is used for all realms. (PW)
#' @param phi_prey The PhiPrey matrix (optional) of dimension no. species x no. 
#'   size bins. If not passed in, it is calculated internally using the 
#'   \code{\link{getPhiPrey}} method. Only used if \code{object} argument is of type 
#'   \code{MizerParams_NPAC}.
#' @param time_range Subset the returned fishing mortalities by time. The time 
#'   range is either a vector of values, a vector of min and max time, or a 
#'   single value. Default is the whole time range. Only used if the 
#'   \code{object} argument is of type \code{MizerSim_NPAC}.
#' @param drop should extra dimensions of length 1 in the output be dropped, 
#'   simplifying the output. Defaults to TRUE.
#' @param ... Other arguments (currently unused).
#' 
#' @note If a \code{MizerParams_NPAC} object is passed in, the method returns a two 
#' dimensional array (predator species x predator size) based on the abundances 
#' also passed in.
#' 
#' If a \code{MizerSim_NPAC} object is passed in, the method returns a three
#' dimensional array (time step x predator species x predator size) with the
#' feeding level calculated at every time step in the simulation.
#' @seealso \code{\link{getPhiPrey}}
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears_NPAC)
#' data(inter)
#' params <- MizerParams_NPAC(NS_species_params_gears_NPAC, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project_NPAC(params, t_max = 20, effort = 0.5)
#' # Get the feeding level at one time step
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' fl <- getFeedingLevel(params,n,n_pp)
#' # Get the feeding level at all saved time steps
#' fl <- getFeedingLevel(sim)
#' # Get the feeding level for time 15 - 20
#' fl <- getFeedingLevel(sim, time_range = c(15,20))
#' }
setGeneric('getFeedingLevel', function(object, n, n_pp, ocean_temp, phi_prey, ...) # PW
	standardGeneric('getFeedingLevel'))

#' getFeedingLevel method for a \code{MizerParams_NPAC} object with already calculated \code{phi_prey} matrix.
#' @rdname getFeedingLevel
setMethod('getFeedingLevel', signature(object='MizerParams_NPAC', n = 'matrix', n_pp='numeric', ocean_temp='numeric', phi_prey='matrix'), # PW
	function(object, n, n_pp, ocean_temp, phi_prey, ...){ # PW
	# Check dims of phi_prey
		if (!all(dim(phi_prey) == c(nrow(object@species_params),length(object@w)))){
			stop("phi_prey argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
		}
		# encountered food = available food * search volume
		encount <- object@search_vol * phi_prey
		# calculate feeding level
		temp_effect <- getTempEffect(object, ocean_temp=ocean_temp) # PW
		f <- temp_effect * (encount/(encount + object@intake_max)) # PW
		f[f>1] <- 1 # Do not allow feeding to exceed intake_max # PW
		return(f)
	}
)

#' getFeedingLevel method for a \code{MizerSim_NPAC} object.
#' @rdname getFeedingLevel
setMethod('getFeedingLevel', signature(object='MizerSim_NPAC', n = 'missing', 
n_pp='missing', ocean_temp='missing', phi_prey='missing'), # PW
	function(object, ocean_temp, time_range=dimnames(object@ocean_temp)$time, drop=FALSE, ...){ # PW
		time_elements <- get_time_elements(object,time_range, slot_name="ocean_temp") # PW
		feed_time <- aaply(which(time_elements), 1, function(x){
			# Necessary as we only want single time step but may only have 1 species which makes using drop impossible
			n <- array(object@n[x,,],dim=dim(object@n)[2:3])
			dimnames(n) <- dimnames(object@n)[2:3]
				phi_prey <- getPhiPrey(object@params, n=n, n_pp=object@n_pp[x,]) # PW
					feed <- getFeedingLevel(object@params, n=n, n_pp = object@n_pp[x,], ocean_temp=object@ocean_temp[x,], phi_prey=phi_prey) # PW
					return(feed)}, .drop=drop)
			return(feed_time)
	}
)

# Predation rate
# Soundtrack: Nick Drake - Pink Moon

#' \code{getPredRate} method for the size based model
#' 
#' Calculates the predation rate of each predator species at size on prey size. 
#' In formulas \deqn{\phi_i(w_p/w) (1-f_i(w)) \gamma_i w^q N_i(w) dw}
#' This method is used by the \code{\link{project_NPAC}} method for performing
#' simulations. In the simulations, it is combined with the interaction matrix
#' (see \code{\link{MizerParams_NPAC}}) to calculate the realised predation mortality
#' (see \code{\link{getM2}}).
#' @param object A \code{MizerParams_NPAC} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param ocean_temp A numeric vectory of the ocean temperature by realm or a
#'   single numeric ocean_temp value which is used for all realms. (PW)
#' @param feeding_level The current feeding level (optional). A matrix of size
#'   no. species x no. size bins. If not supplied, is calculated internally
#'   using the \code{getFeedingLevel()} method.
#' 
#' @return A three dimensional array (predator species x predator size x prey size), 
#'   where the predator size runs over the community size range only and prey size
#'   runs over community plus background spectrum.
#' @export
#' @seealso \code{\link{project_NPAC}}, \code{\link{getM2}}, \code{\link{getFeedingLevel}} and \code{\link{MizerParams_NPAC}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears_NPAC)
#' data(inter)
#' params <- MizerParams_NPAC(NS_species_params_gears_NPAC, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project_NPAC(params, t_max = 20, effort = 0.5)
#' # Get the feeding level at one time step
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' getPredRate(params,n,n_pp)
#' }
setGeneric('getPredRate', function(object, n, n_pp, feeding_level) 
	standardGeneric('getPredRate'))

#' \code{getPredRate} method with \code{feeding_level} argument.
#' @rdname getPredRate
# Called from project_NPAC ->
setMethod('getPredRate', signature(object='MizerParams_NPAC', n = 'matrix', 
n_pp='numeric', feeding_level = 'matrix'), 
	function(object, n, n_pp, feeding_level){ 
		if (!all(dim(feeding_level) == c(nrow(object@species_params),length(object@w)))){
			stop("feeding_level argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
		}
		n_total_in_size_bins <- sweep(n, 2, object@dw, '*', check.margin=FALSE) # N_i(w)dw
		# The next line is a bottle neck
		pred_rate <- sweep(object@pred_kernel,c(1,2),(1-feeding_level)*object@search_vol*n_total_in_size_bins,"*", check.margin=FALSE)
		return(pred_rate)
	}
)

# getM2
# This uses the predation rate which is also used in M2background
# Too much overlap? Inefficient? Same thing is calculated twice

#' getM2 method for the size based model
#'
#' Calculates the total predation mortality \eqn{\mu_{p,i}(w_p)} on each prey
#' species by prey size. This method is used by the \code{\link{project_NPAC}} method
#' for performing simulations.
#' @param object A \code{MizerParams_NPAC} or \code{MizerSim_NPAC} object.
#' @param n A matrix of species abundance (species x size). Only used if
#'   \code{object} argument is of type \code{MizerParams_NPAC}.
#' @param n_pp A vector of the background abundance by size. Only used if
#'   \code{object} argument is of type \code{MizerParams_NPAC}.
#' @param ocean_temp A numeric vector of the ocean temperature by realm or a
#'   single numeric ocean_temp value which is used for all realms. (PW)
#' @param pred_rate An array of predation rates of dimension no. sp x no.
#'   community size bins x no. of size bins in whole spectra (i.e. community +
#'   background, the w_full slot). The array is optional. If it is not provided
#'   it is calculated by the \code{getPredRate()} method.
#' @param time_range Subset the returned fishing mortalities by time. The time
#'   range is either a vector of values, a vector of min and max time, or a
#'   single value. Default is the whole time range. Only used if the
#'   \code{object} argument is of type \code{MizerSim_NPAC}.
#' @param drop Only used when object is of type \code{MizerSim_NPAC}. Should
#'   dimensions of length 1 in the output be dropped, simplifying the output.
#'   Defaults to TRUE
#' @param ... Other arguments (currently unused).
#'
#' @return
#'   If a \code{MizerParams_NPAC} object is passed in, the method returns a two
#'   dimensional array (prey species x prey size) based on the abundances also
#'   passed in. If a \code{MizerSim_NPAC} object is passed in, the method returns a
#'   three dimensional array (time step x prey species x prey size) with the
#'   predation mortality calculated at every time step in the simulation.
#' @seealso \code{\link{getPredRate}} and \code{\link{project_NPAC}}.
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears_NPAC)
#' data(inter)
#' params <- MizerParams_NPAC(NS_species_params_gears_NPAC, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project_NPAC(params, t_max = 20, effort = 0.5)
#' # Get M2 at one time step
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' getM2(params,n,n_pp)
#' # Get M2 at all saved time steps
#' getM2(sim)
#' # Get M2 over the time 15 - 20
#' getM2(sim, time_range = c(15,20))
#' }
setGeneric('getM2', function(object, n, n_pp, ocean_temp, phi_prey, feeding_level, pred_rate, ...) # PW
	standardGeneric('getM2'))

#' \code{getM2} method for \code{MizerParams_NPAC} object with \code{pred_rate} argument.
#' @rdname getM2
setMethod('getM2', signature(object='MizerParams_NPAC', n = 'matrix', 
n_pp='numeric', ocean_temp='numeric', phi_prey='matrix', feeding_level = 'matrix', pred_rate = 'array'), # PW
	function(object,n, n_pp, ocean_temp, phi_prey, feeding_level, pred_rate){ # PW
		if ((!all(dim(pred_rate) == c(nrow(object@species_params),length(object@w),length(object@w_full)))) | (length(dim(pred_rate))!=3)){
			stop("pred_rate argument must have 3 dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),") x no. size bins in community + background (",length(object@w_full),")")
		}
		# get the element numbers that are just species
		idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
		# Interaction is predator x prey so need to transpose so it is prey x pred
		# Sum pred_kernel over predator sizes to give total predation rate of
		# each predator on each prey size
		m2 <- t(object@interaction) %*% colSums(aperm(pred_rate, c(2,1,3)),dims=1)[,idx_sp]
		return(m2)
	}
)

#' \code{getM2} method for \code{MizerSim_NPAC} object.
#' @rdname getM2
setMethod('getM2', signature(object='MizerSim_NPAC', n = 'missing', n_pp='missing', ocean_temp='missing', phi_prey='missing', feeding_level='missing', pred_rate = 'missing'), # PW
		function(object, ocean_temp, time_range=dimnames(object@ocean_temp)$time, drop=TRUE, ...){ # PW
				time_elements <- get_time_elements(object,time_range, slot_name='ocean_temp') # PW
				m2_time <- aaply(which(time_elements), 1, function(x){
			n <- array(object@n[x,,],dim=dim(object@n)[2:3])
			dimnames(n) <- dimnames(object@n)[2:3]
				phi_prey <- getPhiPrey(object@params, n=n, n_pp=object@n_pp[x,]) # PW
				feeding_level <- getFeedingLevel(object@params, n=n, n_pp=object@n_pp[x,], ocean_temp=object@ocean_temp[x,], phi_prey=phi_prey) # PW
				pred_rate <- getPredRate(object@params, n=n, n_pp=object@n_pp[x,], feeding_level=feeding_level) # PW
					m2 <- getM2(object@params, n=n, n_pp = object@n_pp[x,], ocean_temp=object@ocean_temp[x,], phi_prey=phi_prey, feeding_level=feeding_level, pred_rate=pred_rate) # PW
					return(m2)
				}, .drop=drop)
			return(m2_time)
			}
)


#' getM2Background method for the size based model
#'
#' Calculates the predation mortality \eqn{\mu_p(w)} on the background spectrum
#' by prey size. Used by the \code{project_NPAC} method for running size based
#' simulations.
#' @param object A \code{MizerParams_NPAC} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param ocean_temp A numeric vector of the ocean temperature by realm or a 
#'   single numeric ocean_temp value which is used for all realms. (PW)
#' @param pred_rate An array of predation rates of dimension no. sp x no.
#'   community size bins x no. of size bins in whole spectra (i.e. community +
#'   background, the w_full slot). The array is optional. If it is not provided
#'   it is calculated by the \code{getPredRate()} method.
#'
#' @return A vector of predation mortalities by background prey size.
#' @seealso \code{\link{project_NPAC}} and \code{\link{getM2}}.
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears_NPAC)
#' data(inter)
#' params <- MizerParams_NPAC(NS_species_params_gears_NPAC, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project_NPAC(params, t_max = 20, effort = 0.5)
#' # Get M2 of the background spectrum at one time step
#' n <- sim@@n[21,,]
#' n_pp <- sim@@n_pp[21,]
#' getM2Background(params,n,n_pp)
#' }
setGeneric('getM2Background', function(object, n, n_pp, pred_rate)
	standardGeneric('getM2Background'))

#' \code{getM2Background} method with \code{pred_array} argument.
#' @rdname getM2Background
setMethod('getM2Background', signature(object='MizerParams_NPAC', n = 'matrix', 
n_pp='numeric', pred_rate='array'),
	function(object, n, n_pp, pred_rate){
		if ((!all(dim(pred_rate) == c(nrow(object@species_params),length(object@w),length(object@w_full)))) | (length(dim(pred_rate))!=3)){
			stop("pred_rate argument must have 3 dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),") x no. size bins in community + background (",length(object@w_full),")")
		}
		M2background <- colSums(pred_rate,dims=2)
		return(M2background)
	}
)

# getFMortGear
#' Get the fishing mortality by time, gear, species and size
#'
#' Calculates the fishing mortality by gear, species and size at each time step
#' in the \code{effort} argument. Used by the \code{project_NPAC} method to perform
#' simulations.
#' 
#' @param object A \code{MizerParams_NPAC} object or a \code{MizerSim_NPAC} object.
#' @param effort The effort of each fishing gear. Only needed if the object
#'   argument is of class \code{MizerParams_NPAC}. See notes below.
#' @param time_range Subset the returned fishing mortalities by time. The time
#'   range is either a vector of values, a vector of min and max time, or a
#'   single value. Default is the whole time range. Only used if the
#'   \code{object} argument is of type \code{MizerSim_NPAC}.
#' @param ... Other arguments (currently unused).
#' 
#' @return An array. If the effort argument has a time dimension, or a
#'   \code{MizerSim_NPAC} is passed in, the output array has four dimensions (time x
#'   gear x species x size). If the effort argument does not have a time
#'   dimension (i.e. it is a vector or a single numeric), the output array has
#'   three dimensions (gear x species x size).
#' @note Here: fishing mortality = catchability x selectivity x effort.
#' 
#' The \code{effort} argument is only used if a \code{MizerParams_NPAC} object is
#' passed in. The \code{effort} argument can be a two dimensional array (time x
#' gear), a vector of length equal to the number of gears (each gear has a
#' different effort that is constant in time), or a single numeric value (each
#' gear has the same effort that is constant in time). The order of gears in the
#' \code{effort} argument must be the same the same as in the \code{MizerParams_NPAC}
#' object.
#' 
#' If the object argument is of class \code{MizerSim_NPAC} then the effort slot of
#' the \code{MizerSim_NPAC} object is used and the \code{effort} argument is not
#' used.
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears_NPAC)
#' data(inter)
#' params <- MizerParams_NPAC(NS_species_params_gears_NPAC, inter)
#' # Get the fishing mortality when effort is constant
#' # for all gears and time:
#' getFMortGear(params, effort = 1)
#' # Get the fishing mortality when effort is different
#' # between the four gears but constant in time:
#' getFMortGear(params, effort = c(0.5,1,1.5,0.75))
#' # Get the fishing mortality when effort is different
#' # between the four gears and changes with time:
#' effort <- array(NA, dim = c(20,4))
#' effort[,1] <- seq(from=0, to = 1, length=20)
#' effort[,2] <- seq(from=1, to = 0.5, length=20)
#' effort[,3] <- seq(from=1, to = 2, length=20)
#' effort[,4] <- seq(from=2, to = 1, length=20)
#' getFMortGear(params, effort=effort)
#' # Get the fishing mortality using the effort already held in a MizerSim_NPAC object.
#' sim <- project_NPAC(params, t_max = 20, effort = 0.5)
#' getFMortGear(sim)
#' getFMortGear(sim, time_range=c(10,20))
#' }
setGeneric('getFMortGear', function(object, effort, ...)
	standardGeneric('getFMortGear'))

#' \code{getFMortGear} method for \code{MizerParams_NPAC} object with constant effort.
#' @rdname getFMortGear
# Effort is a single value or a numeric vector.
# Effort has no time time dimension
setMethod('getFMortGear', signature(object='MizerParams_NPAC', effort = 'numeric'),
		function(object, effort, ...){
		no_gear <- dim(object@catchability)[1]
		# If a single value, just repeat it for all gears
		if(length(effort) == 1)
			effort <- rep(effort, no_gear)
		if (length(effort) != no_gear)
			stop("Effort must be a single value or a vector as long as the number of gears\n")
		# Streamlined for speed increase - note use of recycling
		out <- object@selectivity
		out[] <- effort * c(object@catchability) * c(object@selectivity)
		return(out)
	}
)

#' \code{getFMortGear} method for \code{MizerParams_NPAC} object with time changing effort.
#' @rdname getFMortGear
# Always returns a 4D array: time x gear x species x size
setMethod('getFMortGear', signature(object='MizerParams_NPAC', effort = 'matrix'),
	function(object, effort, ...){
		no_gear <- dim(object@catchability)[1]
		if (dim(effort)[2] != no_gear)
			stop("Effort array must have a single value or a vector as long as the number of gears for each time step\n")
	# Make the output array - note that we put time as last dimension and then aperm before returning 
	# This is because of the order of the values when we call the other getFMortGear method
	# Fill it up with by calling the other method and passing in each line of the effort matrix
	out <- array(NA, dim=c(dim(object@selectivity), dim(effort)[1]), dimnames= c(dimnames(object@selectivity), list(time = dimnames(effort)[[1]])))
	out[] <- apply(effort, 1, function(x) getFMortGear(object, x))
	out <- aperm(out, c(4,1,2,3))
	return(out)
	}
)

# Returns the fishing mortality: time * gear * species * size
#' \code{getFMortGear} method for \code{MizerSim_NPAC} object.
#' @rdname getFMortGear
setMethod('getFMortGear', signature(object='MizerSim_NPAC', effort='missing'),
	function(object,effort, time_range=dimnames(object@effort)$time, ...){
		time_elements <- get_time_elements(object,time_range, slot_name="effort")
		f_mort_gear <- getFMortGear(object@params, object@effort, ...)
		return(f_mort_gear[time_elements,,,,drop=FALSE])
})

# Total fishing mortality from all gears
# species x size and maybe also by time if effort is time based

#' Get the total fishing mortality from all fishing gears by time, species and
#' size.
#' 
#' Calculates the fishing mortality from all gears by species and size at each
#' time step in the \code{effort} argument.
#' The total fishing mortality is just the sum of the fishing mortalities
#' imposed by each gear.
#' 
#' @param object A \code{MizerParams_NPAC} object or a \code{MizerSim_NPAC} object
#' @param effort The effort of each fishing gear. Only needed if the object
#'   argument is of class \code{MizerParams_NPAC}. See notes below.
#' @param time_range Subset the returned fishing mortalities by time. The time
#'   range is either a vector of values, a vector of min and max time, or a
#'   single value. Default is the whole time range. Only used if the
#'   \code{object} argument is of type \code{MizerSim_NPAC}.
#' @param drop Only used when object is of type \code{MizerSim_NPAC}. Should
#'   dimensions of length 1 be dropped, e.g. if your community only has one
#'   species it might make presentation of results easier. Default is TRUE
#' @param ... Other arguments passed to \code{getFMortGear} method.
#'
#' @return An array. If the effort argument has a time dimension, or object is
#'   of class \code{MizerSim_NPAC}, the output array has three dimensions (time x
#'   species x size). If the effort argument does not have a time dimension, the
#'   output array has two dimensions (species x size).
#' @note Here: fishing mortality = catchability x selectivity x effort.
#'
#' The \code{effort} argument is only used if a \code{MizerParams_NPAC} object is
#' passed in. The \code{effort} argument can be a two dimensional array (time x
#' gear), a vector of length equal to the number of gears (each gear has a
#' different effort that is constant in time), or a single numeric value (each
#' gear has the same effort that is constant in time). The order of gears in the
#' \code{effort} argument must be the same the same as in the \code{MizerParams_NPAC}
#' object.
#'
#' If the object argument is of class \code{MizerSim_NPAC} then the effort slot of the \code{MizerSim_NPAC} object is used and the \code{effort} argument is not used.
#' @export
#' @seealso \code{getFMortGear}, \code{project_NPAC}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears_NPAC)
#' data(inter)
#' params <- MizerParams_NPAC(NS_species_params_gears_NPAC, inter)
#' # Get the total fishing mortality when effort is constant for all gears and time:
#' getFMort(params, effort = 1)
#' # Get the total fishing mortality when effort is different
#' # between the four gears but constant in time:
#' getFMort(params, effort = c(0.5,1,1.5,0.75))
#' # Get the total fishing mortality when effort is different
#' # between the four gears and changes with time:
#' effort <- array(NA, dim = c(20,4))
#' effort[,1] <- seq(from=0, to = 1, length=20)
#' effort[,2] <- seq(from=1, to = 0.5, length=20)
#' effort[,3] <- seq(from=1, to = 2, length=20)
#' effort[,4] <- seq(from=2, to = 1, length=20)
#' getFMort(params, effort=effort)
#' # Get the total fishing mortality using the effort already held in a MizerSim_NPAC
#' object.
#' sim <- project_NPAC(params, t_max = 20, effort = 0.5)
#' getFMort(sim)
#' getFMort(sim, time_range = c(10,20))
#' }
setGeneric('getFMort', function(object, effort, ...)
	standardGeneric('getFMort'))

#' \code{getFMort} method for \code{MizerParams_NPAC} object with constant effort.
#' @rdname getFMort
# Called from project_NPAC -> getZ -> 
setMethod('getFMort', signature(object='MizerParams_NPAC', effort='numeric'),
	function(object, effort, ...){
		fMortGear <- getFMortGear(object, effort, ...)
		fMort <- colSums(fMortGear)
		return(fMort)
})

#' \code{getFMort} method for \code{MizerParams_NPAC} object with time changing effort.
#' @rdname getFMort
setMethod('getFMort', signature(object='MizerParams_NPAC', effort='matrix'),
	function(object, effort, ...){
		fMortGear <- getFMortGear(object, effort, ...)
		fMort <- apply(fMortGear, c(1,3,4), sum)
		return(fMort)
})

#' \code{getFMort} method for \code{MizerSim_NPAC} object.
#' @rdname getFMort
setMethod('getFMort', signature(object='MizerSim_NPAC', effort='missing'),
	function(object, effort, time_range=dimnames(object@effort)$time, drop=TRUE, ...){
		time_elements <- get_time_elements(object,time_range, slot_name="effort")
		fMort <- getFMort(object@params, object@effort, ...)
		return(fMort[time_elements,,,drop=drop])
	}
)


# get total Z
#' getZ method for the size based model
#'
#' Calculates the total mortality \eqn{\mu_i(w)} on each species by size from
#' predation mortality (M2), background mortality (M) and fishing mortality for
#' a single time step.
#' @param object A \code{MizerParams_NPAC} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param ocean_temp A numeric vector of ocean temperature by realm or a 
#'   single numeric ocean_temp value which is used for all realms. (PW)
#' @param effort A numeric vector of the effort by gear or a single numeric
#'   effort value which is used for all gears.
#' @param m2 A two dimensional array of predation mortality (optional). Has
#'   dimensions no. sp x no. size bins in the community. If not supplied is
#'   calculated using the \code{getM2()} method.
#'
#' @return A two dimensional array (prey species x prey size). 
#'
#' @export
#' @seealso \code{\link{getM2}}, \code{\link{getFMort}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears_NPAC)
#' data(inter)
#' params <- MizerParams_NPAC(NS_species_params_gears_NPAC, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project_NPAC(params, t_max = 20, effort = 0.5)
#' # Get the total mortality at a particular time step
#' getZ(params,sim@@n[21,,],sim@@n_pp[21,],effort=0.5)
#' }
setGeneric('getZ', function(object, n, n_pp, effort, m2)
	standardGeneric('getZ'))

#' \code{getZ} method with \code{m2} argument.
#' @rdname getZ
# Called from project_NPAC()
setMethod('getZ', signature(object='MizerParams_NPAC', n = 'matrix', 
n_pp = 'numeric', effort='numeric', m2 = 'matrix'),
	function(object, n, n_pp, effort, m2){
		if (!all(dim(m2) == c(nrow(object@species_params),length(object@w)))){
			stop("m2 argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
		}
		f_mort <- getFMort(object, effort = effort)
		z <- sweep(m2 + f_mort,1,object@species_params$z0,"+", check.margin=FALSE) # not a slow sweep
		return(z)
	}
)


# Energy after metabolism and movement
#' getEReproAndGrowth method for the size based model
#'
#' Calculates the energy available by species and size for reproduction and
#' growth after metabolism and movement have been accounted for. Used by the
#' \code{project_NPAC} method for performing simulations.
#' @param object A \code{MizerParams_NPAC} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param ocean_temp A numeric vector of the ocean temperature by realm or a 
#'   single numeric ocean_temp value which is used for all realms. (PW)
#' @param feeding_level The current feeding level (optional). A matrix of size
#'   no. species x no. size bins. If not supplied, is calculated internally
#'   using the \code{getFeedingLevel()} method.
#'
#' @return A two dimensional array (species x size) 
#' @export
#' @seealso \code{\link{project_NPAC}} and \code{\link{getFeedingLevel}}.
#' @examples
#' \dontrun{
#' data(NS_species_params_gears_NPAC)
#' data(inter)
#' params <- MizerParams_NPAC(NS_species_params_gears_NPAC, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project_NPAC(params, t_max = 20, effort = 0.5)
#' # Get the energy at a particular time step
#' getEReproAndGrowth(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
setGeneric('getEReproAndGrowth', function(object, n, n_pp, ocean_temp, feeding_level) # PW
	standardGeneric('getEReproAndGrowth'))

#' \code{getEReproAndGrowth} method with \code{feeding_level} argument.
#' @rdname getEReproAndGrowth
setMethod('getEReproAndGrowth', signature(object='MizerParams_NPAC', n = 'matrix', 
n_pp = 'numeric', ocean_temp= 'numeric', feeding_level='matrix'), # PW
	function(object, n, n_pp, ocean_temp, feeding_level){ # PW
		if (!all(dim(feeding_level) == c(nrow(object@species_params),length(object@w)))){
stop("feeding_level argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
		}
		# assimilated intake
		e <- sweep(feeding_level * object@intake_max,1,object@species_params$alpha,"*", check.margin=FALSE)
		# Subtract basal metabolism and activity
		temp_effect <- getTempEffect(object, ocean_temp=ocean_temp) # PW
		e <- e - (temp_effect * object@std_metab) - object@activity
		e[e<0] <- 0 # Do not allow negative growth
		return(e)
	}
)

# Energy left for reproduction
# assimilated food intake, less metabolism and activity, split between reproduction and growth

#' getESpawning method for the size based model
#'
#' Calculates the energy available by species and size for reproduction after
#' metabolism and movement have been accounted for.
#' Used by the \code{project_NPAC} method for performing simulations.
#' @param object A \code{MizerParams_NPAC} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param ocean_temp A numeric vector of the ocean temperature by realm or a 
#'   single numeric ocean_temp value which is used for all realms. (PW)
#' @param e The energy available for reproduction and growth (optional). A
#'   matrix of size no. species x no. size bins. If not supplied, is calculated
#'   internally using the \code{getEReproAndGrowth()} method.
#'
#' @return A two dimensional array (prey species x prey size) 
#' @export
#' @seealso \code{\link{project_NPAC}} and \code{\link{getEReproAndGrowth}}.
#' @examples
#' \dontrun{
#' data(NS_species_params_gears_NPAC)
#' data(inter)
#' params <- MizerParams_NPAC(NS_species_params_gears_NPAC, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project_NPAC(params, t_max = 20, effort = 0.5)
#' # Get the energy at a particular time step
#' getESpawning(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
setGeneric('getESpawning', function(object, n, n_pp, e)
	standardGeneric('getESpawning'))

#' \code{getESpawning} method with \code{e} argument.
#' @rdname getESpawning
setMethod('getESpawning', signature(object='MizerParams_NPAC', n = 'matrix', 
n_pp = 'numeric', e = 'matrix'),
	function(object, n, n_pp, e){
		if (!all(dim(e) == c(nrow(object@species_params),length(object@w)))){
			stop("e argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
		}
		e_spawning <- object@psi * e 
		return(e_spawning)
	}
)

#' getEGrowth method for the size based model
#'
#' Calculates the energy \eqn{g_i(w)} available by species and size for growth
#' after metabolism, movement and reproduction have been accounted for. Used by
#' the \code{\link{project_NPAC}} method for performing simulations.
#' @param object A \linkS4class{MizerParams_NPAC} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param ocean_temp A numeric vector of the ocean temperature by realm or a 
#'   single numeric ocean_temp value which is used for all realms. (PW)
#' @param e The energy available for reproduction and growth (optional, although
#'   if specified, e_spawning must also be specified). A matrix of size no.
#'   species x no. size bins. If not supplied, is calculated internally using
#'   the \code{\link{getEReproAndGrowth}} method.
#' @param e_spawning The energy available for spawning (optional, although if
#'   specified, e must also be specified). A matrix of size no. species x no.
#'   size bins. If not supplied, is calculated internally using the
#'   \code{\link{getESpawning}} method.
#' 
#' @return A two dimensional array (prey species x prey size) 
#' @export
#' @seealso \code{\link{project_NPAC}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears_NPAC)
#' data(inter)
#' params <- MizerParams_NPAC(NS_species_params_gears_NPAC, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project_NPAC(params, t_max = 20, effort = 0.5)
#' # Get the energy at a particular time step
#' getEGrowth(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
setGeneric('getEGrowth', function(object, n, n_pp, e_spawning, e)
	standardGeneric('getEGrowth'))

#' \code{getEGrowth} method with \code{e_spawning} and \code{e} arguments.
#' @rdname getEGrowth 
setMethod('getEGrowth', signature(object='MizerParams_NPAC', n = 'matrix', 
n_pp = 'numeric', e_spawning='matrix', e='matrix'),
	function(object, n, n_pp, e_spawning, e){
		if (!all(dim(e_spawning) == c(nrow(object@species_params),length(object@w)))){
			stop("e_spawning argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
		}
		if (!all(dim(e) == c(nrow(object@species_params),length(object@w)))){
			stop("e argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
		}
		# Assimilated intake less activity and metabolism
		# energy for growth is intake - energy for growth
		e_growth <- e - e_spawning
		return(e_growth)
	}
)


#' getRDI method for the size based model
#'
#' Calculates the density independent recruitment (total egg production)
#' \eqn{R_{p,i}} before density dependence, by species. Used by the
#' \code{project_NPAC} method for performing simulations.
#' @param object A \code{MizerParams_NPAC} object.
#' @param n A matrix of species abundance (species x size).
#' @param n_pp A vector of the background abundance by size.
#' @param ocean_temp A numeric vector of the ocean temperature by realm or a 
#'   single numeric ocean_temp value which is used for all realms. (PW)
#' @param e_spawning The energy available for spawning (optional). A matrix of
#'   size no. species x no. size bins. If not supplied, is calculated internally
#'   using the \code{\link{getESpawning}} method.
#' @param sex_ratio Proportion of the population that is female. Default value
#'   is 0.5.
#' @param ... Other arguments (currently unused).
#' 
#' @return A numeric vector the length of the number of species 
#' @export
#' @seealso \code{\link{project_NPAC}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears_NPAC)
#' data(inter)
#' params <- MizerParams_NPAC(NS_species_params_gears_NPAC, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project_NPAC(params, t_max = 20, effort = 0.5)
#' # Get the recruitment at a particular time step
#' getRDI(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
setGeneric('getRDI', function(object, n, n_pp, e_spawning, ...)
	standardGeneric('getRDI'))

#' \code{getRDI} method with \code{e_spawning} argument.
#' @rdname getRDI
setMethod('getRDI', signature(object='MizerParams_NPAC', n = 'matrix', 
n_pp = 'numeric', e_spawning='matrix'),
	function(object, n, n_pp, e_spawning, sex_ratio = 0.5){
		if (!all(dim(e_spawning) == c(nrow(object@species_params),length(object@w)))){
			stop("e_spawning argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
		}
		# Should we put this in the class as part of species_params?
		# Index of the smallest size class for each species
		#w0_idx <- as.vector(tapply(object@species_params$w_min,1:length(object@species_params$w_min),function(w_min,wx) max(which(wx<=w_min)),wx=params@w))
		e_spawning_pop <- (e_spawning*n) %*% object@dw
		rdi <- sex_ratio*(e_spawning_pop * object@species_params$erepro)/object@w[object@species_params$w_min_idx] 
		return(rdi)
	}
)


#' getRDD method for the size based model
#'
#' Calculates the density dependent recruitment (total egg production) \eqn{R_i}
#' for each species. This is the flux entering the smallest size class of each
#' species. The density dependent recruiment is the density independent
#' recruitment after it has been put through the density dependent
#' stock-recruitment relationship function. This method is used by the
#' \code{project_NPAC} method for performing simulations.
#' @param object An \code{MizerParams_NPAC} object
#' @param n A matrix of species abundance (species x size)
#' @param n_pp A vector of the background abundance by size
#' @param ocean_temp A numeric vector of the ocean temperature by realm or a 
#'   single numeric ocean_temp value which is used for all realms. (PW)
#' @param rdi A matrix of density independent recruitment (optional) with
#'   dimensions no. sp x 1. If not specified rdi is calculated internally using
#'   the \code{\link{getRDI}} method.
#' @param sex_ratio Proportion of the population that is female. Default value
#'   is 0.5
#' @param ... Other arguments (currently unused).
#' 
#' @return A numeric vector the length of the number of species. 
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears_NPAC)
#' data(inter)
#' params <- MizerParams_NPAC(NS_species_params_gears_NPAC, inter)
#' # Project with constant fishing effort for all gears for 20 time steps
#' sim <- project_NPAC(params, t_max = 20, effort = 0.5)
#' # Get the energy at a particular time step
#' getRDD(params,sim@@n[21,,],sim@@n_pp[21,])
#' }
setGeneric('getRDD', function(object, n, n_pp, rdi, ...)
	standardGeneric('getRDD'))

#' \code{getRDD} method with \code{rdi} argument.
#' @rdname getRDD
setMethod('getRDD', signature(object='MizerParams_NPAC', n = 'matrix', 
n_pp = 'numeric', rdi='matrix'),
	function(object, n, n_pp, rdi, sex_ratio = 0.5){
		if (!all(dim(rdi) == c(nrow(object@species_params),1))){
			stop("rdi argument must have dimensions: no. species (",nrow(object@species_params),") x 1")
		}
		rdd <- object@srr(rdi = rdi, species_params = object@species_params)
		return(rdd)
})


# get_time_elements
# internal function to get the array element references of the time dimension
# for the time based slots of a MizerSim_NPAC object
# time_range can be character or numeric
# Necessary to include a slot_name argument because the effort and abundance
# slots have different time dimensions
get_time_elements <- function(sim,time_range,slot_name="n"){
	if (!(slot_name %in% c("n","effort","ocean_temp")))
		stop("'slot_name' argument should be 'n', 'effort', 'ocean_temp")
	if (!is(sim,"MizerSim_NPAC"))
		stop("First argument to get_time_elements function must be of class MizerSim_NPAC")
	time_range <- range(as.numeric(time_range))
	# Check that time range is even in object
	sim_time_range <- range(as.numeric(dimnames(slot(sim,slot_name))$time))
	if ((time_range[1] < sim_time_range[1]) | (time_range[2] > sim_time_range[2]))
		stop("Time range is outside the time range of the modell")
	time_elements <- (as.numeric(dimnames(slot(sim,slot_name))$time) >= time_range[1]) & (as.numeric(dimnames(slot(sim,slot_name))$time) <= time_range[2])
	names(time_elements) <- dimnames(slot(sim,slot_name))$time
	return(time_elements)
}