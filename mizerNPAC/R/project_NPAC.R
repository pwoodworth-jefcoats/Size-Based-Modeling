# Project method for the size based modelling package mizer

# Copyright 2012 Finlay Scott and Julia Blanchard. 
# Distributed under the GPL 2 or later 
# Maintainer: Finlay Scott, CEFAS
# Modified to include temperature effect by:
# Phoebe.Woodworth-Jefcoats@noaa.gov
# September 2017
# Expanded to include input plantkon abundances by
# Phoebe.Woodworth-Jefcoats@noaa.gov
# October 2017


# project_NPAC can dispatch with effort & ocean_temp being different classes (missing, numeric,
# array). (PW)

#' project_NPAC method for the size based modelling
#' 
#' Runs the size-based model simulation and projects the size based model
#' through time. \code{project_NPAC} is called using an object of type
#' \linkS4class{MizerParams_NPAC}, an object that contains the effort of the
#' fishing gears through time, and an object that contains the ocean temperature
#' in each realm through time. The method returns an object of type
#' \linkS4class{MizerSim_NPAC} which can then be explored with a range of summary and
#' plotting methods. (PW)
#' 
#' @param object A \code{MizerParams_NPAC} object
#' @param effort The effort of each fishing gear through time. See notes below.
#' @param ocean_temp The ocean temperature in each oceanic realm through time. 
#'   See notes below. (PW)
#' @param n_plank The log10 abundance-at-size for the plankton community. Has the
#'    dimensions of time step (or constant) x size. (PW)
#' @param t_max The maximum time the project_NPACion runs for. The default value is
#'   100. However, this argument is not needed if an array is used for the
#'   \code{effort} argument, in which case this argument is ignored. See notes
#'   below.
#' @param dt Time step of the solver. The default value is 0.1.
#' @param t_save The frequency with which the output is stored. The default
#'   value is 1.
#' @param initial_n The initial populations of the species. See the notes below.
#' @param initial_n_pp The initial population of the background spectrum. It
#'   should be a numeric vector of the same length as the \code{w_full} slot of
#'   the \code{MizerParams_NPAC} argument. By default the \code{cc_pp} slot of the
#'   \code{\link{MizerParams_NPAC}} argument is used.
#' @param ... Currently unused.
#' 
#' @note The \code{effort} argument specifies the level of fishing effort during
#' the simulation. It can be specified in three different ways: \itemize{ \item
#' A single numeric value. This specifies the effort of all fishing gears which
#' is constant through time (i.e. all the gears have the same constant effort). 
#' \item A numerical vector which has the same length as the number of fishing
#' gears. The vector must be named and the names must correspond to the gear
#' names in the \code{MizerParams_NPAC} object. The values in the vector specify the
#' constant fishing effort of each of the fishing gears, i.e. the effort is
#' constant through time but each gear may have a different fishing effort. 
#' \item A numerical array with dimensions time step x gear. This specifies the
#' fishing effort of each gear at each time step. The first dimension, time,
#' must be named numerically and contiguously. The second dimension of the array
#' must be named and the names must correspond to the gear names in the
#' \code{MizerParams_NPAC} argument. }
#' 
#' If effort is specified as an array then the \code{t_max} argument is ignored
#' and the maximum simulation time is the taken from the dimension names.
#' 
#' @note The \code{ocean_temp} argument specifies the ocean temperature during the
#' simulation. It can be specified in three different ways: \itemize{ \item A 
#' single numeric value. This specifies the temperature in all realms which is 
#' constant through time (i.e., all realms have the same constant temperature). 
#' \item A numerical vector which has the same length as the number of realms.
#' The vector must be named and the names must correspond to the realm names in 
#' the \code{MizerParams_NPAC} object. The values in the vectory specify the 
#' constant temperature for each of the realms, i.e., the temperature is 
#' constant through time but each realm may have a different temperature. 
#' \item A numerical array with dimensions time step x realm. This specifies
#' the temperature in each realm at each time step. The first dimenions, time,
#' must be named numerically and contiguously. The second dimension of the array
#' must be named and the names must correspond to the realm names in the 
#' \cide{MizerParams_NPAC} argument. } (PW)
#' 
#' If ocean_temp is specified as an array then the \code{t_max} argument is ignored
#' and the maximum simulation time is taken from the dimension names (PW)	
#' 
#' The \code{initial_n} argument is a matrix with dimensions species x size. The
#' order of species must be the same as in the \code{MizerParams_NPAC} argument. If
#' the initial population is not specified, the argument is set by default by
#' the \code{get_initial_n} function which is set up for a North Sea model.
#' @return An object of type \code{MizerSim_NPAC}.
#' @export
#' @seealso \code{\link{MizerParams_NPAC}}
#' @examples
#' \dontrun{
#' # Data set with different fishing gears
#' data(NS_species_params_gears_NPAC)
#' data(inter)
#' params <- MizerParams_NPAC(NS_species_params_gears_NPAC, inter)
#' # With constant fishing effort for all gears for 20 time steps
#' sim <- project_NPAC(params, t_max = 20, effort = 0.5)
#' # With constant fishing effort which is different for each gear
#' effort <- c(Industrial = 0, Pelagic = 1, Beam = 0.5, Otter = 0.5)
#' sim <- project_NPAC(params, t_max = 20, effort = effort)
#' # With fishing effort that varies through time for each gear
#' gear_names <- c("Industrial","Pelagic","Beam","Otter")
#' times <- seq(from = 1, to = 10, by = 1)
#' effort_array <- array(NA, dim = c(length(times), length(gear_names)),
#' dimnames = list(time = times, gear = gear_names))
#' effort_array[,"Industrial"] <- 0.5
#' effort_array[,"Pelagic"] <- seq(from = 1, to = 2, length = length(times))
#' effort_array[,"Beam"] <- seq(from = 1, to = 0, length = length(times))
#' effort_array[,"Otter"] <- seq(from = 1, to = 0.5, length = length(times))
#' sim <- project_NPAC(params, effort = effort_array)
#' }
setGeneric('project_NPAC', function(object, effort, ocean_temp, n_plank, ...) # PW
	standardGeneric('project_NPAC'))

# No effort is specified - default is to set an effort of 0 (PW)
# No ocean_temp is specified - default is set to a temp of 16.8 C (PW)
# All other arguments passed as ...

#' Project without an effort or ocean_temp argument. (PW)
#' @rdname project_NPAC
setMethod('project_NPAC', signature(object='MizerParams_NPAC', effort='missing', ocean_temp='missing', n_plank='missing'), # PW
	function(object, ...){
		res <- project_NPAC(object, effort=0, ocean_temp=16.8, ...) # PW
		return(res)
})

#' Project with a constant effort and ocean_temp. (PW)
#' @rdname project_NPAC
setMethod('project_NPAC', signature(object='MizerParams_NPAC', effort='numeric', ocean_temp='numeric', n_plank='numeric'), # PW
	function(object, effort, ocean_temp, n_plank, t_max = 100, dt = 0.1, ...){
	#if (!all.equal(t_max %% dt, 0))
		#if (!all((t_max %% dt) == 0))
	if(!all.equal((t_max - floor(t_max / dt) * dt),0))
			stop("t_max must be divisible by dt with no remainder")
		no_gears <- dim(object@catchability)[1]
		if ((length(effort)>1) & (length(effort) != no_gears))
			stop("Effort vector must be the same length as the number of fishing gears\n")
		no_realms <- dim(object@exposure)[1] # PW
		if ((length(ocean_temp)>1) & (length(ocean_temp) != no_realms)) # PW
			stop("Ocean_temp vector must be the same length as the number of realms\n") # PW
	# If more than 1 gear need to check that gear names match
		gear_names <- dimnames(object@catchability)[[1]]
	effort_gear_names <- names(effort)
	if (length(effort) == 1 & is.null(effort_gear_names)){
		effort_gear_names <- gear_names
	}
	if(!all(gear_names %in% effort_gear_names)){
		gear_names_error_message <- paste("Gear names in the MizerParams_NPAC object (", paste(gear_names, collapse=", "), ") do not match those in the effort vector.", sep="")
		stop(gear_names_error_message)
	}
	# If more than 1 realm need to check that realm names match (PS)
		realm_names <- dimnames(object@exposure)[[1]] # PW
	ocean_temp_realm_names <- names(ocean_temp) # PW
	if (length(ocean_temp) == 1 & is.null(ocean_temp_realm_names)){ # PW
		ocean_temp_realm_names <- realm_names # PW
	} # PW
	if(!all(realm_names %in% ocean_temp_realm_names)){ # PW
		realm_names_error_message <- paste("Realm names in the MizerParams_NPAC object (", paste(realm_names, collapse=", "), ") do not match those in the ocean_temp vector.", sep="") # PW
		stop(realm_names_error_message) # PW
	} # PW
		# Set up the effort, ocean_temp, and plankton arrays transposed so we can use the recycling rules (PW)
		plank_sz <- length(object@w_full)# PW
		plank_nm <- names(object@w_full) # PW
	time_dimnames <- signif(seq(from=1,to=t_max,by=dt),3)
		effort_array <- t(array(effort, dim=c(no_gears,length(time_dimnames)), dimnames=list(gear=effort_gear_names,time=time_dimnames)))
		ocean_temp_array <- t(array(ocean_temp, dim=c(no_realms,length(time_dimnames)), dimnames=list(realm=ocean_temp_realm_names,time=time_dimnames))) # PW
		n_plank_array <- t(array(n_plank, dim=c(plank_sz,length(time_dimnames)),dimnames=list(plank_nm,time=time_dimnames))) # PW
		res <- project_NPAC(object,effort_array, ocean_temp_array, n_plank_array, dt=dt, ...) # PW
		return(res)
})

#' Project with time varying effort
#' @rdname project_NPAC
setMethod('project_NPAC', signature(object='MizerParams_NPAC', effort='array', ocean_temp='array', n_plank='array'), # PW
	function(object, effort, ocean_temp, n_plank, t_save=1, dt=0.1, initial_n=get_initial_n(object), initial_n_pp=object@cc_pp, ...){ # PW
		validObject(object)
		# Check that number and names of gears in effort array is same as in MizerParams_NPAC object
		no_gears <- dim(object@catchability)[1]
		if(dim(effort)[2] != no_gears){
			no_gears_error_message <- paste("The number of gears in the effort array (length of the second dimension = ", dim(effort)[2], ") does not equal the number of gears in the MizerParams_NPAC object (", no_gears, ").", sep="")
			stop(no_gears_error_message)
		}
		gear_names <- dimnames(object@catchability)[[1]]
		if(!all(gear_names %in% dimnames(effort)[[2]])){
			gear_names_error_message <- paste("Gear names in the MizerParams_NPAC object (", paste(gear_names, collapse=", "), ") do not match those in the effort array.", sep="")
			stop(gear_names_error_message)
		}
		# Sort effort array to match order in MizerParams_NPAC
		effort <- effort[,gear_names, drop=FALSE]

		# Check that number and names of realms in ocean_temp array is same as in MizerParams_NPAC object (PW)
		no_realms <- dim(object@exposure)[1] # PW
		if(dim(ocean_temp)[2] != no_realms){ # PW
			no_realms_error_message <- paste("The number of realms in the ocean_temp array (length of the second dimension = ", dim(ocean_temp)[2], ") does not equal the number of realms in the MizerParams_NPAC object (", no_realms, ").", sep="") # PW
			stop(no_realms_error_message) # PW
		} # PW
		realm_names <- dimnames(object@exposure)[[1]] # PW
		if(!all(realm_names %in% dimnames(ocean_temp)[[2]])){ # PW
			realm_names_error_message <- paste("Realm names in the MizerParams_NPAC object (", paste(realm_names, collapse=", "), ") do not match those in the ocean_temp array.", sep="") # PW
			stop(realm_names_error_message) # PW
		} # PW
		# Sort ocean_temp array to match order in MizerParams_NPAC (PW)
		ocean_temp <- ocean_temp[,realm_names, drop=FALSE] # PW
		
		# Blow up time dimensions of effort, ocean_temp, and n_plank arrays (PW)
		# i.e. effort might have been passed in using time steps of 1, but actual dt = 0.1, so need to blow up
		# NOTE: This assumes that effort, ocean_temp, and n_plank have the same time step. (PW)
		# NOTE: This fills ocean_temp and n_plank with the time steps from effort (PW)
		if (is.null(dimnames(effort)[[1]])){
			stop("The time dimname of the effort argument must be numeric.")
		}
		if (any(is.na(as.numeric(dimnames(effort)[[1]])))){
			stop("The time dimname of the effort argument must be numeric.")
		}
		if (is.null(dimnames(ocean_temp)[[1]])){ # PW
			stop("The time dimname of the ocean_temp argument must be numeric.") # PW
		} # PW
		if (any(is.na(as.numeric(dimnames(ocean_temp)[[1]])))){ # PW
			stop("The time dimname of the ocean_temp argument must be numeric.") # PW
		} # PW
		if (is.null(dimnames(n_plank)[[1]])){ # PW
			stop("The time dimname of the n_plank argument must be numeric.") # PW
		} # PW
		if (any(is.na(as.numeric(dimnames(n_plank)[[1]])))){ # PW
			stop("The time dimname of the n_plank argument must be numeric.") # PW
		} # PW
		time_effort <- as.numeric(dimnames(effort)[[1]])
		t_max <- time_effort[length(time_effort)]
		# Blow up effort and ocean_temp so that rows are dt spaced (PW)
		time_effort_dt <- seq(from = time_effort[1], to = t_max, by = dt)
		effort_dt <- t(array(NA, dim = c(length(time_effort_dt), dim(effort)[2]), dimnames=list(time = time_effort_dt, dimnames(effort)[[2]])))
		ocean_temp_dt <- t(array(NA, dim = c(length(time_effort_dt), dim(ocean_temp)[2]), dimnames=list(time = time_effort_dt, dimnames(ocean_temp)[[2]]))) # PW
		n_plank_dt <- t(array(NA, dim = c(length(time_effort_dt), dim(n_plank)[2]), dimnames=list(time = time_effort_dt, dimnames(n_plank)[[2]]))) # PW
		for (i in 1:length(time_effort)){
			effort_dt[,time_effort_dt >= time_effort[i]] <- effort[i,]
			ocean_temp_dt[,time_effort_dt >= time_effort[i]] <- ocean_temp[i,] # PW
			n_plank_dt[,time_effort_dt >= time_effort[i]] <- n_plank[i,] # PW
		}
		effort_dt <- t(effort_dt)
		ocean_temp_dt <- t(ocean_temp_dt) # PW
		n_plank_dt <- t(n_plank_dt) # PW

		# Make the MizerSim_NPAC object with the right size
		# We only save every t_save steps
		#if (!all((t_save %% dt) == 0))
		if(!all.equal((t_max - floor(t_max / dt) * dt),0))
			stop("t_save must be divisible by dt with no remainder")
		t_dimnames_index <- as.integer(seq(from = 1+ ((t_save-1) / dt), to = length(time_effort_dt), by = t_save/dt))
		t_dimnames_index <- t_dimnames_index[t_dimnames_index>0]
		t_dimnames <- time_effort_dt[t_dimnames_index]
		sim <- MizerSim_NPAC(object, t_dimnames = t_dimnames) 
		# Fill up the effort array
		sim@effort[] <- effort_dt[t_dimnames_index,]
		sim@ocean_temp[] <- ocean_temp_dt[t_dimnames_index,] # PW
		sim@n_plank[] <- n_plank_dt[t_dimnames_index,] # PW

		# Set initial population
		sim@n[1,,] <- initial_n 
		sim@n_pp[1,] <- initial_n_pp

		# Handy things
		no_sp <- nrow(sim@params@species_params) # number of species
		no_w <- length(sim@params@w) # number of fish size bins
		idx <- 2:no_w
		# If no w_min_idx column in species_params, add one
		if (!("w_min_idx" %in% names(sim@params@species_params)))
			sim@params@species_params$w_min_idx <- 1
		# Hacky shortcut to access the correct element of a 2D array using 1D notation
		# This references the egg size bracket for all species, so for example
		# n[w_minidx_array_ref] = n[,w_min_idx]
		w_min_idx_array_ref <- (sim@params@species_params$w_min_idx-1) * no_sp + (1:no_sp)


		# sex ratio - DO SOMETHING LATER WITH THIS
		sex_ratio <- 0.5

		# Matrices for solver
		A <- matrix(0,nrow=no_sp,ncol=no_w)
		B <- matrix(0,nrow=no_sp,ncol=no_w)
		S <- matrix(0,nrow=no_sp,ncol=no_w)

		# initialise n and nPP
		# We want the first time step only but cannot use drop as there may only be a single species
		n <- array(sim@n[1,,],dim=dim(sim@n)[2:3])
		dimnames(n) <- dimnames(sim@n)[2:3]
		n_pp <- sim@n_pp[1,]
		t_steps <- dim(effort_dt)[1]
		for (i_time in 1:t_steps){
			# Do it piece by piece to save repeatedly calling methods
			# Calculate amount E_{a,i}(w) of available food
			phi_prey <- getPhiPrey(sim@params, n=n, n_pp=n_pp)
			# Calculate amount f_i(w) of food consumed
			feeding_level <- getFeedingLevel(sim@params, n=n, n_pp=n_pp, ocean_temp=ocean_temp_dt[i_time,], phi_prey=phi_prey) # PW
			# Calculate the predation rate
			pred_rate <- getPredRate(sim@params, n=n, n_pp=n_pp, feeding_level=feeding_level)
			# Calculate predation mortality on fish \mu_{p,i}(w)
			m2 <- getM2(sim@params, n=n, n_pp=n_pp, ocean_temp=ocean_temp_dt[i_time,], phi_prey=phi_prey, feeding_level=feeding_level, pred_rate=pred_rate) # PW
			#m2 <- getM2(sim@params, pred_rate=pred_rate)
			# Calculate total mortality \mu_i(w)
			z <- getZ(sim@params, n=n, n_pp=n_pp, effort=effort_dt[i_time,], m2=m2) # PW
			# Calculate predation mortality on the background spectrum
			m2_background <- getM2Background(sim@params, n=n, n_pp=n_pp, pred_rate=pred_rate) # PW
			# Calculate the resources available for reproduction and growth
			e <- getEReproAndGrowth(sim@params, n=n, n_pp=n_pp, ocean_temp=ocean_temp_dt[i_time,], feeding_level=feeding_level) # PW
			# Calculate the resources for reproduction
			e_spawning <- getESpawning(sim@params, n=n, n_pp=n_pp, e=e)
			# Calculate the growth rate g_i(w)
			e_growth <- getEGrowth(sim@params, n=n, n_pp=n_pp,  e_spawning=e_spawning, e=e)
			# R_{p,i}
			rdi <- getRDI(sim@params, n=n, n_pp=n_pp, e_spawning=e_spawning, sex_ratio=sex_ratio) 
			# R_i
			rdd <- getRDD(sim@params, n=n, n_pp=n_pp, rdi=rdi, sex_ratio=sex_ratio)

			# Iterate species one time step forward:
			# See Ken's PDF
			# A_{ij} = - g_i(w_{j-1}) / dw_j dt
			A[,idx] <- sweep(-e_growth[,idx-1,drop=FALSE]*dt, 2, sim@params@dw[idx], "/")
			# B_{ij} = 1 + g_i(w_j) / dw_j dt + \mu_i(w_j) dt
			B[,idx] <- 1 + sweep(e_growth[,idx,drop=FALSE]*dt,2,sim@params@dw[idx],"/") + z[,idx,drop=FALSE]*dt
			# S_{ij} <- N_i(w_j)
			S[,idx] <- n[,idx,drop=FALSE]
			# Boundary condition upstream end (recruitment)
			B[w_min_idx_array_ref] <- 1+e_growth[w_min_idx_array_ref]*dt/sim@params@dw[sim@params@species_params$w_min_idx]+z[w_min_idx_array_ref]*dt
			# Update first size group of n
			n[w_min_idx_array_ref] <- (n[w_min_idx_array_ref] + rdd*dt/sim@params@dw[sim@params@species_params$w_min_idx]) / B[w_min_idx_array_ref]
			# Update n
			for (i in 1:no_sp) # number of species assumed small, so no need to vectorize this loop over species
				for (j in (sim@params@species_params$w_min_idx[i]+1):no_w)
					n[i,j] <- (S[i,j] - A[i,j]*n[i,j-1]) / B[i,j]

			# Dynamics of background spectrum uses a semi-chemostat model (de Roos - ask Ken)
			# We use the exact solution under the assumption of constant mortality during timestep
			tmp <- (sim@params@rr_pp * sim@params@cc_pp / (sim@params@rr_pp + m2_background))
			n_pp <- tmp - (tmp - n_pp) * exp(-(sim@params@rr_pp+m2_background)*dt)
			
			# Reading in the background spectrum at each time step # PW
			n_pp <- n_plank_dt[i_time,] # PW

			# Store results only every t_step steps.
			store <- t_dimnames_index %in% i_time
			if (any(store)){
				sim@n[which(store)+1,,] <- n 
				sim@n_pp[which(store)+1,] <- n_pp
			}
		}
		# and end
		return(sim)
	}
)

#' Calculate initial population abundances for the community populations
#' 
#' This function uses the model parameters and other parameters to calculate 
#' initial population abundances for the community populations. These initial 
#' abundances should be reasonable guesses at the equilibrium values. The 
#' returned population can be passed to the \code{project_NPAC} method.
#' 
#' @param params The model parameters. An object of type \code{MizerParams_NPAC}.
#' @param a A parameter with a default value of 0.35.
#' @param n0_mult Multiplier for the abundance at size 0. Default value is
#' kappa/1000.
#' @export
#' @return A matrix (species x size) of population abundances.
#' @examples
#' \dontrun{
#' data(NS_species_params_gears_NPAC)
#' params <- MizerParams_NPAC(NS_species_params_gears_NPAC)
#' init_n <- get_initial_n(params)
#' }
get_initial_n<- function(params, n0_mult = NULL, a = 0.35){
	if (!is(params,"MizerParams_NPAC"))
		stop("params argument must of type MizerParams_NPAC")
	no_sp <- nrow(params@species_params)
	no_w <- length(params@w)
	initial_n <- array(NA, dim=c(no_sp,no_w))
	dimnames(initial_n) <- dimnames(params@intake_max)
	# N = N0 * Winf^(2*n-q-2+a) * w^(-n-a)
	# Reverse calc n and q from intake_max and search_vol slots (could add get_n as method)
	n <- (log(params@intake_max[,1] / params@species_params$h) / log(params@w[1]))[1]
	q <- (log(params@search_vol[,1] / params@species_params$gamma) / log(params@w[1]))[1]
	# Guessing at a suitable n0 value based on kappa - this was figured out using trial and error and should be updated
	if (is.null(n0_mult)){
		lambda <- 2+q-n
		kappa <- params@cc_pp[1] / (params@w_full[1]^(-lambda))
		n0_mult <- kappa / 1000
	}
	initial_n[] <- unlist(tapply(params@w,1:no_w,function(wx,n0_mult,w_inf,a,n,q)
		n0_mult * w_inf^(2*n-q-2+a) * wx^(-n-a),
		n0_mult=n0_mult, w_inf=params@species_params$w_inf, a=a, n=n, q=q))
	 #set densities at w > w_inf to 0
	initial_n[unlist(tapply(params@w,1:no_w,function(wx,w_inf) w_inf<wx, w_inf=params@species_params$w_inf))] <- 0
	# Also any densities at w < w_min set to 0
	initial_n[unlist(tapply(params@w,1:no_w,function(wx,w_min) w_min>wx, w_min=params@species_params$w_min))] <- 0 
	return(initial_n)
}