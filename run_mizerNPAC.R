# This script will run mizerNPAC and follows the general outline of the examples in the mizer vignette.
# Written by P. Woodworth-Jefcoats, Jan 2018
# The whole package structure is included, but you should run this script from:
# mizerNPAC/R

# Load devtools (because the mizerNPAC package is still in development)
library(devtools)

# Load the package in development 
devtools::load_all()

# Load parameters
params_data <- read.csv("NPAC_species_params.csv")  
params <- MizerParams_NPAC(params_data)

# Load interaction matrix
inter <- read.csv("inter_NPAC.csv", row.names = 1)
inter <- as(inter, "matrix")
params <- MizerParams_NPAC(params_data, interaction = inter)

# Build effort, temperature, and n_pp arrays
# 600 yrs spin-up + 110 year simulation
times <- seq(1391,2100,by=1) 

# Effort
gear_names <- c("Longline")
effort_array <- array(NA, dim = c(length(times), length(gear_names)), dimnames = list(time = times, gear = gear_names))
# For now, all fish are subject to the same level of constant fishing mortality
effort_array[,"Longline"] <- 0 

# Ocean temperature
realm_names <- c("Bigeye","Mahi","BlueShark","Skipjack","Yellowfin","Albacore","StripedMarlin","Wahoo","BigeyeThresherShark","Swordfish","BlueMarlin","ShortfinMakoShark","Lancetfish","Opah") # Because I'm having trouble loading these in with the params - need to fix this
ocean_temp_array <- array(NA, dim = c(length(times), length(realm_names)), dimnames = list(time = times, realm = realm_names))

# n_pp
# values are log10 abundance
# they will need to be transformed (inverse log) and divided by bin width (dw_full)
# LLFG = longline fishing grounds
sizes = names(params@cc_pp)[1:130]
n_pp_array <- array(NA, dim = c(length(times), length(sizes)), dimnames = list(time = times, w = sizes))
plankton_spinup <- read.table("ESM_spinup_AbundAtSize_FullSpectrum_LLFG.dat", header = FALSE, sep = "\t") 
plankton_spinup <- as(plankton_spinup, "matrix")
plankton <- read.table("ESM_annual_dynamic_AbundAtSize_FullSpectrum_LLFG.dat", header = FALSE, sep = "\t")
plankton <- as(plankton, "matrix")

for (t in (times-1390)){
	### If running the dynamic climate change scenario, uncomment lines 45 - 49
	# if (t <= 600){
		# n_pp_array[t,] <- (10^(plankton_spinup[1,]))/params@dw_full
	# } else {
		# n_pp_array[t,] <- (10^(plankton[t-600,]))/params@dw_full
	# }
	### If running the dynamic climate change scenario, comment out line 51
	n_pp_array[t,] <- (10^(plankton_spinup[1,]))/params@dw_full
	
	# Filling the ocean temp array with constant, species-specific 'optimal temps' 
	ocean_temp_array[t,] <- c(23.3,26.6,18.4,27.3,25.3,20.2,24.9,23.2,22.1,24.1,26.5,19.8,24.1,17.9)
}

# Initial n_pp 
# values are log10 abundance
# they will need to be transposed (rows to columns), transformed (inverse log) and divided by bin width (dw_full) 
init_n_pp <- read.table("ESM_InitialNpp_LLFG.dat", header = FALSE, sep = "\t") 
# Run the model forward
sim <- project_NPAC(params, effort = effort_array, dt = 0.1, t_save = 1, ocean_temp = ocean_temp_array, initial_n_pp = (10^t(init_n_pp))/(params@dw_full), n_plank = (n_pp_array))

# View a summary plot of the results
quartz() # avoid potentially overwriting an open plot
plot(sim) # Note: need to up my ggplot2 skills and fix code so that all species plot (right now Yellowfin isn't plotted)





