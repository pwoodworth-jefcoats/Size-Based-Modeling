# The purpose of this script is to run therMizer.
# These scenarios are for the North Pacific, specifically, the footprint of Hawaii's longline fishery for bigeye tuna.
# The code provided here allows for reproduction of the results documented in Woodworth-Jefcoats et al. 2019, in prep.
# Written by Phoebe.Woodworth-Jefcoats@noaa.gov, January 2019

# The full package structure is included, but you should run this script from:
# /R (folder in this repository)

# Load devtools (because the therMizer package is still in development)
library(devtools)

# Load the package in development
devtools::load_all()

# Creating a list of all the CMIP5 model abbreviations
# Their output will provide the plankton densities which determine the background resource as well as ocean temperatures.
models <- c("ESM2G", "ESM2M", "RCC", "ALR", "AMR","ILR", "IMR", "MRI")

# Creating a list of all the fishing effort scenarios
fishing <- c("Fifth", "Halve", "Constant", "Double", "5fold")

# Load parameters
# Initial
params_data <- read.csv("NPAC_species_params.csv")  
params <- therMizerParams(params_data)

# Load interaction matrix
inter <- read.csv("inter_NPAC.csv", row.names = 1)
inter <- as(inter, "matrix")
params <- therMizerParams(params_data, interaction = inter, kappa = 1e12)

# Build effort, temperature, and n_pp arrays
# 600 yrs spin-up + 95 year simulation
times <- seq(1406,2100,by=1) 

# Setting the fishing scenario.  
# Note that a few lines of code could be added to loop through all scenarios. 
f_mort <- fishing[1]

# Setting the climate model to be used.
# Note that a few lines of code could be added to loop through all models.
m <- models[1]

### Effort
gear_names <- c("Longline")
effort_array <- array(NA, dim = c(length(times), length(gear_names)), dimnames = list(time = times, gear = gear_names))
infile_effort <- paste("Effort_", f_mort, sep = "")
effort <- read.table(file = paste(infile_effort,".dat", sep = ""))
effort <- effort[,2]
effort <- as(effort, "matrix")

### Ocean temperature
realm_names <- c("Bigeye","Mahi","BlueShark","Skipjack","Yellowfin","Albacore","StripedMarlin","Wahoo","Swordfish","BlueMarlin","Lancetfish","Opah") # Because I'm having trouble loading these in with the params - need to fix this
# Empty arrays (dynamic and static climate)
ocean_temp_array <- array(NA, dim = c(length(times), length(realm_names)), dimnames = list(time = times, realm = realm_names))
ocean_temp_array_static <- array(NA, dim = c(length(times), length(realm_names)), dimnames = list(time = times, realm = realm_names))
# Read in and format data that will fill arrays
infile_temps <- paste("CMIP5_ocean_temp_array_", m, sep = "")
ocean_temps <- read.table(file = paste(infile_temps,".dat", sep = ""))
ocean_temps <- as(ocean_temps, "matrix")
	
### n_pp
# values are log10 abundance
# they will need to be transformed (inverse log) and divided by bin width (dw_full)
# this is done below when the arrays are filled
sizes = names(params@cc_pp)[1:225]
# Empty arrays (dynamic and static climate)
n_pp_array <- array(NA, dim = c(length(times), length(sizes)), dimnames = list(time = times, w = sizes))
n_pp_array_static <- array(NA, dim = c(length(times), length(sizes)), dimnames = list(time = times, w = sizes))
# Read in and format data that will fill arrays
infile_plankton <- paste("CMIP5_AbundAtSize_FullSpectra_scaled_", m, sep = "")
plankton <- read.table(file = paste(infile_plankton,".dat", sep = ""), header = FALSE, sep = "\t")
plankton <- as(plankton, "matrix")

### Now fill arrays
# We'll do static (no climate change) and dynamic (including climate change) scenarios.
# In the dynamic scenario, the first 600 years are for spin-up.
# During this time, the first input values are repeated at each time step.
# For the input provided, the first value is the 1986 - 2005 mean.
# Note that only 12 of the 14 species included in ocean_temps are used.
for (t in (times-1405)){
	# Dynamic climate scenario
	if (t <= 600){ 
		n_pp_array[t,] <- (10^(plankton[1,]))/params@dw_full
		ocean_temp_array[t,] <- (ocean_temps[1,c(1,2,3,4,5,6,7,8,10,11,13,14)]) 
		effort_array[t,"Longline"] <- effort[1]
	} else {
		n_pp_array[t,] <- (10^(plankton[t-599,]))/params@dw_full # t1 = historical average
		ocean_temp_array[t,] <- (ocean_temps[t-599,c(1,2,3,4,5,6,7,8,10,11,13,14)]) # t1 = historical average
		effort_array[t,"Longline"] <- effort[t-599]
		}
	# Static climate scenario
	n_pp_array_static[t,] <- (10^(plankton[1,]))/params@dw_full
	ocean_temp_array_static[t,] <-(ocean_temps[1,c(1,2,3,4,5,6,7,8,10,11,13,14)])
}
	
# Initial n_pp 
init_n_pp <- n_pp_array[1,]
	
### Run the model forward, doing four simulations:
# Static climate (static)
# Changing temperature only (noPlankton)
# Changing plankton only (noT)
# Full climate change (FullClimate) 
sim_static <- project_therMizer(params, effort = effort_array, dt = 0.1, t_save = 1, ocean_temp = ocean_temp_array_static, initial_n_pp = (init_n_pp), n_plank = (n_pp_array_static))
sim_noPlankton <- project_therMizer(params, effort = effort_array, dt = 0.1, t_save = 1, ocean_temp = ocean_temp_array, initial_n_pp = (init_n_pp), n_plank = (n_pp_array_static))
sim_noT <- project_therMizer(params, effort = effort_array, dt = 0.1, t_save = 1, ocean_temp = ocean_temp_array_static, initial_n_pp = (init_n_pp), n_plank = (n_pp_array))
sim_FullClimate <- project_therMizer(params, effort = effort_array, dt = 0.1, t_save = 1, ocean_temp = ocean_temp_array, initial_n_pp = (init_n_pp), n_plank = (n_pp_array))
	
### Save model runs
outfile_static <- paste("sim_static_",m,"_",f_mort, sep = "")
outfile_noPlankton <- paste("sim_noPlankton_",m,"_",f_mort, sep = "")
outfile_noT <- paste("sim_noT_",m,"_",f_mort, sep = "")
outfile_FullClimate <- paste("sim_FullClimate_",m,"_",f_mort, sep = "")
	
save(sim_static, file = paste(outfile_static,".Rdata", sep = ""), ascii = TRUE)
save(sim_noPlankton, file = paste(outfile_noPlankton,".Rdata", sep = ""), ascii = TRUE)
save(sim_noT, file = paste(outfile_noT,".Rdata", sep = ""), ascii = TRUE)
save(sim_FullClimate, file = paste(outfile_FullClimate,".Rdata", sep = ""), ascii = TRUE)
				
	
# View summary plots of the results 
quartz() # avoid potentially overwriting an open plot
plot(sim_FullClimate) 

