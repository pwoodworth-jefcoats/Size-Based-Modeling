# The purpose of this script is to run therMizer in the North Sea.
# Written by Phoebe.Woodworth-Jefcoats@noaa.gov, June 2020

# The full package structure is included, but you should run this script from:
# therMizer/R

# setwd("/Users/phoebe.woodworth/Documents/Integrated Model/therMizer/R/")
# source("/Users/phoebe.woodworth/Documents/Integrated Model/NorthSea_test/run_therMizer_NorthSea.R")

# Load devtools (because the therMizer package is still in development)
library(devtools)

# Load the package in development
devtools::load_all()

# Load parameters
# Initial
params_data <- read.csv("NS_species_params_gears_temps.csv")  

# Load interaction matrix
inter <- read.csv("inter.csv", row.names = 1)
inter <- as(inter, "matrix")

params <- therMizerParams(params_data, interaction = inter) 

# Build effort, temperature, and n_pp arrays
times <- seq(2050,2059,by=1) 

### Effort
gear_names <- c("Industrial","Pelagic","Beam","Otter")
effort_array <- array(NA, dim = c(length(times), length(gear_names)), dimnames = list(time = times, gear = gear_names))

### Ocean temperature
realm_names <- c("Sprat", "Sandeel", "N.pout", "Herring", "Dab", "Whiting", "Sole", "Gurnard", "Plaice", "Haddock", "Cod", "Saithe") 
ocean_temp_array <- array(NA, dim = c(length(times), length(realm_names)), dimnames = list(time = times, realm = realm_names))

# Read in and format data that will fill arrays
infile_temps <- paste("NorthSea_ocean_temp_array",sep = "")
ocean_temps <- read.table(file = paste(infile_temps,".dat", sep = ""))
ocean_temps <- as(ocean_temps, "matrix")
	
### n_pp
# values are log10 abundance
# they will need to be transformed (inverse log) and divided by bin width (dw_full)
# this is done below when the arrays are filled
all_sizes = names(params@cc_pp)[1:242]
# Empty array
n_pp_array <- array(NA, dim = c(length(times), length(all_sizes)), dimnames = list(time = times, w = all_sizes))

# Read in and format data that will fill arrays
infile_plankton <- paste("NorthSea_AbundAtSize_FullSpectraLong", sep = "")
plankton <- read.table(file = paste(infile_plankton,".dat", sep = ""), header = FALSE, sep = "\t")
plankton <- as(plankton, "matrix")

### Now fill arrays
for (t in (times-2049)){
	n_pp_array[t,] <- (10^(plankton[t+1,]))/params@dw_full[] # t1 = historical average
	ocean_temp_array[t,] <- (ocean_temps[t,]) # t1 = historical average
	effort_array[t, ] <- c(1,1,1,1)
}

	
# Initial n_pp 
init_n_pp <- n_pp_array[1,]
	
### Run the model forward, doing four simulations:
sim_NorthSea <- project_therMizer(params, effort = effort_array, dt = 0.1, t_save = 1, ocean_temp = ocean_temp_array, initial_n_pp = (init_n_pp), n_plank = (n_pp_array))
	
# View summary plots of the results 
quartz() # avoid potentially overwriting an open plot
plot(sim_NorthSea) 

