## Purpose: Reproducible analysis workflow
## Project: High_Severity_Birds
## Author: Zack Steel

## Run Multi-species Occupancy Model
source("Scripts/Analysis/Functions/runstan.R")

runstan(data_path = "Data/modeldata.RData",
        samples = 500, #n sample iterations
        warmup = 500, #n warmup samples
        nc = 4, #n chains
        pc = 4, #number of parallel chains (cores)
        cmdstanr = F,
        outfile = "Models/msom.RData") #path for saving the model (~8GB)

## Run Community models using z-matrix from MSOM
source("Scripts/Analysis/Functions/com_mods.R")
com_mods()