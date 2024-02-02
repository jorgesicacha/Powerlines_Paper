## This is a Script to get stuff done about the power lines paper ###

## Step 1: Call a routine that simulates a regular LGCP
library(sp)
library(raster)
source("src/gen_eco_process.R")

win_seed = c(503.52441,6682.30628) #Somewhere in Norway
win_size = 10 #Size of the square (in kms)

geo_info = readRDS("data/geo_info_sc1.rds")

## Regular area
st_area = geo_info$st_area
res = 0.1 ## Resolution for the ecological covariate
res_gf = 0.1 ## Resolution for the GRF
grf_fam = "Matern"
grf_pars = list(range = 1.2, lambda = 1, sigma2 = 0.85)
eco_fixed_pars = list(beta0 = 2,beta1=-1.1)

# Generating the ecological process data
n_data_sets = 1
eco_process = generate_eco_process(st_area,res,res_gf,st_area,grf_pars,eco_fixed_pars,seed=1)

plot(crop(eco_process$grf.rast,st_area))
points(eco_process$final_pp,pch=19,cex=0.3)

