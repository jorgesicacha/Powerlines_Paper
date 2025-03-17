## Fit PP model

library(INLA)
library(inlabru)
library(sf)
library(terra)

fit_cs_data_model <- function(obs_pp,cs_locs,covs,meshes,grf_priors,samplers){
  
  eco_mesh <- meshes$eco_mesh
  sampling_mesh <- meshes$sampling_mesh
  cmp_text = character()
  eco_formula_text = character()
  #eco_nonlinear_formula_text = character()
  sampling_formula_text = character()
  for(i in names(covs)){
    #print(i)
    for(j in 1:length(covs[[i]])){
    #print(paste0(names(covs[[i]])[j],"_terra"))
    assign(paste0(names(covs[[i]])[j],"_terra"),rast(covs[[i]][[j]]))
    cmp_text = paste0(cmp_text,paste0(names(covs[[i]])[j],"(",paste0(names(covs[[i]])[j],"_terra"),",model='linear')+"))
    if(i=="eco_covs"){eco_formula_text = paste0(eco_formula_text,names(covs[[i]])[j],"+")}
    if(i=="sampling_covs"){sampling_formula_text = paste0(sampling_formula_text,names(covs[[i]])[j],"+")}
  }}
  
  ## -- Defining prior distribution for w1 and w2 -- ##
  eco_spde <- inla.spde2.pcmatern(mesh=eco_mesh,
                                  prior.range=grf_priors$eco_process$prior.range,
                                  prior.sigma=grf_priors$eco_process$prior.sigma)
  
  cs_spde <- inla.spde2.pcmatern(mesh=sampling_mesh,
                                  prior.range=grf_priors$sampling_process$prior.range,
                                  prior.sigma=grf_priors$sampling_process$prior.sigma)
  
  
  cmp <- as.formula(paste0("geometry ~ eco_intercept(1,model='linear')+sampling_intercept(1,model='linear')+",cmp_text,
                           "w1(geometry, model = eco_spde)+w2(geometry, model=cs_spde)"))
  
  lik1_formula <- as.formula(paste0("geometry~eco_intercept+",eco_formula_text,"w1+log(1-exp(-exp(sampling_intercept+",sampling_formula_text,"w2)))"))
  lik1 <- bru_obs("cp",
                  formula = lik1_formula,
                  data = st_as_sf(obs_pp),
                  domain = list(geometry = eco_mesh),
                  samplers = st_as_sf(samplers$eco_process)
  )
  
  lik2_formula <- as.formula(paste0("geometry~sampling_intercept+",sampling_formula_text,"w2"))
  lik2 <- bru_obs("cp",
                  formula = lik2_formula,
                  data = st_as_sf(cs_locs),
                  domain = list(geometry = sampling_mesh),
                  samplers = st_as_sf(samplers$sampling_process)
  )
  #print(cmp)
  #print(lik1$formula)
  #print(lik2$formula)
  mod1 <- bru(cmp,lik1,lik2,options = list(bru_verbose=3,bru_max_iter=50))
  return(mod1)
}