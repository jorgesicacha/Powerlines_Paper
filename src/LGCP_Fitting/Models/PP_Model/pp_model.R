## Fit PP model

library(INLA)
library(inlabru)
library(sf)
library(terra)

fit_pp_model <- function(pp,covs,mesh,grf_priors,samplers){
  
  cmp_text = character()
  formula_text = character()
  for(i in 1:length(covs)){
    assign(paste0(names(covs[i]),"_terra"),rast(covs[[i]]))
      cmp_text = paste0(cmp_text,paste0(names(covs[i]),"(",paste0(names(covs[i]),"_terra"),",model='linear')+"))
      formula_text = paste0(formula_text,names(covs[i]),"+")
  }
  
  ## -- Defining prior distribution for w1 and w2 -- ##
  eco_spde <- inla.spde2.pcmatern(mesh=mesh,
                                  prior.range=grf_priors$prior.range,
                                  prior.sigma=grf_priors$prior.sigma)
  
  
  
  cmp <- as.formula(paste0("geometry ~ eco_intercept(1,model='linear')+",cmp_text,
                           "w1(geometry, model = eco_spde)"))
  
  lik1_formula <- as.formula(paste0("geometry~eco_intercept+",formula_text,"w1"))
  lik1 <- bru_obs("cp",
                  formula = lik1_formula,
                  data = st_as_sf(pp),
                  domain = list(geometry = eco_mesh),
                  samplers = st_as_sf(samplers)
  )

  mod0 <- bru(cmp,lik1,options = list(bru_verbose=3))
  return(mod0)
}