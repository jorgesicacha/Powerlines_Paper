## Fit PP model

library(INLA)
library(inlabru)
library(sf)
library(terra)

fit_ps_data_model <- function(obs_pp,sampling_areas,covs,meshes,grf_priors,samplers){
  eco_mesh <- meshes$eco_mesh
  sampling_mesh <- meshes$sampling_mesh
  cmp_text = character()
  eco_formula_text = character()
  sampling_formula_text = character()
  for(i in names(covs)){
    #print(i)
    for(j in 1:length(covs[[i]])){
      #print(paste0(names(covs[[i]])[j],"_terra"))
      #assign(paste0(names(covs[[i]])[j],"_terra"),rast(covs[[i]][[j]]))
      cmp_text = paste0(cmp_text,paste0(names(covs[[i]])[j],"(",paste0(names(covs[[i]])[j]),",model='linear')+"))
      if(i=="eco_covs"){eco_formula_text = paste0(eco_formula_text,names(covs[[i]])[j],"+")}
      if(i=="sampling_covs"){sampling_formula_text = paste0(sampling_formula_text,names(covs[[i]])[j],"+")}
    }}
  
  ## -- Defining prior distribution for w1 and w2 -- ##
  eco_spde <- inla.spde2.pcmatern(mesh=eco_mesh,
                                  prior.range=grf_priors$eco_process$prior.range,
                                  prior.sigma=grf_priors$eco_process$prior.sigma)
  
  # ps_spde <- inla.spde2.pcmatern(mesh=sampling_mesh,
  #                                prior.range=grf_priors$sampling_process$prior.range,
  #                                prior.sigma=grf_priors$sampling_process$prior.sigma)

  
  cmp <- as.formula(paste0("y0 ~ eco_intercept(1,model='linear')+ps_intercept(1,model='linear')+",cmp_text,
                           "w1(coordinates, model = eco_spde)+w3(coordinates1,copy='w1',fixed=F)"))
  
  eco_cov_terra = terra::rast(eco_cov)
  cmp1 <- geometry ~ eco_intercept(1, model = "linear") +
    eco_cov(eco_cov_terra, model = "linear") +
    w1(geometry, model = eco_spde)
  
  lik1_0 <- bru_obs("cp",
                    formula = geometry ~eco_intercept+eco_cov+w1, #+fun(beta0pref,pref_samp_cov,zeta,w1),
                    data = st_as_sf(obs_pp),
                    domain = list(geometry = eco_mesh),
                    samplers = st_as_sf(samplers$eco_process),
                    options=list(bru_compress_cp=FALSE)
                    )
  
  #mod_test <- bru(cmp1,lik1_0,options=list(bru_verbose=3))
  for(i in names(covs)){
    for(j in 1:length(covs[[i]])){
      lik1_0$data[[names(covs[[i]])[j]]] <- raster::extract(covs[[i]][[j]],lik1_0$data)
    }
  }
  A_ps_data_point_area <- matrix(0,nrow=nrow(lik1_0$data),ncol=length(sampling_areas))

  whichpoly <- over(as_Spatial(lik1_0$data),sampling_areas)$id
  
  for(i in 1:length(whichpoly)){
    if(!is.na(whichpoly[i])){
      A_ps_data_point_area[i,whichpoly[i]] <- 1
    }
  }
  
  A_ps_data_point_area = Matrix(A_ps_data_point_area)
  eco_meshpoints <- eco_mesh$loc[,1:2]
  points.eco_mesh <- SpatialPoints(coords=cbind(eco_meshpoints[,1],eco_meshpoints[,2]),
                                   proj4string =crs(samplers$eco_process))
  
  whichpoly <- over(points.eco_mesh,sampling_areas)
  whichnapoly <- which(is.na(whichpoly$id))
  mesh.intersect <- eco_mesh$loc[-whichnapoly,1:2]
  mi_1 <- whichpoly$id[-whichnapoly]
  tmi1 = table(mi_1)
  print(length(tmi1))
  print(tmi1)
  mesh.intersect.id <- mi_1
  A_area <- INLA::inla.spde.make.A(mesh = eco_mesh, loc = mesh.intersect,
                                   block = mesh.intersect.id, block.rescale = "none")
  
  A_area <- as.matrix(A_area)
  A_area <- t(apply(A_area,1,function(x){
    x[x>0]  =1
    x = x/sum(x)
    return(x)
  }))
  A_area = Matrix(A_area)
  print(dim(A_area))
  
  
  exdatal1 <- list(y0 = lik1_0$response_data$BRU_response_cp,
                   eco_intercept = rep(1,nrow(lik1_0$data)),
                   ps_intercept=as.vector(as.matrix(A_ps_data_point_area)%*%rep(1,length(sampling_areas))),
                   coordinates = st_coordinates(lik1_0$data)[,c("X","Y")],
                   coordinates1 = eco_mesh$loc[,1:2],
                   x=(st_coordinates(lik1_0$data)[,"X"]),
                   y=(st_coordinates(lik1_0$data)[,"Y"]))
  
  for(i in names(covs)){
    for(j in 1:length(covs[[i]])){
      if(i == "eco_covs"){
        print(names(covs[[i]])[j])
        exdatal1[[names(covs[[i]])[j]]] <- lik1_0$data[[names(covs[[i]])[j]]]
      }
      else{
        print(names(covs[[i]])[j])
        exdatal1[[names(covs[[i]])[j]]] <- as.vector(as.matrix(A_ps_data_point_area)%*%as.vector(as.matrix(A_area)%*%extract(covs[[i]][[j]],points.eco_mesh)))
      }
    }
  }
  lik1_formula_text <- paste0("y0~as.vector(-1+eco_intercept+",eco_formula_text,
         "w1+log(ilogit(ps_intercept+",sampling_formula_text,"wmat%*%w3)))")
  #print(substr(lik1_formula_text, nchar(lik1_formula_text)-3, nchar(lik1_formula_text)-3))
  if(substr(lik1_formula_text, nchar(lik1_formula_text)-3, nchar(lik1_formula_text)-3)=="+"){
    lik1_formula_text <- paste0(substr(lik1_formula_text, 1, nchar(lik1_formula_text) - 4), substr(lik1_formula_text, nchar(lik1_formula_text)-2, nchar(lik1_formula_text)))
  }
  #print(lik1_formula_text)
  #print("*****")
  lik1_formula <- as.formula(lik1_formula_text)
  print(lik1_formula_text)
  lik1_formula <- as.formula("y0~as.vector(eco_intercept+eco_cov + w1 + log(ilogit(ps_intercept+ps_cov)))")
  lik1 <- bru_obs(family="poisson",
                  formula = lik1_formula,
                  data = exdatal1,
                  domain = list(coordinates = eco_mesh),
                  samplers = st_as_sf(samplers$eco_process),
  )  
  lik1$E <- lik1_0$E
  
  # cmp2 <- y0 ~ eco_intercept(1, model = "linear") + ps_intercept(1, model = "linear") + 
  #   eco_cov(eco_cov, model = "linear") + ps_cov(ps_cov, model = "linear") + 
  #   w1(coordinates1, model = eco_spde)
  ilogit <- function(y){plogis(y)}
  #mod_test <- bru(cmp2,lik1,options=list(bru_verbose=3))
  exdata <- list(y1 = sampling_areas$which_lines,
                 ps_intercept=rep(1,length(sampling_areas)),
                 coordinates1 = sampling_mesh$loc[,1:2])
  #print(class(exdata$coordinates1))
  
  for(j in 1:length(covs[["sampling_covs"]])){
      exdata[[names(covs[["sampling_covs"]])[j]]] <- as.vector(as.matrix(A_area)%*%extract(covs[["sampling_covs"]][[j]],points.eco_mesh))
  }
  
  lik2_formula <- as.formula(paste0("y1~ps_intercept+",sampling_formula_text,"as.matrix(A_area)%*%w3"))
  lik2_formula <- as.formula("y1 ~ ps_intercept + ps_cov")
  
  lik2 <- like("binomial",
               formula = lik2_formula,
               data = exdata,
               domain = list(coordinates = eco_mesh),
               samplers = st_as_sf(samplers$sampling_process))

  print(cmp)
  print(lik1$formula)
  print(lik2$formula)

  wmat <- as.matrix(A_ps_data_point_area)%*%as.matrix(A_area)
  #cmp3 <- y0 ~ ps_intercept(1, model = "linear") + ps_cov(ps_cov, model = "linear") +
  # w3(coordinates1, model=eco_spde)
  #mod2 <- bru(cmp3,lik2,options = list(bru_verbose=3))
  mod2 <- bru(cmp,lik1,lik2,options=list(bru_verbose=3,bru_max_iter=50,bru_method=list(rel_tol=0.2)))
  return(mod2)
}