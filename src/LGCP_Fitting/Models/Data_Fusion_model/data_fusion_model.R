library(INLA)
library(inlabru)
library(sf)
library(terra)

fit_data_fusion_model <- function(cs_obs_pp,cs_locs,ps_obs_pp,ps_sampling_areas,covs,meshes,grf_priors,samplers){
  eco_mesh <- meshes$eco_mesh
  cs_sampling_mesh <- meshes$cs_sampling_mesh
  ps_sampling_mesh <- meshes$ps_sampling_mesh
  
  cmp_text = character()
  eco_formula_text = character()
  cs_sampling_formula_text = character()
  ps_sampling_formula_text = character()
  
  covs_in_cmp <- list()
  for(i in names(covs)){
    #print(i)
    for(j in 1:length(covs[[i]])){
      #print(paste0(names(covs[[i]])[j],"_terra"))
      #assign(paste0(names(covs[[i]])[j],"_terra"),rast(covs[[i]][[j]]))
      if(!names(covs[[i]])[j]%in%covs_in_cmp){
        cmp_text = paste0(cmp_text,paste0(names(covs[[i]])[j],"(",paste0(names(covs[[i]])[j]),",model='linear')+"))
        covs_in_cmp = c(covs_in_cmp,names(covs[[i]])[j])
      }
      if(i=="eco_covs"){eco_formula_text = paste0(eco_formula_text,names(covs[[i]])[j],"+")}
      if(i=="cs_sampling_covs"){cs_sampling_formula_text = paste0(cs_sampling_formula_text,names(covs[[i]])[j],"+")}
      if(i=="ps_sampling_covs"){ps_sampling_formula_text = paste0(ps_sampling_formula_text,names(covs[[i]])[j],"+")}
      }}
  
  ## -- Defining prior distribution for w1 and w2 -- ##
  eco_spde <- inla.spde2.pcmatern(mesh=eco_mesh,
                                  prior.range=grf_priors$eco_process$prior.range,
                                  prior.sigma=grf_priors$eco_process$prior.sigma)
  
  cs_spde <- inla.spde2.pcmatern(mesh=cs_sampling_mesh,
                                  prior.range=grf_priors$cs_sampling_process$prior.range,
                                  prior.sigma=grf_priors$cs_sampling_process$prior.sigma)
  
  
  # ps_spde <- inla.spde2.pcmatern(mesh=sampling_mesh,
  #                                prior.range=grf_priors$sampling_process$prior.range,
  #                                prior.sigma=grf_priors$sampling_process$prior.sigma)
  
  
  cmp <- as.formula(paste0("y1 ~ eco_intercept(1,model='linear')+ps_sampling_intercept(1,model='linear')+cs_sampling_intercept(1,model='linear')+",cmp_text,
                           "w1(coordinates, model = eco_spde)+w2(coordinates, model=cs_spde)+w3(coordinates1,copy='w1',fixed=F)"))
  
  lik1_0 <- bru_obs("cp",
                    formula = geometry ~ eco_intercept+eco_cov+w1,
                    data = st_as_sf(cs_obs_pp),
                    domain = list(geometry = eco_mesh),
                    samplers = st_as_sf(samplers$eco_process),
                    options=list(bru_compress_cp=FALSE)
  )

  for(i in names(covs)){
    for(j in 1:length(covs[[i]])){
      lik1_0$data[[names(covs[[i]])[j]]] <- raster::extract(covs[[i]][[j]],lik1_0$data)
    }
  }
  # cmp_test <- geometry ~ eco_intercept(1, model = "linear") + ps_sampling_intercept(1, 
  #                                                                                   model = "linear") + cs_sampling_intercept(1, model = "linear") + 
  #   eco_cov(eco_cov_terra, model = "linear") + cs_sampling_cov(cs_sampling_cov, 
  #                                                        model = "linear") + ps_sampling_cov(ps_sampling_cov, model = "linear") + 
  #   w1(geometry, model = eco_spde) + w2(geometry, model = cs_spde) + 
  #   w3(coordinates2, copy = "w1", fixed = F)
  #mod_test <- bru(cmp_test,lik1_0,options=list(bru_verbose=3))
  
  exdatal1 <- list(y1 = lik1_0$response_data$BRU_response_cp,
                   eco_intercept = rep(1,nrow(lik1_0$data)),
                   cs_sampling_intercept = rep(1,nrow(lik1_0$data)),
                   coordinates = st_coordinates(lik1_0$data)[,c("X","Y")],
                   #coordinates1 = eco_mesh$loc[,1:2],
                   x=(st_coordinates(lik1_0$data)[,"X"]),
                   y=(st_coordinates(lik1_0$data)[,"Y"]))

  for(i in names(covs)){
    for(j in 1:length(covs[[i]])){
        print(names(covs[[i]])[j])
        exdatal1[[names(covs[[i]])[j]]] <- lik1_0$data[[names(covs[[i]])[j]]]
    }
  }
  
  lik1_formula_text <- paste0("y1~eco_intercept+",eco_formula_text,
                              "w1+log(1-exp(-exp(cs_sampling_intercept+",cs_sampling_formula_text,"w2)))")
  #print(lik1_formula_text)
  #print("*****")
  lik1_formula <- as.formula(lik1_formula_text)
  
  A_cs_sampling_locs <- INLA::inla.spde.make.A(mesh = cs_sampling_mesh, loc = eco_mesh$loc[,1:2])
  
  lik1 <- bru_obs(family="poisson",
                  formula = lik1_formula,
                  data = exdatal1,
                  #domain = list(coordinates = eco_mesh),
                  samplers = st_as_sf(samplers$eco_process),
  )  
  lik1$E <- lik1_0$E
  #mod_test1 = bru(cmpa,lik1,options=list(bru_verbose=3))
  ####
  
  lik2_0 <- bru_obs("cp",
                    formula = geometry ~cs_sampling_intercept+cs_sampling_cov+w2,
                    data = st_as_sf(cs_locs),
                    domain = list(geometry = cs_sampling_mesh),
                    samplers = st_as_sf(samplers$cs_sampling_process),
                    options=list(bru_compress_cp=FALSE)
  )
  
  for(i in names(covs)){
    for(j in 1:length(covs[[i]])){
      lik2_0$data[[names(covs[[i]])[j]]] <- raster::extract(covs[[i]][[j]],lik2_0$data)
    }
  }
  #mod_test <- bru(cmp_test,lik2_0,options=list(bru_verbose=3))
  exdatal2 <- list(y1 = lik2_0$response_data$BRU_response_cp,
                   cs_sampling_intercept = rep(1,nrow(lik2_0$data)),
                   coordinates = st_coordinates(lik2_0$data)[,c("X","Y")],
                   x=(st_coordinates(lik2_0$data)[,"X"]),
                   y=(st_coordinates(lik2_0$data)[,"Y"]))
  
  for(i in names(covs)){
    for(j in 1:length(covs[[i]])){
      print(names(covs[[i]])[j])
      exdatal2[[names(covs[[i]])[j]]] <- lik2_0$data[[names(covs[[i]])[j]]]
    }
  }
  
  lik2_formula_text <- paste0("y1~ cs_sampling_intercept+",cs_sampling_formula_text,"w2")
  
  lik2_formula <- as.formula(lik2_formula_text)
  
  lik2 <- bru_obs(family="poisson",
                  formula = lik2_formula,
                  data = exdatal2,
                  #domain = list(coordinates = cs_sampling_mesh),
                  samplers = st_as_sf(samplers$cs_sampling_process),
  )  
  lik2$E <- lik2_0$E
  
  #mod_test2 <- bru(cmp,lik2,options=list(bru_verbose=3))
  #mod_test3 <- bru(cmp,lik1,lik2,options=list(bru_verbose=3,bru_max_iter=50))
  ############
  
  lik3_0 <- bru_obs("cp",
                    formula = geometry ~ eco_intercept+eco_cov+w1, #+fun(beta0pref,pref_samp_cov,zeta,w1),
                    data = st_as_sf(ps_obs_pp),
                    domain = list(geometry = eco_mesh),
                    samplers = st_as_sf(samplers$eco_process),
                    options=list(bru_compress_cp=FALSE)
  )
  
  for(i in names(covs)){
    for(j in 1:length(covs[[i]])){
      lik3_0$data[[names(covs[[i]])[j]]] <- raster::extract(covs[[i]][[j]],lik3_0$data)
    }
  }
  A_ps_data_point_area <- matrix(0,nrow=nrow(lik3_0$data),ncol=length(ps_sampling_areas))
  
  whichpoly <- over(as_Spatial(lik3_0$data),ps_sampling_areas)$id
  
  for(i in 1:length(whichpoly)){
    if(!is.na(whichpoly[i])){
      A_ps_data_point_area[i,whichpoly[i]] <- 1
    }
  }
  
  A_ps_data_point_area = Matrix(A_ps_data_point_area)
  eco_meshpoints <- eco_mesh$loc[,1:2]
  points.eco_mesh <- SpatialPoints(coords=cbind(eco_meshpoints[,1],eco_meshpoints[,2]),
                                   proj4string =crs(samplers$eco_process))
  
  whichpoly <- over(points.eco_mesh,ps_sampling_areas)
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
  
  
  exdatal3 <- list(y3 = lik3_0$response_data$BRU_response_cp,
                   eco_intercept = rep(1,nrow(lik3_0$data)),
                   ps_sampling_intercept=as.vector(as.matrix(A_ps_data_point_area)%*%rep(1,length(ps_sampling_areas))),
                   coordinates = st_coordinates(lik3_0$data)[,c("X","Y")],
                   coordinates1 = eco_mesh$loc[,1:2],
                   x=(st_coordinates(lik3_0$data)[,"X"]),
                   y=(st_coordinates(lik3_0$data)[,"Y"]))
  
  for(i in names(covs)){
    for(j in 1:length(covs[[i]])){
      if(i == "eco_covs"){
        print(names(covs[[i]])[j])
        exdatal3[[names(covs[[i]])[j]]] <- lik3_0$data[[names(covs[[i]])[j]]]
      }
      else{
        print(names(covs[[i]])[j])
        exdatal3[[names(covs[[i]])[j]]] <- as.vector(as.matrix(A_ps_data_point_area)%*%as.vector(as.matrix(A_area)%*%extract(covs[[i]][[j]],points.eco_mesh)))
      }
    }
  }
  lik3_formula_text <- paste0("y3~as.vector(-1+eco_intercept+",eco_formula_text,
                              "w1+log(ilogit(ps_sampling_intercept+",ps_sampling_formula_text,"wmat%*%w3)))")
  #print(substr(lik1_formula_text, nchar(lik1_formula_text)-3, nchar(lik1_formula_text)-3))
  #if(substr(lik1_formula_text, nchar(lik1_formula_text)-3, nchar(lik1_formula_text)-3)=="+"){
  #  lik1_formula_text <- paste0(substr(lik1_formula_text, 1, nchar(lik1_formula_text) - 4), substr(lik1_formula_text, nchar(lik1_formula_text)-2, nchar(lik1_formula_text)))
  #}
  #print(lik1_formula_text)
  #print("*****")
  lik3_formula <- as.formula(lik3_formula_text)
  lik3_formula <- as.formula("y3~as.vector(-1+eco_intercept+eco_cov+w1+log(ilogit(ps_sampling_intercept+ps_sampling_cov)))")
  lik3 <- bru_obs(family="poisson",
                  formula = lik3_formula,
                  data = exdatal3,
                  #domain = list(coordinates = eco_mesh),
                  samplers = st_as_sf(samplers$eco_process),
  )  
  lik3$E <- lik3_0$E
  
  # cmp2 <- y0 ~ eco_intercept(1, model = "linear") + ps_intercept(1, model = "linear") + 
  #   eco_cov(eco_cov, model = "linear") + ps_cov(ps_cov, model = "linear") + 
  #   w1(coordinates1, model = eco_spde)
  ilogit <- function(y){plogis(y)}
  #mod_test <- bru(cmp2,lik1,options=list(bru_verbose=3))
  exdatal4 <- list(y4 = sampling_areas$which_lines,
                 ps_intercept=rep(1,length(sampling_areas)),
                 coordinates1 = eco_mesh$loc[,1:2])
  #print(class(exdata$coordinates1))
  
  for(j in 1:length(covs[["ps_sampling_covs"]])){
    exdatal4[[names(covs[["ps_sampling_covs"]])[j]]] <- as.vector(as.matrix(A_area)%*%extract(covs[["ps_sampling_covs"]][[j]],points.eco_mesh))
  }
  
  lik4_formula <- as.formula(paste0("y4~ps_intercept+",ps_sampling_formula_text,"as.matrix(A_area)%*%w3"))
  lik4_formula <- as.formula(paste0("y4~ps_sampling_intercept+ps_sampling_cov"))
  lik4 <- like("binomial",
               formula = lik4_formula,
               data = exdatal4,
               #domain = list(coordinates = eco_mesh),
               samplers = st_as_sf(samplers$ps_sampling_process))
  
  print(cmp)
  print(lik1$formula)
  print(lik2$formula)
  
  wmat <- as.matrix(A_ps_data_point_area)%*%as.matrix(A_area)
  cmp3 <- y0 ~ ps_intercept(1, model = "linear") + ps_sampling_cov(ps_sampling_cov, model = "linear") +
   w3(coordinates1, model=eco_spde)
  #mod2 <- bru(cmp3,lik4,options = list(bru_verbose=3))
  mod3 <- bru(cmp,lik1,lik2,lik3,lik4,options=list(bru_verbose=3,bru_max_iter=50,bru_method=list(rel_tol=0.2)))
  return(mod3)
}