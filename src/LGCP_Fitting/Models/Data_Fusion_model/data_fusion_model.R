library(INLA)
library(inlabru)
library(sf)
library(terra)
library(sp)

source("spde-book-functions.R")

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
  
  # Dual mesh definition
  dmesh <- book.mesh.dual(eco_mesh)
  crs(dmesh) <- eco_mesh$crs
  dmesh_sf <- st_as_sf(dmesh)
  dmesh_sf$mesh_id <- seq_len(nrow(dmesh_sf))
  
  # Compute intersections: only those that intersect
  intersections <- st_intersection(dmesh_sf, samplers$eco_process)
  
  # Calculate area of each intersected piece
  intersections$area <- st_area(intersections)
  
  # Summarize total intersection area per mesh cell
  areas_df <- aggregate(area ~ mesh_id, data = intersections, FUN = sum)
  
  # Initialize area vector and fill with computed areas (zero for non-intersecting)
  w <- numeric(nrow(dmesh_sf))
  w[areas_df$mesh_id] <- as.numeric(areas_df$area)
  
  nv <- eco_mesh$n
  n.c1 <- nrow(cs_obs_pp)
  # Define data: marks = 1 for observed points, 0 for mesh cells
  Y <- c(rep(1, n.c1), rep(0, nv))
  E <- c(rep(0, n.c1), w)
  
  # Putting together obs and mesh points
  sf_from_mat <- st_as_sf(as.data.frame(eco_mesh$loc[,1:2]), coords = c(1, 2), crs = st_crs(cs_obs_pp))
  sf_from_obs <- st_as_sf(as.data.frame(st_coordinates(cs_obs_pp)), coords = c(1, 2), crs = st_crs(cs_obs_pp))
  cs_full_sf <- rbind(sf_from_obs, sf_from_mat)
  
  exdatal1 <- list(y = Y,
                   E = E,
                   coordinates = st_coordinates(cs_full_sf)
  )
  
  for(i in names(covs)){
    for(j in 1:length(covs[[i]])){
        exdatal1[[names(covs[[i]])[j]]] <- raster::extract(covs[[i]][[j]],cs_full_sf)[,2]
    }
  }
  
  lik1_formula_text <- paste0("y1~eco_intercept+",eco_formula_text,
                              "w1+log(1-exp(-exp(cs_sampling_intercept+",cs_sampling_formula_text,"w2)))")
  lik1_formula <- as.formula(lik1_formula_text)
  
  lik1 <- bru_obs(family="poisson",
                  formula = lik1_formula,
                  data = exdatal1,
                  samplers = st_as_sf(samplers$eco_process),
  )  
  lik1$E <- lik1_0$E
  
  dmesh_sampling <- book.mesh.dual(cs_sampling_mesh)
  crs(dmesh_sampling) <- cs_sampling_mesh$crs
  dmesh_sampling_sf <- st_as_sf(dmesh_sampling)
  dmesh_sampling_sf$mesh_id <- seq_len(nrow(dmesh_sampling_sf))
  
  # Compute intersections: only those that intersect
  intersections <- st_intersection(dmesh_sampling_sf, samplers$cs_sampling_process)
  
  # Calculate area of each intersected piece
  intersections$area <- st_area(intersections)
  
  # Summarize total intersection area per mesh cell
  areas_df <- aggregate(area ~ mesh_id, data = intersections, FUN = sum)
  
  # Initialize area vector and fill with computed areas (zero for non-intersecting)
  w_cs <- numeric(nrow(dmesh_sampling_sf))
  w_cs[areas_df$mesh_id] <- as.numeric(areas_df$area)
  
  cs_nv <- cs_sampling_mesh$n
  n.c2 <- nrow(cs_locs)
  # Define data: marks = 1 for observed points, 0 for mesh cells
  Y_cs <- c(rep(1, n.c2), rep(0, cs_nv))
  E_cs <- c(rep(0, n.c2), w_cs)
  
  # Putting together obs and mesh points
  sf_from_cs_mat <- st_as_sf(as.data.frame(cs_sampling_mesh$loc[,1:2]), coords = c(1, 2), crs = st_crs(cs_obs_pp))
  sf_from_cs_obs <- st_as_sf(as.data.frame(st_coordinates(cs_locs)), coords = c(1, 2), crs = st_crs(cs_obs_pp))
  cs_full_sf <- rbind(sf_from_cs_obs, sf_from_cs_mat)
  
  exdatal2 <- list(y_cs = Y_cs,
                   E = E_cs,
                   coordinates = st_coordinates(cs_full_sf)
  )
  
  for(i in names(covs)){
    for(j in 1:length(covs[[i]])){
      exdatal2[[names(covs[[i]])[j]]] <- raster::extract(covs[[i]][[j]],cs_full_sf)[,2]
    }
  }
  
  lik2_formula_text <- paste0("y_cs~eco_intercept+",eco_formula_text,
                              "w1+log(1-exp(-exp(cs_sampling_intercept+",cs_sampling_formula_text,"w2)))")
  lik2_formula <- as.formula(lik2_formula_text)
  
  lik2 <- bru_obs(family="poisson",
                  formula = lik2_formula,
                  data = exdatal2,
                  samplers = st_as_sf(samplers$cs_sampling_process),
  )  
  lik2$E <- lik2$E
  
  n.c3 <- nrow(ps_obs_pp)
  # Define data: marks = 1 for observed points, 0 for mesh cells
  Y_ps <- c(rep(1, n.c3), rep(0, nv))
  E <- c(rep(0, n.c3), w)
  
  # Putting together obs and mesh points
  sf_from_ps_obs <- st_as_sf(as.data.frame(st_coordinates(ps_obs_pp)), coords = c(1, 2), crs = st_crs(ps_obs_pp))
  ps_full_sf <- rbind(sf_from_ps_obs, sf_from_mat)
  
  # Building A-matrices
  n_poly <- nrow(ps_sampling_areas)
  ps_sampling_areas$map_id <- seq_len(n_poly)
  
  A_ps_data_point_area <- matrix(0,nrow=n.c3 + nv,ncol=n_poly)
  
  
  whichpoly <- st_join(ps_full_sf, ps_sampling_areas, join = st_intersects, left = T, largest = T)$map_id
  
  
  for(i in 1:length(whichpoly)){
    if(!is.na(whichpoly[i])){
      A_ps_data_point_area[i,whichpoly[i]] <- 1
    }
  }
  
  A_ps_data_point_area = Matrix(A_ps_data_point_area)
  
  whichpoly_mesh <- whichpoly[-(1:n.c3),]
  whichnapoly <- which(is.na(whichpoly_mesh$map_id))
  mesh.intersect <- sf_from_mat[-whichnapoly,]
  mi_1 <- whichpoly_mesh$map_id[-whichnapoly]
  mesh.intersect.id <- mi_1
  A_area <- INLA::inla.spde.make.A(mesh = eco_mesh, loc = as.matrix(mesh.intersect),
                                   block = mesh.intersect.id, block.rescale = "sum")
  
  
  exdatal3 <- list(y_ps = Y_ps,
                   E = E,
                   coordinates = st_coordinates(ps_full_sf)
  )
  
  
  for(i in names(covs)){
    for(j in 1:length(covs[[i]])){
      if(i == "eco_covs"){
        print(names(covs[[i]])[j])
        exdatal3[[names(covs[[i]])[j]]] <- raster::extract(covs[[i]][[j]],ps_full_sf)[,2]
      }
      else{
        print(names(covs[[i]])[j])
        exdatal3[[names(covs[[i]])[j]]] <- as.vector(as.matrix(A_ps_data_point_area)%*%as.vector(as.matrix(A_area)%*%extract(covs[[i]][[j]],sf_from_mat)))
      }
    }
  }
  
  lik3_formula_text <- paste0("y_ps~as.vector(-1+eco_intercept+",eco_formula_text,
                              "w1+log(ilogit(ps_intercept+",sampling_formula_text,")))")
  
  if(substr(lik3_formula_text, nchar(lik3_formula_text)-3, nchar(lik3_formula_text)-3)=="+"){
    lik3_formula_text <- paste0(substr(lik3_formula_text, 1, nchar(lik3_formula_text) - 4), substr(lik3_formula_text, nchar(lik3_formula_text)-2, nchar(lik3_formula_text)))
  }
  
  lik3_formula <- as.formula(lik3_formula_text)
  
  lik3 <- bru_obs(family="poisson",
                  formula = lik3_formula,
                  data = exdatal3,
                  domain = list(coordinates = eco_mesh),
                  samplers = st_as_sf(samplers$eco_process),
  )  
  lik3$E <- exdatal3$E
  
  if(!("map_id"%in%names(obs_pp))){
    sampled_polys = unique(st_join(
      ps_obs_pp,
      ps_sampling_areas["map_id"],       # only bring in the ID column
      join = st_within,     # keeps only points inside polygons
      left = TRUE           # keeps all points; non-matches get NA
    )$map_id)
    
  }else{
    sampled_polys <- unique(ps_obs_pp$map_id)
  }
  
  ps_sampling_areas$which_lines <- 0
  ps_sampling_areas$which_lines[sampled_polys] <- 1
  exdata <- list(y4 = ps_sampling_areas$which_lines,
                 ps_intercept=rep(1,nrow(ps_sampling_areas))
                 )
  
  
  for(j in 1:length(covs[["ps_sampling_covs"]])){
    ext_vals = extract(covs[["ps_sampling_covs"]][[j]],sf_from_mat)[,2]
    exdata[[names(covs[["ps_sampling_covs"]])[j]]] <- as.vector(as.matrix(A_area)%*%ext_vals)
  }
  
  lik4_formula <- as.formula(paste0("y4~ps_intercept+",ps_sampling_formula_text))
  
  lik4 <- like("binomial",
               formula = lik4_formula,
               data = exdata,
               samplers = st_as_sf(samplers$ps_sampling_process))
  
  
  
  ## -- Defining prior distribution for w1 and w2 -- ##
  eco_spde <- inla.spde2.pcmatern(mesh=eco_mesh,
                                  prior.range=grf_priors$eco_process$prior.range,
                                  prior.sigma=grf_priors$eco_process$prior.sigma)
  
  cs_spde <- inla.spde2.pcmatern(mesh=cs_sampling_mesh,
                                  prior.range=grf_priors$cs_sampling_process$prior.range,
                                  prior.sigma=grf_priors$cs_sampling_process$prior.sigma)
  
  
  cmp <- as.formula(paste0("y1 ~ eco_intercept(1,model='linear')+ps_sampling_intercept(1,model='linear')+cs_sampling_intercept(1,model='linear')+",cmp_text,
                           "w1(coordinates, model = eco_spde)+w2(coordinates, model=cs_spde)"))
  
  
  mod3 <- bru(cmp,lik1,lik2,lik3,lik4,options=list(bru_verbose=3,bru_max_iter=10,bru_method=list(rel_tol=0.2)))
  return(mod3)
}