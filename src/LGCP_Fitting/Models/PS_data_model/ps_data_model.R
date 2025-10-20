## Fit PP model

library(INLA)
library(inlabru)
library(sf)
library(terra)
library(sp)

source("spde-book-functions.R")

fit_ps_data_model <- function(obs_pp,sampling_areas,covs,meshes,grf_priors,samplers){
  eco_mesh <- meshes$eco_mesh
  sampling_mesh <- meshes$sampling_mesh
  cmp_text = character()
  eco_formula_text = character()
  sampling_formula_text = character()
  for(i in names(covs)){
    for(j in 1:length(covs[[i]])){
      cmp_text = paste0(cmp_text,paste0(names(covs[[i]])[j],"(",paste0(names(covs[[i]])[j]),",model='linear')+"))
      if(i=="eco_covs"){eco_formula_text = paste0(eco_formula_text,names(covs[[i]])[j],"+")}
      if(i=="sampling_covs"){sampling_formula_text = paste0(sampling_formula_text,names(covs[[i]])[j],"+")}
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
  n.c1 <- nrow(obs_pp)
  # Define data: marks = 1 for observed points, 0 for mesh cells
  Y <- c(rep(1, n.c1), rep(0, nv))
  E <- c(rep(0, n.c1), w)
  
  # Putting together obs and mesh points
  sf_from_mat <- st_as_sf(as.data.frame(eco_mesh$loc[,1:2]), coords = c(1, 2), crs = st_crs(obs_pp))
  sf_from_obs <- st_as_sf(as.data.frame(st_coordinates(obs_pp)), coords = c(1, 2), crs = st_crs(obs_pp))
  full_sf <- rbind(sf_from_obs, sf_from_mat)
  
  # Building A-matrices
  n_poly <- nrow(sampling_areas)
  sampling_areas$map_id <- seq_len(n_poly)
  
  A_ps_data_point_area <- matrix(0,nrow=n.c1 + nv,ncol=n_poly)
  
  
  whichpoly <- st_join(full_sf, sampling_areas, join = st_intersects, left = T, largest = T)$map_id

  
  for(i in 1:length(whichpoly)){
    if(!is.na(whichpoly[i])){
      A_ps_data_point_area[i,whichpoly[i]] <- 1
    }
  }
  
  A_ps_data_point_area = Matrix(A_ps_data_point_area)
  
  whichpoly_mesh <- whichpoly[-(1:n.c1),]
  whichnapoly <- which(is.na(whichpoly_mesh$map_id))
  mesh.intersect <- sf_from_mat[-whichnapoly,]
  mi_1 <- whichpoly_mesh$map_id[-whichnapoly]
  mesh.intersect.id <- mi_1
  A_area <- INLA::inla.spde.make.A(mesh = eco_mesh, loc = as.matrix(mesh.intersect),
                                   block = mesh.intersect.id, block.rescale = "sum")
  
  
  exdatal1 <- list(y = Y,
                   E = E,
                   coordinates = st_coordinates(full_sf)
                   )
  
  
  for(i in names(covs)){
    for(j in 1:length(covs[[i]])){
      if(i == "eco_covs"){
        print(names(covs[[i]])[j])
        exdatal1[[names(covs[[i]])[j]]] <- raster::extract(covs[[i]][[j]],full_sf)[,2]
      }
      else{
        print(names(covs[[i]])[j])
        exdatal1[[names(covs[[i]])[j]]] <- as.vector(as.matrix(A_ps_data_point_area)%*%as.vector(as.matrix(A_area)%*%extract(covs[[i]][[j]],sf_from_mat)))
      }
    }
  }
  
  lik1_formula_text <- paste0("y0~as.vector(-1+eco_intercept+",eco_formula_text,
                              "w1+log(ilogit(ps_intercept+",sampling_formula_text,")))")
  
  if(substr(lik1_formula_text, nchar(lik1_formula_text)-3, nchar(lik1_formula_text)-3)=="+"){
    lik1_formula_text <- paste0(substr(lik1_formula_text, 1, nchar(lik1_formula_text) - 4), substr(lik1_formula_text, nchar(lik1_formula_text)-2, nchar(lik1_formula_text)))
  }

  lik1_formula <- as.formula(lik1_formula_text)

  #lik1_formula <- as.formula("y0~as.vector(eco_intercept+eco_cov + w1 + log(ilogit(ps_intercept+ps_cov+ps_eco_cov)))")
  #lik1_formula_alt <-as.formula("y0~as.vector(eco_intercept+eco_cov + w1 + test_offset)")
  lik1 <- bru_obs(family="poisson",
                  formula = lik1_formula,
                  data = exdatal1,
                  domain = list(coordinates = eco_mesh),
                  samplers = st_as_sf(samplers$eco_process),
  )  
  lik1$E <- lik1_0$E
  
  if(!("map_id"%in%names(obs_pp))){
    sampled_polys = unique(st_join(
      obs_pp,
      sampling_areas["map_id"],       # only bring in the ID column
      join = st_within,     # keeps only points inside polygons
      left = TRUE           # keeps all points; non-matches get NA
    )$map_id)
    
  }else{
    sampled_polys <- unique(obs_pp$map_id)
  }
  
  sampling_areas$which_lines <- 0
  sampling_areas$which_lines[sampled_polys] <- 1
  exdata <- list(y1 = sampling_areas$which_lines,
                 ps_intercept=rep(1,nrow(sampling_areas)),
                 coordinates1 = sampling_mesh$loc[,1:2])

  
  for(j in 1:length(covs[["sampling_covs"]])){
    ext_vals = extract(covs[["sampling_covs"]][[j]],sf_from_mat)[,2]
    exdata[[names(covs[["sampling_covs"]])[j]]] <- as.vector(as.matrix(A_area)%*%ext_vals)
  }
  
  lik2_formula <- as.formula(paste0("y1~ps_intercept+",sampling_formula_text))

  lik2 <- like("binomial",
               formula = lik2_formula,
               data = exdata,
               domain = list(coordinates = eco_mesh),
               samplers = st_as_sf(samplers$sampling_process))
  
  
  

  ## -- Defining prior distribution for w1
  eco_spde <- inla.spde2.pcmatern(mesh=eco_mesh,
                                  prior.range=grf_priors$eco_process$prior.range,
                                  prior.sigma=grf_priors$eco_process$prior.sigma)
  
  cmp <- as.formula(paste0("y0 ~ eco_intercept(1,model='linear')+ps_intercept(1,model='linear')+",cmp_text,
                           "w1(coordinates, model = eco_spde)"))  
  
  
  mod2 <- bru(cmp,lik1,lik2,options=list(bru_verbose=3,bru_max_iter=50,bru_method=list(rel_tol=0.2)))
  return(mod2)
  
}
