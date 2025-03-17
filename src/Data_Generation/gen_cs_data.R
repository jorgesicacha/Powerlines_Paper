## Generating Citizen Science data
library(sf)
library(doParallel)
source("src/Data_Generation/studyarea_gen.R")
## The roads wil be generated following the same idea of the powerlines
#mapview(roads$st_area)+mapview(roads1)

generate_dist2roads_rast = function(roads,st_area,res){
  window = extent(st_area)
  st_outw = data.frame(xmin = window@xmin - 10,
                       xmax = window@xmax + 10,
                       ymin = window@ymin - 10,
                       ymax = window@ymax + 10)
  xgrid = seq(st_outw$xmin + 0.5*res,st_outw$xmax-0.5*res,res)
  ygrid = seq(st_outw$ymin + 0.5*res,st_outw$ymax-0.5*res,res)
  gridlocs = expand.grid(xgrid,ygrid)
  gridlocs_sp = SpatialPointsDataFrame(coords = gridlocs ,data=data.frame(id=1:nrow(gridlocs)),proj4string = crs(st_area))
  dist_grid <- st_as_sf(gridlocs_sp)
  dists_0 <- st_distance(dist_grid,st_as_sf(roads1))
  dist2road <- apply(dists_0,1,min)
  gridlocs_sp$dist = dist2road
  r <- raster(gridlocs_sp)
  r1<-disaggregate(r, fact=res(r)/c(res,res))
  dist.rast = rasterize(gridlocs_sp@coords,r1,gridlocs_sp$dist, fun=mean,na.rm=T)
  return(dist.rast)
}

generate_data_cs = function(eco_process_idx,thinning_pars,thinning_cov,gridlocs_gf_sp1,st_area,res_gf,window){
  eco_process_data <<- eco_process$data[[eco_process_idx]]
  print("I think I've produced eco_process_data")
  print(class(eco_process_data))
  seed = eco_process_data$seed
  orig_lgcp = eco_process_data$final_pp
  z0 = matrix(rnorm(nrow(A_m)),ncol=1)
  z = thinning_pars$fixed_effects$gamma0 + thinning_pars$fixed_effects$gamma1*extract(thinning_cov,gridlocs_gf_sp1) + A_m%*%z0
  zgf = A_m%*%z0
  print("Done with GRF")
  dft = data.frame(z =rep(NA,length(gridlocs_gf_sp1)),
                   zgf=rep(NA,length(gridlocs_gf_sp1)))
  dft$z = z
  dft$zgf = zgf
  grf_sp = SpatialPointsDataFrame(gridlocs_gf_sp1,data=dft)
  r <- raster(grf_sp)
  r1<-disaggregate(r, fact=res(r)/c(res_gf,res_gf))
  grf.rast = rasterize(grf_sp@coords,r1,grf_sp$z, fun=mean,na.rm=T)
  grf_sp$thinprob = 1-exp(-exp(grf_sp$z))
  thinprob.rast = 1-exp(-exp(grf.rast))
  #plot(thinprob.rast)
  #lines(roads1,col="red",lwd=2)
  
  ## Two remaining tasks: a) Generate auxiliary data. b) Thin the original PP
  
  ## Generate auxiliary data
  
  grf_sp$lambda = exp(grf_sp$z)
  grf_sp_vals = grf_sp[complete.cases(grf_sp@data),]
  lambda.rast = rasterize(grf_sp@coords,r1,grf_sp$lambda, fun=mean,na.rm=T)
  print (" Done with lambda raster")
  #mapView(lambda.rast)
  ### Simulate as an NHPP since conditional on lambda, a LGCP is a NHPP
  lambda_max = max(grf_sp$lambda,na.rm = T)
  N_max = rpois(1,lambda_max*gArea(st_area))
  ## Randomly localize the Nmax points in the buffers
  pts = list()
  sp_locs = SpatialPoints(coords = matrix(c(runif(N_max,window@xmin,window@xmax),
                                            runif(N_max,window@ymin,window@ymax)),ncol=2),
                          proj4string=crs(st_area))
  pts_tmp = raster::intersect(sp_locs,st_area)
  S_max = pts_tmp[1:N_max]  
  
  ## Finally, the points of S_max will be thinned according to their intensity
  lambda_S_max = extract(lambda.rast,S_max)
  ## Estimating lambda if missing: Nearest cell
  lambda_S_max[is.na(lambda_S_max)] <-  grf_sp_vals$lambda[sapply(which(is.na(lambda_S_max)),function(x){
    which.min(spDistsN1(grf_sp_vals@coords,cbind(S_max@coords[x,1],S_max@coords[x,2])))
  })]
  
  thinprobs = lambda_S_max/lambda_max
  S_max = SpatialPointsDataFrame(S_max,data=data.frame(ret=sapply(thinprobs,function(x){rbinom(1,1,x)})))
  cs_sampling_pp = S_max[S_max$ret==1,]
  
  ## Thin the original pp
  orig_lgcp$thinprob = extract(thinprob.rast,orig_lgcp)
  orig_lgcp$thinprob[is.na(orig_lgcp$thinprob)] <-  grf_sp$thinprob[sapply(which(is.na(orig_lgcp$thinprob)),function(x){
    which.min(spDistsN1(grf_sp@coords,cbind(orig_lgcp@coords[x,1],orig_lgcp@coords[x,2])))
  })]
  thinned_pp = SpatialPointsDataFrame(orig_lgcp@coords,data=data.frame(ret=sapply(orig_lgcp$thinprob,function(x){rbinom(1,1,x)})))
  thinned_pp = thinned_pp[which(thinned_pp$ret==1),]
  
  return(list(grf_sp=grf_sp,thinprob.rast=thinprob.rast,cs_sampling_pp=cs_sampling_pp,thinned_pp=thinned_pp))
}


generate_cs_data = function(eco_process,thinning_pars,thinning_cov,st_area){
  ## Begin by generating the thinning probability raster and auxiliary data
  ## It'll be assumed as a LGCP as well
  eco_process <<- eco_process
  thinning_cov <<- thinning_cov
  window = extent(st_area)
  st_outw = data.frame(xmin = window@xmin - 10,
                       xmax = window@xmax + 10,
                       ymin = window@ymin - 10,
                       ymax = window@ymax + 10)
  ## Now I have to simulate a Gaussian RF, the resolution has to be high
  xgrid_gf = seq(st_outw$xmin + 0.5*res_gf,st_outw$xmax-0.5*res_gf,res_gf)
  ygrid_gf = seq(st_outw$ymin + 0.5*res_gf,st_outw$ymax-0.5*res_gf,res_gf)
  gridlocs_gf_sp = SpatialPoints(expand.grid(xgrid_gf,ygrid_gf),proj4string = crs(st_area))
  
  ## Op2
  idx_int <<- rgeos::gIntersects(gridlocs_gf_sp,st_area,byid = TRUE)
  gridlocs_gf_sp1 <<- raster::intersect(gridlocs_gf_sp,st_area)
  ### the distance between the points in the grid
  grid_gf_dist = spDists(gridlocs_gf_sp1)
  gf_cov = matern_cov_fun(grid_gf_dist,thinning_pars$grf_pars$range,
                          thinning_pars$grf_pars$lambda,thinning_pars$grf_pars$sigma2)
  sp_gf_cov = as.spam(gf_cov)
  ## Cholesky decomposition of gf_cov
  A_m <<- as.matrix(spam::chol(sp_gf_cov))
  window1 <<- window
  
  cl <- makeCluster(8)
  
  
  clusterEvalQ(cl, {
   library(pbapply)
    library(raster)
    library(sp)
    library(rgeos)
  })
  
  
  
  # Export the necessary variables to the cluster
  clusterExport(cl, list("generate_data_cs","seeds","A_m","thinning_pars","thinning_cov","gridlocs_gf_sp1","st_area","window1","res_gf","eco_process","window1"))
  
  library(pbapply)
  print("Starting parallel")
  data = pblapply(1:length(eco_process$data),function(x){
    tmp <- generate_data_cs(1,thinning_pars,thinning_cov,gridlocs_gf_sp1,st_area,res_gf,window1)
    return(tmp)
  },cl=cl)
  
  stopCluster(cl)
  
  return(tmp)
}