## Generating the GRF and the LGCP for the dead bird locations ##
library(spam)
library(raster)
library(doParallel)
library(sf)

matern_cov_fun = function(dist,range,lambda=1,sigma2){ ##Following definition in Blangiardo
  kappa = sqrt(8)/range
  ifelse(dist>0,(sigma2/(gamma(lambda)*(2^(lambda-1))))*((kappa*dist)^lambda)*besselK(kappa*dist,lambda),sigma2)
}

generate_data <- function(seed,eco_fixed_pars,cov1.rast,gridlocs_gf_sp1,std_area,window,res_gf,pl_sp,proj_crs){
  set.seed(seed)
  z0 <- matrix(rnorm(nrow(A_m)),ncol=1)
  z = eco_fixed_pars$beta0 + eco_fixed_pars$beta1*extract(cov1.rast,gridlocs_gf_sp1) + A_m%*%z0
  zgf = A_m%*%z0

  dft = data.frame(z =rep(NA,length(gridlocs_gf_sp1)),
                   zgf=rep(NA,length(gridlocs_gf_sp1)))
  dft$z = z
  dft$zgf = zgf

  grf_sp = SpatialPointsDataFrame(gridlocs_gf_sp1,data=dft)
  r <- raster(grf_sp)
  r1<-disaggregate(r, fact=res(r)/c(res_gf,res_gf))
  grf.rast = rasterize(grf_sp@coords,r1,grf_sp$z, fun=mean,na.rm=T)

  ## Idea of a "baseline" intensity
  grf_sp$lambda = exp(grf_sp$z)
  grf_sp_vals = grf_sp[complete.cases(grf_sp@data),]
  ## Crop the raster and try again
  #grf.rast.crop = crop(grf.rast,st_area)
  grf_sp_lambda = raster::intersect(grf_sp,std_area)
  r <- raster(grf_sp_lambda)
  r1<-disaggregate(r, fact=res(r)/c(res_gf,res_gf))
  lambda.rast = rasterize(grf_sp_lambda@coords,r1,as.numeric(grf_sp_lambda$lambda), fun=mean,na.rm=T)

  lambda_max = max(as.numeric(grf_sp_lambda$lambda,na.rm = T))
  N_max = rpois(1,lambda_max*st_area(st_as_sf(pl_sp)))
  ## Randomly localize the Nmax points in the buffers
  pts = list()
  sp_locs = SpatialPoints(coords = matrix(c(runif(1.2*N_max*st_area(st_as_sf(std_area))/st_area(st_as_sf(pl_sp)),window@xmin,window@xmax),
                                            runif(1.2*N_max*st_area(st_as_sf(std_area))/st_area(st_as_sf(pl_sp)),window@ymin,window@ymax)),ncol=2))
  crs(sp_locs)<- crs(pl_sp)
  pts_tmp = raster::intersect(sp_locs,pl_sp)
  S_max = pts_tmp[1:N_max]

  # Finally, the points of S_max will be thinned according to their intensity
  lambda_S_max = extract(lambda.rast,S_max)
  ## Estimating lambda if missing: Nearest cell
  lambda_S_max[is.na(lambda_S_max)] <-  grf_sp_vals$lambda[sapply(which(is.na(lambda_S_max)),function(x){
    which.min(spDistsN1(grf_sp_vals@coords,cbind(S_max@coords[x,1],S_max@coords[x,2])))
  })]

  thinprobs = lambda_S_max/lambda_max
  S_max = SpatialPointsDataFrame(S_max,data=data.frame(ret=sapply(thinprobs,function(x){rbinom(1,1,x)})))
  final_pp = S_max[which(S_max$ret==1),]

  return(list(grf_sp=grf_sp,final_pp=final_pp,lambda.rast=lambda.rast,
              grf.rast=grf.rast,seed=seed))
  
}


generate_eco_process = function(nsims,std_area,res,res_gf,pl_sp,grf_pars,eco_fixed_pars){
print("Starting here")
window <- extent(std_area)
st_outw = data.frame(xmin = window@xmin - 10,
                     xmax = window@xmax + 10,
                     ymin = window@ymin - 10,
                     ymax = window@ymax + 10)


st_outw_inner = data.frame(xmin = window@xmin - 1,
                     xmax = window@xmax + 1,
                     ymin = window@ymin - 1,
                     ymax = window@ymax + 1)

inner_window = SpatialPolygons(list(Polygons(list(Polygon(matrix(c(st_outw_inner$xmin,st_outw_inner$ymin,
                                                                 st_outw_inner$xmax,st_outw_inner$ymin,
                                                                 st_outw_inner$xmax,st_outw_inner$ymax,
                                                                 st_outw_inner$xmin,st_outw_inner$ymax,
                                                                 st_outw_inner$xmin,st_outw_inner$ymin
                                                                 ),ncol=2,byrow=T))),ID=1)),proj4string = crs(std_area))


xgrid = seq(st_outw$xmin + 0.5*res,st_outw$xmax-0.5*res,res)
ygrid = seq(st_outw$ymin + 0.5*res,st_outw$ymax-0.5*res,res)
gridlocs = expand.grid(xgrid,ygrid)
gridcov = apply(gridlocs,1,function(x){cos((x[2]-x[1]))})
proj_crs <<- "+proj=utm +zone=32 +ellps=GRS80 +units=km +no_defs"
eco_cov.sp = SpatialPointsDataFrame(gridlocs,data=data.frame(eco_cov=gridcov),proj4string = crs(std_area))
r <- raster(eco_cov.sp)
r1<-disaggregate(r, fact=res(r)/c(res,res))
cov1.rast <<- rasterize(eco_cov.sp@coords,r1,eco_cov.sp$eco_cov, fun=mean,na.rm=T)
print("Covariate generated")

## Now I have to simulate a Gaussian RF, the resolution has to be high
xgrid_gf = seq(st_outw$xmin + 0.5*res_gf,st_outw$xmax-0.5*res_gf,res_gf)
ygrid_gf = seq(st_outw$ymin + 0.5*res_gf,st_outw$ymax-0.5*res_gf,res_gf)
gridlocs_gf_sp = SpatialPoints(expand.grid(xgrid_gf,ygrid_gf),proj4string = crs(std_area))

## Op2
#idx_int = rgeos::gIntersects(gridlocs_gf_sp,inner_window,byid = TRUE)

idx_int <- st_intersects(st_as_sf(gridlocs_gf_sp), st_as_sf(inner_window), sparse = TRUE)
idx_int <- !is.na(idx_int)

gridlocs_gf_sp1 <<- raster::intersect(gridlocs_gf_sp,inner_window)
### the distance between the points in the grid
grid_gf_dist = spDists(gridlocs_gf_sp1)
gf_cov = matern_cov_fun(grid_gf_dist,grf_pars$range,grf_pars$lambda,grf_pars$sigma2)

sp_gf_cov = as.spam(gf_cov)
# Cholesky decomposition of gf_cov
A_m <<- as.matrix(spam::chol(sp_gf_cov))
window1 <<- window
print("window assigned")
print("OK?")
#print(class(A_m))
cl <- makeCluster(2)

clusterEvalQ(cl, {
library(pbapply)
library(raster)
library(sp)
library(sf)
#library(rgeos)
})
print("Issue here1?")
seeds <<- sample(1e6,nsims,replace = F)
std_area <<- std_area
# Export the necessary variables to the cluster
clusterExport(cl, list("generate_data","seeds","A_m","eco_fixed_pars","cov1.rast","gridlocs_gf_sp1","std_area","window1","seeds","res_gf","pl_sp","proj_crs"))
print("Issue here2?")

#model <- RMmatern(nu=grf_pars$lambda,var=grf_pars$sigma2, scale=1/grf_pars$range) 
## define the locations:
#simu <- RFsimulate(model, x=gridlocs_gf_sp1@coords)
#z = eco_fixed_pars$beta0 + eco_fixed_pars$beta1*extract(cov1.rast,gridlocs_gf_sp1) + simu$variable1
#zgf = simu$variable1
library(pbapply)
print("OK until here")

data = pblapply(seeds,function(x){

  tmp <- generate_data(x,eco_fixed_pars,cov1.rast,gridlocs_gf_sp1,std_area,window1,res_gf,pl_sp,proj_crs)
  return(tmp)
},cl=cl)

stopCluster(cl)

print("LGCP generated")
return(list(data=data,cov1.rast=cov1.rast))
}
