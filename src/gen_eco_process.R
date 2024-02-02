## Generating the GRF and the LGCP for the dead bird locations ##
library(spam)
library(rgeos)
matern_cov_fun = function(dist,range,lambda=1,sigma2){ ##Following definition in Blangiardo
  kappa = sqrt(8)/range
  ifelse(dist>0,(sigma2/(gamma(lambda)*(2^(lambda-1))))*((kappa*dist)^lambda)*besselK(kappa*dist,lambda),sigma2)
}

generate_eco_process = function(st_area,res,res_gf,pl_sp,grf_pars,eco_fixed_pars,seed){
set.seed(seed)
window = extent(st_area)
st_outw = data.frame(xmin = window@xmin - 10,
                     xmax = window@xmax + 10,
                     ymin = window@ymin - 10,
                     ymax = window@ymax + 10)

xgrid = seq(st_outw$xmin + 0.5*res,st_outw$xmax-0.5*res,res)
ygrid = seq(st_outw$ymin + 0.5*res,st_outw$ymax-0.5*res,res)
gridlocs = expand.grid(xgrid,ygrid)
gridcov = apply(gridlocs,1,function(x){cos((x[2]-x[1]))})
proj_crs = "+proj=utm +zone=32v +ellps=GRS80 +units=km +no_defs"
eco_cov.sp = SpatialPointsDataFrame(gridlocs,data=data.frame(eco_cov=gridcov),proj4string = crs(proj_crs))
r <- raster(eco_cov.sp)
r1<-disaggregate(r, fact=res(r)/c(res,res))
cov1.rast = rasterize(eco_cov.sp@coords,r1,eco_cov.sp$eco_cov, fun=mean,na.rm=T)
print("Covariate generated")
## Now I have to simulate a Gaussian RF, the resolution has to be high
xgrid_gf = seq(st_outw$xmin + 0.5*res_gf,st_outw$xmax-0.5*res_gf,res_gf)
ygrid_gf = seq(st_outw$ymin + 0.5*res_gf,st_outw$ymax-0.5*res_gf,res_gf)
gridlocs_gf_sp = SpatialPoints(expand.grid(xgrid_gf,ygrid_gf),proj4string = crs(proj_crs))

## Op2
idx_int = rgeos::gIntersects(gridlocs_gf_sp,pl_sp,byid = TRUE)
gridlocs_gf_sp1 = raster::intersect(gridlocs_gf_sp,pl_sp)
### the distance between the points in the grid
grid_gf_dist = spDists(gridlocs_gf_sp1)
gf_cov = matern_cov_fun(grid_gf_dist,grf_pars$range,grf_pars$lambda,grf_pars$sigma2)
sp_gf_cov = as.spam(gf_cov)
## Cholesky decomposition of gf_cov
A = as.matrix(spam::chol(sp_gf_cov))
z0 = matrix(rnorm(nrow(A)),ncol=1)
z = eco_fixed_pars$beta0 + eco_fixed_pars$beta1*extract(cov1.rast,gridlocs_gf_sp1) + A%*%z0
zgf = A%*%z0
print("Done with GRF")
dft = data.frame(z =rep(NA,length(gridlocs_gf_sp)),
                 zgf=rep(NA,length(gridlocs_gf_sp)))
dft$z[idx_int] = z
dft$zgf[idx_int] = zgf
grf_sp = SpatialPointsDataFrame(gridlocs_gf_sp,data=dft)
r <- raster(grf_sp)
r1<-disaggregate(r, fact=res(r)/c(res_gf,res_gf))
grf.rast = rasterize(grf_sp@coords,r1,grf_sp$z, fun=mean,na.rm=T)
## GRF generated
#mapView(grf.rast)+mapview(gridlocs_gf_sp1,cex=0.1)
print(" Done with GRF raster")
## Idea of a "baseline" intensity
grf_sp$lambda = exp(grf_sp$z)
grf_sp_vals = grf_sp[complete.cases(grf_sp@data),]
lambda.rast = rasterize(grf_sp@coords,r1,grf_sp$lambda, fun=mean,na.rm=T)
print (" Done with lambda raster")
#mapView(lambda.rast)
### Simulate as an NHPP since conditional on lambda, a LGCP is a NHPP
lambda_max = max(grf_sp$lambda,na.rm = T)
N_max = rpois(1,lambda_max*gArea(pl_sp))
## Randomly localize the Nmax points in the buffers
pts = list()
sp_locs = SpatialPoints(coords = matrix(c(runif(1.2*N_max*gArea(st_area)/gArea(pl_sp),window@xmin,window@xmax),
                                            runif(1.2*N_max*gArea(st_area)/gArea(pl_sp),window@ymin,window@ymax)),ncol=2),
                              proj4string=crs(proj_crs))
pts_tmp = raster::intersect(sp_locs,pl_sp)
S_max = pts_tmp[1:N_max]  

## Finally, the points of S_max will be thinned according to their intensity
lambda_S_max = extract(lambda.rast,S_max)
## Estimating lambda if missing: Nearest cell
lambda_S_max[is.na(lambda_S_max)] <-  grf_sp_vals$lambda[sapply(which(is.na(lambda_S_max)),function(x){
  which.min(spDistsN1(grf_sp_vals@coords,cbind(S_max@coords[x,1],S_max@coords[x,2])))
})]

thinprobs = lambda_S_max/lambda_max
S_max = SpatialPointsDataFrame(S_max,data=data.frame(ret=sapply(thinprobs,function(x){rbinom(1,1,x)})))
final_pp = S_max[S_max$ret==1,]
print("LGCP generated")
return(list(grf_sp=grf_sp,final_pp=final_pp,cov1.rast=cov1.rast,lambda.rast=lambda.rast,
            grf.rast=grf.rast))
}