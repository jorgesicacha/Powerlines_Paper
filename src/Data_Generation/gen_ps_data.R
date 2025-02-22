## Generating professional surveys data
library(deldir)
library(sp)
#library(rgeos)
## First I need to convert the existing powerlines into non-overlapping polygons

## This is a bit difficult though, but it is based on a voronoi tesselation

extract_intersections = function(geo_info,dens_size,grid_win_size){
  pl_sp = geo_info$buffered_pls
  window = extent(pl_sp)
  st_outw = data.frame(xmin = window@xmin,
                       xmax = window@xmax,
                       ymin = window@ymin,
                       ymax = window@ymax)
  
  ## Now, I'll generate a grid to get few points that determine the tesellation
  pl_xseq = seq(st_outw$xmin-.5*grid_win_size,st_outw$xmax+.5*grid_win_size,grid_win_size)
  pl_yseq = seq(st_outw$ymin-.5*grid_win_size,st_outw$ymax+.5*grid_win_size,grid_win_size)
  xy=expand.grid(pl_xseq,pl_yseq)
  pl_grid = SpatialPointsDataFrame(xy,
                                   data=data.frame(ID=1:nrow(xy)),
                                   proj4string = crs(pl_sp))
  gridded(pl_grid) = TRUE
  grid <- as(pl_grid, "SpatialPolygons")
  print("Fine grid generated!")
  pts_0 = st_segmentize(st_as_sf(geo_info$powerlines),dens_size)
  pts_0_1 = lapply(1:nrow(pts_0),function(x){
    ll = lapply(1:length(pts_0[x,]$geometry[[1]]),function(y){
      #print(x$geometry[[1]][[y]])
      #print((y))
      df_temp = data.frame(x=pts_0[x,]$geometry[[1]][[y]][,1],
                           y=pts_0[x,]$geometry[[1]][[y]][,2],
                           segm = y
      )
      return(df_temp)
    })
    return(do.call('rbind',ll))
  })
  pts_0_1 = do.call('rbind',pts_0_1)
  pts_sp = SpatialPointsDataFrame(pts_0_1[,c('x','y')],
                                  data=data.frame(segm=pts_0_1$segm),
                                  proj4string = crs(geo_info$buffered_pls))
  
  oo = over(pts_sp,grid)
  print("Done with over pts_sp/grid")
  tab_oo = table(oo)
  idx_check = tab_oo[which(tab_oo>1)]
  
  grid_centers_inters = xy[as.numeric(names(idx_check)),] ## These are the proposed intersection centroids
  ## Refinement
  test=spDists(SpatialPoints(grid_centers_inters,proj4string = crs(pl_sp)))
  diag(test)=999
  check = which(test<.01,arr.ind = T)
  rm_check = apply(check,1,function(x){x[1]<x[2]})
  check1 = ifelse(length(which(rm_check))==1,t(as.matrix(check[which(rm_check),])),
                  check[which(rm_check),])
  final_intersections = SpatialPointsDataFrame(grid_centers_inters[-check1,],
                                               data=data.frame(ID=1:(nrow(test)-length(which(rm_check)))),
                                               proj4string = crs(pl_sp))
  return(final_intersections)
}

gen_prof_survey_sampling_areas = function(geo_info,dens_size,dens_size_vor,pt_buff_size,grid_win_size){
  powerlines = geo_info$powerlines
  buffers = geo_info$buffered_pls
  ## We begin by getting the intersection of the existing lines
  final_intersections = extract_intersections(geo_info,dens_size =dens_size,grid_win_size =grid_win_size)
  print("Intersections found!")
  ## Then, we buffer the intersection points
  buff_intersections = as(st_buffer(st_as_sf(final_intersections),dist = pt_buff_size
  ),"Spatial")
  ## I remove the segments within the buffered point of the line object
  #tt_diff = gDifference(powerlines,buff_intersections)
  #tt_diff <- st_difference(st_as_sf(powerlines), st_as_sf(buff_intersections[1,]))
  powerlines_sf <- st_as_sf(powerlines)
  buff_sf <- st_as_sf(buff_intersections)
  
  # Sequential difference
  tt_diff <- powerlines_sf
  for (i in 1:nrow(buff_sf)) {
    tt_diff <- st_difference(tt_diff, buff_sf[i, ])
  }
  tt_diff <- tt_diff[, c("ID", "geometry")]
  ## Then, I densify the resulting object and create the voronoi tesselation
  pts_0 = st_segmentize(st_as_sf(tt_diff),dens_size_vor)
  pts_0_1 = lapply(1:nrow(pts_0),function(x){
    print(x)
    ll = lapply(1:length(pts_0[x,]$geometry[[1]]),function(y){
      df_temp = data.frame(x=pts_0[x,]$geometry[[1]][[y]][,1],
                           y=pts_0[x,]$geometry[[1]][[y]][,2],
                           segm = y
      )
      return(df_temp)
    })
    return(do.call('rbind',ll))
  })
  pts_0_1 = do.call('rbind',pts_0_1)
  pts_sp = SpatialPointsDataFrame(pts_0_1[,c('x','y')],
                                  data=data.frame(segm=pts_0_1$segm),
                                  proj4string = crs(pl_sp))
  
  ## Now, it's time for voronoi Tesselation
  teselation <- deldir(pts_sp@coords[,1], pts_sp@coords[,2])
  tiles <- tile.list(teselation)
  polys = vector(mode='list', length=length(tiles))
  require(sp)
  for (i in seq(along=polys)) {
    pcrds = cbind(tiles[[i]]$x, tiles[[i]]$y)
    pcrds = rbind(pcrds, pcrds[1,])
    polys[[i]] = Polygons(list(Polygon(pcrds)), ID=as.character(i))
  }
  SP = SpatialPolygons(polys,proj4string=crs(pl_sp))
  voronoi = SpatialPolygonsDataFrame(SP, data=data.frame(x=pts_sp@coords[,1], 
                                                         y=pts_sp@coords[,2],
                                                         segm = pts_sp$segm,
                                                         row.names=sapply(slot(SP, 'polygons'), 
                                                                          function(x) slot(x, 'ID')))
  )
  ## Finally, we intersect this tesselation with the buffer we had :)
  vor_int = raster::intersect(voronoi,pl_sp)
  
  library(dplyr)
  
  vor_int_f <- st_as_sf(vor_int) %>%
    group_by(segm) %>%
    summarise(geometry = st_union(geometry))
  
  vor_int_f <- do.call("rbind",lapply(1:nrow(vor_int_f),function(x){st_cast(vor_int_f[x,],"POLYGON")}))
  vor_int_f$segm <- seq(1,nrow(vor_int_f))
  
  #vor_int_f =gUnaryUnion(vor_int,id=vor_int$segm)
  vor_int_f <- as(vor_int_f, "Spatial")
  
  vor_int_f$id = 1:length(vor_int_f)
  print("Professional sampling areas identified")
  return(vor_int_f)  
}

gen_prof_survey_data = function(mesh,std_area,eco_final_pp,eco_grf_sp,prof_survey_sampling_areas,prof_sampling_pars,seed,res){
  ## It all starts with a covariate
  true_eco_pp = eco_final_pp
  set.seed(seed)
  window = extent(std_area)
  st_outw = data.frame(xmin = window@xmin - 10,
                       xmax = window@xmax + 10,
                       ymin = window@ymin - 10,
                       ymax = window@ymax + 10)
  
  xgrid = seq(st_outw$xmin + 0.5*res,st_outw$xmax-0.5*res,res)
  ygrid = seq(st_outw$ymin + 0.5*res,st_outw$ymax-0.5*res,res)
  gridlocs = expand.grid(xgrid,ygrid)
  gridcov = apply(gridlocs,1,function(x){cos((x[2]+x[1]))})
  prof_samp_cov.sp = SpatialPointsDataFrame(gridlocs,data=data.frame(prof_samp_cov=gridcov))
  crs(prof_samp_cov.sp)= crs(pl_sp)
  r <- raster(prof_samp_cov.sp)
  r1<-disaggregate(r, fact=res(r)/c(res,res))
  ps_cov.rast = rasterize(prof_samp_cov.sp@coords,r1,prof_samp_cov.sp$prof_samp_cov, fun=mean,na.rm=T)
  print("Professional Sampling Covariate generated")
  ## The base for getting all we need
  grid_spdf = eco_grf_sp
  grid_spdf$fixed_int = rep(1,nrow(grid_spdf))
  grid_spdf$fixed_cov = extract(ps_cov.rast,grid_spdf)

  ## Now that information is in place, I need to generate the parameter phi at each area
  ##I need to get the same I have in grid_spdf, but for the mesh points, I'll do it based on a NN approach

  meshpoints <- mesh$loc[,1:2]
  points.mesh <- SpatialPoints(coords=cbind(meshpoints[,1],meshpoints[,2]),
                               proj4string =crs(st_area))

  NN_idx = sapply(1:mesh$n,function(x){
    #idx = which.min(spDistsN1(grid_spdf,points.mesh[x]))
    idx = spDistsN1(grid_spdf,points.mesh[x])
    return(idx)
  })

  NN_idx_1 = apply(NN_idx,2,order,decreasing=F)
  ## First choice
  first_nn_idx = apply(NN_idx,2,which.min)

  first_mesh_grid_df = data.frame(z=grid_spdf$z[first_nn_idx],
                                  zgf=grid_spdf$zgf[first_nn_idx],
                                  lambda=grid_spdf$lambda[first_nn_idx],
                                  fixed_int = grid_spdf$fixed_int[first_nn_idx],
                                  fixed_cov = grid_spdf$fixed_cov[first_nn_idx])

  vars_0 = apply(first_mesh_grid_df,2,function(x){
    length(which(is.na(x)))
  })
  check_vars = names(vars_0[which(vars_0>0)])

  for(x in check_vars){
    tmp = which(is.na(first_mesh_grid_df[x]))
    for(y in tmp){
      ii=0
      jj=1
      while(ii==0){
        tmp_val=grid_spdf@data[x][NN_idx_1[jj,y],1]
        if(!is.na(tmp_val)){
          first_mesh_grid_df[x][y,1]=tmp_val
          ii=ii+1
        }
        jj=jj+1
      }
    }
  }

  print("Information at mesh calculated")
  points.mesh <- SpatialPointsDataFrame(points.mesh,data=first_mesh_grid_df)

  mi_0 <- over(points.mesh,prof_survey_sampling_areas)

  mesh.intersect <- points.mesh[which(!is.na(mi_0$id)),]
  mi_1 <- mi_0$id[which(!is.na(mi_0$id))]
  tmi1 = table(mi_1)

  ## This is a validation step
  if(length(which(is.na(tmi1)))==0){print("Approximation based on mesh is OK!")}
  else{print("Check the resolution of the mesh!")}
  #print(head(mi_1))
  mesh.intersect.id <- mi_1
  print("a1")
  print(mesh$n)
  print(class(mesh.intersect))
  print(mesh.intersect.id)
  A.area <- INLA::inla.spde.make.A(mesh = mesh, loc = mesh.intersect,
                                   block = mesh.intersect.id, block.rescale = "none")
  A.area <- as.matrix(A.area)
  A.area <- t(apply(A.area,1,function(x){
    x[x>0]  =1
    x = x/sum(x)
    return(x)
  }))
  A.area = Matrix(A.area)

  print("a2")
  logit <- function(y){log(y/(1-y))}
  ilogit <- function(y){plogis(y)}
  print("a")
  prof_survey_sampling_areas$pref_linpreds <- logit(as.numeric(A.area%*%ilogit(prof_sampling_pars$gamma_0*points.mesh$fixed_int+
                                                                                 prof_sampling_pars$gamma_1*points.mesh$fixed_cov+
                                                                                 prof_sampling_pars$pref_samp_par*points.mesh$zgf)))

  prof_survey_sampling_areas$samp_prob <- as.numeric(A.area%*%ilogit(prof_sampling_pars$gamma_0*points.mesh$fixed_int+
                                                                       prof_sampling_pars$gamma_1*points.mesh$fixed_cov+
                                                                       prof_sampling_pars$pref_samp_par*points.mesh$zgf))
  print("b")
  prof_survey_sampling_areas$ps_cov <- as.numeric(A.area%*%(points.mesh$fixed_cov))
  prof_survey_sampling_areas$zgf <- as.numeric(A.area%*%(points.mesh$zgf))
  # prof_survey_sampling_areas$ps_cov <- as.numeric(A.area%*%ilogit(points.mesh$fixed_cov))
  # prof_survey_sampling_areas$zgf <- as.numeric(A.area%*%ilogit(points.mesh$zgf))
  
  prof_survey_sampling_areas$which_lines <- sapply(prof_survey_sampling_areas$samp_prob,function(x){rbinom(1,1,x)})

  prof_sf <- st_as_sf(prof_survey_sampling_areas)
  survey_sf <- st_as_sf(prof_surveys_data)
  # Plot using ggplot2
  ggplot() +
    geom_sf(data = prof_sf, aes(fill = as.factor(which_lines)), color = "black") +  # Polygons
    geom_sf(data = survey_sf, aes(color = as.factor(which_lines)), size = 2) +     # Points
    scale_fill_brewer(palette = "Set1", name = "Polygon Lines") +                  # Polygon legend
    scale_color_brewer(palette = "Dark2", name = "Point Lines") +                 # Point legend
    theme_minimal() +
    labs(title = "Survey Data with Sampling Areas",
         x = "Easting (km)", y = "Northing (km)")
  
  prof_surveys_data = raster::intersect(true_eco_pp,prof_survey_sampling_areas)
  prof_surveys_data = prof_surveys_data[which(prof_surveys_data$which_lines==1),]
  
  
  print("c")
  return(list(ps_cov.rast=ps_cov.rast,
              prof_survey_sampling_areas=prof_survey_sampling_areas,
              prof_surveys_data=prof_surveys_data))
}




## Just testing out all
# grid_spdf$pref_linpred = prof_sampling_pars$gamma_0*grid_spdf$fixed_int+ prof_sampling_pars$gamma_1*grid_spdf$fixed_cov+ prof_sampling_pars$pref_samp_par*grid_spdf$zgf
# grid_spdf$samp_prob = ilogit(prof_sampling_pars$gamma_0*grid_spdf$fixed_int+ prof_sampling_pars$gamma_1*grid_spdf$fixed_cov+ prof_sampling_pars$pref_samp_par*grid_spdf$zgf)
# 
# r <- raster(grid_spdf)
# r1<-disaggregate(r, fact=res(r)/c(res,res))
# ps_linpred.rast = rasterize(grid_spdf@coords,r1,grid_spdf$pref_linpred, fun=mean,na.rm=T)
# ps_sampprob.rast = rasterize(grid_spdf@coords,r1,grid_spdf$samp_prob, fun=mean,na.rm=T)
# 
# ##
# spplot(prof_survey_sampling_areas['samp_prob'],zcol=terrain.colors(256))
# 
# gen_prof_surveys_data = function(geo_info,){
#   
#   
#   
#   
# }
# 
# class(pl_sp@polygons[[1]]@Polygons[[1]])
# 
# ## Getting the extent of the powerlines
# 
# ## First I'll make a dense grid of points so that I can determine the intersections of 
# ## my powerlines
# window = extent(pl_sp)
# st_outw = data.frame(xmin = window@xmin,
#                      xmax = window@xmax,
#                      ymin = window@ymin,
#                      ymax = window@ymax)
# 
# ## Now, I'll generate a grid to get few points that determine the tesellation
# pl_xseq = seq(st_outw$xmin-.5*.007,st_outw$xmax+.5*.007,.007)
# pl_yseq = seq(st_outw$ymin-.5*.007,st_outw$ymax+.5*.007,.007)
# xy=expand.grid(pl_xseq,pl_yseq)
# pl_grid = SpatialPointsDataFrame(xy,
#                                  data=data.frame(ID=1:nrow(xy)),
#                                  proj4string = crs(pl_sp))
# gridded(pl_grid) = TRUE
# grid <- as(pl_grid, "SpatialPolygons")
# 
# pts_0 = st_segmentize(st_as_sf(geo_info$powerlines),.01)
# pts_0_1 = lapply(1:nrow(pts_0),function(x){
# print(x)
# ll = lapply(1:length(pts_0[x,]$geometry[[1]]),function(y){
#   #print(x$geometry[[1]][[y]])
#   #print((y))
#   df_temp = data.frame(x=pts_0[x,]$geometry[[1]][[y]][,1],
#                        y=pts_0[x,]$geometry[[1]][[y]][,2],
#                        segm = y
#                        )
# return(df_temp)
# })
# return(do.call('rbind',ll))
# })
# pts_0_1 = do.call('rbind',pts_0_1)
# pts_sp = SpatialPointsDataFrame(pts_0_1[,c('x','y')],
#                                 data=data.frame(segm=pts_0_1$segm),
#                                 proj4string = crs(pl_sp))
# 
# oo = over(pts_sp,grid)
# tab_oo = table(oo)
# idx_check = tab_oo[which(tab_oo>1)]
# 
# grid_centers_inters = xy[as.numeric(names(idx_check)),] ## These are the proposed intersection centroids
# ## Refinement
# test=spDists(SpatialPoints(grid_centers_inters,proj4string = crs(pl_sp)))
# diag(test)=999
# check = which(test<.01,arr.ind = T)
# rm_check = apply(check,1,function(x){x[1]<x[2]})
# check1 = ifelse(length(which(rm_check))==1,t(as.matrix(check[which(rm_check),])),
#                check[which(rm_check),])
# final_intersections = SpatialPointsDataFrame(grid_centers_inters[-check1,],
#                                              data=data.frame(ID=1:(nrow(test)-length(which(rm_check)))),
#                                     proj4string = crs(pl_sp))
# 
# pts_1 = st_segmentize(st_as_sf(geo_info$powerlines),3)
# pts_1_1 = lapply(1:nrow(pts_1),function(x){
#   print(x)
#   ll = lapply(1:length(pts_1[x,]$geometry[[1]]),function(y){
#     #print(x$geometry[[1]][[y]])
#     #print((y))
#     df_temp = data.frame(x=pts_1[x,]$geometry[[1]][[y]][,1],
#                          y=pts_1[x,]$geometry[[1]][[y]][,2]
#     )
#     return(df_temp)
#   })
#   return(do.call('rbind',ll))
# })
# pts_1_1 = SpatialPoints(do.call('rbind',pts_1_1),
#                                  proj4string = crs(pl_sp))
# ## Refining centroids
# rm_idxs = lapply(1:nrow(final_intersections),function(x){
#   tmp_dist = spDistsN1(pts_1_1,final_intersections[x,])
#   which_rm = which(tmp_dist<.1)
#   return(which_rm)
# })
# rm_idxs = do.call('c',rm_idxs)
# rm_idxs = unique(rm_idxs)
# 
# voronoi_pts = rbind(pts_1_1[-rm_idxs],final_intersections)
# 
# plot(pl_sp)
# points(voronoi_pts,pch=19,cex=.5,col="red")
# 
# teselation <- deldir(voronoi_pts@coords[,1], voronoi_pts@coords[,2])
# tiles <- tile.list(teselation)
# polys = vector(mode='list', length=length(tiles))
# require(sp)
# for (i in seq(along=polys)) {
#   pcrds = cbind(tiles[[i]]$x, tiles[[i]]$y)
#   pcrds = rbind(pcrds, pcrds[1,])
#   polys[[i]] = Polygons(list(Polygon(pcrds)), ID=as.character(i))
# }
# SP = SpatialPolygons(polys,proj4string=crs(pl_sp))
# voronoi = SpatialPolygonsDataFrame(SP, data=data.frame(x=voronoi_pts@coords[,1], 
#                                                        y=voronoi_pts@coords[,2], row.names=sapply(slot(SP, 'polygons'), 
#                                                                                            function(x) slot(x, 'ID')))
#                                    )
# 
# rr = raster::intersect(voronoi,pl_sp)
# plot(rr)
# points(voronoi_pts,col="red",pch=19,cex=.5)
# over(rr,voronoi_pts)
# plot(pl_sp)
# for(i in 1:length(rr)){
#   plot(rr[i,],col=i,add=T)
# }
# plot(voronoi,add=T)
