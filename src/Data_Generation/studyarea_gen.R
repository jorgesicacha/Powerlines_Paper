## Devicing a simple example for data fusion and power lines ""
library(sp)
library(raster)
library(mapview)
### The study area ... a regular area ##
gen_area_pls = function(win_seed,win_size,n_centers,n_powerlines_per_center,buff_size,mode){
  
  proj_crs = "+proj=utm +zone=32 +ellps=GRS80 +units=km +no_defs"
  x0 = win_seed[1]
  y0 = win_seed[2]
  x1 = x0 + win_size
  y1 = y0 + win_size
  coordsmat <- matrix(c(x0,y0,x1,y0,x1,y1,x0,y1,x0,y0),ncol=2,byrow=T)
  aa <- SpatialPolygons(list(Polygons(list(Polygon(coordsmat)),ID=1)),proj4string = raster::crs(proj_crs))
  #mapView(aa)
  #plot(aa)
  
  ### A simple task: designing a powerline network within the study area 
  #n_centers = 4
  centers = list()
  for(i in 1:n_centers){
    centers[[i]] = data.frame(x = runif(1,x0,x1),
                              y=runif(1,y0,y1))
  }
  
  ## Generate lines in 8 directions
  
  #n_powerlines_per_center = 4
  
  lines_per_center = function(center,n_powerlines_per_center){
    ref_df = data.frame(dist=c(0,runif(2*n_powerlines_per_center,0,5)),
                        angle=c(0,seq(0,((2*n_powerlines_per_center)-1)*pi/n_powerlines_per_center,pi/n_powerlines_per_center)),
                        line=c(0,rep(1:n_powerlines_per_center,2)))
    #center = data.frame(x=5,y=5)
    
    lines_nodes =do.call('rbind',apply(ref_df,1,function(x){
      pars = data.frame(dist=x[1],angle=x[2])
      a = center$x + pars$dist*cos(pars$angle)
      b = center$y + pars$dist*sin(pars$angle)
      
      #return(data.frame(x=a,y=b))  
      #}))
      if(!(a>=x0 & a<=x1) & !(b>=y0 & b<=y1)){
        if(a>x0 & b>y0){
          corr_dists = c(a-x1,b-y1)
          tops = c(x1,y1)
          d_corr = (tops[which.max(corr_dists)]-center[which.max(corr_dists)][[1]])/sin(pars$angle)
          a = center$x + d_corr*cos(pars$angle)
          b = center$y + d_corr*sin(pars$angle)
        }
        else{
          if(a>x0 & b<y0){
            corr_dists = c(a-x1,y0-b)
            d_corr = ifelse(corr_dists[1]==max(corr_dists),
                            (x1-center$x)/sin(-pars$angle),
                            (center$y-y0)/sin(-pars$angle))
            a = center$x + d_corr*cos(pars$angle)
            b = center$y + d_corr*sin(pars$angle)        
          }
          else{
            if(a<x0 & b>y0){
              corr_dists = c(x0-a,b-y1)
              d_corr = ifelse(corr_dists[2]==max(corr_dists),
                              (y1-center$y)/sin(pars$angle),(center$x-x0)/sin(pars$angle))
              a = center$x + d_corr*cos(pars$angle)
              b = center$y + d_corr*sin(pars$angle)  
            }
            else{
              if(a<x0 & b<y0){
                corr_dists = c(x0-a,y0-b)
                tops=c(x0,y0)
                d_corr = (center[which.max(corr_dists)][[1]] - tops[which.max(corr_dists)])/sin(-pars$angle)
                a = center$x + d_corr*cos(pars$angle)
                b = center$y + d_corr*sin(pars$angle)            
              }
            }}
        }
      }
      else{
        if(!(a>=x0 & a<=x1)){
          if(a>x0){
            d_corr = ifelse(sin(pars$angle)==0,x1-center$x,(x1-center$x)/sin(-pars$angle))
            a = center$x + d_corr*cos(pars$angle)
            b = center$y + d_corr*sin(pars$angle)
          }
          else{
            d_corr = ifelse(round(sin(pars$angle),10)==0,center$x-x0,ifelse(sin(pars$angle)>0,(center$x-x0)/sin(pars$angle),(center$x-x0)/sin(-pars$angle)))
            a = center$x + d_corr*cos(pars$angle)
            b = center$y + d_corr*sin(pars$angle)   
          }
        }
        else{
          if(!(b>=y0 & b<=y1)){
            if(b>y0){
              d_corr = (y1-center$y)/sin(pars$angle)
              a = center$x + d_corr*cos(pars$angle)
              b = center$y + d_corr*sin(pars$angle)      
            }
            else{
              d_corr = (center$y-y0)/sin(-pars$angle)
              a = center$x + d_corr*cos(pars$angle)
              b = center$y + d_corr*sin(pars$angle)    
            }
          }
        }
      }
      return(data.frame(x=a,y=b))
    }))
    ref_df = cbind(ref_df,lines_nodes)
    return(ref_df)
  }
  
  LinesList = lapply(1:n_centers,function(i){
    a1 = lines_per_center(centers[[i]],n_powerlines_per_center)
    ID = as.character(i)
    Lines = lapply(1:n_powerlines_per_center,function(j){Line(a1[which(a1$line==j),c("x","y")])})
    return(Lines(Lines,ID))
  })
  
  powerlines = SpatialLines(LinesList,proj4string = crs(proj_crs))
  #mapView(aa,alpha.regions=0.1) + mapView(powerlines)
  
  ## Creating the buffers
  library(raster)
  if(mode =="powerlines"){
    buff_pl = buffer(powerlines,buff_size)
  }
  if(mode =="roads"){
    buff_pl = powerlines  
  }
  #mapView(aa,alpha.regions=0.1) + mapView(powerlines) + mapview(buff_pl)
  return(list(st_area = aa,powerlines = powerlines,buffered_pls = buff_pl))
}
