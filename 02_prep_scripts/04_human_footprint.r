###################################################################### #
# Calculating human footprint metrics that can be summed over species  #
###################################################################### #

# set options ----
calc.delta <- FALSE
calc.even <- FALSE
calc.quant <- TRUE

# load packages ----

library(tidyverse)
library(sf)
library(rWCVP)
library(cli)
library(raster)
library(stars)

# prevent spherical geometry errors
sf_use_s2(FALSE)

# load tdwg data ----

# tdwg regions
tdwg <- rWCVPdata::wgsrpd3 


#extra teeny buffer
eps = sqrt(.Machine$double.eps)*2

#crop to clean up
tdwg <- tdwg %>% 
  st_crop(xmin=-180+eps, xmax=180-eps, ymin=-90+eps, ymax=83.62361) %>% 
  st_transform(crs="ESRI:54009")


# load HFP data and aggregate from raw (if files do not alredy exist) ----

# 1993 footprint
if(file.exists("01_raw_data/human_footprint/aggregated/hfp_1993.tif")){
  hfp1 <- raster("01_raw_data/human_footprint/aggregated/hfp_1993.tif")
} else {
  
  hfp1 <- raster("01_raw_data/human_footprint/HFP1993.tif") %>% 
    aggregate(fact=5)
  writeRaster(hfp1, "01_raw_data/human_footprint/aggregated/hfp_1993.tif")
}

#2009 footprint
if(file.exists("01_raw_data/human_footprint/aggregated/hfp_2009.tif")){
  hfp2 <- raster("01_raw_data/human_footprint/aggregated/hfp_2009.tif")
} else {
  hfp2 <- raster("01_raw_data/human_footprint/HFP2009.tif")%>% 
    aggregate(fact=5)
  writeRaster(hfp2, "01_raw_data/human_footprint/aggregated/hfp_2009.tif")
}

#calculate change in hfp
if(file.exists("01_raw_data/human_footprint/aggregated/hfp_change.tif")){
  hfp_change <- raster("01_raw_data/human_footprint/aggregated/hfp_change.tif")
} else {
  hfp_change <- hfp2-hfp1
  writeRaster(hfp_change, "01_raw_data/human_footprint/aggregated/hfp_change.tif")
}

############################################################################  #

# Change in HFP between 1993-2009 ----
if(calc.delta==TRUE){
# _bin the raster ----
breaks <- c(-30,-5,-1,1,5,30)
hfp_change_bins <- cut(hfp_change,breaks=breaks)

#convert to sf and merge all polygons with same value
if(file.exists("01_raw_data/human_footprint/sf/hfp_delta.rds")){
  hfp_delta <- readRDS("01_raw_data/human_footprint/sf/hfp_delta.rds")
} else { 
hfp_delta <- st_as_sf(st_as_stars(hfp_change_bins), as_points = FALSE, merge=TRUE)
hfp_delta <- aggregate(hfp_delta, list(hfp_delta$layer), function(x) x[1])%>% 
  st_transform(crs="ESRI:54009")
saveRDS(hfp_delta, "01_raw_data/human_footprint/sf/hfp_delta.rds")
}


# _plots for sense-checking of bins ----

  delta_pal <- rev(RColorBrewer::brewer.pal(5, "RdBu"))
  
  ggplot()+
    geom_sf(data=hfp_delta, aes(fill=as.factor(layer)), colour="transparent")+
    geom_sf(data=tdwg, fill="transparent", colour="gray20")+
    scale_fill_manual(values=delta_pal, 
                      labels=c("Strong decrease",
                               "Weak decrease",
                               "No change",
                               "Weak increase",
                               "Strong increase"))
  ggsave("05_figures/hfp_delta.svg", width=12, height=8)
  ggsave("05_figures/hfp_delta.png", width=23, height=14.3)
  
#   _calculating composition by finding intersections ----

  # set up blank matrix to fill  
  tdwg_composition <- matrix(nrow=nrow(tdwg), ncol=nrow(hfp_delta), data = 0)
  
  # progress bar setup
  cli_progress_bar(
    total = nrow(tdwg),
    format = "{'HFP delta'} | {pb_bar} {pb_percent}"
  )
  
  # for each tdwg region, fill in the matrix with intersection areas 
  for (i in 1:nrow(tdwg_composition)){
    
    cli_progress_update()
    tdwgi <- st_make_valid(tdwg[i,])
    
    intersects <- st_intersects(tdwgi ,
                                hfp_delta,
                                sparse = FALSE)
    intersecting_polygons <- st_intersection(hfp_delta, tdwgi)
    tdwg_composition[i,which(intersects==TRUE)] <- as.numeric(units::set_units(st_area(intersecting_polygons),"km^2"))
  }
  
  tdwg_composition <- tdwg_composition %>% 
    as.data.frame() %>% 
    `colnames<-`(paste0("hfp_delta_class",1:5)) %>% 
    mutate(area_code_l3 = tdwg$LEVEL3_COD)
  
  write_csv(tdwg_composition, "04_output/tdwg_hfp_delta_composition_km.csv")
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
}
  # HFP 2009 (even sized bins) ----
if(calc.even==TRUE){  
  # _bin the raster ----
breaks <- c(-1,10,20,30,40,51)
hfp2_bins_even <- cut(hfp2, breaks=breaks)
  
  #convert to sf and merge all polygons with same value
if(file.exists("01_raw_data/human_footprint/sf/hfp_2009_even.rds")){
  hfp_2009_even <- readRDS("01_raw_data/human_footprint/sf/hfp_2009_even.rds")
} else { 
hfp_2009_even <- st_as_sf(st_as_stars(hfp2_bins_even), as_points = FALSE, merge=TRUE) 
hfp_2009_even <- aggregate(hfp_2009_even, list(hfp_2009_even$layer), function(x) x[1]) %>% 
  st_transform(crs="ESRI:54009")
saveRDS(hfp_2009_even, "01_raw_data/human_footprint/sf/hfp_2009_even.rds")
}

  
  # _plots for sense-checking of bins ----
  
hfp_pal <- rev(RColorBrewer::brewer.pal(5, "RdYlGn"))
  ggplot()+
    geom_sf(data=hfp_2009_even, aes(fill=as.factor(layer)), colour="transparent")+
    geom_sf(data=tdwg, fill="transparent", colour="gray20")+
    scale_fill_manual(values=hfp_pal)
  ggsave("05_figures/hfp_2009_even.svg", width=12, height=8)
  ggsave("05_figures/hfp_2009_even.png", width=23, height=14.3)
  

  #   _calculating composition by finding intersections ----
  
  # set up blank matrix to fill  
  tdwg_composition <- matrix(nrow=nrow(tdwg), ncol=nrow(hfp_2009_even), data = 0)
  
  # progress bar setup
  cli_progress_bar(
    total = nrow(tdwg),
    format = "{'HFP even'} | {pb_bar} {pb_percent}"
  )
  
  # for each tdwg region, fill in the matrix with intersection areas 
  for (i in 1:nrow(tdwg_composition)){
    
    cli_progress_update()
    tdwgi <- st_make_valid(tdwg[i,])
    
    intersects <- st_intersects(tdwgi ,
                                hfp_2009_even,
                                sparse = FALSE)
    intersecting_polygons <- st_intersection(hfp_2009_even, tdwgi)
    tdwg_composition[i,which(intersects==TRUE)] <- as.numeric(units::set_units(st_area(intersecting_polygons),"km^2"))
  }
  
  tdwg_composition <- tdwg_composition %>% 
    as.data.frame() %>% 
    `colnames<-`(paste0("hfp_2009_class",1:10)) %>% 
    mutate(area_code_l3 = tdwg$LEVEL3_COD)
  
  write_csv(tdwg_composition, "04_output/tdwg_hfp_2009_even_composition_km.csv")
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
}
# HFP 2009 (quantile-based bins) ----
if(calc.quant==TRUE){
  # Note: Because >10% of the raster = 0, 10 equal quantiles are not possible. 
  # Instead, we have binned the values at the 20th, 30th, 40th, 50th, 60th,
  # 70th, 80th, 90th and 95th percentiles. 
  
  # _bin the raster ----
  # 10 quantiles at 10% intervals except for bottom 20% and top 5% 
  breaks <- quantile(values(hfp2), c(0,0.2,0.3,0.4,0.5, 0.6,0.7,0.8,0.9,0.95, 1), na.rm=TRUE)
  breaks[1] <- -1
  hfp2_bins_quant10 <- cut(hfp2, breaks=unique(breaks))
  
  #convert to sf and merge all polygons with same value
  if(file.exists("01_raw_data/human_footprint/sf/hfp_2009_quant10.rds")){
    hfp_2009_quant10 <- readRDS("01_raw_data/human_footprint/sf/hfp_2009_quant10.rds")
  } else { 
  hfp_2009_quant10 <- st_as_sf(st_as_stars(hfp2_bins_quant10), as_points = FALSE, merge=TRUE)
  hfp_2009_quant10 <- aggregate(hfp_2009_quant10, list(hfp_2009_quant10$layer), function(x) x[1])%>% 
    st_transform(crs="ESRI:54009")
  saveRDS(hfp_2009_quant10, "01_raw_data/human_footprint/sf/hfp_2009_quant10.rds")
  }
  
  
  # _plots for sense-checking of bins ----
  
  hfp_pal <- rev(RColorBrewer::brewer.pal(10, "RdYlGn"))
    
    ggplot()+
      geom_sf(data=hfp_2009_quant10, aes(fill=as.factor(layer)), colour="transparent")+
      geom_sf(data=tdwg, fill="transparent", colour="gray20")+
      scale_fill_manual(values=hfp_pal)
    ggsave("05_figures/hfp_2009_quant10.svg", width=12, height=8)
    ggsave("05_figures/hfp_2009_quant10.png", width=23, height=14.3)
  
  
  #   _calculating composition by finding intersections ----
  
  # set up blank matrix to fill  
  tdwg_composition <- matrix(nrow=nrow(tdwg), ncol=nrow(hfp_2009_quant10), data = 0)
  
  # progress bar setup
  cli_progress_bar(
    total = nrow(tdwg),
    format = "{'HFP quantile'} | {pb_bar} {pb_percent}"
  )
  
  # for each tdwg region, fill in the matrix with intersection areas 
  for (i in 1:nrow(tdwg_composition)){
    
    cli_progress_update()
    tdwgi <- st_make_valid(tdwg[i,])
    
    intersects <- st_intersects(tdwgi ,
                                hfp_2009_quant10,
                                sparse = FALSE)
    intersecting_polygons <- st_intersection(hfp_2009_quant10, tdwgi)
    tdwg_composition[i,which(intersects==TRUE)] <- as.numeric(units::set_units(st_area(intersecting_polygons),"km^2"))
  }
  
  tdwg_composition <- tdwg_composition %>% 
    as.data.frame() %>% 
    `colnames<-`(paste0("hfp_2009_class",1:10)) %>% 
    mutate(area_code_l3 = tdwg$LEVEL3_COD)
  
  write_csv(tdwg_composition, "04_output/tdwg_hfp_2009_quant10_composition_km.csv")
}
  
  



