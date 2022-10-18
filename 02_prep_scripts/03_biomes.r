################################################################################ #
# Calculating TDWG predictor values that can be averaged or summed over species  #
################################################################################ #

# load packages ----

library(tidyverse)
library(sf)
library(rWCVP)
library(cli)
library(raster)
library(stars)

# prevent spherical geometry errors
sf_use_s2(FALSE)

# load data ----

# tdwg regions
tdwg <- rWCVPdata::wgsrpd3 


#extra teeny buffer
eps = sqrt(.Machine$double.eps)*2

#crop to clean up
tdwg <- tdwg %>% 
  st_crop(xmin=-180+eps, xmax=180-eps, ymin=-90+eps, ymax=83.62361) %>% 
  st_make_valid()


# WWF ECOREGIONS/BIOMES  ####
  
  # load data
  wwf <- st_read("01_raw_data/ecoregions/wwf_terr_ecos.shp")
  
  key <- read_csv("01_raw_data/ecoregions/ecoregions_key.csv") 
  key$biome_name <-  factor(key$biome_name, levels=key$biome_name)
  
  # aggregate to biome
  wwf_agg_bio <- aggregate(wwf,by=list(wwf$BIOME), FUN = "first") %>%
    left_join(key)
  
  #fn to erase
  st_erase = function(x, y) st_difference(x, st_union(st_combine(y)))
  #erase biome from tdwg
  # tdwgnoeco <- st_erase(tdwg, wwf_agg_bio) # keeps throwing an error

  
  
  
  # plot the two datasets ----
<<<<<<< HEAD
  # ggplot()+
  #   geom_sf(data=wwf_agg_bio, colour="transparent", aes(fill=biome_name))+
  #   geom_sf(data=tdwg, colour="gray20", fill="transparent")+
  #   scale_fill_manual(values=wwf_agg_bio$fill)+
  #   coord_sf(expand=FALSE, crs="ESRI:54009")+
  #   guides(fill=guide_legend("Biome"))
  # 
  # ggsave("05_figures/tdwg_biomes.svg", width=12, height=8)
  # ggsave("05_figures/tdwg_biomes.png", width=23, height=14.3)
=======
  ggplot()+
    geom_sf(data=wwf_agg_bio, colour="transparent", aes(fill=biome_name))+
    geom_sf(data=tdwg, colour="gray20", fill="transparent")+
    scale_fill_manual(values=wwf_agg_bio$fill)+
    coord_sf(expand=FALSE, crs="ESRI:54009")+
    theme_void()+
    guides(fill=guide_legend("Biome"))
  
  ggsave("05_figures/tdwg_biomes.svg")
  ggsave("05_figures/tdwg_biomes.png")
>>>>>>> e2233f77e239498a0f3152c638a49e4ee97ebf87
  
  # calculate areas of each biome per tdwg region ----
  
  # set up matrix
  tdwg_composition <- matrix(nrow=nrow(tdwg), ncol=nrow(wwf_agg_bio), data = 0)
  
  #transform to mollweide projection for area
  tdwg <- tdwg %>% 
    st_transform(crs="ESRI:54009")
  
  wwf_agg_bio <- wwf_agg_bio %>% 
    st_transform(crs="ESRI:54009")
  
  
  # # fill in the matrix with intersection areas for each tdwg region
  # 
  # cli_progress_bar(
  #   total = nrow(tdwg),
  #   format = "Biomes | {pb_bar} {pb_percent}"
  # )
  # 
  # for (i in 1:nrow(tdwg_composition)){
  #   cli_progress_update()
  #   tdwgi <- st_make_valid(tdwg[i,])
  #   
  #   intersects <- st_intersects(tdwgi ,
  #                               wwf_agg_bio,
  #                               sparse = FALSE)
  #   intersecting_polygons <- st_intersection(wwf_agg_bio, tdwgi)
  #   tdwg_composition[i,which(intersects==TRUE)] <- as.numeric(units::set_units(st_area(intersecting_polygons),"km^2"))
  # }
  # 
  # 
  # 
  # 
  # tdwg_composition <- tdwg_composition %>% 
  #   as.data.frame() %>% 
  #   `colnames<-`(wwf_agg_bio$biome_name) %>% 
  #   mutate(area_code_l3 = tdwg$LEVEL3_COD)
  # 
  # write_csv(tdwg_composition, "04_output/tdwg_biome_composition_km.csv")
  
  
  
  #######################
  
  cli_progress_bar(
    total = nrow(tdwg),
    format = "Biomes | {pb_bar} {pb_percent}"
  )
  
  for (i in 49:nrow(tdwg_composition)){
    cli_progress_update()
    tdwgi <- st_make_valid(tdwg[i,])
    
    intersects <- st_intersects(tdwgi ,
                                wwf_agg_bio,
                                sparse = FALSE)
    noeco <- st_erase(tdwgi, wwf_agg_bio)
    if(length(st_area(noeco))>0){
    
    ggplot()+
      geom_sf(data=noeco, colour="transparent", fill="red")+
      geom_sf(data=tdwgi, colour="gray20", fill="transparent")+
      coord_sf(expand=FALSE, crs="ESRI:54009")
    ggsave(filename=paste0("noecoplots/noeco_",tdwgi$LEVEL3_COD,".png"))
    
    }
  
  }
