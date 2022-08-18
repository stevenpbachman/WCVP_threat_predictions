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
tdwg <- rWCVPdata::wgsprd3 


#extra teeny buffer
eps = sqrt(.Machine$double.eps)*2

#crop to clean up
tdwg <- tdwg %>% 
  st_crop(xmin=-180+eps, xmax=180-eps, ymin=-90+eps, ymax=83.62361)


# WWF ECOREGIONS/BIOMES  ####

if(!file.exists("04_output/tdwg_biome_composition_km.csv")){
  
  # load data
  wwf <- st_read("01_raw_data/ecoregions/wwf_terr_ecos.shp")
  
  key <- read_csv("01_raw_data/ecoregions/ecoregions_key.csv") 
  key$biome_name <-  factor(key$biome_name, levels=key$biome_name)
  
  # aggregate to biome
  wwf_agg_bio <- aggregate(wwf,by=list(wwf$BIOME), FUN = "first") %>%
    left_join(key)
  
  
  # plot the two datasets ----
  ggplot()+
    geom_sf(data=wwf_agg_bio, colour="transparent", aes(fill=biome_name))+
    geom_sf(data=tdwg, colour="gray20", fill="transparent")+
    scale_fill_manual(values=wwf_agg_bio$fill)+
    coord_sf(expand=FALSE, crs="ESRI:54009")+
    theme_void()+
    guides(fill=guide_legend("Biome"))
  
  #ggsave("03_figures/tdwg_biomes.svg")
  #ggsave("03_figures/tdwg_biomes.png")
  
  # calculate areas of each biome per tdwg region ----
  
  # calculate binary intersections
  tdwg_composition <- st_intersects(tdwg ,
                                    wwf_agg_bio,
                                    sparse = FALSE)
  
  #transform to mollweide projection for area
  tdwg <- tdwg %>% 
    st_transform(crs="ESRI:54009")
  
  wwf_agg_bio <- wwf_agg_bio %>% 
    st_transform(crs="ESRI:54009")
  
  
  # fill in the matrix with intersection areas for each tdwg region
  
  cli_progress_bar(
    total = nrow(tdwg),
    format = "Calculating intersections | {pb_bar} {pb_percent}"
  )
  
  for (i in 1:nrow(tdwg_composition)){
    cli_progress_update()
    tdwgi <- st_make_valid(tdwg[i,])
    intersecting_biomes <- st_intersection(wwf_agg_bio, tdwgi)
    tdwg_composition[i,which(tdwg_composition[i,]==TRUE)] <- st_area(intersecting_biomes)
  }
  
  # scale areas to sq km (from sq m) and round to nearest
  tdwg_composition_km <- as.data.frame(round(tdwg_composition/10e6)) %>% 
    `colnames<-`(key$biome_name) %>% 
    `rownames<-`(tdwg$LEVEL3_COD)
  
  write_csv(tdwg_composition_km, "04_output/tdwg_biome_composition_km.csv")
}