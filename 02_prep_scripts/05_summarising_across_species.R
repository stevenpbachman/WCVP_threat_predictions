###################################################################### #
# Summing biome and HFP over species                                   #
###################################################################### #


# load packages ----
library(tidyverse)
library(rWCVP)
library(sf)

# get WCVP list of accepted species and distributions ----

#load data
dist <- rWCVPdata::wcvp_distributions
names <- rWCVPdata::wcvp_names

tdwg_names <- rWCVPdata::wgsrpd3 %>% 
  st_drop_geometry() 


tdwg <- rWCVPdata::wgsrpd3 %>% 
  st_transform(crs="ESRI:54009") %>% 
  mutate(area=as.numeric(units::set_units(st_area(.),"km^2"))) %>% 
  st_drop_geometry() %>% 
  select(area_code_l3=LEVEL3_COD,
         area)

tdwg_1 <- rWCVPdata::wgsrpd3 %>% 
  st_transform(crs="ESRI:54009") %>% 
  mutate(area=as.numeric(st_area(.))/10e6) %>% 
  st_drop_geometry() %>% 
  select(area_code_l3=LEVEL3_COD,
         area)

biomes <- read_csv("04_output/tdwg_biome_composition_km.csv") 
# biomes$area_code_l3 <- tdwg_names$LEVEL3_COD[order(tdwg_names$LEVEL3_COD)]
biomes<- biomes %>% left_join(tdwg)

# Manually adding things without any ecoregion
#Mozambique Channel Islands and Nauru are Tropical moist
biomes[which(biomes$area_code_l3 %in% c("MCI","NRU")),1] <- 1
#Marcus Island/Minamitori Shima is Tropical dry
biomes[which(biomes$area_code_l3 %in% c("MCS")),2] <- 1
#Selvagens Islands are Mediterranean 
biomes[which(biomes$area_code_l3 %in% c("SEL")),12] <- 1

biomes[,1:16] <- biomes[,1:16]/rowSums(biomes[,1:16])*biomes$area

# checks to make sure it worked as expected
# check <- biomes
# check[,1:16] <- check[,1:16]/check$area
# check$totalprop <- rowSums(check[,1:16])

hfp_2009 <- read_csv("04_output/tdwg_hfp_2009_quant10_composition_km.csv") %>% left_join(tdwg)
hfp_2009[,1:10] <- hfp_2009[,1:10]/rowSums(hfp_2009[,1:10])*hfp_2009$area

#Manually adding things with no footprint data
manual <- read_csv("01_raw_data/manually_estimated_hfp2009.csv") %>% left_join(tdwg)
manual[,4:13] <- manual[,4:13]/rowSums(manual[,4:13])*manual$area
manual <- manual %>%
  select(colnames(hfp_2009))
hfp_2009 <- hfp_2009 %>% 
  filter(!area_code_l3 %in% manual$area_code_l3) %>% 
  rbind(manual)

# colnames(hfp_2009) <- paste0("hfp_2009_class",1:10)
# hfp_2009$area_code_l3 <- tdwg_names$LEVEL3_COD[order(tdwg_names$LEVEL3_COD)]

# checks to make sure it worked as expected
check <- hfp_2009
check[,1:10] <- check[,1:10]/check$area
check$totalprop <- rowSums(check[,1:10])
highestclass10 <- check$area_code_l3[order(check$hfp_2009_class10)]
check %>% 
  pivot_longer(1:10, names_to = "class") %>%
  mutate(area_code_l3=factor(area_code_l3, levels=highestclass10),
         class=factor(class, levels=colnames(hfp_2009)[1:10])) %>% 
  ggplot(aes(x=area_code_l3, y=value, fill=class)) +
  geom_bar(position = position_stack(), stat = "identity", width = 1)+
  scale_fill_viridis_d()+
  ggpubr::theme_pubclean()+
  coord_cartesian(expand=FALSE)

hfp_delta <- read_csv("04_output/tdwg_hfp_delta_composition_km.csv")%>% left_join(tdwg)
# colnames(hfp_delta) <- paste0("hfp_delta_class",1:5)
# hfp_delta$area_code_l3 <- tdwg_names$LEVEL3_COD[order(tdwg_names$LEVEL3_COD)]
hfp_delta$totalprop <- rowSums(hfp_delta[,1:5]/hfp_delta$area)
#manually add no change to places that are not on the hfp map
hfp_delta$hfp_delta_class3[which(hfp_delta$totalprop==0)] <- 1

hfp_delta <- hfp_delta %>% select(-totalprop)
hfp_delta[,1:5] <- hfp_delta[,1:5]/rowSums(hfp_delta[,1:5])*hfp_delta$area

# filter to accepted species
df <- names %>% 
  filter(taxon_status=="Accepted",
         taxon_rank == "Species",
         genus_hybrid == "",
         species_hybrid == "",
         infraspecies == "") #should be redundant but better to be safe

code_corrections = c("SAO"="SOA",
                     "SRI"="SRL",
                     "QLS"="QLD")

# join ids to distributions
accepted_distributions <-
  dist %>%
  filter(extinct + introduced + location_doubtful == 0) %>%
  select(area_code_l3, plant_name_id) %>%
  inner_join(
    df,
    by="plant_name_id"
  )

# clean up code names
accepted_distributions <-
  accepted_distributions %>%
  mutate(area_code_l3=str_to_upper(area_code_l3),
         area_code_l3=recode(area_code_l3, !!! code_corrections)) %>%
  filter(area_code_l3 %in% tdwg_names$LEVEL3_COD)

# remove infraspecies ids and deduplicate by species and region
accepted_distributions <-
  accepted_distributions %>%
  select(area_code_l3, taxon_name, plant_name_id) %>%
  group_by(plant_name_id, area_code_l3) %>%
  filter(row_number() == 1) %>%
  ungroup()

no_dist <- df[which(!df$plant_name_id %in% accepted_distributions$plant_name_id),]

wcvp_preds <- accepted_distributions %>% 
  left_join(tdwg_1) %>% 
  left_join(biomes %>% select(-area)) %>% 
  left_join(hfp_2009 %>% select(-area)) %>% 
  left_join(hfp_delta %>% select(-area)) %>% 
  group_by(plant_name_id, taxon_name) %>% 
  summarise(across(where(is.numeric), list(sum))) %>% 
  `colnames<-`(gsub("_1","", colnames(.)))

  wcvp_preds_percent <- wcvp_preds 
  wcvp_preds_percent[,4:ncol(wcvp_preds_percent)] <- 
    wcvp_preds_percent[,4:ncol(wcvp_preds_percent)]/wcvp_preds_percent$area*10 #because the decimals got messed up somewhere?

write_csv(wcvp_preds_percent, "tdwg-derived-predictors_percent-of-range.csv")
write_csv(wcvp_preds, "tdwg-derived-predictors_rawareas.csv")

