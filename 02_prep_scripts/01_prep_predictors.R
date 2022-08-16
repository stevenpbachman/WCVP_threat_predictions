
# load in and wrangle the predictor data

library(dplyr)
library(readr)
library(tidyr)
library(glue)
library(tidyverse)
library(scales)

# load in the new phyt special version of WCVP
wcvp_names <- read_delim("01_raw_data/wcvp_names.txt",
                               "|", 
                               escape_double = FALSE, 
                               trim_ws = TRUE)

# filter on the accepted names only
wcvp_acc = wcvp_names %>%
  filter(taxon_rank == "Species",
         taxon_status == "Accepted")

# get total accepted
tot_acc <- nrow(wcvp_acc)

# thin out the table
wcvp_thin <- wcvp_acc %>%
  select(plant_name_id,taxon_name,family,genus, lifeform_description, climate_description)

# link orders and higher groups
APGIV_families_orders_IUCN_crosswalk <- read.csv("01_raw_data/APGIV_families_orders_IUCN_crosswalk.csv")
apg_fams <- APGIV_families_orders_IUCN_crosswalk %>%
  select(!IUCN.Family)

wcvp_thin <-  left_join(wcvp_thin, apg_fams, by = c("family" = "APG.IV.Family"))

# get total number of accepted angiosperms
angio = wcvp_thin %>%
  filter(complete.cases(HigherGroups))

glue("There are", "{nrow(angio)}", "accepted angiosperm species", .sep=" ")

# set factors
angio$climate_description <- as.factor(angio$climate_description)
angio$APG.IV.Order  <- as.factor(angio$APG.IV.Order )
angio$family <- as.factor(angio$family)
angio$genus <- as.factor(angio$genus)
angio$lifeform_description <- as.factor(angio$lifeform_description)

str(angio)

#### Pred 1 - 2 family and genus ####
# done

# check for missing values
missing_family <- angio %>%
  filter(!complete.cases(family))

missing_genus <- angio %>%
  filter(!complete.cases(genus))

missing_lifeform <- angio %>%
  filter(!complete.cases(lifeform_description))

missing_climate <- angio %>%
  filter(!complete.cases(climate_description))

glue("There are {nrow(missing_family)} missing families", 
     "There are {nrow(missing_genus)} missing genera", 
     "There are {nrow(missing_lifeform)} missing lifeforms", 
     "There are {nrow(missing_climate)} missing climate descriptions", 
     .sep="\n ")

# pred 1 = family
# pred 2 = genus

#### Pred 3 - dominant life form (Humphreys) ####
life_forms <- angio %>%
  select(plant_name_id, lifeform_description)

# check the unique life form values
life_forms_unique <- data.frame(lifeform = unique(life_forms$lifeform_description))

# import the mapping used by Humphreys
life_form_mapping <-  read_csv("01_raw_data/life_form_mapping.csv")

# join mapping to wcvp_thin
angio <- left_join(angio, life_form_mapping, by = "lifeform_description")

angio$humphreys_lifeform <- as.factor(angio$humphreys_lifeform)

# pred 3 = humphreys_lifeform

#### Pred 4 - count of lifeforms ####
#humphreys_unique <- data.frame(lifeform = unique(life_form_mapping$humphreys_lifeform))

# need to think about this - what about 'sometimes' and 'somewhat'?

# try first separation with comma?
#life_forms_tidy = life_forms_unique %>%
#  separate(lifeform, into = c("first", "second"), sep = ", ")


#### Pred 5- count of TDWG regions per species ####
# load in the new phyt special version of WCVP distributions
wcvp_geog <- read.table("01_raw_data/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8") 

# thin out the files and filter out introduced, extinct and location doubtfull
geog_thin <- wcvp_geog %>%
  filter(introduced == "0" & extinct == "0", location_doubtful == "0"  ) %>%
  select(plant_name_id, area_code_l3)

# get the count of TDWG regions per species
count_geog <- geog_thin %>%
  group_by(plant_name_id) %>%
  summarise(L3_count = n_distinct(area_code_l3))

# join mapping to wcvp_thin
angio <- left_join(angio, count_geog, by = "plant_name_id")

#### remove NAs ####
#for final list to try model on 
angio_pred_v1 <-  angio %>%
  filter(complete.cases(family),
         complete.cases(climate_description),
         complete.cases(HigherGroups),         
         complete.cases(humphreys_lifeform),
         complete.cases(L3_count))

glue("We have {nrow(angio_pred_v1)} angiosperm species to model", "({percent_format(accuracy=0.1)(nrow(angio_pred_v1)/nrow(angio))})",
     .sep=" ")


#### link RL assessments ####
rl <- read.csv("01_raw_data/redlistJul2022_wcvpNewPhyt.csv")

# now join threat status to wcvp 
angio_pred_v1 <- left_join(angio_pred_v1, rl, by = c("plant_name_id" = "accepted_plant_name_id"))

str(angio_pred_v1)
# save it   
#write_csv(angio_pred_v1, paste0("03_analysis_scripts/angio_pred_v1.csv"))


  



  
  
  