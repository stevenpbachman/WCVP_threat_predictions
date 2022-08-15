
# load in and wrangle the predictor data

library(dplyr)
library(readr)
library(tidyr)

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

wcvp_thin %>%
  group_by(HigherGroups)%>%
  count()

unique(wcvp_thin$HigherGroups)

#### Pred 1 - 2 family and genus ####


# check for missing values
missing_family <- wcvp_thin %>%
  filter(!complete.cases(family))
nrow(missing_family)

missing_genus <- wcvp_thin %>%
  filter(!complete.cases(genus))
nrow(missing_genus)

missing_lifeform <- wcvp_thin %>%
  filter(!complete.cases(lifeform_description))
nrow(missing_lifeform)

missing_climate <- wcvp_thin %>%
  filter(!complete.cases(climate_description))
nrow(missing_climate)

# pred 1 = family
# pred 2 = genus

#### Pred 3 - dominant life form (Humphreys) ####
life_forms <- wcvp_thin %>%
  select(plant_name_id, lifeform_description)

# check the unique life form values
life_forms_unique <- data.frame(lifeform = unique(life_forms$lifeform_description))

# import the mapping used by Humphreys
life_form_mapping <-  read_csv("01_raw_data/life_form_mapping.csv")

# join mapping to wcvp_thin
wcvp_thin <- left_join(wcvp_thin, life_form_mapping, by = "lifeform_description")
  
# pred 3 = humphreys_lifeform

#### Pred 4 - count of lifeforms ####
humphreys_unique <- data.frame(lifeform = unique(life_form_mapping$humphreys_lifeform))

# first row selects the columns you want
# second row does the separation by comma
life_forms_tidy = life_forms_unique %>%
  separate(lifeform, into = c("first", "second"), sep = " or ")

# try first separation with comma?




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
wcvp_thin <- left_join(wcvp_thin, count_geog, by = "plant_name_id")

#### remove NAs ####
#for final list to try model on 
wcvp_pred_v1 <-  wcvp_thin %>%
  filter(complete.cases(family),
         complete.cases(climate_description),
         complete.cases(HigherGroups),         
         complete.cases(humphreys_lifeform),
         complete.cases(count))

#### link RL assessments ####
redlistJul2022_wcvpNewPhyt <- read.csv("01_raw_data/redlistJul2022_wcvpNewPhyt.csv")

rl <- redlistJul2022_wcvpNewPhyt

#cats <- rl %>%
#  group_by(category) %>%
 # count()

rl <- rl %>%
  mutate(threat_stat = case_when(
    category == "CR" ~ "Threatened",
    category == "EN" ~ "Threatened",
    category == "VU" ~ "Threatened",
    category == "NT" ~ "Not_Threatened",
    category == "LC" ~ "Not_Threatened",
    category == "LR/cd" ~ "Threatened",
    category == "LR/nt" ~ "Threatened",
    category == "LR/lc" ~ "Not_Threatened",
    category == "EW" ~ "Threatened",
    category == "EX" ~ "Threatened",
    category == "DD" ~ "Data_deficient"
    )) %>%
  select(scientific_name, category, threat_stat, accepted_plant_name_id)


  
# now join threat status to wcvp 
wcvp_pred_v1 <- left_join(wcvp_pred_v1, rl, by = c("plant_name_id" = "accepted_plant_name_id"))


  
  
  