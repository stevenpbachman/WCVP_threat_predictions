# load in and wrangle the predictor data

# libraries ----
library(tidyverse) # data manipulation libaries
library(glue)      # string interpolation

# load data ----
# load in the new phyt special version of WCVP
wcvp_names <- read_delim("01_raw_data/wcvp_names.txt",
                         delim="|", escape_double=FALSE, 
                         trim_ws = TRUE)

wcvp_dist <- 
  read_delim("01_raw_data/wcvp_distribution.txt",
             delim="|", escape_double=FALSE, trim_ws=TRUE) |>
  filter(introduced + extinct + location_doubtful == 0) |>
  filter(! is.na(area_code_l3))

phylo_vectors <- 
  read_csv("01_raw_data/eigenvectors_selected_BS_252_0.2.csv") |>
  rename("genus"="...1") |>
  rename_with(~str_replace(.x, "c", "pvr"), starts_with("c"))

biome_names <- c(
  "biome_trop.mf"="Tropical and subtropical moist broadleaf forests",
  "biome_trop.df"="Tropical and subtropical dry broadleaf forests",
  "biome_trop.cf"="Tropical and subtropical coniferous forests",
  "biome_temp.bf"="Temperate broadleaf and mixed forests",
  "biome_temp.cf"="Temperate coniferous forests",
  "biome_boreal"="Boreal forests/taiga",
  "biome_trop.grass"="Tropical and subtropical grasslands, savannas, and shrublands",
  "biome_temp.grass"="Temperate grasslands, savannas, and shrublands",
  "biome_flood.grass"="Flooded grasslands and savannas",
  "biome_mont.grass"="Montane grasslands and shrublands",
  "biome_tundra"="Tundra",
  "biome_med.f"="Mediterranean forests, woodlands, and shrub",
  "biome_deserts"="Deserts and xeric shrublands",
  "biome_mangrove"="Mangrove",
  "biome_lake"="Lake",
  "biome_rock.ice"="Rock and Ice"
)

tdwg_vars <- 
  read_csv("01_raw_data/tdwg-derived-predictors_percent-of-range.csv") |>
  rename(!!! biome_names) |>
  mutate(across(-c(plant_name_id, taxon_name, area), ~.x/100))

apg_families <- 
  read_csv("01_raw_data/APGIV_families_orders_IUCN_crosswalk.csv") |>
  select("higher_groups"="HigherGroups", "order"="APG IV Order",
         "family"="APG IV Family")

#redlist <- read_csv("01_raw_data/redlistJul2022_wcvpNewPhyt.csv") |>
redlist <- read_csv("01_raw_data/redlist_2022_2_all.csv") |>
  add_count(accepted_plant_name_id) |>
  #filter(n == 1 | scientific_name == match_name & match_status == "Accepted") |>
  filter(n == 1 | scientificName == match_name & accepted_taxon_status == "Accepted") |> # get rid of dups
  select(-n)  |>
  rename(category = redlistCategory)

# remove species duplicated on Red List - accepted and synonym
# redlist <- redlist |>
#   dplyr::filter(internalTaxonId != 38347) |>
#   dplyr::filter(scientificName != "Aster chimanimaniensis")
# dplyr::filter(scientificName != "Henckelia smitinandii")

wcvp_YoD <- read_csv("01_raw_data/wcvp_YoD.csv") |>
  distinct(plant_name_id, year)

# generate species list ----
# filter on the accepted names only
species_list <-  
  wcvp_names |>
  filter(taxon_rank == "Species",
         taxon_status == "Accepted",
         !is.na(taxon_name),
         is.na(genus_hybrid),
         is.na(species_hybrid)) |>
  select(plant_name_id, taxon_name, family, genus, lifeform_description, climate_description) |>
  left_join(apg_families, by="family")
            #, relationship="many-to-one")

# get total number of accepted angiosperms
species_list <- filter(species_list, !is.na(higher_groups))
species_list <- filter(species_list, plant_name_id %in% wcvp_dist$plant_name_id)

cat(glue::glue("There are", "{nrow(species_list)}", "accepted angiosperm species", .sep=" "))

cat(str(species_list))

# check missing values ----
cat("Missing values from WCVP:")
species_list |>
  select(family, genus, lifeform_description, climate_description) |>
  summarise(across(everything(), ~sum(is.na(.x))))

# standardise lifeform ----

# import the mapping used by Humphreys
life_form_mapping <-  read_csv("01_raw_data/life_form_mapping.csv")

# join mapping to wcvp_thin
predictors <- left_join(species_list, life_form_mapping, by="lifeform_description")
#, relationship="many-to-one")

# count regions per species ----

region_counts <- 
  wcvp_dist |>
  filter(introduced + extinct + location_doubtful == 0) |>
  group_by(plant_name_id) |>
  summarise(
    L3_count=n_distinct(area_code_l3)
  )

predictors <- inner_join(predictors, region_counts, by="plant_name_id")
#, relationship="one-to-one")

# join phylovectors ----

predictors <-
  predictors |>
  left_join(
    phylo_vectors,
    by="genus")
#, relationship="many-to-one")

# join tdwg-based predictors ----

predictors <-
  predictors |>
  left_join(
    tdwg_vars,
    by=c("plant_name_id", "taxon_name"))
    #,relationship="one-to-one")

# link RL assessments ----
predictors <- 
  predictors |>
  left_join(
    redlist |> select(category, accepted_plant_name_id), 
    by=c("plant_name_id"="accepted_plant_name_id"))
#, relationship="one-to-one")

#Ilex pseudoumbelliformis - 861317-az

# link year of description ----
predictors <- 
  predictors |>
  left_join(
    wcvp_YoD |> select(plant_name_id, year), 
    by="plant_name_id")
#,relationship="one-to-one")

# check missing values again ----
cat("Missing predictor values:")
missing_preds <- predictors |>
  summarise(across(everything(), ~sum(is.na(.x)))) |>
  pivot_longer(everything()) |>
  filter(value > 0) |>
  arrange(desc(value))

write_csv(missing_preds, "05_results/review1/missing_predictor_values.csv")

missing <- 
  predictors |>
  summarise(across(everything(), ~sum(is.na(.x)))) |>
  pivot_longer(everything()) |>
  filter(value > 0) |>
  arrange(desc(value))

cat(missing$name)

cat(glue::glue("{sum(!is.na(predictors$category))} assessed species"))

# save to file ----
now <- format(Sys.time(), "%Y%m%d-%H%M%S")
name <- paste("predictors-angiosperm", now, sep="-")
write_csv(predictors, file.path("output", paste0(name, ".csv")))
