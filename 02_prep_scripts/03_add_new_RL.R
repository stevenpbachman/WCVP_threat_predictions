#03_add_new_RL


# reviewer asked to re-run model on updated RL 2022-2
# download new RL - saved in 01_raw_data
# run a new RL exact match on the special names. use only phontic function for the fuzzy match

library(tidyverse)
library(rWCVP)

apg_families <- 
  read_csv("01_raw_data/APGIV_families_orders_IUCN_crosswalk.csv") |>
  select("higher_groups"="HigherGroups", "order"="APG IV Order",
         "family"="APG IV Family")

# assessment csv to get year and category
assessrl_2022_1 <- read_csv("01_raw_data/assessments_v2022_1.csv") %>%
  select(scientificName, internalTaxonId, redlistCategory, redlistCriteria, yearPublished, assessmentDate)

assessrl_2022_2 <- read_csv("01_raw_data/assessments_v2022_2.csv") %>%
  select(internalTaxonId, redlistCategory, redlistCriteria, yearPublished, assessmentDate)

# load the special version of WCVP
names <- read_delim("G:/Users/sb42kg/WCVP_threat_predictions/01_raw_data/wcvp_names.txt",
                    delim="|", escape_double=FALSE, 
                    trim_ws = TRUE)

#names$plant_name_id <- as.factor(names$plant_name_id)

# original RL data - name matched to WCVP
rl_2022_1 <- read_csv("01_raw_data/redlistJul2022_wcvpNewPhyt.csv") |>
  add_count(accepted_plant_name_id) |>
  filter(n == 1 | scientific_name == match_name & match_status == "Accepted") |>
  select(-n)

# new RL 2022-2 not matched yet
rl_2022_2 <- read.csv(here("G:/Users/sb42kg/WCVP_threat_predictions/01_raw_data/taxonomy_v2022_2.csv"))
rl_2022_2$familyName <- tolower(rl_2022_2$family)
rl_2022_2$familyName <- R.utils::capitalize(rl_2022_2$family)

# filter out the non-vascular from the new RL
rl_2022_2_vasc <- dplyr::left_join(rl_2022_2, apg_families, by = c("familyName" = "family")) |>
  dplyr::filter(complete.cases(order))

# get the species that are new in 2022-2
rl_diffs <- left_join(rl_2022_2_vasc, rl_2022_1, by = c("scientificName" = "scientific_name")) |>
  filter(!complete.cases(match_name)) |>
  select( ,1:15) |>
  rename(authority = authority.x)

# 3240 names in 2022-2 that are not in 2022-1

#### 1. exact matches first ####
# now exact match on the diffs - new RL taxa
rl_diffs_em <- wcvp_match_exact(rl_diffs, wcvp_names = names, name_col = "scientificName", 
                                 author_col = "authority", id_col = "internalTaxonId")

# write down and manually sort the multi matches
write_csv(rl_diffs_em, "rl_diffs_em.csv")
# then import back in
rl_diffs_em <- read_csv("rl_diffs_em.csv") %>%
  filter(complete.cases(match_type))

# for the exact matches - link to accepted name
accepted_matches <- rl_diffs_em %>%
  left_join(names, by=c("wcvp_accepted_id"="plant_name_id")) %>%
  mutate(keep=case_when(
    taxon_status == "Accepted" & (wcvp_status != "Synonym" | wcvp_homotypic) ~
      "Matched to an accepted name",
    TRUE ~ "Not matched to an accepted name"
  ))

count(accepted_matches, keep)

# just need to join back to the assessment csv to get year and category
accepted_matches <- left_join(accepted_matches, assessrl_2022_2, by = "internalTaxonId" )

#final exact matches
final_em <- 
  accepted_matches %>%
  filter(keep == "Matched to an accepted name") %>%
  select(internalTaxonId,scientificName, authority, redlistCategory, redlistCriteria, yearPublished, assessmentDate,
         match_name=wcvp_name, match_status=wcvp_status,
         accepted_plant_name_id=wcvp_accepted_id, ipni_id,
         accepted_taxon_name=taxon_name, accepted_taxon_authors=taxon_authors)

#### 2. now the fuzzy matches ####
# import the diffs in and keep the non-exact matches
rl_diffs_fm <- read_csv("rl_diffs_em.csv") %>%
  filter(!complete.cases(match_type)) %>%
  select(internalTaxonId,scientificName, authority,familyName,genusName,speciesName,authority)

# run the phonetic match only (fuzzy currently not working due to 'special' format issue)
rl_diffs_pm <- phonetic_match(rl_diffs_fm, names, name_col = "scientificName")

rl_diffs_pm <- rl_diffs_pm %>%
  filter(str_detect(match_type, "Fuzzy")) %>%
  mutate(
    keep = case_when( #set up a keep column
      match_similarity < 0.9 ~ NA_real_, # fill with blank for dissimilar names
      #wcvp_author_edit_distance == 0 ~ 1, # fill with 1 if authors identical
      match_edit_distance <= 1 ~ 1 # fill with 1 if only one letter different
    )
  )

#how many did this resolve?
table(rl_diffs_pm$keep, useNA = "always")

# write down and double check the multi matches and any rogue phonetic matches
write_csv(rl_diffs_pm, "rl_diffs_pm.csv")
# then import back in
rl_diffs_pm <- read_csv("rl_diffs_pm.csv") 

# for the phonetic matches - link to accepted name
accepted_p_matches <- rl_diffs_pm %>%
  left_join(names, by=c("wcvp_accepted_id"="plant_name_id")) %>%
  mutate(keep=case_when(
    taxon_status == "Accepted" & (wcvp_status != "Synonym" | wcvp_homotypic) ~
      "Matched to an accepted name",
    TRUE ~ "Not matched to an accepted name"
  ))

count(accepted_p_matches, keep)

# just need to join back to the assessment csv to get year and category
accepted_p_matches <- left_join(accepted_p_matches, assessrl_2022_2, by = "internalTaxonId" )

#final exact matches
final_pm <- 
  accepted_p_matches %>%
  filter(keep == "Matched to an accepted name") %>%
  select(internalTaxonId,scientificName, authority, redlistCategory, redlistCriteria, yearPublished, assessmentDate,
         match_name=wcvp_name, match_status=wcvp_status,
         accepted_plant_name_id=wcvp_accepted_id, ipni_id,
         accepted_taxon_name=taxon_name, accepted_taxon_authors=taxon_authors)

#### 3. merge the matches to get the new species between RL 2022-1 and 2022-2
redlist_2022_2_extras <- bind_rows(final_pm, final_em) %>%
  rename(accepted_taxon_status = match_status) 

write_csv(redlist_2022_2_extras, "01_raw_data/redlist_2022-2_extras_only.csv")

# now merge the extras with the original to get the full 2022-2 list for the modelling
# add extra columns
redlist_2022_1_all <- left_join(rl_2022_1, assessrl_2022_1, by = c("scientific_name" = "scientificName" )) %>%
  select(-category, -match_status) %>%
  rename(scientificName = scientific_name)

# final merge - the original 2022-1 and the cleaned/matched extras
redlist_2022_2_all <- bind_rows(redlist_2022_2_extras, redlist_2022_1_all)
# note that duplicates might have been introduced so check and remove
redlist_2022_2_all <- redlist_2022_2_all %>%
  distinct(accepted_plant_name_id, .keep_all = TRUE)
# write to file
write_csv(redlist_2022_2_all, "01_raw_data/redlist_2022_2_all.csv")

# now the version with just the 'recent' assessments <= 10 yeas
redlist_2022_2_last_ten_years <- redlist_2022_2_all %>%
  filter(yearPublished >= 2012)
# note that duplicates might have been introduced so check and remove
redlist_2022_2_last_ten_years <- redlist_2022_2_last_ten_years %>%
  distinct(accepted_plant_name_id, .keep_all = TRUE)
# and write to file
write_csv(redlist_2022_2_last_ten_years, "01_raw_data/redlist_2022_2_last_ten_years.csv")



###################
# import back in
redlist_2022_2_all <- read_csv("01_raw_data/redlist_2022_2_all.csv") 

