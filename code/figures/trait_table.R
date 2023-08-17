library(tidyverse)


# read in species mean traits
traits <- readRDS("data/cleaned/trait_data.rds") %>%
  select(species = species_wfo, p50, rd_max)


# add observation counts
p50_n <- readRDS("data/cleaned/p50_database.rds") %>%
  count(species_wfo) %>%
  select(species = 1, n_p50 = 2)

rd_n <- readRDS("data/cleaned/rooting_data.rds") %>%
  count(species_wfo) %>%
  select(species = 1, n_rd = 2)


# add family
fams <- read_csv("data/phylo/phylo_species_list.csv") %>%
  select(species, family)


# add leaf phenology
pheno <- read_csv("data/traits/leaf_phenology.csv")


out <- left_join(traits, fams) %>% 
  left_join(pheno) %>% 
  left_join(p50_n) %>% 
  left_join(rd_n) %>%
  mutate_at(vars(n_p50, n_rd), ~paste0("(", ., ")")) %>%
  select(species, family, leaf_phenology, p50, n_p50, rd_max, n_rd)

write_csv(out, "results/tables/trait_table.csv")
