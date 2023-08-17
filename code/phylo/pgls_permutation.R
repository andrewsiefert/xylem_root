library(dplyr)
library(phytools)
library(nlme)
library(doParallel)
library(ggplot2)


# read in traits and tree
species <- readRDS("data/cleaned/trait_data.rds") %>% 
  select(species = species_wfo)

p50 <- readRDS("data/cleaned/p50_database.rds") %>% 
  select(species = species_wfo, p50 = p50) %>%
  semi_join(species)
rd <- readRDS("data/cleaned/rooting_data.rds") %>% 
  select(species = species_wfo, rd = dr) %>%
  semi_join(species)


tree <- read.tree("data/phylo/p50_rd_phylo.tre")




nrep <- 1e4


fits <- list()

for(i in 1:nrep) {
  
  trait <- species %>%
    left_join(p50 %>% group_by(species) %>% sample_n(1) %>% ungroup()) %>%
    left_join(rd %>% group_by(species) %>% sample_n(1) %>% ungroup()) %>%
    mutate(p50 = scale(sqrt(abs(p50))),
           rd = scale(log(rd))) %>%
    as.data.frame()
  
  tree2 <- tree
  
  # to ensure correct order of tips and traits
  rownames(trait) <- trait$species %>% str_replace_all(" ", "_")       ###add row names
  p.dist.mat <- cophenetic(tree2)         ###to ensure correct order
  trait <- trait[row.names(p.dist.mat),] ###to ensure correct order
  tree2$tip.label <- str_replace_all(tree2$tip.label, "_", " ")
  rownames(trait) <- tree2$tip.label
  
  
  # test phylogenetic correlation
  spp <- trait$species
  corLambda <- corPagel(value=1, phy=tree2, form=~spp)
  
  fits[[i]] <- gls(rd ~ p50, data = trait, correlation = corLambda)  
  
}
 
saveRDS(fits, "results/phylo/pgls_permutations.rds")

rm(fits)
gc()



# Distributions -----------------------------------------------------------


fits <- readRDS("results/phylo/pgls_permutations.rds")

slopes <- sapply(fits, function(i) -coefficients(i)[2])
pvals <- sapply(fits, function(i) summary(i)[[19]][2,4])

qplot(slopes) + theme_bw() + labs(x = "Slope") + scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
ggsave("results/figures/pgls_permutation_slopes.jpg")

qplot(pvals) + theme_bw() + labs(x = "P-value") + scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
ggsave("results/figures/pgls_permutation_pvals.jpg")
