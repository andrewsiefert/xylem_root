library(tidyverse)
library(phytools)
library(nlme)


# read in traits and tree
trait <- readRDS("data/cleaned/trait_data.rds") %>% 
  select(species= species_wfo, p50, rd = rd_max) %>%
  as.data.frame()
tree <- read.tree("data/phylo/p50_rd_phylo.tre")

# to ensure correct order of tips and traits
rownames(trait) <- trait$species %>% str_replace_all(" ", "_")       ###add row names
p.dist.mat <- cophenetic(tree)         ###to ensure correct order
trait <- trait[row.names(p.dist.mat),] ###to ensure correct order
tree$tip.label <- str_replace_all(tree$tip.label, "_", " ")
rownames(trait) <- tree$tip.label


# test phylogenetic correlation
spp <- trait$species
corLambda <- corPagel(value=1, phy=tree, form=~spp)

m <- gls(scale(log(rd)) ~ scale(sqrt(abs(p50))), data = trait, correlation = corLambda)  
summary(m)


# test phylogenetic signal
p50 <- setNames(sqrt(abs(trait$p50)), trait$species)
rd <- setNames(log(trait$rd), trait$species)

p50_sig <- phylosig(tree, p50, method = "lambda", test = T)
rd_sig <- phylosig(tree, rd, method = "lambda", test = T)

