library(tidyverse)
library(mgcv)

trait <- readRDS("results/models/trait_model_slim.rds")
no_trait <- readRDS("results/models/no_trait_model_slim.rds")

lrt <- anova(no_trait, trait, test = "Chisq")

aic_tab <- AIC(no_trait, trait) %>% 
  rownames_to_column() %>% 
  rename(model = 1) %>%
  arrange(AIC) %>%
  mutate(AICdif = AIC - AIC[1])

lrt <- lst(lrt, aic_tab)
saveRDS(lrt, "results/model_comparison.rds")
