library(tidyverse)
library(mgcv)

source("code/transformers.R")


d <- readRDS("data/cleaned/p50_rd_data_ecoregion.rds")

ecoregion <- distinct(d, ecoregion) %>%
  mutate(eco = factor(ecoregion))

species <- distinct(d, species) %>%
  mutate(spp = factor(species))

d <- d %>% 
  left_join(ecoregion) %>% 
  left_join(species) %>% 
  select(pres, arid, wtd, precip_s, p50, rd = rd_max, spp, eco) %>%
  na.omit() %>%
  mutate(arid = transform(arid, "arid_sqrt"),
         arid2 = arid^2,
         wtd = transform(wtd, "wtd_log"),
         wtd2 = wtd^2,
         ps = transform(precip_s, "ps_log"), 
         ps2 = ps^2,
         p50 = transform(p50, "P50_sqrt"),
         p502 = p50^2,
         rd = transform(rd, "rd_log"),
         rd2 = rd^2) 

gc()

start <- Sys.time()

m <- bam(pres ~ arid*wtd*ps*p50*rd + arid2 + wtd2 + ps2 + p502 + rd2 + s(spp, bs = "re") + s(eco, bs = "re"), 
         family = "binomial", data = d, discrete = T)

end <- Sys.time()
print(end-start)

saveRDS(m, "results/models/trait_model.rds")


m_slim <- m
m_slim$model <- NULL
m_slim$linear.predictors <- NULL
m_slim$fitted.values <- NULL
m_slim$wt <- NULL
m_slim$y <- NULL
m_slim$prior.weights <- NULL
m_slim$offset <- NULL
m_slim$weights <- NULL
m_slim$residuals <- NULL

saveRDS(m_slim, "results/models/trait_model_slim.rds")
