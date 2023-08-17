library(tidyverse)
library(lme4)



# Read in and prepare RSIP rooting depth and texture data --------------------------------

rsip <- read_csv("data/traits/rsip_w_texture.csv") %>%
  rename_all(tolower) %>%
  rename(texture = 3)

rd <- rsip %>% 
  rename_all(tolower) %>%
  rename(texture = 3) %>%
  mutate(binomial = str_extract(species, "[A-Za-z]+ [A-Za-z-]+"), 
         texture = factor(texture, levels = c("Wl", "F", "M", "C"))) %>%
  select(binomial, texture, dr) %>%
  filter(dr > 0.1)

# get accepted species names
spp_lookup <- read_csv("data/traits/rd_species_lookup.csv") %>% select(-flag)
rd <- left_join(rd, spp_lookup) %>% na.omit()

# keep species in 
keep_spp <- readRDS("data/cleaned/trait_data.rds") %>% distinct(species_wfo)

rd <- inner_join(rd, keep_spp) %>%
  mutate(log_rd = log(dr))


# Fit model ---------------------------------------------------------------

fit <- lmer(log_rd ~ texture + (1|species_wfo), data = rd)
summary(fit)




# Get expected mean rooting depth by soil texture -------------------------

# function to get predicted rooting depth for each soil texture level
preds <- function(x) predict(x, newdata = data.frame(texture = unique(rd$texture)), re.form = NA)

# calculate parametric bootstrap means and confidence intervals
soil_boot <- bootMer(fit, preds, nsim = 1e3, parallel = "multicore", ncpus = 10)
rd_pred <- exp(soil_boot$t)
df <- tibble(texture = c("Coarse", "Medium", "Fine", "Waterlogged") %>% fct_relevel("Fine", "Medium"),
             est = apply(rd_pred, 2, mean),
             se = apply(rd_pred, 2, sd),
             lo = est-se,
             hi = est+se)


# Plot soil texture effects -----------------------------------------------

ggplot(df, aes(x = texture, y = est, ymin = lo, ymax = hi)) + 
  geom_pointrange(size = 0.6, linewidth = 1, color = 4) + 
  labs(x = "Soil texture", y = "Rooting depth (m)") +
  theme_bw() +
  theme(aspect.ratio = 1)

ggsave("results/figures/soil_texture_effects.jpg")
