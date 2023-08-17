library(tidyverse)
library(mgcv)


source("code/transformers.R")


# Read in model and data --------------------------------------------------

m <- readRDS("results/models/trait_model_slim.rds")

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

trait <- readRDS("data/cleaned/trait_data_sPlot.rds") %>%
  dplyr::select(species, p50, rd_max) %>%
  na.omit()

plot <- readRDS("data/cleaned/plot_clim_data.rds")

gc()

# Prepare stuff for plotting ----------------------------------------------

post <- rmvn(2000, coef(m)[1:37], vcov(m)[1:37,1:37]) %>%
  as_tibble() %>%
  setNames(names(coef(m))[1:37]) %>%
  janitor::clean_names()

coef <- coef(m)
eco <- coef[str_detect(names(coef), "eco")]
spp <- coef[str_detect(names(coef), "spp")]

eco_int <- weighted.mean(eco, table(d$eco))
spp_int <- weighted.mean(spp, table(d$spp))
re_int <- eco_int + spp_int

arid_r <- quantile(d$arid, c(0.01, 0.99))
arid_seq <- seq(arid_r[1], arid_r[2], length.out = 50)
arid_levs <- quantile(plot$arid, c(0.05, 0.95)) %>% transform("arid_sqrt")

wtd_r <- quantile(d$wtd, c(0.01, 0.99))
wtd_seq <- seq(wtd_r[1], wtd_r[2], length.out = 50)
wtd_levs <- quantile(plot$wtd, c(0.05, 0.95)) %>% transform("wtd_log")

ps_r <- quantile(d$ps, c(0.01, 0.99))
ps_seq <- seq(ps_r[1], ps_r[2], length.out = 50)
ps_levs <- quantile(plot$precip_s, c(0.05, 0.95)) %>% transform("ps_log")

p50_seq <- quantile(d$p50, seq(0.01, 0.99, length.out = 50))
p50_levs <- quantile(trait$p50, c(0.05, 0.95)) %>% transform("P50_sqrt")

p50_labs <- c(-2, -4, -6, -8)
p50_breaks <- transform(p50_labs, "P50_sqrt")

rd_seq <- quantile(d$rd, seq(0.01, 0.99, length.out = 50))
rd_levs <- quantile(trait$rd_max, c(0.05, 0.95)) %>% transform("rd_log")

label_p50 <- function(x, bt = TRUE) {
  if(bt) x <- backtransform(x, "P50_sqrt")
  lab <- x %>% round(1) %>% paste("P50 =", ., "MPa") %>% fct_reorder(x)
  return(lab)
}

label_rd <- function(x, bt = TRUE) {
  if(bt) x <- backtransform(x, "rd_log")
  lab <- x %>% round(1) %>% paste("Rooting depth =", ., "m") %>% fct_reorder(-x)
  return(lab)
}

label_arid <- function(x, bt = TRUE) {
  if(bt) x <- backtransform(x, "arid_sqrt")
  lab <- x %>% round(1) %>% paste("Aridity Index =", .) %>% fct_reorder(x)
  lab <- ifelse(x == min(x), paste(lab, "(arid)"), paste(lab, "(humid)"))
  return(lab)
}

label_wtd <- function(x, bt = TRUE) {
  if(bt) x <- backtransform(x, "wtd_log")
  lab <- x %>% round(1) %>% paste("WTD =", ., "m") %>% fct_reorder(-x)
  return(lab)
}

label_ps <- function(x, bt = TRUE) {
  if(bt) x <- backtransform(x, "ps_log")
  lab <- x %>% round(1) %>% paste("Precip seasonality =", .) %>% fct_reorder(-x)
  return(lab)
}

# axis labels
p50_lab <- bquote(P[50]~(MPa))
rd_lab <- "Rooting depth (m)"
wtd_lab <- "Water table depth (m)"
ps_lab <- "Precipitation seasonality (CV)"

gc()


# Single trait effects -----------------------------------------------------------

## P50 ----

grid <- expand.grid(p50 = p50_seq,
                    arid = arid_levs,
                    wtd = wtd_levs, 
                    ps = ps_levs) %>%
  as_tibble() %>%
  mutate(p502 = p50^2, 
         y = 1)

X <- model.matrix(y ~ arid*wtd*ps*p50 + p502, data = grid) %>% 
  as_tibble() %>% 
  janitor::clean_names() %>%
  as.matrix()

beta <- t(post[colnames(X)])

pred <- plogis(X %*% beta + re_int)

df <- grid %>%
  mutate(arid = label_arid(arid),
         wtd = label_wtd(wtd),
         ps = ifelse(ps == min(ps), "Low seasonality (CV = 11.6)", "High seasonality (CV = 49.1)") %>% fct_reorder(ps),
         est = apply(pred, 1, mean),
         lo = apply(pred, 1, quantile, 0.05),
         hi = apply(pred, 1, quantile, 0.95))

labs <- df %>%
  mutate(p50 = max(p50), est = max(hi)*1.05) %>%
  distinct(arid, wtd, p50, est) %>%
  mutate(panel = c("C", "D", "A", "B"))

p50 <- ggplot(df, aes(x = p50, y = est)) + 
  geom_path(aes(color = ps)) + 
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = ps), alpha = 0.2) +
  geom_text(data = labs, aes(label = panel), fontface = "bold", size = 4) +
  scale_x_reverse(breaks = p50_breaks, labels = p50_labs) +
  scale_y_continuous(expand = expansion(mult = 0.07)) +
  scale_color_manual(values = c("#56B4E9", "#E69F00")) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00")) +
  facet_grid(wtd~arid) +
  labs(x = p50_lab, 
       y = "Occurrence probability", 
       color = "", fill = "") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = "bottom", 
        aspect.ratio = 0.8)



## Rooting depth ----

grid <- expand.grid(rd = rd_seq,
                    arid = arid_levs,
                    wtd = wtd_levs, 
                    ps = ps_levs) %>%
  as_tibble() %>%
  mutate(rd2 = rd^2, 
         y = 1)

X <- model.matrix(y ~ arid*wtd*ps*rd + rd2, data = grid) %>% 
  as_tibble() %>% 
  janitor::clean_names() %>%
  as.matrix()

beta <- t(post[colnames(X)])

pred <- plogis(X %*% beta + re_int)

df <- grid %>%
  mutate(rd = backtransform(rd, "rd_log"),
         arid = label_arid(arid),
         wtd = label_wtd(wtd),
         ps = ifelse(ps == min(ps), "Low seasonality (CV = 11.6)", "High seasonality (CV = 49.1)") %>% fct_reorder(ps),
         est = apply(pred, 1, mean),
         lo = apply(pred, 1, quantile, 0.05),
         hi = apply(pred, 1, quantile, 0.95))

labs <- df %>%
  mutate(rd = min(rd), est = max(hi)*1.05) %>%
  distinct(arid, wtd, rd, est) %>%
  mutate(panel = c("H", "G", "E", "F"))

rd <- ggplot(df, aes(x = rd, y = est)) + 
  geom_path(aes(color = ps)) + 
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = ps), alpha = 0.2) +
  geom_text(data = labs, aes(label = panel), fontface = "bold", size = 4.5) +
  facet_grid(wtd~arid) +
  scale_x_log10() +
  scale_y_continuous(expand = expansion(mult = 0.07)) +
  scale_color_manual(values = c("#56B4E9", "#E69F00")) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00")) +
  labs(x = rd_lab, 
       y = "Occurrence probability", 
       color = "", fill = "") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = "bottom", 
        aspect.ratio = 0.8)


# Combine figures

ggpubr::ggarrange(p50, rd, nrow = 2, common.legend = T, legend = "bottom")

ggsave("results/figures/single_trait_effects.jpg", height = 9, width = 5.5)
