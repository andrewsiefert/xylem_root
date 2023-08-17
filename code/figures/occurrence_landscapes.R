library(tidyverse)
library(mgcv)
library(aws)
library(terra)
library(ks)
library(modelr)

source("code/transformers.R")



# Setup -------------------------------------------------------------------

# read in model
m <- readRDS("results/models/trait_model_slim.rds")

# read in data
d <- readRDS("data/cleaned/p50_rd_data_ecoregion.rds")

# convert ecoregion and species to factors
ecoregion <- distinct(d, ecoregion) %>%
  mutate(eco = factor(ecoregion))

species <- distinct(d, species) %>%
  mutate(spp = factor(species))

# prepare data for modeling
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

# read in trait data
trait <- readRDS("data/cleaned/trait_data_sPlot.rds") %>%
  dplyr::select(species, p50, rd_max) %>%
  na.omit()

# read in plot climate data
plot <- readRDS("data/cleaned/plot_clim_data.rds")


# sample from posterior distribution
post <- rmvn(2000, coef(m)[1:37], vcov(m)[1:37,1:37]) %>%
  as_tibble() %>%
  setNames(names(coef(m))[1:37]) %>%
  janitor::clean_names()

# extract ecoregion and species random effects
coef <- coef(m)
eco <- coef[str_detect(names(coef), "eco")]
spp <- coef[str_detect(names(coef), "spp")]

# calculate weighted mean ecoregion and species random effects
eco_int <- weighted.mean(eco, table(d$eco))
spp_int <- weighted.mean(spp, table(d$spp))
re_int <- eco_int + spp_int

# get aridity values for plotting
arid_r <- quantile(d$arid, c(0.01, 0.99))
arid_seq <- seq(arid_r[1], arid_r[2], length.out = 50)
arid_levs <- quantile(plot$arid, c(0.05, 0.95)) %>% transform("arid_sqrt")

# get water table depth values for plotting
wtd_r <- quantile(d$wtd, c(0.01, 0.99))
wtd_seq <- seq(wtd_r[1], wtd_r[2], length.out = 50)
wtd_levs <- quantile(plot$wtd, c(0.05, 0.95)) %>% transform("wtd_log")

# get precipitation seasonality values for plotting
ps_r <- quantile(d$ps, c(0.01, 0.99))
ps_seq <- seq(ps_r[1], ps_r[2], length.out = 50)
ps_levs <- quantile(plot$precip_s, c(0.05, 0.95)) %>% transform("ps_log")

# get P50 values for plotting
p50_seq <- quantile(d$p50, seq(0.01, 0.99, length.out = 50))
p50_levs <- quantile(trait$p50, c(0.05, 0.95)) %>% transform("P50_sqrt")

p50_labs <- c(-2, -4, -6, -8)
p50_breaks <- transform(p50_labs, "P50_sqrt")

# get rooting depth values for plotting
rd_seq <- quantile(d$rd, seq(0.01, 0.99, length.out = 50))
rd_levs <- quantile(trait$rd_max, c(0.05, 0.95)) %>% transform("rd_log")

# functions to create plot labels
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

# colors for masks in 2D plots
t_col <- function(color, percent = 50) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3], max = 255, 
               alpha = (100 - percent) * 255 / 100)
  invisible(t.col)
}

clear <- t_col("white", 100)
mid <- t_col("white", 25)


# Plot trait effects -----------------------------------------------------------

# create convex hull of trait values
tr <- d %>% 
  dplyr::select(p50, rd) %>% 
  distinct() %>%
  rename(p50 = p50) %>%
  mutate(p50 = sqrt(abs(backtransform(p50, "P50_sqrt"))),
         rd = log(backtransform(rd, "rd_log")))

hull_pts <- chull(tr)
hull <- tr[c(hull_pts, hull_pts[1]),] %>% 
  as_tibble()

p <- sp::Polygon(hull)
ps <- sp::Polygons(list(p),1)
sps <- sp::SpatialPolygons(list(ps))

# create mask to highlight trait convex hull
mask <- expand.grid(x = modelr::seq_range(tr$p50, 200), 
                    y = modelr::seq_range(tr$rd, 200)) %>%
  mutate(value = 1) %>%
  raster::rasterFromXYZ() %>%
  raster::crop(raster::extent(sps)) %>%
  raster::mask(sps, inverse = T) %>%
  as.data.frame(xy = T) %>%
  as_tibble() %>%
  mutate(x = x, y = exp(y))

# get predictions for trait-environment combinations
grid <- expand.grid(arid = arid_levs, 
                    wtd = wtd_levs,
                    ps = ps_levs,
                    p50 = modelr::seq_range(d$p50, 50, expand = 0.01), 
                    rd = modelr::seq_range(transform(exp(tr$rd), "rd_log"), 50, expand = 0.01)) %>%
  as_tibble() %>%
  mutate(arid2 = arid^2,
         wtd2 = wtd^2,
         ps2 = ps^2,
         p502 = p50^2,
         rd2 = rd^2, 
         y = 1)

X <- model.matrix(y ~ arid*wtd*ps*p50*rd + arid2 + wtd2 + ps2 + p502 + rd2, data = grid) %>% 
  as_tibble() %>% 
  janitor::clean_names() %>%
  as.matrix()

beta <- t(post[colnames(X)])

pred <- plogis(X %*% beta + re_int) %>% rowMeans()

# prepare dataframe for plotting
df <- grid %>%
  mutate(p50 = backtransform(p50, "P50_sqrt") %>% abs() %>% sqrt(),
         rd = backtransform(rd, "rd_log"),
         arid = label_arid(arid),
         wtd = label_wtd(wtd),
         ps = label_ps(ps),
         pres = pred)

p50_labs <- -seq(0.5, 4, by = 0.5)^2
p50_breaks <- sqrt(-p50_labs)


## Low seasonality ----
low <- df %>% 
  filter(ps == "Precip seasonality = 11.6") %>%
  ggplot(aes(x = p50, y = rd)) + 
  geom_raster(aes(fill = pres %>% pmax(0.0001) %>% pmin(0.3))) + 
  geom_contour(aes(z = pres), color = "white", breaks = seq(0, 1, 0.02)) +
  scale_fill_viridis_c(trans = "log", breaks = c(0.001, 0.01, 0.1)) +   
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        aspect.ratio = 1) +
  labs(title = "Low seasonality (CV = 12)", x = p50_lab, y = rd_lab, fill = "Probability of\noccurrence") +
  ggnewscale::new_scale_fill() +
  geom_raster(data = mask, aes(x = x, y = y, fill = value), show.legend = F) +
  scale_fill_gradient2(low = 1, mid = mid, high = 2, na.value = clear) + 
  facet_grid(wtd~arid) +
  geom_point(data = tr %>% mutate(rd = exp(rd)), alpha = 0.2) +
  scale_x_continuous(expand=c(0,0), limits = c(max(tr$p50), min(tr$p50)),
                     oob = scales::squish_infinite, labels = p50_labs, breaks = p50_breaks) +
  scale_y_log10(expand=c(0,0), limits = range(df$rd), oob = scales::squish_infinite) 

## High seasonality ----
high <- df %>% 
  filter(ps == "Precip seasonality = 49.1") %>%
  ggplot(aes(x = p50, y = rd)) + 
  geom_raster(aes(fill = pres %>% pmax(0.0001) %>% pmin(0.3))) + 
  geom_contour(aes(z = pres), color = "white", breaks = seq(0, 1, 0.02)) +
  scale_fill_viridis_c(trans = "log", breaks = c(0.001, 0.01, 0.1)) +   
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        aspect.ratio = 1) +
  labs(title = "High seasonality (CV = 49)", x = p50_lab, y = rd_lab, fill = "Probability of\noccurrence") +
  ggnewscale::new_scale_fill() +
  geom_raster(data = mask, aes(x = x, y = y, fill = value), show.legend = F) +
  scale_fill_gradient2(low = 1, mid = mid, high = 2, na.value = clear) + 
  facet_grid(wtd~arid) +
  geom_point(data = tr %>% mutate(rd = exp(rd)), alpha = 0.2) +
  scale_x_continuous(expand=c(0,0), limits = c(max(tr$p50), min(tr$p50)),
                     oob = scales::squish_infinite, labels = p50_labs, breaks = p50_breaks) +
  scale_y_log10(expand=c(0,0), limits = range(df$rd), oob = scales::squish_infinite) 

## Combine----
ggpubr::ggarrange(low, high, nrow = 2, common.legend = T, legend = "right")
ggsave("results/figures/p50_rd_effects.jpg", height = 10, width = 6)



# Environmental responses -------------------------------------------------


p50_levs2 <- quantile(trait$p50, c(0.1, 0.9)) %>% transform("P50_sqrt")
rd_levs2 <- quantile(trait$rd_max, c(0.1, 0.9)) %>% transform("rd_log")

env <- distinct(d, arid, wtd)

n <- 500
bw <- Hpi(cbind(env$arid, env$wtd))*4

de <- kde(cbind(env$arid, env$wtd), H = bw, gridsize = n)

r <- terra::rast(t(de$estimate)[n:1,])
ext(r) <- c(min(env$arid), max(env$arid), min(env$wtd), max(env$wtd))

cont <- as.contour(r, nlevels = 100)

poly <- as.polygons(cont[cont$level == cont$level[2]])
poly2 <- as.polygons(cont[cont$level == cont$level[1]])

r2 <- rast(ncols=200, nrows=200, extent = ext(poly2), vals = 1)
mask <- rasterize(poly, r2, fun = sum) %>%
  as.data.frame(xy = T, na.rm = F) %>%
  as_tibble() %>%
  rename(arid = 1, wtd = 2, mask = 3)


grid <- expand.grid(p50 = p50_levs2,
                    rd = rd_levs2,
                    arid = seq_range(mask$arid, 50),
                    wtd = seq_range(mask$wtd, 50),
                    ps = ps_levs) %>%
  na.omit() %>%
  #filter(arid < 3.5) %>%
  as_tibble() %>%
  mutate(arid2 = arid^2,
         wtd2 = wtd^2, 
         p502 = p50^2,
         rd2 = rd^2,
         ps2 = ps^2,
         y = 1)

X <- model.matrix(y ~ arid*wtd*ps*p50*rd + arid2 + wtd2 + ps2 + p502 + rd2, data = grid) %>% 
  as_tibble() %>% 
  janitor::clean_names() %>%
  as.matrix()

beta <- t(post[colnames(X)])

pred <- plogis(X %*% beta + re_int) %>% rowMeans()

df <- grid %>%
  mutate(p50_b = label_p50(p50),
         rd_b = label_rd(rd),
         arid = backtransform(arid, "arid_sqrt"),
         wtd = backtransform(wtd, "wtd_log") + 0.01,
         ps = label_ps(ps),
         pres = pred)

mask <- mask %>%
  mutate(arid = backtransform(arid, "arid_sqrt"),
         wtd = backtransform(wtd, "wtd_log") + 0.01)

## Low seasonality ---- 
low <- df %>%
  filter(ps == "Precip seasonality = 11.6") %>%
  ggplot(aes(x = arid, y = wtd)) + 
  geom_raster(aes(fill = pres %>% pmin(0.11))) +
  geom_contour(aes(z = pres), color = "white", breaks = seq(0, 1, 0.01)) +
  scale_fill_viridis_c(trans = "sqrt", breaks = c(0.01, 0.05, 0.1)) +   
  labs(title = "Low seasonality (CV = 12)", 
       x = "Aridity Index (P:PET)",
       y = "Water table depth (m)",
       fill = "Probability of\noccurrence") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = mask, aes(x = arid, y = wtd, fill = mask), show.legend = F) +
  scale_fill_gradient2(low = clear, mid = clear, high = clear, na.value = mid) + 
  facet_grid(rd_b~p50_b) +
  scale_x_sqrt(expand = c(0,0), limits = range(df$arid), oob = scales::squish_infinite) +
  scale_y_log10(expand = c(0,0), limits = range(df$wtd), oob = scales::squish_infinite,) +
  theme_bw()

## High seasonality ----
high <- df %>%
  filter(ps == "Precip seasonality = 49.1") %>%
  ggplot(aes(x = arid, y = wtd)) + 
  geom_raster(aes(fill = pres %>% pmin(0.11))) +
  geom_contour(aes(z = pres), color = "white", breaks = seq(0, 1, 0.01)) +
  scale_fill_viridis_c(trans = "sqrt", breaks = c(0.01, 0.05, 0.1)) +   
  labs(title = "High seasonality (CV = 49)", 
       x = "Aridity Index (P:PET)",
       y = "Water table depth (m)",
       fill = "Probability of\noccurrence") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = mask, aes(x = arid, y = wtd, fill = mask), show.legend = F) +
  scale_fill_gradient2(low = clear, mid = clear, high = clear, na.value = mid) + 
  facet_grid(rd_b~p50_b) +
  scale_x_sqrt(expand = c(0,0), limits = range(df$arid), oob = scales::squish_infinite) +
  scale_y_log10(expand = c(0,0), limits = range(df$wtd), oob = scales::squish_infinite,) +
  theme_bw()


## Combine ----
ggpubr::ggarrange(low, high, nrow = 2, common.legend = T, legend = "right")
ggsave("results/figures/arid_wtd_responses.jpg", height = 10, width = 6)

