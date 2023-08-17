library(tidyverse)
library(ggExtra)

set.seed(123)


traits <- readRDS("data/cleaned/trait_data.rds")
pheno <- read_csv("data/traits/leaf_phenology.csv")

d <- traits %>% 
  select(species = species_wfo, p50, rd_max) %>% 
  left_join(pheno) %>%
  distinct() %>% 
  na.omit() %>%
  mutate(p50_t = sqrt(-p50), 
         pheno = ifelse(leaf_phenology == "d", "Deciduous", "Evergreen"))

p50_labs <- -seq(0.5, 4, by = 0.5)^2
p50_breaks <- sqrt(-p50_labs)

# sample labeled species using bins and thresholds
ss <- d %>%
  sample_n(n()) %>%
  mutate(p50_bin = ntile(p50_t, 7), 
         rd_bin = ntile(log(rd_max), 7)) %>%
  group_by(p50_bin, rd_bin) %>%
  mutate(bin_n = 1:n()) %>%
  ungroup() %>%
  filter(bin_n == 1 | p50 > -1.5 | p50 < -6.5 | rd_max < 0.9 | rd_max > 12)

# quantiles
q <- expand.grid(p50_t = sqrt(c(1.5, 6.9)), 
                 rd_max = c(0.6, 12.9)) %>%
  mutate(label = c("Vulnerable confronters",
                   "Resistant confronters",
                   "Vulnerable avoiders",
                   "Resistant avoiders"), 
         pheno = "Deciduous")


p <- d %>% 
  ggplot(aes(x = p50_t, y = rd_max, shape = pheno)) + 
  geom_hline(yintercept = median(d$rd_max), color = "gray") +
  geom_vline(xintercept = median(d$p50_t), color = "gray") +
  geom_point(alpha = 0.25, aes(color = pheno)) + 
  scale_x_reverse(breaks = p50_breaks, labels = p50_labs) +
  scale_y_log10(limits = range(d$rd_max)) + 
  #geom_text(aes(label = species), fontface = 3, size = 3) +
  labs(y = "Rooting depth (m)") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = bquote(P[50]~(MPa)), shape = NULL, color = NULL) +
  ggrepel::geom_text_repel(data = ss, aes(label = species, color = pheno), size = 3, max.overlaps = 5, 
                           fontface = 3, min.segment.length = 0.2, box.padding = 0.1, show.legend = F) +
  geom_point(data = ss, aes(color = pheno)) +
  geom_point(data = q, size = 4, color = "gray40", show.legend = F) +
  scale_color_hue(l = 30) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12)) 

p


svg("results/figures/p50_rd_space_pheno.svg", height = 9.5, width = 9)
ggMarginal(p, size = 12, groupFill = T, color = NA)
dev.off()

