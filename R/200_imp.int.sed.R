## 200_imp.int.sed.R ##

## generate plots from gradistat (G2Sd) outputs

# Set up ####
### load packages ####
ld_pkgs <- c("tidyverse","tictoc","ggthemes")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)

tictoc::tic.clearlog();tic("SET UP")

### load metadata ####
source("R/00_meta_setMeta.R")

# load data
df <- readRDS(file="output/models/grad_all_raw.Rdat")

# Produce plots ####

## Sediment type ####
df$stat$fowa %>% 
# df$sedim$texture %>% 
  ggplot(.,
         aes(
           x = year,
           # y = texture,
           y = skewness.fw.phi,
           shape = shore,
           colour = shore,
         )
         ) +
  geom_jitter() +
  facet_grid(.~zone1)
