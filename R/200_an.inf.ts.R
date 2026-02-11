# 200_an.inf.ts.R ####
### Multivariate analysis of intertidal invertebrate time series data

# Set up ####
## load packages ####
ld_pkgs <- c("tidyverse","ggplot2","vegan","ggdendro",
             "mvabund",
             "tictoc"
)
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)

tictoc::tic.clearlog();tic("SET UP")
print("Setting up")
## set metadata ####
source("R/00_meta_setMeta.R")

## load data ####
source("R/100_imp.inf0.R")
toc(log=TRUE)

names(dfw0)
dfw0 %>% 
  # convert to long
  pivot_longer(cols = -c(year,transect,shore,rep,
                         zone1,mesh,core.area_m2,Flag)) %>% 
  ## drop Fragments
  mutate(tmp = case_when(
    is.na(Flag)~"none",
    Flag == "#frag"~"#frag",
    Flag == "#agg" ~ "#agg",
    Flag == "#juv" ~"juv",
    TRUE ~ NA
  )) %>% select(-Flag) %>% 
  filter(tmp !="#frag") %>% 
  # calculate mean by station
  ### drop rep
  select(-rep) %>% 
  group_by(across(-c(value))) %>% 
  summarise(value = mean(value),.groups = "drop") %>% 
  ungroup() -> dfl

## plot Scolelepis 
dfl %>% 
  filter(value!=0) %>% 
  filter(grepl("^Scolelepis",name)) %>% 
  ggplot(.,
         aes(
           x=year,
           y=value,
         )
         )+
  geom_point()+
  geom_smooth()+
  facet_grid(shore~zone1)
