# inf_gllvm.R ####

# Set up ####
## load packages ####
ld_pkgs <- c("tidyverse","ggplot2","vegan","ggdendro",#data vis
             "dendextend",#data vis
             "ggtext",#data vis
             "ggpp",
             "mvabund",
             "Hmsc",
             "ggpubr",
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
  ungroup() %>% 
  ##OPTIONAL: DROP WASH SAMPLES
  filter(zone1 != "Wash") %>% 
  ##remove zeroes
  filter(value != 0) %>% 
  #widen
  pivot_wider(names_from = name,values_from = value,
              values_fill = 0)-> dfw
names(dfw)

##to do: generate function to remove any taxa which appear in very few samples
##and run gllvm models to test role of zone1 + (1|shore) on abundance