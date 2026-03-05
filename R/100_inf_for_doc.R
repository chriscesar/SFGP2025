# 100_inf_for_doc.R ####
### wrangle infaunal data for doc appendices

## load abundance & biomass time series data
source("R/100_imp.inf0.R")

### load packages
ld_pkgs <- c("tidyverse")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)

dfw0 %>% select(.,-Flag) %>% 
  filter(year == cur.yr) %>% 
  select(-mesh,-core.area_m2) %>% 
  #pivot longer, remove zeroes and rewiden
  pivot_longer(cols = -c(year,transect,shore,rep,zone1)) %>% 
  filter(value != 0) %>% 
  # reduce duplicate taxa
  group_by(across(-value)) %>% 
  summarise(value = sum(value),.groups = "drop") %>% 
  #rewiden
  pivot_wider(names_from = name, values_from = value,values_fill = 0) %>% 
  # pivot longer & calc mean across reps
  pivot_longer(cols = -c(year,transect,shore,rep,zone1)) %>% 
  select(-rep) %>% 
  group_by(across(-value)) %>% 
  summarise(value = mean(value),.groups = "drop") %>% 
  filter(value !=0) %>% 
  #rewiden for export
  pivot_wider(names_from = name,values_from = value,values_fill = 0) -> dfw
write.csv(dfw,
          file = "output/inf_tax.csv",row.names = FALSE)
