# 100_imp_int_inv_lab.R ####
## Import data from lab returns and process for inclusion into master data

# Set up ####
### load packages ####
ld_pkgs <- c("tidyverse","tictoc")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)

tictoc::tic.clearlog();tictoc::tic("SET UP & LOAD DATA")

### load metadata ####
source("R/00_meta_setMeta.R")

# load data
## abundance
df0_abund <- readxl::read_xlsx(
  paste0(
    file.path(
      dirname(fol)),
    # "/2025/LabReturns/Int_Wash_forR.xlsx"),
    "/2025/LabReturns/Int_Linc_forR.xlsx"),
  sheet = "Abundance")

## biomass
df0_biomass <- readxl::read_xlsx(
  paste0(
    file.path(
      dirname(fol)),
    # "/2025/LabReturns/Int_Wash_forR.xlsx"),
    "/2025/LabReturns/Int_Linc_forR.xlsx"),
  sheet = "Biomass")

tictoc::toc(log=TRUE)

# convert to long and remove zero values ####
## abundance
str(df0_abund)
df0_abund %>% 
  ##lengthen data
  tidyr::pivot_longer(cols = -c(Transect,Shore,Rep,code,method),
               values_transform = as.character) %>% 
  ## remove zero values
  dplyr::filter(value !="0") -> df_abund_l

write.csv(
  df_abund_l,
  file = paste0(
    file.path(
      dirname(fol)),
    "/2025/LabReturns/df_abund_l.csv"
    ))

## biomass
df0_biomass %>% 
  ##lengthen data
  tidyr::pivot_longer(cols = -c(Transect,Shore,Rep,code,method),
                      values_transform = as.character) %>% 
  ## remove zero values
  dplyr::filter(value !="0") -> df_biomass_l

write.csv(
  df_biomass_l,
  file = paste0(
    file.path(
      dirname(fol)),
    "/2025/LabReturns/df_biomass_l.csv"
  ))
