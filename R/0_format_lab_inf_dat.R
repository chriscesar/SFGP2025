# 0_format_lab_inf_dat.R #
# rejig lab return data for addition to master data

# lab return sheet with Abundanc AND Biomass have been transposed
# Aphia ID's as column headers (format: paste0("aph_",AphiaID)),
# with sample names as row headers
# there were duplicated column names; e.g., accounting for 'whole' and 'fragments'
# of taxa. 'P' values of such taxa have been removed and Biomass values have been
# summed. Each taxon should have only a single column.
# 'Empty' Aphia ID's have been deleted.
# Samples with no taxa given 'AFAUNAL' tag.
# Occurrences of "-" in sample name removed, e.g. "T1N Low A - 0.5 mm".  This causes
# issues with opening in Excel after dplyr::separate() function

# Set up ####
## load packages ####
ld_pkgs <- c("tidyverse","readxl")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)

## set metadata ####
source("R/00_meta_setMeta.R")

## load data ####
# as_tibble(read_xlsx(paste0(projfol,"2024/data/tmpFiles/2024_inf_for_formatting.xlsx"),
#                     sheet=1)) %>% 
#   mutate(across(everything(), as.character)) %>% 
#   pivot_longer(., cols=-Sample,
#                names_to = "AphiaID",
#                values_to = "Value",
#                #names_repair = "unique"
#                ) %>% 
#   filter(., !is.na(Value)) %>% 
#   separate(., Sample, sep=" ",
#            into=c("Col1",
#                   "Col2",
#                   "Col3",
#                   "Col4",
#                   "Col5",
#                   "Col6",
#                   "Col7",
#                   "Col8"
#                   )) %>% 
#   write.csv(., file = paste0(projfol,"2024/data/tmpFiles/tmplong.csv"),
#             row.names = FALSE)
######

as_tibble(read_xlsx(paste0(projfol,"2024/data/tmpFiles/inf_tmp.xlsx"),
                    sheet=1)) %>% 
  mutate(across(everything(), as.character)) %>%
  pivot_longer(., cols=-Sample,
               names_to = "AphiaID",
               values_to = "Value",
               names_repair = "minimal"
  ) %>% 
  filter(., !is.na(Value)) %>% 
  separate(., Sample, sep=" ",
           into=c("Col1",
                  "Col2",
                  "Col3",
                  "Col4",
                  "Col5",
                  "Col6",
                  "Col7",
                  "Col8"
           )) -> ee 

ee$AphiaUSE <- gsub("\\.\\.\\.\\d+$", "", ee$AphiaID)

write.csv(ee, file = paste0(projfol,"2024/data/tmpFiles/tmplong.csv"),
          row.names=FALSE)
