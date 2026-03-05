## 100_imp.epi0.R ##
## Import epifaunal & Crangon data from current year and combine with historic

# Set up ####
### load packages ####
ld_pkgs <- c("tidyverse","readxl")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)

### load metadata ####
source("R/00_meta_setMeta.R")

# copy & load data ####

file_name <- "epi.ts.long.xlsx"

# Set the source and destination file paths
source_file <- paste0(fol,file_name)
destination_file <- paste0("data/",file_name)

# Check if the file exists in the data folder
if (!file.exists(destination_file)) {
  # If not, copy the file from the source folder to the data folder
  file.copy(source_file, destination_file)
  cat("File copied successfully.\n")
} else {
  # If the file already exists in the data folder, do nothing
  cat("File already exists in the data folder.\n")
}
rm(file_name,source_file)

df0 <- as_tibble(read_xlsx(destination_file,
                           sheet = "epi.ts.long",
                           guess_max = 50000))#blanks at top of data means read in as logical
df_cra_ts <- as_tibble(read_xlsx(destination_file,
                                 sheet = "df.cra.ts",
                                 guess_max = 50000#blanks at top of data means read in as logical
)
)

df0 %>% 
  ## remove taxa flagged as 'Remove'
  filter(!grepl("^Rem", Flag)) %>% 
  #filter(., Kingdom == "Animalia") %>% 
  mutate(.,abundUSE = as.numeric(ifelse(abund == "x",0,abund))) %>% #View()
  dplyr::select(., -c(taxonRecorded,Kingdom:Comments)) %>% 
  group_by(across(-abundUSE)) %>% 
  summarise(abundUSE = sum(abundUSE), .groups = "drop") %>% 
  pivot_wider(
    .,
    names_from = taxonUse,
    values_from = abundUSE,
    values_fill = list(abundUSE = 0)
  ) %>%
  ungroup() %>% 
  relocate("AFAUNAL", .after = Flag)-> dfw

## export data for reporting
df0 %>% 
  filter(year == cur.yr) %>% 
  # remove superfluous columns
  select(-c(zone1, taxonRecorded, Comments,year)) %>% 
  mutate(code=paste0(transect,"_",depth,"_",mon)) %>%
  relocate(code) %>% 
  select(-c(transect,depth,mon,Flag,DataSource)) %>% 
  # group and sum across lab and vessel
  group_by(across(-abund)) %>% 
  summarise(abund = sum(abund),.groups = "drop") %>% ungroup() %>% 
  pivot_wider(names_from = code,
              values_from = abund,
              values_fill = 0) %>% 
  write.csv(.,file = "output/epi_cur_taxa_report.csv",row.names = FALSE)
  
# TIDY UP ####
# rm(list = ls(pattern = "^df"))
rm(list = ls(pattern = "^cb"))
rm(cur.yr,destination_file,file_name,fol,gisfol,perm,ppi,source_file)

detach(package:readxl, unload=TRUE)
detach(package:tidyverse, unload=TRUE)
