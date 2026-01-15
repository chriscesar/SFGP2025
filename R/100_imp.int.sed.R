## 100_imp.int.sed.R ##
## Import intertidal sediment data for later assessment

# Set up ####
### load packages ####
ld_pkgs <- c("tidyverse","readxl","tictoc")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)

tictoc::tic.clearlog();tic("SET UP")

### load metadata ####
source("R/00_meta_setMeta.R")

ab <- c("T1N","T1","T1S")
inside <- c("T4","T7","T8","T11","T12")
inside2 <- "T13"
bel <- c("T15","T21","T22","T23","T24","T25","T26")
wash <- c("WA1","WA2","WA3","WA4","WA5","WA6")

df_sed <- readxl::read_xlsx(paste0(fol,"sed.data.ALL.USE.xlsx"), sheet = "AllDat") %>% 
  filter(., DetUse != "Remove: metadata") %>%  # drop unneeded rows
  mutate(., year = lubridate::year(SAMP_SAMPLE_DATE))
df_sed$zone1 <- ifelse(
  df_sed$Transect %in% ab,"Above",
  ifelse(
    df_sed$Transect %in% inside,"Inside",
    ifelse(
      df_sed$Transect %in% inside2,"Inside2",
      ifelse(
        df_sed$Transect %in% bel,"Below",
        ifelse(
          df_sed$Transect %in% wash,"Wash",NA
        )))))

df_sed$zone1 <- factor(df_sed$zone1, levels = c("Above","Inside","Inside2","Below","Wash"))
df_sed$Shore <- factor(df_sed$Shore, levels = c("Upper","Mid","Low","Surf"))
df_sed$Transect <- factor(df_sed$Transect,
                          levels = c(
                            "T1N", "T1","T1S",
                            "T4","T7", "T8", "T11", "T12",
                            "T13",
                            "T15", "T21", "T22", "T23","T24", "T25","T26",
                            "WA1","WA2","WA3","WA4","WA5","WA6"
                          ))

# load and append older data ####
df_sed_old <- read.csv(file = paste0(fol,"sed.psa.hi.ts.csv"))

df_sed_bulk <- readxl::read_xlsx(paste0(fol,"sed.psa.bulkWIP_use.xlsx"),
                                 sheet="sed.bulk.ts.out")
df_sed_bulk$transect <- factor(df_sed_bulk$transect, levels=c(
  "T1N","T1","T1S","T4","T11","T7", "T8", "T12", "T13", "T15", "T17",
  "T20","T21", "T22", "T23", "T24", "T25", "T26","WA1","WA2","WA3","WA4","WA5","WA6"
))

df_sed_bulk$shore <- factor(df_sed_bulk$shore,levels=c("Upper","Mid","Low","Surf"))
df_sed_bulk$zone1 <- factor(df_sed_bulk$zone1, levels=c("Above","Inside","Inside2","Below","Wash"))

names(df_sed)
names(df_sed_old)

#rm(ab,inside,inside2,bel,wash)
toc(log=TRUE)

##
tic("prep data for gradistat")
# prep data for gradistat ####
df_sed %>% 
  filter(.,str_detect(DetUse, 'phi')) %>% 
  # filter(.,year == cur.yr) %>%
  # filter(str_starts(Transect, "WA")) %>% #names(.)
  mutate("code" = paste0(year,".",Transect,".",Shore,".",method)) %>% #View(.)
  dplyr::select(.,code,sediment_um,MEAS_RESULT) %>% 
  # pivot_wider(.,names_from = code,values_from = MEAS_RESULT
  #             ) -> sed_grad
  tidyr::pivot_wider(

    names_from  = code,
    values_from = MEAS_RESULT,
    values_fn   = mean
  )->sed_grad0

## convert to DF and create row names ####

sed_grad0 <- as.data.frame(sed_grad0)
row.names(sed_grad0) <- sed_grad0$sediment_um
sed_grad0 <- sed_grad0[,-1]
toc(log=TRUE)

# Run Gradistat & save outputs ####
tic("Run Gradistat & generate factors")
grad_out <- granstat(sed_grad0)

## Index values ####
grad_out$index <- grad_out$index %>% 
  separate(
    col = samples,
    into = c("year", "transect", "shore", "method"),
    sep = "\\.",     # period is a regex special char, so escape it
    remove = FALSE,  # keep original 'samples' if you want
    fill = "right",  # if some rows are missing parts, fill with NA on the right
    extra = "merge"  # if there are extra pieces, merge them into the last column
    ) %>% 
  mutate(year = as.integer(year)) %>% 
  # assign zones
  mutate(zone1 = case_when(
    transect %in% ab ~ "Above",
    transect %in% bel ~ "Below",
    transect %in% inside ~ "Inside",
    transect %in% inside2 ~ "Inside2",
    transect %in% wash ~ "Wash",
  )) %>% 
  mutate(shore = factor(shore, levels = c("Upper","Mid","Low","Surf")),
         transect = factor(transect,levels = c(
           "T1N", "T1","T1S",
           "T4","T7", "T8", "T11", "T12",
           "T13",
           "T15", "T21", "T22", "T23","T24", "T25","T26",
           "WA1","WA2","WA3","WA4","WA5","WA6"
         )),
         zone1 = factor(zone1, levels = c("Above","Inside","Inside2",
                                          "Below","Wash"))) %>% 
  relocate(zone1, .after = "method")

## Stat values ####
### arith ####
grad_out$stat$arith <- grad_out$stat$arith %>% 
  separate(
    col = samples,
    into = c("year", "transect", "shore", "method"),
    sep = "\\.",     # period is a regex special char, so escape it
    remove = FALSE,  # keep original 'samples' if you want
    fill = "right",  # if some rows are missing parts, fill with NA on the right
    extra = "merge"  # if there are extra pieces, merge them into the last column
  ) %>% 
  mutate(year = as.integer(year)) %>% 
  # assign zones
  mutate(zone1 = case_when(
    transect %in% ab ~ "Above",
    transect %in% bel ~ "Below",
    transect %in% inside ~ "Inside",
    transect %in% inside2 ~ "Inside2",
    transect %in% wash ~ "Wash",
  )) %>% 
  mutate(shore = factor(shore, levels = c("Upper","Mid","Low","Surf")),
         transect = factor(transect,levels = c(
           "T1N", "T1","T1S",
           "T4","T7", "T8", "T11", "T12",
           "T13",
           "T15", "T21", "T22", "T23","T24", "T25","T26",
           "WA1","WA2","WA3","WA4","WA5","WA6"
         )),
         zone1 = factor(zone1, levels = c("Above","Inside","Inside2",
                                          "Below","Wash"))) %>% 
  relocate(zone1, .after = "method")

### geom ####
grad_out$stat$geom <- grad_out$stat$geom %>% 
  separate(
    col = samples,
    into = c("year", "transect", "shore", "method"),
    sep = "\\.",     # period is a regex special char, so escape it
    remove = FALSE,  # keep original 'samples' if you want
    fill = "right",  # if some rows are missing parts, fill with NA on the right
    extra = "merge"  # if there are extra pieces, merge them into the last column
  ) %>% 
  mutate(year = as.integer(year)) %>% 
  # assign zones
  mutate(zone1 = case_when(
    transect %in% ab ~ "Above",
    transect %in% bel ~ "Below",
    transect %in% inside ~ "Inside",
    transect %in% inside2 ~ "Inside2",
    transect %in% wash ~ "Wash",
  )) %>% 
  mutate(shore = factor(shore, levels = c("Upper","Mid","Low","Surf")),
         transect = factor(transect,levels = c(
           "T1N", "T1","T1S",
           "T4","T7", "T8", "T11", "T12",
           "T13",
           "T15", "T21", "T22", "T23","T24", "T25","T26",
           "WA1","WA2","WA3","WA4","WA5","WA6"
         )),
         zone1 = factor(zone1, levels = c("Above","Inside","Inside2",
                                          "Below","Wash"))) %>% 
  relocate(zone1, .after = "method")

### Folk Ward ####
grad_out$stat$fowa <- grad_out$stat$fowa %>% 
  separate(
    col = samples,
    into = c("year", "transect", "shore", "method"),
    sep = "\\.",     # period is a regex special char, so escape it
    remove = FALSE,  # keep original 'samples' if you want
    fill = "right",  # if some rows are missing parts, fill with NA on the right
    extra = "merge"  # if there are extra pieces, merge them into the last column
  ) %>% 
  mutate(year = as.integer(year)) %>% 
  # assign zones
  mutate(zone1 = case_when(
    transect %in% ab ~ "Above",
    transect %in% bel ~ "Below",
    transect %in% inside ~ "Inside",
    transect %in% inside2 ~ "Inside2",
    transect %in% wash ~ "Wash",
  )) %>% 
  mutate(shore = factor(shore, levels = c("Upper","Mid","Low","Surf")),
         transect = factor(transect,levels = c(
           "T1N", "T1","T1S",
           "T4","T7", "T8", "T11", "T12",
           "T13",
           "T15", "T21", "T22", "T23","T24", "T25","T26",
           "WA1","WA2","WA3","WA4","WA5","WA6"
         )),
         zone1 = factor(zone1, levels = c("Above","Inside","Inside2",
                                          "Below","Wash"))) %>% 
  relocate(zone1, .after = "method")

## Stat values ####
### texture ####
grad_out$sedim$texture <- grad_out$sedim$texture %>% 
  separate(
    col = samples,
    into = c("year", "transect", "shore", "method"),
    sep = "\\.",     # period is a regex special char, so escape it
    remove = FALSE,  # keep original 'samples' if you want
    fill = "right",  # if some rows are missing parts, fill with NA on the right
    extra = "merge"  # if there are extra pieces, merge them into the last column
  ) %>% 
  mutate(year = as.integer(year)) %>% 
  # assign zones
  mutate(zone1 = case_when(
    transect %in% ab ~ "Above",
    transect %in% bel ~ "Below",
    transect %in% inside ~ "Inside",
    transect %in% inside2 ~ "Inside2",
    transect %in% wash ~ "Wash",
  )) %>% 
  mutate(shore = factor(shore, levels = c("Upper","Mid","Low","Surf")),
         transect = factor(transect,levels = c(
           "T1N", "T1","T1S",
           "T4","T7", "T8", "T11", "T12",
           "T13",
           "T15", "T21", "T22", "T23","T24", "T25","T26",
           "WA1","WA2","WA3","WA4","WA5","WA6"
         )),
         zone1 = factor(zone1, levels = c("Above","Inside","Inside2",
                                          "Below","Wash"))) %>% 
  mutate(texture = factor(texture, levels = c(
    "Mud",
    "Slightly Gravelly Mud",
    "Sandy Mud",
    "Slightly Gravelly Sandy Mud",
    "Gravelly Mud",
    "Muddy Sand",
    "Slightly Gravelly Muddy Sand",
    "Gravelly Muddy Sand",
    "Sand",
    "Slightly Gravelly Sand",
    "Gravelly Sand",
    "Muddy Gravel",
    "Muddy Sandy Gravel",
    "Sandy Gravel",
    "Gravel"
    ))) %>%
  relocate(zone1, .after = "method")

### description ####
grad_out$sedim$descript <- grad_out$sedim$descript %>% 
  separate(
    col = samples,
    into = c("year", "transect", "shore", "method"),
    sep = "\\.",     # period is a regex special char, so escape it
    remove = FALSE,  # keep original 'samples' if you want
    fill = "right",  # if some rows are missing parts, fill with NA on the right
    extra = "merge"  # if there are extra pieces, merge them into the last column
  ) %>% 
  mutate(year = as.integer(year)) %>% 
  # assign zones
  mutate(zone1 = case_when(
    transect %in% ab ~ "Above",
    transect %in% bel ~ "Below",
    transect %in% inside ~ "Inside",
    transect %in% inside2 ~ "Inside2",
    transect %in% wash ~ "Wash",
  )) %>% 
  mutate(shore = factor(shore, levels = c("Upper","Mid","Low","Surf")),
         transect = factor(transect,levels = c(
           "T1N", "T1","T1S",
           "T4","T7", "T8", "T11", "T12",
           "T13",
           "T15", "T21", "T22", "T23","T24", "T25","T26",
           "WA1","WA2","WA3","WA4","WA5","WA6"
         )),
         zone1 = factor(zone1, levels = c("Above","Inside","Inside2",
                                          "Below","Wash"))) %>% 
  relocate(zone1, .after = "method")

## save ####
saveRDS(grad_out, file="output/models/grad_all_raw.Rdat")
toc(log=TRUE)

unlist(tictoc::tic.log())

# tidy up ####