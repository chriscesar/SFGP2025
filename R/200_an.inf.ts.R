# 200_an.inf.ts.R ####
### Multivariate analysis of intertidal invertebrate time series data

# Set up ####
## load packages ####
ld_pkgs <- c("tidyverse","ggplot2","vegan","ggdendro",
             "mvabund","ggtext","ggh4x",
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

# load taxon info
dftx <- readxl::read_xlsx("data/inf_ts_longRAW_USE.xlsx",
                          sheet = "taxa")

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

# Taxa of interest ####
# Polychaetes (e.g., Scolelepis spp.) — often most negatively affected and slowest to recover. Wooldridge et al. 2016
# Haustoriid amphipods and donacid clams—can rebound faster; useful for detecting early recovery. Wooldridge et al. 2016
# Peracarids / isopods (Eurydice) and talitrids—sediment/organic matter sensitive; respond to wrack and grain size. [Staudt et al. 2021]



## plot Scolelepis 
dfl %>% 
  filter(value!=0) %>% 
  filter(grepl("^Bathyporeia",name)) %>% 
  mutate(value_m3 = value/core.area_m2) %>% 
  mutate(zone1 = factor(zone1,
                        levels = c("Above","Inside",
                                   "Inside2","Below",
                                   "Wash"))) %>% 
  mutate(shore = factor(shore,
                        levels=c("Mid","Low"))) %>% 
  ggplot(.,
         aes(
           x=year,
           y=log(value_m3+1),
         )
         )+
  geom_jitter()+
  geom_smooth(method="gam")+
  facet_grid(shore~zone1)

dfl %>%
  # remove AFAUNAL
  filter(name != "AFAUNAL") %>% 
  #remove presence only
  filter(value >=0) %>% 
  ## calculate abundance per m2
  mutate(value = value/core.area_m2) %>% 
  select(-core.area_m2, -tmp) %>% 
  select(-transect,-mesh) %>% 
  mutate(zone1 = factor(zone1,
                        levels = c("Above","Inside",
                                   "Inside2","Below",
                                   "Wash"))) %>% 
  mutate(shore = factor(shore,
                        levels=c("Mid","Low"))) %>% 
  group_by(across(-value)) %>% 
  summarise(value = mean(value),.groups = "drop") %>% ungroup() %>% 
  ggplot(.,
         aes(
           x = year,
           y = log10(value+1),
         )
         )+
  geom_point(aes(
    group = name
  ),
  alpha = 0.2) +
  facet_grid(shore ~ zone1)+
  geom_smooth(method="gam")
  
# append taxon MTG for plotting

dftx <- dftx %>% select(ScientificName_accepted:MTG) %>% distinct()

dfl %>% 
  left_join(.,dftx,
            by = c("name" = "ScientificName_accepted"))->dfltx

dfltx %>% 
  select(year,transect,shore,zone1,core.area_m2,MTG,value) %>% 
  mutate(value_m3 = value/core.area_m2) %>% 
  select(-core.area_m2) %>% 
  filter(value_m3 >=0) %>% 
  filter(MTG != "AFAUNAL") %>% 
  select(-value) %>% 
  group_by(across(-value_m3)) %>% 
  summarise(value_m3 = sum(value_m3),.groups = "drop")->dfltmp

dfltmp %>% 
  group_by(MTG) %>%
  summarise(any_nonzero = any(value_m3 != 0), .groups = "drop") %>%
  filter(any_nonzero) %>%
  pull(MTG) -> kp

dfltmp %>% 
  filter(MTG %in% kp) %>% 
  select(-transect) %>% 
  group_by(across(-value_m3)) %>% 
  summarise(value_m3 = mean(value_m3),.groups = "drop") %>% 
  mutate(zone1 = factor(zone1,
                        levels = c("Above","Inside",
                                   "Inside2","Below",
                                   "Wash"))) %>% 
  mutate(shore = factor(shore,
                        levels=c("Mid","Low"))) %>% 
  mutate(MTG2 = sub("_.*", "", MTG)) ->dfltmp_MTG2

dfltmp_MTG2 %>% 
  ggplot(.,
         aes(
           x = year,
           y = log10(value_m3+1),
         ))+
  geom_line(aes(
    group = MTG2
  ),
  alpha = 0.2
  )+
  geom_jitter(aes(
    group = MTG2,
    col=MTG2
    ),
    alpha = 0.2)+
  geom_smooth(method = "gam",
              aes(
                group = MTG2,
                col=MTG2,
              ),
              se=FALSE
              )+
  facet_grid(shore ~ zone1)

####
# Extract polychaetes from dfl data and plot with lines for individual taxa



dfltx %>% 
  mutate(zone1 = factor(zone1,
                        levels = c("Above","Inside",
                                   "Inside2","Below",
                                   "Wash"))) %>% 
  mutate(shore = factor(shore,levels=c("Mid","Low"))) %>% 
  filter(MTG == "Mollusc_Bivalve") %>% 
  mutate(value_m2 = value/core.area_m2) %>% 
  # ## sums by Family
  # select(year,zone1,transect,shore,Family,value_m2) %>%
  ## sums by Class
  select(year,zone1,transect,shore,Class,value_m2) %>% 
  group_by(across(-value_m2)) %>% 
  summarise(value_m2 = sum(value_m2),.groups = "drop") %>% ungroup() %>% 
  # mean by shore_zone
  select(-transect) %>% 
  group_by(across(-value_m2)) %>% 
  summarise(value_m2 = sum(value_m2),.groups = "drop") %>% ungroup() %>%
  ggplot(.,aes(
    x = year,
    y = log10(value_m2+1),
    )
  )+
  geom_line(aes(
    group = Class,
    colour = Class,
    ),
    show.legend = FALSE)+
  facet_grid(shore ~ zone1)+
  geom_smooth(method="gam",
              aes(group = Class,
                  colour = Class,),
              se=FALSE, show.legend = FALSE
              )+
  geom_jitter(aes(group = Class,
             colour = Class,),show.legend = FALSE)



dfltx %>% 
  mutate(zone1 = factor(zone1,
                        levels = c("Above","Inside",
                                   "Inside2","Below",
                                   "Wash"))) %>% 
  mutate(shore = factor(shore,levels=c("Mid","Low"))) %>% 
  filter(MTG != "AFAUNAL") %>% 
  mutate(MTG_label = case_when(
    MTG == "Arthropod_Amphipod"~ "Amphipods",
    MTG == "Mollusc_Bivalve"~ "Bivalve Molluscs",
    MTG == "Arthropod_Pycnogonid"~ "Pycnogonids",
    MTG == "Bryozoan"~ "Bryozoans",
    MTG == "Annelid_Polychaete"~ "Polychaetes",
    MTG == "Ascidian"~ "Ascidians",
    MTG == "Arthropod_Barnacle"~ "Barnacles",
    MTG == "Hydrozoan"~ "Hydrozoan",
    MTG == "Arthropod_Shrimp_Crab"~ "Shrimps & Crabs",
    MTG == "Nemertea"~ "Nemerteans",
    MTG == "Arthropod_Hexapod"~ "Hexapods",
    MTG == "Arthropod_Copepod"~ "Copepods",
    MTG == "Mollusc_Gastropod"~ "Gastropod Molluscs",
    MTG == "Algae_brown"~ "Algae",
    MTG == "Anemone"~ "Anemones",
    MTG == "Annelid_Polychaete" ~ "Polychaetes",
    MTG == "Annelid_Oligochaete"~ "Oligochaetes",
    MTG == "Arthropod_IsopodMysid"~ "Isopods & Mysids",
    MTG == "Flatworm"~ "Flatworms",
    MTG == "Mollusc_Chiton"~ "Chitons",
    MTG == "Arthropod_Ostracod"~ "Ostracods",
    MTG == "Nematoda"~ "Nematodes",
    MTG == "Porifera"~ "Porifera",
    TRUE ~ NA
    )
    )->df_temp

tic("Generate plots")
for(i in unique(df_temp$MTG)){
  print(i)
  data_tmp <- df_temp %>% filter(MTG == i)
  
  zones <- unique(data_tmp$zone1)
  bg_cols <- cbPaletteFill[c(1:4,7)]
  bg_cols <- setNames(bg_cols[seq_along(zones)], zones)
  
  strip_elems <- lapply(zones, function(z)
    element_rect(fill = bg_cols[[z]], color = "black", linewidth = 1)
  )

  data_tmp %>% 
    mutate(value_m2 = value/core.area_m2) %>% 
    # ## sums by Family
    # select(year,zone1,transect,shore,Family,value_m2) %>%
    ## sums by Class
    select(year,zone1,transect,shore,MTG,MTG_label,value_m2) %>% 
    group_by(across(-value_m2)) %>% 
    summarise(value_m2 = sum(value_m2),.groups = "drop") %>% ungroup() %>% 
    # mean by shore_zone
    select(-transect) %>% 
    group_by(across(-value_m2)) %>% 
    summarise(value_m2 = sum(value_m2),.groups = "drop") %>% ungroup() %>%
    ggplot(.,aes(
      x = year,
      y = log10(value_m2+1),
    )
    )+
    facet_grid2(shore ~ zone1, #ncol = 2,
                strip = strip_themed(background_x = strip_elems)) +
    # geom_smooth(method="gam",
    #             aes(group = MTG,
    #                 colour = MTG,),
    #             se=FALSE, show.legend = FALSE
    # )+
    labs(
      title = unique(data_tmp$MTG_label),
      y = expression(bold(Log[10(n+1)]~"abundance (individuals"~m^-2*")"))
      )+
    geom_jitter(aes(group = MTG,
                    colour = MTG,),show.legend = FALSE, col=1)+
    theme(legend.title=element_blank(),legend.direction="horizontal",
          legend.position = "inside",
          legend.position.inside = c(.625, 1.4),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face=2, size=12),
          strip.text = element_text(face=2,size = 16,),
          axis.title = element_text(face=2),
          axis.text = element_text(face=2),
          axis.text.x = element_text(angle=270, vjust=0.5, size = 14,face=2),
          # plot.subtitle = element_markdown(face=2),
          plot.title = element_markdown(face=2,size = 18),
    ) -> pl
  
  png(file = paste0("figs/inf.ts_",vegan::make.cepnames(i),".png"),
      width = 16 * ppi, height = 8 * ppi, res = ppi)
  print(pl)
  dev.off();
  rm(pl,data_tmp)
  print(paste0("Saved image for ",i))
  }
toc(log=TRUE)
