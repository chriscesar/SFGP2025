# mapfun.R ####
## create site map for intertidal points

### load packages
ld_pkgs <- c("tidyverse","ggplot2","sf","rgdal","maps",
             "ggpubr", "ggspatial","ggrepel")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
source("R/00_datfol.R")
source("R/00_meta_setMeta.R")

ppi <- 300
cbPalette <- c( ### colourblind-friendly chart colour palette
  # comments reflect Red/Green/Blue equivalents
  "#0072B2", #000, 114, 178
  "#e79f00", #231, 159, 0
  "#009E73", #000, 158, 115
  "#9ad0f3", #154, 208, 243
  "#CC79A7", #204, 121, 167
  "#DEAAC6", #222, 170, 198 - The Wash
  "#000000", #0, 0, 0
  "#D55E00", #213, 94, 0
  "#F0E442"  #240, 228, 66
)

### GIS folder
# load data and convert to lat/long ####
df0 <- as_tibble(openxlsx::read.xlsx(paste0(projfol,"GIS/2025/IntertidalStations2025.xlsx"),
                                     sheet = 1))
df0$Zone <- factor(df0$Zone,
                   levels=c("Above","Inside","Inside2","Below","Wash"))

# GGPLOT STATIC MAPS ######
#### attempt at ggplot2 version
base_0 <- st_read(paste0(gisfol,"shapes/Coast_polygon.shp"))
basedf <- fortify(base_0)
towns_pt_0 <- st_read(paste0(gisfol,"shapes/","Urban_pt_TwnCty.shp"))
towns_pt_df <- as_tibble(fortify(towns_pt_0))
towns_pt_df$East <- st_coordinates(towns_pt_0)[,"X"]
towns_pt_df$North <- st_coordinates(towns_pt_0)[,"Y"]
towns_area <- st_read(paste0(gisfol,"shapes/","Urban_areas_250k.shp"))

## load designated shapes ####
# sac #
sh_sac0 <- st_read(paste0(gisfol,"shapes/","sac_10k.shp"))
sh_sac <- fortify(sh_sac0)

# spa #
sh_spa0 <- st_read(paste0(gisfol,"shapes/","Spa_10k.shp"))
sh_spa <- fortify(sh_spa0)

# Ramsar #
sh_rams0 <- st_read(paste0(gisfol,"shapes/","Ramsar_10k.shp"))
sh_rams <- fortify(sh_rams0)

png(file = "figs/siteMap.png",
    width=8*ppi, height=12*ppi, res=ppi)
ggplot()+
  geom_sf(data = base_0, fill = "#B4D79D") +
  geom_sf(data = towns_area[towns_area$DESCRIPTIO == "Large Urban Area polygon",],
          #fill="darkgrey")+
          fill="#CCCCCC")+
  geom_point(data = df0,
             aes(x = Eastings,
                 y = Northings,
                 fill = Zone#,
                 #shape = Tag
             ),
             colour = 1,
             pch = 21,
             size = 4, 
             inherit.aes = FALSE,
             show.legend = FALSE)+
  geom_text_repel(data = df0,
                  segment.colour="grey",
                  nudge_x = 0.1,
                  # point.padding = 0.5,
                  # point.padding = 4.5,
                  point.padding = .5,
                  force_pull = 8,
                  aes(x = Eastings,
                      y = Northings,
                      label = Transect), 
                  fontface = "bold",
                  force = 0.5,
                  inherit.aes = FALSE,
                  seed = pi,
                  bg.colour="white",
                  direction = "x") +
  geom_text_repel(data = towns_pt_df[towns_pt_df$NAME == "SKEGNESS" |
                                       towns_pt_df$NAME == "MABLETHORPE",],
                  aes(x = East, y = North, label = NAME), 
                  inherit.aes = FALSE,
                  hjust = 1,
                  vjust = 1,
                  nudge_x = -1000,
                  segment.colour=NA,
                  bg.colour="white",
                  fontface = "bold.italic")+
  coord_sf(xlim = c(534750, 563000),
           ylim = c(344500, 389070))+
  scale_fill_manual(values = cbPalette)+
  # scale_shape_manual(values=c(3,21))+
  labs(title = "Location of intertidal transects surveyed as part of the Saltfleet to Gibraltar\nPoint Strategy 2025")+
  ggthemes::theme_few()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  ggspatial::annotation_scale(location="br", width_hint=0.5)+
  ggspatial::annotation_north_arrow(which_north = "grid", location = "bl",
                                    pad_x = unit(0.0, "cm"),
                                    pad_y = unit(0.5, "cm"),
                                    height = unit(2, "cm"),
                                    width = unit(2, "cm"),
                                    style = north_arrow_fancy_orienteering
  )
dev.off()

# to do: ####
## add WA proposed samples
## produce version with WFD WBs & protected sites ###
### not looking great... ##
ggplot()+
  geom_sf(data = base_0, fill = "#B4D79D") +
  geom_sf(data = towns_area[towns_area$DESCRIPTIO == "Large Urban Area polygon",],
          #fill="darkgrey")+
          fill="#CCCCCC")+
  geom_sf(data=sh_sac, fill = "lightblue",alpha=0.5)+
  ggplot2::geom_sf_text(data=sh_sac,aes(label=NAME))+
  geom_text_repel(data = towns_pt_df[towns_pt_df$NAME == "SKEGNESS" |
                                       towns_pt_df$NAME == "MABLETHORPE",],
                  aes(x = East, y = North, label = NAME), 
                  inherit.aes = FALSE,
                  hjust = 1,
                  vjust = 1,
                  nudge_x = -1000,
                  segment.colour=NA,
                  bg.colour="white",
                  fontface = "bold.italic")+
  coord_sf(xlim = c(534750, 563000),
           ylim = c(344500, 389070))+
  scale_fill_manual(values = cbPalette)+
  # scale_shape_manual(values=c(3,21))+
  # labs(title = "Location of intertidal transects surveyed as part of the Saltfleet to Gibraltar\nPoint Strategy 2024")+
  ggthemes::theme_few()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  ggspatial::annotation_scale(location="br", width_hint=0.5)+
  ggspatial::annotation_north_arrow(which_north = "grid", location = "bl",
                                    pad_x = unit(0.0, "cm"),
                                    pad_y = unit(0.5, "cm"),
                                    height = unit(2, "cm"),
                                    width = unit(2, "cm"),
                                    style = north_arrow_fancy_orienteering
  )

# Offshore transects ####
## read survey log ####
df_log <- readxl::read_xlsx(paste0(fol,"vesselLogs.xlsx"),
                            sheet = "ts_trawls_for_Arc") %>% 
  dplyr::filter(Year == cur.yr) %>% 
  dplyr::mutate(Survey = factor(Survey,levels = c("Sep","Oct")))

png(file = "figs/siteMap_epifaunal.png",
    width=8*ppi, height=12*ppi, res=ppi)
ggplot()+
  geom_sf(data = base_0, fill = "#B4D79D") +
  geom_sf(data = towns_area[towns_area$DESCRIPTIO == "Large Urban Area polygon",],
          #fill="darkgrey")+
          fill="#CCCCCC")+
  geom_segment(
    data = df_log,
    aes(x = ShootEast, y = ShootNorth,
        xend = HaulEast, yend = HaulNorth,
        colour = Survey
        ),
    size = 1,
    lineend = "round"
  ) +
  geom_text_repel(data = df_log %>% dplyr::filter(Site=="C" & Survey=="Sep"),
                  aes(
                    x = ShootEast,
                    y= ShootNorth,
                    label = Transect
                  ),
                  segment.colour = NA,
                  fontface = "bold",
                  force = 0.5,
                  inherit.aes = FALSE,
                  point.padding = .5,
                  force_pull = 8,
                  seed = pi,
                  size=5,
                  bg.colour="white",
                  direction = "x",
                  nudge_x = 1500,
                  nudge_y = 500
                  )+
  geom_text_repel(data = towns_pt_df[towns_pt_df$NAME == "SKEGNESS" |
                                       towns_pt_df$NAME == "MABLETHORPE",],
                  aes(x = East, y = North, label = NAME), 
                  inherit.aes = FALSE,
                  hjust = 1,
                  vjust = 1,
                  nudge_x = -1000,
                  segment.colour=NA,
                  bg.colour="white",
                  fontface = "bold.italic")+
  coord_sf(xlim = c(540750, 569000),
           ylim = c(347000, 391570))+
  scale_colour_manual(values = c(2,1))+
  # scale_fill_manual(values = cbPalette)+
  # scale_shape_manual(values=c(3,21))+
  labs(title = "Location of epifaunal trawls surveyed as part of the Saltfleet to Gibraltar\nPoint Strategy 2025")+
  ggthemes::theme_few()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face=2,size = 12),
        )+
  ggspatial::annotation_scale(location="br", width_hint=0.5)+
  ggspatial::annotation_north_arrow(which_north = "grid", location = "bl",
                                    pad_x = unit(0.0, "cm"),
                                    pad_y = unit(0.5, "cm"),
                                    height = unit(2, "cm"),
                                    width = unit(2, "cm"),
                                    style = north_arrow_fancy_orienteering
  )
dev.off()

# Tidy up ####
rm(list = ls(pattern = "^base"))
rm(list = ls(pattern = "^town"))
rm(list = ls(pattern = "^df"))
rm(list = ls(pattern = "^sh"))
rm(list = ls(pattern = "^cb"))
rm(fol, gisfol,ppi,transects_df,cur.yr,perm,projfol,sum_zero)

detach(package:ggrepel, unload=TRUE)
detach(package:ggspatial, unload=TRUE)
detach(package:ggpubr, unload=TRUE)
detach(package:maps, unload=TRUE)
detach(package:rgdal, unload=TRUE)
detach(package:sf, unload=TRUE)
detach(package:ggplot2, unload=TRUE)
detach(package:tidyverse, unload=TRUE)
