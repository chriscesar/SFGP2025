# 200_an_inf_TRT.R ####
# generate plots for traits ####
# Set up ####
### load packages ####
ld_pkgs <- c("tidyverse","tictoc","ggh4x")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)
### load metadata ####
source("R/00_meta_setMeta.R")

df0 <- as_tibble(readxl::read_xlsx("output/inf.int.traitsUSE.xlsx", sheet="inf.int.traits.ts_OUT")) %>% 
  mutate(shore = factor(shore, levels = c("Mid","Low")),
         zone1 = factor(zone1, levels = c("Above", "Inside", "Inside2","Below","Wash"))) %>% 
  select(.,-zsort,-sshort)

# trtloop <- data.frame(var = c("b_","ed_","f_","ld_","lh_","m_","mob_","sp_","sr_")) %>% 
#   mutate(trait_lable = case_when(
#     var == "b_" ~"Bioturbation",
#     var == "ed_" ~"Egg Development",
#     var == "f_" ~"Feeding Mode",
#     var == "l_" ~"Lifespan (years)",
#     var == "ld_" ~"Larval development",
#     var == "lh_" ~"Living Habit",
#     var == "m_" ~"Morphology",
#     var == "mob_" ~"Mobility",
#     var == "sp_" ~"Sediment Position",
#     var == "sr_" ~"Maximum Size (mm)",
#     TRUE ~ NA
#   ))

tic("Format data for plotting")
df0 %>% 
  ## remove label columns
  select(.,-c(BIOTURB,EGG_DEVT,FEEDING_MODE,LIFESPAN,LARV_DEVT,LIVING_HAB,
  MORPH,MOBIL,SED_POSN,MAX_SIZE,SPECIES_TYPOLOGY,SPECIES_TYPOLOGY_2)) %>% 
  pivot_longer(.,-c(year,shore,zone1,mesh),
               names_to = "trt"
                ) %>% 
  mutate(Trait = case_when(
    startsWith(trt, "b_") ~ "Bioturbation",
    startsWith(trt,"ed_") ~"Egg Development",
    startsWith(trt,"f_") ~"Feeding Mode",
    startsWith(trt, "l_") ~"Lifespan (years)",
    startsWith(trt, "ld_") ~"Larval development",
    startsWith(trt, "lh_") ~"Living Habit",
    startsWith(trt, "m_") ~"Morphology",
    startsWith(trt, "mob_") ~"Mobility",
    startsWith(trt, "sp_") ~"Sediment Position",
    startsWith(trt, "sr_") ~"Maximum Size (mm)",
    TRUE ~ NA_character_
  )) %>% 
  mutate(modality = sub("^[^_]*_", "", trt)) %>% 
  mutate(trt_lb = vegan::make.cepnames(Trait)
         ) %>% 
  mutate(modality = factor(modality, levels = c(
    #Bioturbation
    "None","Surface_deposition","Diffusive_mixing",	"Downward_conveyer","Upward_conveyor",
    #Egg Devt
    "Asexual",	"Sexual_benthic",	"Sexual_pelagic","Sexual_brooded",
    #Feeding
    "Suspension","Surface_deposit","Subsurface_deposit","Predator","Scavenger", "Parasite",
    #Livespan
    "Less_than_1","1_to_3","3_to_10","More_than_10",
    #LarvalDispersal
    "Pelagic_planktotrophic",	"Pelagic_lecithotrophic","Benthic_direct",
    #LivingHabit
    "Free_living","Burrow_dwelling","Tube_dwelling","Attached_to_substratum","Crevice_hole_under_stones","Epi_endo_biotic",
    #Morphology
    "Crustose","Cushion","Exoskeleton","Soft","Stalked","Tunic",
    #Mobility
    "Sessile","Burrower","Crawl_creep_climb","Swim",
    #SedimentPosn
    "Surface","Shallow_infauna_0_to_5cm","Mid_depth_infauna_5_to_10cm","Deep_infauna_more_than_10cm",
    #maxSize
    "Less_than_10","11_to_20","21_to_100","101_to_200","201_to_500","More_than_500"
    ))
  ) -> df_trt
toc(log=TRUE)

tic("Generate plots")
# Build a strip background list that matches your facet order
zones <- unique(df_trt$zone1)
bg_cols <- cbPaletteFill[c(1:4,7)]
bg_cols <- setNames(bg_cols[seq_along(zones)], zones)

strip_elems <- lapply(zones, function(z)
  element_rect(fill = bg_cols[[z]], color = "black", linewidth = 1)
  )

for(i in 1:length(unique(df_trt$Trait))){
trt_tmp <- unique(df_trt$Trait)[i]
dftmp <- droplevels(df_trt %>% filter(Trait == trt_tmp))
print(unique(dftmp$Trait))
ggplot(dftmp,
       aes(
         x = year,
         y = value,
         colour = modality,
         shape = modality,
       )) +
  geom_line(size=1.5,alpha=0.7,
            aes(linetype = modality)
            )+
  # geom_point(size=3)+
  # facet_grid(shore ~ zone1)+
  facet_grid2(shore ~ zone1, #ncol = 2,
              strip = strip_themed(background_x = strip_elems)) +
  labs(title = trt_tmp)+
  scale_colour_manual(values = cbPalette3)+
  scale_y_continuous(breaks = c(0,1))+
  scale_linetype_manual(values = rep(c(1,3,2),4))+
  theme(
    plot.title = element_text(face = 2,size=16),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(face= 2),
    axis.text.x = element_text(face= 2,size=14),
    strip.text = element_text(face=2, size = 14),
    legend.title = element_blank(),
    legend.text = element_text(face=2, size = 14),
    strip.background.y = element_rect(color = "black",fill = "grey95",
                                      linewidth = 1)
  ) -> pl
  png(
    file = paste0("figs/inf.Trt.ts.",unique(dftmp$trt_lb),".png"),
    width = 15 * ppi,
    height = 10 * ppi,
    res = ppi
    )
    print(pl)
    dev.off()
    rm(trt_tmp,dftmp)
}
toc(log=TRUE)
