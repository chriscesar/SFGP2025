# 200_an.inf_Wash.R ####
### Display of Wash assemblages

# Set up ####
## load packages ####
ld_pkgs <- c("tidyverse","ggplot2","vegan","ggtext",
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

# format Wash data ####
# sum duplicated taxa ####
df %>% #names()
  filter(!is.na(count)) %>%
  filter(.,zone1 == "Wash") %>% #names()
  ## remove zero values
  filter(.,count !=0) %>% #names()
  ## convert -9999 to 1
  mutate(count = ifelse(count < 0, #replace <0 values with 1
                        1,
                        count)) %>%
  dplyr::select(.,
                -taxonReported,
                -zone2.1,-zone2.2,
                -yr.trn,-yr.trn.sh,-yr.trn.sh.meth,-yr.trn.sh.meth.rep,
                -Kingdom,
                -Phylum,
                -Class,
                -Order,
                -Family,
                -Genus,
                -Species,
                -Comment,
                -MTG,
                -MxTx,
                -Flag,
  ) %>% 
  group_by(across(c(!count)
  )) %>% 
  summarise(.,count=sum(count), .groups = "drop") %>% ungroup() %>% 
  ## remove superfluous cols
  ### widen and fill gaps with 0:
  pivot_wider(.,
              names_from=taxonUSE,
              values_from=count,
              values_fill=list(count = 0)
  ) -> dfwash_all


# plots ####
dfwash_all %>% 
  ## lengthen
  pivot_longer(.,cols = -c("year","transect","shore","rep","zone1","mesh",
                           "core.area_m2"),
               names_to = "Taxon",values_to="Abundance") %>% 
  ## calculate mean across reps
  dplyr::select(.,-rep) %>% 
  group_by(across(!"Abundance")) %>% 
  summarise(., Abundance = sum(Abundance),.groups = "drop") %>% ungroup() %>% 
  filter(.,mesh=="1.0mm") -> dfwash_mean

Taxon <- unique(dfwash_mean$Taxon)
xxx <- vegan::make.cepnames(Taxon)

xx <- data.frame("Taxon" = Taxon, "Taxon_lab" = xxx)

dfwash_mean <- left_join(dfwash_mean, xx, by = "Taxon")
dfwash_mean$shore <- factor(dfwash_mean$shore, levels=c("Mid","Low"))

ntax <- length(unique(dfwash_mean$Taxon))
minyr <- min(dfwash_mean$year)

dfwash_mean <- dfwash_mean %>% filter(Abundance!=0)

## add 'empty' values for plotting ####
dfwash_mean <- bind_rows(dfwash_mean,tibble(year=c(2020,2022,2023),
                                            shore=c("Low","Low","Mid"),
                                            Taxon_lab=c("Ceraedul","Ceraedul","Ceraedul"),
                                            Taxon = rep("Cerastoderma edule",3),
                                            Abundance = rep(0,3)
                                            )
                         )

dfwash_mean$shore <- factor(dfwash_mean$shore, levels=c("Mid","Low"))

df0 %>% select(taxonUSE,MTG,MxTx) %>% distinct()->df_names

df_names %>% 
  mutate(MxTx2 = case_when(
    MxTx == "Annelid" ~ "Annelid",
    MxTx == "Mollusc" ~ "Mollusc",
    MxTx == "Arthropod" ~ "Arthropod",
    MxTx %in% c("Nematoda","Nemertea", "Bryozoan","Hydrozoan","Flatworm",
                "Ascidian","Animalia","Porifera","Anemone") ~ "Other",
    TRUE ~ NA
  )) -> df_names

dfwash_mean %>% 
  left_join(.,df_names,by=c("Taxon" = "taxonUSE")) -> dfwash_mean_tx

dfwash_mean_tx$MxTx <- factor(dfwash_mean_tx$MxTx)

ord <- order(dfwash_mean_tx$MxTx2, dfwash_mean_tx$Abundance, dfwash_mean_tx$Taxon_lab)
levs <- unique(dfwash_mean_tx$Taxon_lab[ord])
dfwash_mean_tx$Taxon_lab <- factor(dfwash_mean_tx$Taxon_lab, levels = levs)

############################
############################
# TO DO ####
## arrange taxon names by MxTx for plotting
############################
############################

(dfwash_mean_tx %>% 
  filter(MxTx!="Algae") %>% 
  #filter(Abundance != 0) %>% 
  # filter(!is.na(shore)) %>% 
  ggplot(.,
         aes(
           y = as.factor(year),
           x = factor(Taxon_lab,levels=)
         )) +
  geom_vline(xintercept = seq(from = 1, to = ntax,by=1), colour="grey",lty=2)+
  geom_point(aes(size = log(Abundance+1),
                 # colour = log(Abundance+1)
                 colour=MxTx2,
                 ),
             # show.legend = FALSE
             )+
  facet_wrap(.~shore)+
  scale_x_discrete(limits=rev)+
  guides(size = "none",
         colour = guide_legend(override.aes = list(size=8)))+
  labs(title = paste0(
    "Taxa recorded in intertidal transects in The Wash since ",
    minyr),
    subtitle = "Point colours indicate the major taxonomic group. Point sizes indicate the relative abundances",
    caption = "Taxon names abbreviated using the vegan::makecepnames() function in R
    Note that prior to 2024, only a single tranect was visited. Three transects were visited in 2024, six transects were visited in 2025")+
  theme(
    axis.text = element_text(face=2),
    plot.title = element_text(face=2),
    strip.text = element_text(face=2),
    axis.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(face=2),
    
  )+
  coord_flip()->pl)
ggsave(plot = pl, filename = "figs/WashTaxa2.png",
       width = 16,height = 9,units = "in");rm(pl)

# tidy up ####
rm(list=ls(pattern = "^df"))
rm(list=ls(pattern = "^cb"))
rm(xxx,xx,cur.yr,fol,gisfol,minyr,ntax,perm,ppi,projfol,source_file,Taxon,sum_zero)

detach("package:tidyverse", unload=TRUE)
detach("package:vegan", unload=TRUE)
detach("package:ggplot2", unload=TRUE)
detach("package:tictoc", unload=TRUE)
