# 200_an.inf_v2.R ####
### Multivariate analysis of intertidal invertebrate data

# Set up ####
## load packages ####
ld_pkgs <- c("tidyverse","ggplot2","vegan",
             "mvabund",
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

# Extract current year's data ####
tic("Extract current year's data")
dfw0 %>% 
  dplyr::select(-Flag) %>% 
  dplyr::filter(.,year == cur.yr) %>% 
  dplyr::filter(.,zone1 != "Wash") %>% 
  ## remove 'empty' columns
  dplyr::select(where(~ !is.numeric(.)| !all(.==0))) %>% 
  ## calculate means across replicates
  pivot_longer(.,cols = -c(year, transect, shore,,rep,zone1, mesh, core.area_m2),
               names_to = "Taxon", values_to = "Abundance") %>%
  ## replace -9999 with 1
  mutate(Abundance = if_else(Abundance < 0, 1, Abundance)) %>% 
  dplyr::select(.,-rep) %>% 
  group_by(across(c(!Abundance)
  )) %>% 
  summarise(.,Abundance=mean(Abundance), .groups = "drop") %>% ungroup() %>% 
  pivot_wider(.,names_from = Taxon, values_from = Abundance) %>% 
  dplyr::select(.,-AFAUNAL) -> df.cur

## isolate taxon data ####
df.cur %>% dplyr::select(.,
                         -c(year,transect,
                            shore,zone1,mesh,
                            core.area_m2)) -> df.tx
toc(log=TRUE)

# MVABUND MODELS ####
tic("MVABUND MODELS")
## version 1 ####
df.cur$zone_model <- factor(df.cur$zone1, levels = c("Inside","Above",
                                                     "Inside2","Below"))
# m1 <- mvabund::manyglm(mvabund::mvabund(df.tx)~df.cur$zone_model,
#                        family = "negative.binomial")
# m1.summary <- summary(m1, nBoot = 9999)
# anova_m1 <- mvabund::anova.manyglm(m1,p.uni = "adjusted")
# 
# saveRDS(m1, file = "output/models/mvabund_mod1_use.Rdat")
# saveRDS(m1.summary, file = "output/models/mvabund_mod1_summary_use.Rdat")
# saveRDS(anova_m1, file = "output/models/mvabund_mod1_anova_use.Rdat")

m1 <- readRDS("output/models/mvabund_mod1_use.Rdat")
m1.summary <- readRDS("output/models/mvabund_mod1_summary_use.Rdat")
anova_m1 <- readRDS("output/models/mvabund_mod1_anova_use.Rdat")

## version 2 ####
m2 <- mvabund::manyglm(mvabund::mvabund(df.tx)~df.cur$zone_model*df.cur$shore,
                       family = "negative.binomial")
# m2.summary <- summary(m2, nBoot = 9999)
# anova_m2 <- mvabund::anova.manyglm(m2,p.uni = "adjusted")
tic("anova_m2_shrink");anova_m2_shrink <- mvabund::anova.manyglm(m2,
                                                                p.uni = "adjusted",
                                                                cor.type = "shrink");toc(log=TRUE)
# #anova_pw_m2 <- anova.manyglm(m2, pairwise.comp = df.cur$zone_model*df.cur$shore)
# 
# saveRDS(m2, file = "output/models/mvabund_mod2_use.Rdat")
# saveRDS(m2.summary, file = "output/models/mvabund_mod2_summary_use.Rdat")
# saveRDS(anova_m2, file = "output/models/mvabund_mod2_anova_use.Rdat")
# saveRDS(anova_m2_shrink, file = "output/models/mvabund_mod2_anova_shrink_use.Rdat")

# pairwise comparisons
m2 <- readRDS("output/models/mvabund_mod2_use.Rdat")
m2.summary <- readRDS("output/models/mvabund_mod2_summary_use.Rdat")
anova_m2 <- readRDS("output/models/mvabund_mod2_anova_use.Rdat")
anova_m2_shrink <- readRDS("output/models/mvabund_mod2_anova_shrink_use.Rdat")

#Create an Interaction Factor
df.cur$zone_shore <- interaction(df.cur$zone_model, df.cur$shore)
anova_pw_m2 <- anova.manyglm(m2, pairwise.comp = df.cur$zone_shore)

toc(log=TRUE)

uni <- t(anova_m2$uni.test)
uni_shrink <- t(anova_m2_shrink$uni.test)
#interaction effects
sig_interaction <- rownames(uni)[uni[ , "df.cur$zone_model:df.cur$shore"] < 0.05]
sig_interaction
sig_interaction_shrink <- rownames(uni_shrink)[uni_shrink[ , "df.cur$zone_model:df.cur$shore"] < 0.05]
sig_interaction_shrink

#main effects
sig_zone  <- rownames(uni)[uni[ , "df.cur$zone_model"] < 0.05]
sig_shore <- rownames(uni)[uni[ , "df.cur$shore"] < 0.05]

sig_zone_shrink  <- rownames(uni_shrink)[uni_shrink[ , "df.cur$zone_model"] < 0.05]
sig_shore_shrink <- rownames(uni_shrink)[uni_shrink[ , "df.cur$shore"] < 0.05]


coef_mat <- coef(m2)
head(coef_mat)

species_effects <- data.frame(
  species = rownames(uni),
  p_zone = uni[, "df.cur$zone_model"],
  p_shore = uni[, "df.cur$shore"],
  p_interaction = uni[, "df.cur$zone_model:df.cur$shore"]#,
  # coef_zone = t(coef_mat[, "(Intercept)"]), # or your selected contrast
  # coef_shore = coef_mat[, "df.cur$shoreMid"]
  ) 

species_effects <- as_tibble(cbind(species_effects,t(coef_mat)))

species_effects %>% 
  filter(p_interaction<0.05) %>% View()

## version 3 ####
m3 <- mvabund::manyglm(mvabund::mvabund(df.tx)~df.cur$zone_shore,
                       family = "negative.binomial")
# m3.summary <- summary(m3, nBoot = 9999)
# anova_m3 <- mvabund::anova.manyglm(m3,p.uni = "adjusted")
tic("anova_m2_shrink");anova_m3_shrink <- mvabund::anova.manyglm(m3,
                                                                 p.uni = "adjusted",
                                                                 cor.type = "shrink");toc(log=TRUE)
# #anova_pw_m2 <- anova.manyglm(m2, pairwise.comp = df.cur$zone_model*df.cur$shore)
# 
saveRDS(m3, file = "output/models/mvabund_mod3_use.Rdat")
saveRDS(m3.summary, file = "output/models/mvabund_mod3_summary_use.Rdat")
saveRDS(anova_m3, file = "output/models/mvabund_mod3_anova_use.Rdat")
saveRDS(anova_m3_shrink, file = "output/models/mvabund_mod3_anova_shrink_use.Rdat")

# pairwise comparisons
m3 <- readRDS("output/models/mvabund_mod3_use.Rdat")
m3.summary <- readRDS("output/models/mvabund_mod3_summary_use.Rdat")
anova_m3 <- readRDS("output/models/mvabund_mod3_anova_use.Rdat")
anova_m3_shrink <- readRDS("output/models/mvabund_mod3_anova_shrink_use.Rdat")

## use 'shrink' version
# Bathyporeia.pelagica	Conopeum.reticulum	Copepoda	Gastrosaccus.spinifer	Haustorius.arenarius	Nematoda	Nemertea	Nephtys.cirrosa	Scolelepis..Scolelepis..squamata
# pull out taxa with significant values
as.data.frame(t(anova_m3_shrink$uni.p)) %>%
  # remove empty column
  select(-"(Intercept)") %>% 
  mutate(tax= row.names(.)) %>% 
  as_tibble()-> df_pvals

df_pvals %>%
  rename(zone_shore = `df.cur$zone_shore`) %>% 
  #retain only significant taxa
  filter(zone_shore <0.051) %>% 
  ## remove "."
  mutate(across(where(is.character), ~ gsub("\\.", " ", .))
         ) -> df_pvals_sig

sigtx <- as.vector(df_pvals_sig$tax)

### plot mod3 ####
png(file = "figs/inf.current.yr.sig_shrink_.png",
    width=12*ppi, height=8*ppi, res=ppi)
set.seed(42423);df.cur %>% 
  pivot_longer(cols = -c(year,transect,shore,zone1,mesh,core.area_m2,zone_model,
                         zoneplot,zone_shore), 
               names_to = "Taxon", 
               values_to = "Mean_abundance") %>% 
  filter(Taxon %in% sigtx) %>%
  mutate(zone1 = factor(zone1, levels = c("Above","Inside",
                                          "Inside2","Below"))) %>% 
  droplevels(.) %>% 
  # mutate(Taxon = vegan::make.cepnames(Taxon)) %>% 
  ggplot(aes(
    # y = Taxon,
    y = zone1,
    x = log(Mean_abundance + 1),
    # x = Mean_abundance,
    group = zoneplot
  )) +
  geom_hline(yintercept = seq(0.5, (4+0.5), by = 1), 
             color = "gray", linetype = "dashed") +  # Add grid lines between taxa
  geom_point(
    position = position_jitter(width = 0.01, height=0.4,seed = pi),
    alpha=0.5,size=4,
    aes(
      col = zoneplot,
      group = zoneplot),
    show.legend = FALSE
  ) +
  facet_grid(shore~Taxon)+
  scale_y_discrete(limits=rev)+
  scale_color_manual(values=cbPaletteTxt)+
  labs(
    title = "Infaunal abundances of fauna significantly differing in intertidal assemblages",
    subtitle = paste0("Monitored as part of the ",cur.yr, " SGPBM"),
    y="",
    x="log mean abundance (n+1)",
    caption=
      paste0("<b>Taxa displated showed significantly different abundances between different levels of nourishment zone:shore height variables<br>",
             "Taxon correlation structure estimated using the Ledoit-Wolf shrinkage estimator")
    )+
  theme(
    strip.text = element_text(face=2,size = 12),
    axis.text.y = element_text(face=2,size = 12),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face=2),
    plot.title = element_text(face=2),
    plot.subtitle = element_text(face=2),
    plot.caption = ggtext::element_markdown(),
  )
dev.off()


## expt: include shore as random effect ####
# Use geepack::geeglm() to specify a correlation structure for shore
# The idea is to account for the repeated measures within shore, treating
# it as a random effect.
# Fit the model using the manyglm() function
# Convert species abundance data into mvabund object
abundance_data <- mvabund(df.cur %>% dplyr::select(.,-c(year,transect,shore,,
                                                        zone1,mesh,core.area_m2,
                                                        zone_model,zoneplot)))
# Fit manyglm with GEE correlation structure for shore
mod01 <- manyglm(abundance_data ~ zone1, data = df.cur, family = "negative.binomial",
                 corStr = corCompSymm(form = ~ 1 | shore))  # Compound symmetry for repeated measures
summary(mod01)

# PLOT ####
tic("Plot")
nums <- ncol(df.tx)

df.cur$zoneplot <- factor(df.cur$zone1, levels=c("Above","Inside","Inside2","Below"))
df.cur$shore <- factor(df.cur$shore, levels=c("Mid","Low"))

png(file = "figs/inf.current.yr.png",
    width=12*ppi, height=10*ppi, res=ppi)
df.cur %>% 
  pivot_longer(cols = -c(year,transect,shore,zone1,mesh,core.area_m2,zone_model,
                         zoneplot,zone_shore), 
               names_to = "Taxon", 
               values_to = "Mean_abundance") %>% 
  mutate(Taxon = vegan::make.cepnames(Taxon)) %>% 
  ggplot(aes(
    y = Taxon,
    x = log(Mean_abundance + 1),
    # x = Mean_abundance,
    group = zoneplot
  )) +
  geom_hline(yintercept = seq(0.5, (nums+0.5), by = 1), 
             color = "gray", linetype = "dashed") +  # Add grid lines between taxa
  geom_point(position = position_jitter(width = 0.01, height=0.4,seed = pi),
             alpha=0.5,size=2,
             aes(
               #shape = zoneplot,
               col = zoneplot,
               group = zoneplot),
             show.legend = FALSE) +
  facet_grid(shore~zoneplot)+
  scale_y_discrete(limits=rev)+
  scale_color_manual(values=cbPaletteTxt)+
  labs(
    title = "Infaunal abundances of intertidal assembalges",
    subtitle = paste0("Monitored as part of the ",cur.yr, " SGPBM"),
    y="",
    x="log mean abundance (n+1)"
  )+
  theme(
    strip.text = element_text(face=2,size = 12),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face=2),
    plot.title = element_text(face=2),
    plot.subtitle = element_text(face=2),
  )
dev.off()
toc(log=TRUE)

## plot significant taxa ####
png(file = "figs/inf.current.yr.sig.png",
    width=12*ppi, height=10*ppi, res=ppi)
set.seed(42423);df.cur %>% 
  pivot_longer(cols = -c(year,transect,shore,zone1,mesh,core.area_m2,zone_model,
                         zoneplot,zone_shore), 
               names_to = "Taxon", 
               values_to = "Mean_abundance") %>% 
  filter(Taxon %in% c("Scolelepis (Scolelepis) squamata",
                      "Nematoda",
                      "Haustorius arenarius")) %>% 
  mutate(zone1 = factor(zone1, levels = c("Above","Inside",
                                          "Inside2","Below"))) %>% 
  droplevels(.) %>% 
  # mutate(Taxon = vegan::make.cepnames(Taxon)) %>% 
  ggplot(aes(
    # y = Taxon,
    y = zone1,
    x = log(Mean_abundance + 1),
    # x = Mean_abundance,
    group = zoneplot
  )) +
  geom_hline(yintercept = seq(0.5, (4+0.5), by = 1), 
             color = "gray", linetype = "dashed") +  # Add grid lines between taxa
  geom_point(
    position = position_jitter(width = 0.01, height=0.4,seed = pi),
    alpha=0.5,size=4,
    aes(
    col = zoneplot,
    group = zoneplot),
    show.legend = FALSE
    ) +
  facet_grid(shore~Taxon)+
  scale_y_discrete(limits=rev)+
  scale_color_manual(values=cbPaletteTxt)+
  labs(
    title = "Infaunal abundances of fauna significantly differing in intertidal assemblages",
    subtitle = paste0("Monitored as part of the ",cur.yr, " SGPBM"),
    y="",
    x="log mean abundance (n+1)",
    caption="Taxa displated showed significantly different abundances between"
  )+
  theme(
    strip.text = element_text(face=2,size = 12),
    axis.text.y = element_text(face=2,size = 12),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face=2),
    plot.title = element_text(face=2),
    plot.subtitle = element_text(face=2),
  )
dev.off()


unlist(tictoc::tic.log())

# DEVT: Visualisation of trends ####
df$shore <- factor(df$shore, levels =c("Mid","Low"))
df$zone1 <- factor(df$zone1,levels = c("Above","Inside","Inside2","Below","Wash"))
df %>% 
  ## remove presence-only
  filter(.,count>0) %>% 
  # remove Wash samples
  filter(., zone1 !="Wash") %>% #names()
  filter(., mesh=="1.0mm") %>% 
  ## drop cols, calculate means, widen, add zeroes, lengthen
  dplyr::select(.,year,transect, shore, zone1, mesh,taxonUSE, count) %>%
  group_by(across(!count)) %>% 
  summarise(count = mean(count)) %>% 
  pivot_wider(.,names_from = taxonUSE, values_from = count,values_fill = 0) %>% 
  pivot_longer(.,cols=-c(year,transect,shore,zone1,mesh),
               names_to = "Taxon", values_to = "Count"
  ) %>%
  ggplot(., aes(y = log(Count+1), x = as.factor(year),
                colour=zone1))+
  facet_grid(shore~zone1)+
  geom_smooth(method="loess", 
              #aes(group=taxonUSE),
              se=FALSE,alpha=0.8)+
  scale_x_discrete(limits=rev)+
  geom_point(aes(group=Taxon),alpha=0.5,
             position=position_jitter(width = .25,height=0.05),
             show.legend = FALSE)+
  coord_flip()+
  scale_colour_manual(values = cbPalette)+
  #coord_fixed(ylim(0,NA))+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(face=2),
        axis.text = element_text(face=2),
        strip.text = element_text(face=2),
        legend.position = "none")

############

df.cur %>% 
  pivot_longer(cols = -c(year,transect,shore,zone1,mesh,core.area_m2,zone_model,
                         zoneplot,zone_shore), 
               names_to = "Taxon", 
               values_to = "Mean_abundance") %>% 
  filter(Mean_abundance !=0) %>% 
  ggplot(aes(
    y = Taxon,
    x = log(Mean_abundance + 1),
    # x = Mean_abundance,
    group = zoneplot
  )) +
  geom_hline(yintercept = seq(0.5, (nums+0.5), by = 1), 
             color = "gray", linetype = "dashed") +  # Add grid lines between taxa
  geom_point(position = position_jitter(width = 0.01, height=0.4,seed = pi),
             alpha=0.5,size=2,
             aes(
               shape = zoneplot,
               col = zoneplot,
               group = zoneplot),
             show.legend = FALSE) +
  # facet_grid(shore~zoneplot)+
  scale_y_discrete(limits=rev)+
  scale_color_manual(values=cbPaletteTxt)+
  labs(
    title = "Infaunal abundances of intertidal assembalges",
    subtitle = paste0("Monitored as part of the ",cur.yr, " SGPBM"),
    y="",
    x="log mean abundance (n+1)"
  )+
  theme(
    strip.text = element_text(face=2,size = 12),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face=2),
    plot.title = element_text(face=2),
    plot.subtitle = element_text(face=2),
  )
  

# tidy up ####
rm(list=ls(pattern = "^df"));rm(list=ls(pattern = "^m"))
rm(list=ls(pattern = "^ano"));rm(list=ls(pattern = "^cb"))
rm(cur.yr,fol,gisfol,nums,perm,ppi,sum_zero,projfol,source_file)

detach("package:tidyverse", unload=TRUE)
detach("package:vegan", unload=TRUE)
detach("package:mvabund", unload=TRUE)
detach("package:tictoc", unload=TRUE)
