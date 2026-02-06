## 200_imp.int.sed.2.R ##

## generate plots & compare values from gradistat (G2Sd) outputs

# Set up ####
### load packages ####
ld_pkgs <- c("tidyverse","tictoc","ggthemes","ggh4x",
             "ggridges","lme4","lmerTest")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)

tictoc::tic.clearlog();tic("SET UP")

### load metadata ####
source("R/00_meta_setMeta.R")

# load data
df <- readRDS(file="output/models/grad_all_raw.Rdat")
toc(log=TRUE)

zones <- unique(df$index$zone1)
bg_cols <- c("#E894BC","#238ACD","#23B888","#FFB935","#B5EAFE")
bg_cols <- setNames(bg_cols[seq_along(zones)], zones)

strip_elems <- lapply(zones, function(z)
  element_rect(fill = bg_cols[[z]], color = "black", linewidth = 1)
  )

df$index %>% #names()
  ungroup() %>% 
  dplyr::filter(method=="5cm") %>% 
  ## keep only D10,D50, and D90 & lengthen
  dplyr::select(c(year,shore,zone1,D10,D50,D90)) %>% 
  tidyr::pivot_longer(.,cols = c(D10,D50,D90)
                        ) %>% 
  ## convert to phi
  mutate(value_phi = -log(value/1000,base = 2)) %>% 
  ## group & summarise
  ggplot(.,
         aes(
           x = year,
           y = value_phi,
           colour = shore,
           shape = shore,
         )) +
  geom_hline(yintercept = 0, lty=2, colour = "lightgrey")+
  geom_jitter(width = 0.25,
              size = 4, alpha = 0.7
              )+
  scale_shape_manual(values=c(15:18))+
  scale_fill_manual(values=cbPaletteshr) +
  scale_color_manual(values=cbPaletteshr) +
  geom_smooth(
    method="gam",
    se=FALSE,
    )+
  ggh4x::facet_grid2(
    rows = vars(name),
    cols = vars(zone1),
    strip = strip_themed(
      background_x = strip_elems)
  )+
  # facet_grid(name~zone1) +
  theme(
    plot.title = element_text(face=2),
    strip.text = element_text(face=2),
    axis.text.x = element_text(face=2,angle=270,vjust=0.5),
    axis.text.y = element_text(face=2),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    # legend.direction = "horizontal",
    # legend.position = "inside",
    legend.position.inside = c(0.8,0.1),
    legend.title = element_blank(),
    legend.text = element_text(face = 2)
  )

# Dn MODELS ####
## reorganise factors
df$index$zone1 <- factor(df$index$zone1,
                         levels = c(
                           "Inside",
                           "Above", "Inside2","Below","Wash"
                         ))

## generate data for current year only & convert to phi values
df_cur <- df$index %>% filter(year==cur.yr) %>% 
  filter(method == "5cm") %>% 
  filter(zone1 !="Wash")
df_cur$D10_phi <- -log(df_cur$D10/1000,base = 2)
df_cur$D50_phi <- -log(df_cur$D50/1000,base = 2)
df_cur$D90_phi <- -log(df_cur$D90/1000,base = 2)

## D10 phi ####
anova(mod2 <- lmer(D10_phi ~ zone1 + (1|shore),
                   data = df_cur,
                   REML=TRUE))
summary(mod2)
d <- as.data.frame(ls_means(mod2, test.effs = "Group",pairwise = TRUE))
d[d$`Pr(>|t|)`<0.051,]
# sjPlot::plot_model(mod2,show.values=TRUE, show.p=TRUE)
visreg::visreg(mod2)
rm(mod2,d)

## D50_phi ####
anova(mod2 <- lmer(D50_phi ~ zone1 + (1|shore),
                   data = df_cur,
                   REML=TRUE))
summary(mod2)
d <- as.data.frame(ls_means(mod2, test.effs = "Group",pairwise = TRUE))
d[d$`Pr(>|t|)`<0.051,]
# sjPlot::plot_model(mod2,show.values=TRUE, show.p=TRUE)
visreg::visreg(mod2)
rm(mod2,d)

## D90_phi ####
anova(mod2 <- lmer(D90_phi ~ zone1 + (1|shore),
                   data = df_cur,
                   REML=TRUE))
summary(mod2)
d <- as.data.frame(ls_means(mod2, test.effs = "Group",pairwise = TRUE))
d[d$`Pr(>|t|)`<0.051,]
# sjPlot::plot_model(mod2,show.values=TRUE, show.p=TRUE)
visreg::visreg(mod2)
rm(mod2,d)

# phi values, etc ####
df$stat$fowa %>% 
  filter(year==cur.yr) %>% 
  filter(method == "5cm") %>% 
  filter(zone1 != "Wash") -> df_stat_cur

names(df_stat_cur)

df_stat_cur$zone1 <- factor(df_stat_cur$zone1, levels = c(
  "Inside","Above","Inside2","Below"
  ))

## phi ####
anova(mod2 <- lmer(mean.fw.phi ~ zone1 + (1|shore),
                   data = df_stat_cur,
                   REML=TRUE))
summary(mod2)
d <- as.data.frame(ls_means(mod2, test.effs = "Group",pairwise = TRUE))
d[d$`Pr(>|t|)`<0.051,]
# sjPlot::plot_model(mod2,show.values=TRUE, show.p=TRUE)
visreg::visreg(mod2)
rm(mod2,d)

## sorting ####
anova(mod2 <- lmer(sd.fw.phi ~ zone1 + (1|shore),
                   data = df_stat_cur,
                   REML=TRUE))
summary(mod2)
d <- as.data.frame(ls_means(mod2, test.effs = "Group",pairwise = TRUE))
d[d$`Pr(>|t|)`<0.051,]
# sjPlot::plot_model(mod2,show.values=TRUE, show.p=TRUE)
visreg::visreg(mod2)
rm(mod2,d)

## Skew ####
df_stat_cur[which.min(df_stat_cur$skewness.fw.phi),]
df_stat_cur[which.max(df_stat_cur$skewness.fw.phi),]

anova(mod2 <- lmer(skewness.fw.phi ~ zone1 + (1|shore),
                   data = df_stat_cur,
                   REML=TRUE))
summary(mod2)
d <- as.data.frame(ls_means(mod2, test.effs = "Group",pairwise = TRUE))
d[d$`Pr(>|t|)`<0.051,]
# sjPlot::plot_model(mod2,show.values=TRUE, show.p=TRUE)
visreg::visreg(mod2)
rm(mod2,d)

## Kurtosis ####
df_stat_cur[which.min(df_stat_cur$kurtosis.fw.phi),]
df_stat_cur[which.max(df_stat_cur$kurtosis.fw.phi),]

anova(mod2 <- lmer(kurtosis.fw.phi ~ zone1 + (1|shore),
                   data = df_stat_cur,
                   REML=TRUE))
summary(mod2)
d <- as.data.frame(ls_means(mod2, test.effs = "Group",pairwise = TRUE))
d[d$`Pr(>|t|)`<0.051,]
# sjPlot::plot_model(mod2,show.values=TRUE, show.p=TRUE)
visreg::visreg(mod2)
rm(mod2,d)

# PLOT SEDIMENT TS ####
df$stat$fowa %>% 
  ungroup() %>% 
  dplyr::filter(method=="5cm") %>% 
  dplyr::filter(year==cur.yr) %>% 
  ggplot(.,
         aes(
           x = year,
           # y = sd.fw.phi,
           # y = mean.fw.phi,
           y = kurtosis.fw.phi,
           colour = shore,
           shape = shore,
         )) +
  ##kurtosis lines:
  geom_hline(yintercept = c(0.67,0.9,1.11,1.5,3))+
  # geom_hline(yintercept = 0, lty=2, colour = "lightgrey")+
  geom_jitter(width = 0.25,
              size = 4, alpha = 0.7
  ) +
  scale_shape_manual(values=c(15:18))+
  scale_fill_manual(values=cbPaletteshr) +
  scale_color_manual(values=cbPaletteshr) +
  geom_smooth(
    method="gam",
    se=FALSE,
  )+
  # ggh4x::facet_grid2(
  #   rows = vars(name),
  #   cols = vars(zone1),
  #   strip = strip_themed(
  #     background_x = strip_elems)
  # )+
  # facet_grid(name~zone1) +
  facet_wrap(.~zone1)+
  theme(
    plot.title = element_text(face=2),
    strip.text = element_text(face=2),
    axis.text.x = element_text(face=2,angle=270,vjust=0.5),
    axis.text.y = element_text(face=2),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    # legend.direction = "horizontal",
    # legend.position = "inside",
    legend.position.inside = c(0.8,0.1),
    legend.title = element_blank(),
    legend.text = element_text(face = 2)
  )

# PLOT CURRENT YEAR STATS ####
## Kurtosis ####
df$stat$fowa %>% 
  ungroup() %>% 
  dplyr::filter(method=="5cm") %>% 
  dplyr::filter(year==cur.yr) %>% 
  ## assign kurtosis label
  dplyr::mutate(kurt_lab = case_when(
    kurtosis.fw.phi <.67 ~ "Very platykurtic",
    kurtosis.fw.phi >.67 & kurtosis.fw.phi <.9 ~"Platykurtic",
    kurtosis.fw.phi >.9 & kurtosis.fw.phi <1.11 ~"Mesokurtic",
    kurtosis.fw.phi >1.11 & kurtosis.fw.phi <1.5 ~"Leptokurtic",
    kurtosis.fw.phi >1.5 & kurtosis.fw.phi <3 ~"Very leptokurtic",
    kurtosis.fw.phi >3 ~"Extremely leptokurtic",
    TRUE~"NA"
  )) %>% 
  dplyr::mutate(kurt_lab = factor(kurt_lab, levels = c(
    "Very platykurtic",
    "Platykurtic",
    "Mesokurtic",
    "Leptokurtic",
    "Very leptokurtic",
    "Extremely leptokurtic"
  ))) %>% 
  ggplot(.,
         aes(
           x = zone1,
           y = kurt_lab,
           colour = shore,
           shape = shore,
         )) +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5),lty=2,col="grey")+
  geom_hline(yintercept = c(1.5,2.5,3.5,4.5),lty=2,col="grey")+
  geom_jitter(width = 0.25,
              size = 4, alpha = 0.7
  ) +
  labs(
    title = "Kurtosis"
  )+
  theme(
    plot.title = element_text(face=2),
    strip.text = element_text(face=2),
    axis.text.x = element_text(face=2,vjust=0.5),
    axis.text.y = element_text(face=2),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    # legend.direction = "horizontal",
    # legend.position = "inside",
    legend.position.inside = c(0.8,0.1),
    legend.title = element_blank(),
    legend.text = element_text(face = 2)
  )

## Sorting ####
df$stat$fowa %>% 
  ungroup() %>% 
  dplyr::filter(method=="5cm") %>% 
  dplyr::filter(year==cur.yr) %>% 
  ## assign kurtosis label
  dplyr::mutate(sort_lab = case_when(
    sd.fw.phi <1.27 ~ "Very well sorted",
    sd.fw.phi >1.27 & sd.fw.phi <1.41 ~"Well sorted",
    sd.fw.phi >1.41 & sd.fw.phi <1.62 ~"Moderately well sorted",
    sd.fw.phi >1.62 & sd.fw.phi <2 ~"Moderately sorted",
    sd.fw.phi >2 & sd.fw.phi <4 ~"Poorly sorted",
    sd.fw.phi >4 & sd.fw.phi <16 ~"Very poorly sorted",
    sd.fw.phi >16 ~"Extremely poorly sorted",
    TRUE~"NA"
  )) %>% View()
  dplyr::mutate(kurt_lab = factor(kurt_lab, levels = c(
    "Very platykurtic",
    "Platykurtic",
    "Mesokurtic",
    "Leptokurtic",
    "Very leptokurtic",
    "Extremely leptokurtic"
  ))) %>% 
  ggplot(.,
         aes(
           x = zone1,
           y = kurt_lab,
           colour = shore,
           shape = shore,
         )) +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5),lty=2,col="grey")+
  geom_hline(yintercept = c(1.5,2.5,3.5,4.5),lty=2,col="grey")+
  geom_jitter(width = 0.25,
              size = 4, alpha = 0.7
  ) +
  labs(
    title = "Kurtosis"
  )+
  theme(
    plot.title = element_text(face=2),
    strip.text = element_text(face=2),
    axis.text.x = element_text(face=2,vjust=0.5),
    axis.text.y = element_text(face=2),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    # legend.direction = "horizontal",
    # legend.position = "inside",
    legend.position.inside = c(0.8,0.1),
    legend.title = element_blank(),
    legend.text = element_text(face = 2)
  )
