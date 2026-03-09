## 200_imp.int.sed.version2_USE.R ##

## generate plots & compare values from gradistat (G2Sd) outputs

# Set up ####
### load packages ####
ld_pkgs <- c("tidyverse","tictoc","ggthemes","ggh4x","patchwork",
             "ggridges","lme4","lmerTest","emmeans","multcomp",
             "multcompView")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)

tictoc::tic.clearlog();tic("SET UP")

### load metadata ####
source("R/00_meta_setMeta.R")

# load data
df0 <- readxl::read_xlsx(paste0(fol,"sed.data.ALL.USE.xlsx"), sheet = "GradistatOut")
# this version is the data generated in Gradistat
toc(log=TRUE)

# define factors
df0$transect <- factor(df0$transect,
                       levels = c(
                         "T1N", "T1", "T1S", "T4", "T7", "T8", "T11", "T12",
                         "T13", "T15", "T17","T20", "T21", "T22", "T23", "T24",
                         "T25", "T26", "WA1", "WA2", "WA3", "WA4", "WA5", "WA6"
                         ))
df0$shore <- factor(df0$shore, levels = c("Upper","Mid","Low","Surf"))
df0$zone1 <- factor(df0$zone1, levels = c("Above", "Inside", "Inside2","Below","Wash"))

df0$sedName <- factor(df0$`SEDIMENT NAME`,levels = rev(c(
  "Sandy Very Coarse Gravel",
  "Sandy Coarse Gravel",
  "Sandy Medium Gravel",
  "Fine Gravel",
  "Sandy Fine Gravel",
  "Fine Silty Sandy Fine Gravel",
  "Very Fine Gravel",
  "Medium Silty Sandy Very Fine Gravel",
  "Sandy Very Fine Gravel",
  "Slightly Very Fine Gravelly Very Coarse Sand",
  "Very Fine Gravelly Very Coarse Sand",
  "Slightly Medium Gravelly Coarse Sand",
  "Slightly Fine Gravelly Coarse Sand",
  "Slightly Very Fine Gravelly Coarse Sand",
  "Very Fine Gravelly Coarse Sand",
  "Very Coarse Gravelly Medium Sand",
  "Slightly Coarse Gravelly Medium Sand",
  "Coarse Gravelly Medium Sand",
  "Medium Gravelly Medium Sand",
  "Slightly Medium Gravelly Medium Sand",
  "Slightly Fine Gravelly Medium Sand",
  "Fine Gravelly Medium Sand",
  "Slightly Very Fine Gravelly Medium Sand",
  "Slightly Very Fine Gravelly Medium Silty Fine Sand",
  "Very Fine Gravelly Medium Sand",
  "Well Sorted Medium Sand",
  "Moderately Well Sorted Medium Sand",
  "Slightly Very Fine Gravelly Fine Silty Medium Sand",
  "Very Coarse Silty Fine Sand",
  "Coarse Gravelly Fine Sand",
  "Slightly Coarse Gravelly Fine Sand",
  "Slightly Coarse Gravelly Fine Silty Fine Sand",
  "Medium Gravelly Fine Sand",
  "Medium Gravelly Fine Silty Fine Sand",
  "Slightly Medium Gravelly Very Coarse Silty Fine Sand",
  "Slightly Medium Gravelly Medium Silty Fine Sand",
  "Slightly Medium Gravelly Fine Sand",
  "Slightly Medium Gravelly Fine Silty Fine Sand",
  "Slightly Fine Gravelly Fine Sand",
  "Slightly Fine Gravelly Fine Silty Fine Sand",
  "Slightly Very Fine Gravelly Coarse Silty Fine Sand",
  "Slightly Very Fine Gravelly Very Coarse Silty Fine Sand",
  "Slightly Very Fine Gravelly Fine Sand",
  "Slightly Very Fine Gravelly Fine Silty Fine Sand",
  "Fine Gravelly Fine Sand",
  "Fine Gravelly Fine Silty Fine Sand",
  "Slightly Fine Gravelly Very Coarse Silty Fine Sand",
  "Slightly Fine Gravelly Medium Silty Fine Sand",
  "Very Fine Gravelly Fine Sand",
  "Very Fine Gravelly Very Coarse Silty Fine Sand",
  "Very Fine Gravelly Fine Silty Fine Sand",
  "Moderately Well Sorted Fine Sand",
  "Well Sorted Fine Sand",
  "Moderately Well Sorted Very Fine Sand",
  "Slightly Medium Gravelly Very Fine Sand",
  "Slightly Coarse Gravelly Very Coarse Silty Very Fine Sand",
  "Slightly Medium Gravelly Very Coarse Silty Very Fine Sand",
  "Slightly Fine Gravelly Very Coarse Silty Very Fine Sand",
  "Slightly Fine Gravelly Medium Silty Very Fine Sand",
  "Slightly Medium Gravelly Fine Silty Very Fine Sand",
  "Slightly Very Fine Gravelly Very Fine Sand",
  "Slightly Very Fine Gravelly Very Coarse Silty Very Fine Sand",
  "Fine Sandy Very Coarse Silt",
  "Slightly Very Fine Gravelly Fine Sandy Very Coarse Silt",
  "Slightly Very Fine Gravelly Very Fine Sandy Very Coarse Silt",
  "Slightly Very Fine Gravelly Coarse Silt",
  "Slightly Very Fine Gravelly Fine Sandy Medium Silt",
  "Slightly Very Fine Gravelly Very Fine Sandy Medium Silt",
  "Slightly Medium Gravelly Fine Sandy Fine Silt",
  "Slightly Very Fine Gravelly Fine Sandy Fine Silt",
  "Fine Silt")))

# plots ####
zones <- unique(df0$zone1)
bg_cols <- c("#B5EAFE","#FFB935","#23B888","#238ACD","#E894BC"
  
  # "#E894BC","#238ACD","#23B888","#FFB935","#B5EAFE"
  )
bg_cols <- setNames(bg_cols[seq_along(zones)], zones)

strip_elems <- lapply(zones, function(z)
  element_rect(fill = bg_cols[[z]], color = "black", linewidth = 1)
  )

## Sediment TYPE ####
cbPaletteshr_border <- c("#C44845","#735473","#695723","#242323")
## create broad "type" variable & assign colours
df0$SedType <- str_extract(df0$sedName, "\\S+$")

SedCols <- ifelse(df0$SedType == "Gravel","green",
                  ifelse(df0$SedType == "Sand","blue",
                         ifelse(df0$SedType == "Silt","sienna",NA
                         )))
tf <- df0$method=="5cm"

SedCols <- SedCols[tf]

df0 <- df0 %>% mutate(sediment_col = case_when(
  SedType == "Gravel" ~"darkblue",
  SedType == "Sand" ~"orange",
  SedType == "Silt" ~"brown",
  TRUE ~ NA
))

png(
  file = "figs/sed_types_ts.png",
  width = 15 * ppi,
  height = 10 * ppi,
  res = ppi
)
set.seed(222)
df0 %>% 
  filter(method=="5cm") %>% 
  ggplot(.,
       aes(
         x = Year,
         y = sedName
         )
       )+
  geom_hline(yintercept = c(7.5,58.5),
             col = "grey", linetype = "dashed")+
  geom_jitter(aes(fill = shore,
                  shape = shore), 
              position = position_jitter(width = 0.25,
                                         height = 0.25),
              size = 7, alpha = 0.5)+
  # scale_shape_manual(values=c(15:18))+
  scale_shape_manual(values=c(21:24))+
  # scale_shape_manual(values=c(1,2,3,5))+
  scale_fill_manual(values=cbPaletteshr) +
  scale_color_manual(values=cbPaletteshr_border)+
  xlab("") + ylab("")+
  # facet_wrap(.~zone1)+
  ggh4x::facet_grid2(.~zone1,
    strip = strip_themed(
      background_x = strip_elems)
  )+
  # geom_hline(yintercept = c(seq(from = 0.5, to = 22.5, by = 1)),
  #            col = "grey", linetype = "dashed")+
  # geom_vline(xintercept = c(0.5,1.5,2.5,3.5,4.5),
  #            col = "grey", linetype = "dashed")+
  labs(title = paste0("Sediment types recorded since 2009"))+
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(face=2),
    axis.text.y = element_text(size=8,face=2),
    # axis.text.y = ggtext::element_markdown(size=8,face=2,colour = SedCols),
    legend.text = element_text(size = 10,face=2),
    plot.title = element_text(face=2),
    strip.text = element_text(face=2),
    )
dev.off()

## SORTING ####
png(
  file = "figs/sed_sort_cur.png",
  width = 15 * ppi,
  height = 10 * ppi,
  res = ppi
)
set.seed(14236);df0 %>% 
  mutate(`Folk and Ward (description) SORTING` = factor(df0$`Folk and Ward (description) SORTING`,
                                                        levels=c(
                                                          "Very Poorly Sorted",
                                                          "Poorly Sorted",
                                                          "Moderately Sorted",
                                                          "Moderately Well Sorted",
                                                          "Well Sorted"
                                                        ))) %>% 
  filter(Year == cur.yr) %>% 
  filter(method=="5cm") %>% 
ggplot(.,
       aes(
         x = 0,
         y = `Folk and Ward (description) SORTING`,
         shape = shore,
         colour = shore
       ))+
  geom_hline(yintercept = c(1.5:5.5),colour = "grey",linetype=2)+
  # geom_vline(xintercept = c(1.5:4.5),colour = "grey",linetype=2)+
  geom_jitter(aes(fill = shore,
                  shape = shore), 
              position = position_jitter(width = 0.05,
                                         height = 0.25),
              size = 7, alpha = 0.5)+
  scale_shape_manual(values=c(21:24))+
  # scale_shape_manual(values=c(1,2,3,5))+
  scale_fill_manual(values=cbPaletteshr) +
  scale_color_manual(values=cbPaletteshr_border)+
  facet_wrap(.~zone1,scales = "free_x")+
  labs(
    title = paste0("Sediment sorting classifications recorded in ",cur.yr))+
  ggh4x::facet_grid2(.~zone1,
                     strip = strip_themed(
                       background_x = strip_elems)
  )+
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=12,face=2),
    legend.text = element_text(size = 12,face=2),
    plot.title = element_text(face=2),
    axis.title = element_blank(),
    strip.text = element_text(face=2),
  )
dev.off()

# sorting_TS ####
df_sort <- df0 %>% 
  filter(method == "5cm") %>% 
  filter(zone1 != "Wash") %>% 
  filter(Year == cur.yr) %>% 
  mutate(shore = factor(shore, levels = c("Surf","Low","Mid","Upper")))

summary(mod1 <- lmer(`Folk and Ward (phi) SORTING` ~ shore*zone1 +(1|transect),
           data=df_sort))
anova(mod1)

# Pairwise post‑hoc comparisons
emmeans(mod1, pairwise ~ shore)
emmeans(mod1, pairwise ~ zone1)

# Pairwise comparisons of shore heights within zones
height_within_mgmt <- emmeans(mod1, 
                              pairwise ~ shore | zone1)

height_within_mgmt
plot(height_within_mgmt, comparison = TRUE)+
  theme(axis.title.y = element_blank(),
        axis.text = element_text(face=2)
  )


png(
  file = "figs/sed_sort_ts.png",
  width = 15 * ppi,
  height = 10 * ppi,
  res = ppi
)
set.seed(14236);df0 %>% 
  mutate(`Folk and Ward (description) SORTING` = factor(df0$`Folk and Ward (description) SORTING`,
                                                        levels=c(
                                                          "Very Poorly Sorted",
                                                          "Poorly Sorted",
                                                          "Moderately Sorted",
                                                          "Moderately Well Sorted",
                                                          "Well Sorted"
                                                        ))) %>%
  filter(method=="5cm") %>% 
  ggplot(.,
         aes(
           x = Year,
           y = `Folk and Ward (description) SORTING`,
           # shape = shore,
           # colour = zone1
         ))+
  geom_hline(yintercept = c(1.5:5.5),colour = "grey",linetype=2)+
  geom_jitter(aes(fill = zone1,
                  shape = shore), 
              position = position_jitter(width = 0.05,
                                         height = 0.25),
              size = 7, alpha = 0.5,
              show.legend = FALSE
              )+
  scale_shape_manual(values=c(21:24))+
  scale_fill_manual(values=cbPalette[c(1:4,7)]) +
  scale_color_manual(values=cbPalette[c(1:4,7)])+
  labs(
    title = paste0("Sediment sorting classifications recorded since ",min(df0$Year)))+
  ggh4x::facet_grid2(shore~zone1,
                     strip = strip_themed(
                       background_x = strip_elems)
  )+
  theme(
    legend.title = element_blank(),
    # axis.text = element_text(size = 12),
    # axis.text.x = element_blank(),
    axis.text = element_text(size=12,face=2),
    legend.text = element_text(size = 12,face=2),
    plot.title = element_text(face=2),
    axis.title = element_blank(),
    strip.text = element_text(face=2,size=12),
  )
dev.off()

# skew_TS ####
png(
  file = "figs/sed_skew_ts.png",
  width = 15 * ppi,
  height = 10 * ppi,
  res = ppi
)
set.seed(14236);df0 %>% 
  mutate(`Folk and Ward (description) SKEWNESS` = factor(df0$`Folk and Ward (description) SKEWNESS`,
                                                        levels=c(
                                                          "Very Coarse Skewed",
                                                          "Coarse Skewed",
                                                          "Symmetrical",
                                                          "Fine Skewed",
                                                          "Very Fine Skewed"
                                                        ))) %>%
  filter(method=="5cm") %>% 
  ggplot(.,
         aes(
           x = Year,
           y = `Folk and Ward (description) SKEWNESS`,
           # shape = shore,
           # colour = zone1
         ))+
  geom_hline(yintercept = c(1.5:5.5),colour = "grey",linetype=2)+
  geom_jitter(aes(fill = zone1,
                  shape = shore), 
              position = position_jitter(width = 0.05,
                                         height = 0.25),
              size = 7, alpha = 0.5,
              show.legend = FALSE
  )+
  scale_shape_manual(values=c(21:24))+
  scale_fill_manual(values=cbPalette[c(1:4,7)]) +
  scale_color_manual(values=cbPalette[c(1:4,7)])+
  labs(
    title = paste0("Sediment skewness classifications recorded since ",min(df0$Year)))+
  ggh4x::facet_grid2(shore~zone1,
                     strip = strip_themed(
                       background_x = strip_elems)
  )+
  theme(
    legend.title = element_blank(),
    # axis.text = element_text(size = 12),
    # axis.text.x = element_blank(),
    axis.text = element_text(size=12,face=2),
    legend.text = element_text(size = 12,face=2),
    plot.title = element_text(face=2),
    axis.title = element_blank(),
    strip.text = element_text(face=2,size=12),
  )
dev.off()

# kurt_TS ####
png(
  file = "figs/sed_kurt_ts.png",
  width = 15 * ppi,
  height = 10 * ppi,
  res = ppi
)
set.seed(14236);df0 %>% 
  mutate(`Folk and Ward (description) KURTOSIS` = factor(df0$`Folk and Ward (description) KURTOSIS`,
                                                         levels=c(
                                                           "Extremely Leptokurtic",
                                                           "Very Leptokurtic",
                                                           "Leptokurtic",
                                                           "Mesokurtic",
                                                           "Platykurtic",
                                                           "Very Platykurtic"
                                                           ))) %>%
  filter(method=="5cm") %>% 
  ggplot(.,
         aes(
           x = Year,
           y = `Folk and Ward (description) KURTOSIS`,
           # shape = shore,
           # colour = zone1
         ))+
  geom_hline(yintercept = c(1.5:5.5),colour = "grey",linetype=2)+
  geom_jitter(aes(fill = zone1,
                  shape = shore), 
              position = position_jitter(width = 0.05,
                                         height = 0.25),
              size = 7, alpha = 0.5,
              show.legend = FALSE
  )+
  scale_shape_manual(values=c(21:24))+
  scale_fill_manual(values=cbPalette[c(1:4,7)]) +
  scale_color_manual(values=cbPalette[c(1:4,7)])+
  labs(
    title = paste0("Sediment kurtosis classifications recorded since ",min(df0$Year)))+
  ggh4x::facet_grid2(shore~zone1,
                     strip = strip_themed(
                       background_x = strip_elems)
  )+
  theme(
    legend.title = element_blank(),
    # axis.text = element_text(size = 12),
    # axis.text.x = element_blank(),
    axis.text = element_text(size=12,face=2),
    legend.text = element_text(size = 12,face=2),
    plot.title = element_text(face=2),
    axis.title = element_blank(),
    strip.text = element_text(face=2,size=12),
  )
dev.off()


# phi_TS ####
png(
  file = "figs/sed_phi_ts.png",
  width = 15 * ppi,
  height = 10 * ppi,
  res = ppi
)
set.seed(14236);df0 %>% 
  filter(method=="5cm") %>% 
  ggplot(.,
         aes(
           x = Year,
           y = `Folk and Ward (phi) MEAN`,
           # shape = shore,
           # colour = zone1
         ))+
  geom_hline(yintercept = c(-1,4),colour = "grey",linetype=2)+
  geom_jitter(aes(fill = zone1,
                  shape = shore), 
              position = position_jitter(width = 0.05,
                                         height = 0.25),
              size = 7, alpha = 0.5,
              show.legend = FALSE
  )+
  scale_shape_manual(values=c(21:24))+
  scale_fill_manual(values=cbPalette[c(1:4,7)]) +
  scale_color_manual(values=cbPalette[c(1:4,7)])+
  geom_smooth(method = "gam")+
  labs(
    title = paste0("Mean sediment phi values recorded since ",min(df0$Year)))+
  ggh4x::facet_grid2(shore~zone1,
                     strip = strip_themed(
                       background_x = strip_elems)
  )+
  theme(
    legend.title = element_blank(),
    # axis.text = element_text(size = 12),
    # axis.text.x = element_blank(),
    axis.text = element_text(size=12,face=2),
    legend.text = element_text(size = 12,face=2),
    plot.title = element_text(face=2),
    axis.title = element_blank(),
    strip.text = element_text(face=2,size=12),
  )->pl_phi
print(pl_phi)
dev.off()

# phi_Dx ####
png(
  # file = "figs/sed_phi_ts.png",
  width = 15 * ppi,
  height = 10 * ppi,
  res = ppi
)
set.seed(14236);df0 %>% 
  filter(method=="5cm") %>% 
  select(Code:zone1,"D10 (f)":"D90 (f)") %>% 
  rename(D10 = "D10 (f)",
         D50 = "D50 (f)",
         D90 = "D90 (f)"
  ) %>% 
  pivot_longer(cols = D10:D90) %>% 
  ggplot(.,
         aes(
           x = Year,
           y = value,
           ))+
  geom_hline(yintercept = c(-1,4),colour = "grey",linetype=2)+
  geom_jitter(aes(fill = name,
                  shape = name), 
              position = position_jitter(width = 0.05,
                                         height = 0.25),
              size = 4, alpha = 0.5,
              show.legend = FALSE
  )+
  scale_shape_manual(values=c(21,22,24))+
  scale_fill_manual(values=cbPalette[c(1:4,7)]) +
  scale_color_manual(values=cbPalette[c(1:4,7)])+
  geom_smooth(method = "gam",se=FALSE,
              aes(colour = name)
              )+
  labs(
    title = paste0("Mean sediment phi values recorded since ",min(df0$Year)))+
  ggh4x::facet_grid2(shore~zone1,
                     strip = strip_themed(
                       background_x = strip_elems)
  )+
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size=12,face=2),
    legend.text = element_text(size = 12,face=2),
    plot.title = element_text(face=2),
    axis.title = element_blank(),
    strip.text = element_text(face=2,size=12),
  )->pl_Dx
print(pl_Dx)
dev.off()

# phi_TS & Dx####
png(
  file = "figs/sed_phiDx_ts.png",
  width = 15 * ppi,
  height = 13 * ppi,
  res = ppi
)
set.seed(14236);df0 %>% 
  filter(method=="5cm") %>% 
  ggplot(.,
         aes(
           x = Year,
           y = `Folk and Ward (phi) MEAN`,
           # shape = shore,
           # colour = zone1
         ))+
  geom_hline(yintercept = c(-1,4),colour = "grey",linetype=2)+
  geom_jitter(aes(fill = zone1,
                  shape = shore), 
              position = position_jitter(width = 0.05,
                                         height = 0.25),
              size = 7, alpha = 0.5,
              show.legend = FALSE
  )+
  scale_shape_manual(values=c(21:24))+
  scale_fill_manual(values=cbPalette[c(1:4,7)]) +
  scale_color_manual(values=cbPalette[c(1:4,7)])+
  geom_smooth(method = "gam")+
  # labs(
    # title = paste0("Mean sediment phi values recorded since ",min(df0$Year)))+
  ggh4x::facet_grid2(shore~zone1,
                     strip = strip_themed(
                       background_x = strip_elems)
  )+
  theme(
    legend.title = element_blank(),
    # axis.text = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.text = element_text(size=12,face=2),
    legend.text = element_text(size = 12,face=2),
    plot.title = element_text(face=2),
    axis.title = element_blank(),
    strip.text = element_text(face=2,size=14),
  )->pl_phi

set.seed(14236);df0 %>% 
  filter(method=="5cm") %>% 
  select(Code:zone1,"D10 (f)":"D90 (f)") %>% 
  rename(D10 = "D10 (f)",
         D50 = "D50 (f)",
         D90 = "D90 (f)"
  ) %>% 
  pivot_longer(cols = D10:D90) %>% 
  ggplot(.,
         aes(
           x = Year,
           y = value,
         ))+
  geom_hline(yintercept = c(-1,4),colour = "grey",linetype=2)+
  geom_smooth(method = "gam",se=FALSE,
              aes(colour = name,lty=name)
  )+
  geom_jitter(aes(fill = name,
                  shape = name), 
              position = position_jitter(width = 0.05,
                                         height = 0.25),
              size = 2, alpha = 0.5,
              show.legend = FALSE
  )+
  scale_shape_manual(values=c(21,22,24))+
  scale_fill_manual(values=cbPalette[c(1:4,7)]) +
  scale_color_manual(values=cbPalette[c(1:4,7)])+
  ggh4x::facet_grid2(shore~zone1,
                     strip = strip_themed(
                       background_x = strip_elems)
  )+
  labs(caption = "Lines represent generalised additive model predictor")+
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size=14,face=2),
    legend.text = element_text(size = 14,face=2),
    plot.title = element_text(face=2),
    plot.caption = element_text(face=2, size = 12),
    axis.title = element_blank(),
    strip.text = element_text(face=2,size=14),
    strip.text.x = element_blank()
  )->pl_Dx

pl_out <- pl_phi/pl_Dx & patchwork::plot_annotation(tag_levels="A")
print(pl_out)
dev.off()

