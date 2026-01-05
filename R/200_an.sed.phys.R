# 200_an.sed.phys.R ####
### import & analyse cone penetrometer and beach angle data

# Set up ####
### load packages ####
ld_pkgs <- c("tidyverse","ggthemes","lmerTest","effects","tictoc",
             "sjPlot","visreg","mgcv","gratia","patchwork")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)

### load metadata ####
source("R/00_meta_setMeta.R")

## import & format data ####
df <- as_tibble(read.csv(paste0(fol,"sed.phys.ts.csv"),
                         stringsAsFactors = FALSE))

##order factors
df$transect <- factor(df$transect,
                      levels =
                        c("T1N","T1","T1S","T4", "T7", "T8",
                          "T11", "T12", "T13", "T15", "T17",
                          "T24","T25","T26",
                          "T20", "T21","T22", "T23","WA1","WA2","WA3","WA4",
                          "WA5","WA6"))
df$zone1 <- factor(df$zone1,
                   levels = c("Above", "Inside","Inside2","Below","Wash"))

df$shore <- factor(df$shore,
                   levels=c("Upper","Mid","Low","Surf", "ALL"))

cur_dat0 <- droplevels(subset(df, year == cur.yr))

cur_dat <- droplevels(cur_dat0 %>% filter(., shore != "ALL"))

##=============###
# Wave type   ####
##=============###
subset(cur_dat, type == "WaveClass") %>% ## pull out angle data
  View()

##=============###
# Beach Angle  ####
##=============###
ang_cu <- cur_dat %>% filter(., type == "angle") %>% ## pull out angle data
  mutate(., value = as.numeric(value))

png(file = paste0("figs/sed.",cur.yr,".sed.mor.angle.png"),
    width=12*ppi, height=6*ppi, res=ppi)
set.seed(pi); ggplot(data = ang_cu, aes(y = value, x = zone1,
                                        fill=zone1))+
  geom_boxplot(outlier.colour=NA)+
  geom_jitter(aes(shape = zone1),
              width = 0.3, height = 0.05,size = 3,stroke = 1,
              show.legend = FALSE)+
  scale_shape_manual(values = c(21:(length(unique(ang_cu$zone1))+21)))+
  labs(title = paste0("Beach slopes recorded in the ",cur.yr," monitoring programme"))+
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom",#c(.5,1.02),
        legend.direction = "horizontal",
        legend.key.size = unit(1,"cm"),#change legend size!
        legend.title=element_blank(),
        legend.text=element_text(size = 12,
                                 face = "bold"),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(0,0,0,0),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 12),
        strip.text = element_text(face="bold"),
        strip.background = element_rect(color = "black",fill = "grey95", size = 1),
  )+
  xlab("") + ylab("Beach angle (°)")+
  scale_fill_manual(values=cbPaletteFill[c(1:4,7)])+
  scale_x_discrete(breaks=NULL)+
  facet_wrap(~shore)+
  guides(fill=guide_legend(nrow=1))
dev.off()

# png(file = paste0("output/figs/sed.",cur.yr,".sed.mor.angle_v2.png"),
#     width=12*ppi, height=6*ppi, res=ppi)
# set.seed(pi); df %>% 
#   filter(!(zone1 %in% c("Above", "Inside", "Inside2", "Below") & year == 2025)) %>% 
#   filter(., type == "angle") %>% 
#   mutate(., value = as.numeric(value)) %>% 
#   ggplot(data = ., aes(y = value, x = zone1,
#                        fill=zone1))+
#   geom_boxplot(outlier.colour=NA)+
#   geom_jitter(aes(shape = zone1),
#               width = 0.3, height = 0.05,size = 3,stroke = 1,
#               show.legend = FALSE)+
#   scale_shape_manual(values = c(21:(length(unique(ang_cu$zone1))+21)))+
#   labs(title = paste0("Beach slopes recorded since 2008"),
#        subtitle = "2024 data excluded from non-Wash stations")+
#   theme(plot.title = element_text(face = "bold"),
#         legend.position = "bottom",#c(.5,1.02),
#         legend.direction = "horizontal",
#         legend.key.size = unit(1,"cm"),#change legend size!
#         legend.title=element_blank(),
#         legend.text=element_text(size = 12,
#                                  face = "bold"),
#         legend.box.margin=margin(-10,-10,-10,-10),
#         legend.margin=margin(0,0,0,0),
#         axis.text.y = element_text(size = 12),
#         axis.title.y = element_text(size = 14),
#         strip.text.x = element_text(size = 12),
#         strip.text = element_text(face="bold")
#   )+
#   xlab("") + ylab("Beach angle (°)")+
#   scale_fill_manual(values=cbPaletteFill[c(1:4,7)])+
#   scale_x_discrete(breaks=NULL)+
#   facet_wrap(~shore)+
#   guides(fill=guide_legend(nrow=1))
# dev.off()

df %>% 
  filter(!(zone1 %in% c("Above", "Inside", "Inside2", "Below") & year == cur.yr)) %>% 
  filter(., type == "angle") %>% 
  mutate(., value = as.numeric(value)) %>% 
  summarise(mn = mean(value),
            sd = sd(value))

## Get summary ####
## mean by zone & shore: current year
(mean_23 <- ang_cu %>%
   group_by(zone1, shore) %>%
   summarise("mean.angle"=mean(value), sd(value))) %>% 
  arrange(.,by = mean.angle)

## mean by zone: current year
(mean_23_zone <- ang_cu %>%
    group_by(zone1) %>%
    summarise("sdangle"=mean(value), sd(value)))

## mean by zone & shore: all years
(mean_all <- subset(df, type == "angle") %>%
    group_by(zone1, shore) %>%
    summarise("mean.angle"=mean(as.numeric(value), na.rm=TRUE),
              "sd.angle"=sd(as.numeric(value),na.rm=TRUE)))


### mean by station: current year
(mn_st_cur <- subset(df, type == "angle") %>%
    filter(.,year == cur.yr) %>% 
    group_by(transect, shore) %>% 
    summarise("mean.angle"=mean(as.numeric(value),na.rm=TRUE),
              "sd.angle"=sd(as.numeric(value),na.rm=TRUE),
              .groups = "drop") %>% 
    arrange(.,by=mean.angle))
head(mn_st_cur);tail(mn_st_cur)

print(mean_23); print(mean_all)

df_tm <- df
df_tm$shore <- factor(df_tm$shore, levels=c("Upper","Mid","Low"))
png(file = "figs/sed.ts.mor.angle.png",
    width=17*ppi, height=10*ppi, res=ppi)
df_tm %>% 
  # dplyr::filter(type == "angle") %>% 
  # dplyr::filter(year !=2024) %>%
  dplyr::filter(type == "angle",
                !(year == 2024 & zone1 != "Wash")) %>% #View()
  ggplot(.,
       aes(y = as.numeric(value), x = year, fill = zone1))+
  # geom_hline(yintercept = seq(0,12,by=2),colour="lightgrey",linetype=2)+
  geom_vline(xintercept = seq(2008,cur.yr,by=1),colour="lightgrey",linetype=2)+
  geom_boxplot(aes(group=year),outlier.shape = NA)+
  geom_jitter(width = 0.1, height = 0,alpha=0.3)+
  # geom_smooth(method = "loess", colour = "red", span = .9)+
  geom_smooth(method="gam")+
  facet_grid(shore~zone1)+
  scale_colour_manual(name = "", values=cbPalette)+
  scale_fill_manual(name = "", values=cbPaletteFill[c(1:4,8)])+
  scale_x_continuous(breaks = seq(2008, cur.yr, by = 1))+
  xlab("Year") + ylab("Angle")+
  theme(plot.title = element_text(face = "bold"),
        legend.position="none",
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(angle = 270,vjust=.25, face = 2, size=12),
        strip.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        strip.text = element_text(face="bold"),
        axis.title = element_text(face="bold")
  )+
  labs(
    title = paste0("Beach slopes recorded since 2008 as part of the SFGPBM programme"),
    subtitle = "Non-Wash data from 2024 contained methodological errors and are excluded.\nBlue lines indicate generalised additive model.",
        )
  dev.off()

  
## Statistical comparisons ####
tic("statistical comparison: Angle")
ang_cu_lnc <- droplevels(ang_cu %>% filter(., zone1 != "Wash"))

##re-level to compare all with Inside
ang_cu_lnc$zone1 <- factor(ang_cu_lnc$zone1,levels=c("Inside","Above","Inside2","Below"))

ang_cur_mod1 <- lmerTest::lmer(value ~ zone1 + (1|shore) , data = ang_cu_lnc,REML=FALSE)
ang_cur_mod2 <- lmerTest::lmer(value ~ zone1 + (transect|shore) , data = ang_cu_lnc,REML=FALSE)
anova(ang_cur_mod3 <- lmerTest::lmer(value ~ zone1 + (1|transect)+(1|shore) , data = ang_cu_lnc,REML=FALSE))
anova(ang_cur_mod3.1 <- lmerTest::lmer(value ~ zone1 + (shore|transect) , data = ang_cu_lnc,REML=TRUE))
anova(ang_cur_mod4 <- lmerTest::lmer(value ~ zone1*shore + (1|transect) , data = ang_cu_lnc,REML=FALSE))

# model checking: which is best? Lowest AIC = 'better' fit
anova(ang_cur_mod1, ang_cur_mod2, ang_cur_mod3,ang_cur_mod3.1, ang_cur_mod4)###compare model fits
AIC(ang_cur_mod1,ang_cur_mod2,ang_cur_mod3,ang_cur_mod3.1,ang_cur_mod4) %>% 
  as.data.frame() %>%
  arrange(AIC) %>% print()
# ang_cur_mod1 has lowest AIC
plot(performance::compare_performance(ang_cur_mod1,
                                      ang_cur_mod2,
                                      ang_cur_mod3,
                                      ang_cur_mod3.1,
                                      ang_cur_mod4,
                                      rank = TRUE))
# visreg::visreg(ang_cur_mod1)
visreg::visreg(ang_cur_mod3.1)

# go with model ang_cur_mod1
anova(ang_cur_mod1)
summary(ang_cur_mod1)

lmerTest::ls_means(ang_cur_mod1, test.effs = "Group",pairwise = TRUE)
lmerTest::ls_means(ang_cur_mod3.1, test.effs = "Group",pairwise = TRUE)
sjPlot::plot_model(ang_cur_mod1,show.values=TRUE, show.p=TRUE)
visreg::visreg(ang_cur_mod1)
performance::check_predictions(ang_cur_mod3.1)
toc(log = TRUE)

### WASH ONLY: compare transects ####
ang_cu %>% filter(., zone1=="Wash") ->ang_cu_wa
summary(m_ang_wa <- lm(value ~ shore, data=ang_cu_wa))

## time series ####
png(file = "figs/sed.ang.ts.png",
    width=12*ppi, height=6*ppi, res=ppi)
df %>% 
  # filter(.,type=="angle") %>% 
  dplyr::filter(type == "angle",
                !(year == 2024 & zone1 != "Wash")) %>% #View()
  #filter(., zone1 != "Wash") %>%
  droplevels(.) %>% 
  mutate(value = as.numeric(value)) %>% 
  filter(., value >= 0) %>% 
  dplyr::select(.,c(year,shore,type,value,zone1)) %>% 
  ggplot(data=., aes(x=year, y=value, fill=zone1))+
  # ylim(0,12)+
  geom_vline(xintercept = seq(2008,cur.yr,by=1),linetype=2, colour="lightgrey")+
  # geom_hline(yintercept = seq(0,12,by=2),linetype=2, colour="lightgrey")+
  geom_boxplot(outlier.shape = NA,
               aes(group=as.factor(year)), show.legend = FALSE)+
  geom_jitter(show.legend = FALSE, alpha=0.25)+
  facet_grid(shore~zone1)+
  scale_fill_manual(values=cbPaletteFill[c(1:4,7)])+
  scale_y_continuous(breaks = seq(0,12,by = 4))+
  geom_smooth(
    show.legend = FALSE, col="red",
    method = "gam"
    # method = "loess",span=0.9,
  )+
  labs(
    title="Beach slope values recorded since 2008",
    subtitle = "Non-Wash data from 2024 contained methodological errors and are excluded.\nBlue lines indicate generalised additive model.",
    y="Beach angle (°)"
    )+
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, face = 2),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        strip.text = element_text(face="bold"),
        strip.background = element_rect(color = "black",fill = "grey95", size = 1),
        )
dev.off()

###=============###
# Compaction ####
###=============###
com_cu <- cur_dat %>% filter(., type == "cone") %>% 
  mutate(value = as.numeric(value))

(com_cu_sum <- com_cu %>% 
    group_by(.,zone1) %>% 
    summarise(mean_cu=mean(value,na.rm=TRUE),sd.cu=sd(value,na.rm=TRUE)) %>% 
    ungroup())
com_cu_st <- com_cu %>% 
  filter(.,year==cur.yr) %>% 
  group_by(transect, shore) %>% 
  summarise("mean.cone"=mean(as.numeric(value),na.rm=TRUE),
            "sd.cone"=sd(as.numeric(value),na.rm=TRUE),
            .groups = "drop")
com_cu_st[which.min(com_cu_st$mean.cone),]
com_cu_st[which.max(com_cu_st$mean.cone),]

### boxplot ###
png(file = paste0("figs/sed.",cur.yr,".sed.mor.pen.png"),
    width=12*ppi, height=6*ppi, res=ppi)
set.seed(pi); com_cu %>% 
  # filter(., zone1 != "Wash") %>%
  ggplot(data = ., aes(y = value, x = zone1,
                       fill=zone1))+
  geom_boxplot(outlier.colour=NA)+
  geom_jitter(aes(shape = zone1),
              width = 0.3, height = 0.05,size = 3,stroke = 1,
              show.legend = FALSE)+
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom",#c(.5,1.02),
        legend.direction = "horizontal",
        # legend.box = "horizontal",
        legend.key.size = unit(1,"cm"),#change legend size!
        legend.title=element_blank(),
        legend.text=element_text(size = 12,
                                 face = "bold"),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.margin=margin(0,0,0,0),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 12),
        strip.text = element_text(face="bold"),
        strip.background = element_rect(color = "black",fill = "grey95", size = 1),
  )+
  xlab("") + ylab("Cone index")+
  scale_fill_manual(values=cbPaletteFill[c(1:4,8)])+
  scale_x_discrete(breaks=NULL)+
  scale_shape_manual(values = c(21:(length(unique(com_cu$zone1))+21)))+
  facet_wrap(~shore)+
  guides(fill=guide_legend(nrow=1))+
  labs(title = paste0("Sediment compaction recorded in the ",cur.yr," monitoring programme"),
       subtitle = "Higher values indicate more compacted sediments")
dev.off()

## Get summary ####
#require(dplyr)
(mean_cur <- com_cu %>%
   group_by(zone1, shore) %>%
   summarise(mean(value,na.rm=TRUE), sd(value,na.rm=TRUE)))

(mean_cur_zone <- com_cu %>%
    group_by(zone1) %>%
    summarise(mean(value,na.rm=TRUE), sd(value,na.rm=TRUE)))

(mean_all <- subset(df, type == "cone") %>%
    group_by(zone1, shore) %>%
    summarise(mean(as.numeric(value),na.rm=TRUE),
              sd(as.numeric(value),na.rm=TRUE)))

(mean_cur_tran <- com_cu %>%
    group_by(transect,shore) %>%
    summarise(mean(value,na.rm=TRUE), sd(value,na.rm=TRUE)))

### mean by station
mn_st_cur <- subset(df, type == "cone") %>%
  group_by(transect, shore) %>% 
  summarise(mean(as.numeric(value),na.rm=TRUE),
            sd(as.numeric(value),na.rm=TRUE),
            .groups = "drop")

## Time series ####
df_tm <- df
df_tm$shore <- factor(df_tm$shore, levels=c("Upper","Mid","Low"))

png(file = "figs/sed.ts.mor.pen.png",
    width=12*ppi, height=6*ppi, res=ppi)
df %>% 
  filter(.,type=="cone") %>% 
  # filter(., zone1 != "Wash") %>%
  droplevels(.) %>% 
  mutate(value = as.numeric(value)) %>% 
  filter(., value >= 0) %>% 
  dplyr::select(.,c(year,shore,type,value,zone1)) %>% 
  ggplot(.,
         aes(y = as.numeric(value), x = year, fill = zone1))+
  geom_vline(xintercept = seq(2008,cur.yr,by=1),linetype=2, colour="lightgrey")+
  geom_boxplot(aes(group=year),outlier.shape = NA)+
  geom_jitter(width = 0.1, height = 0,alpha=0.3)+
  # geom_smooth(method = "loess", colour = "red", span = .9)+
  geom_smooth(method = "gam", colour = "red", #span = .9
  )+
  facet_grid(shore~zone1)+
  scale_colour_manual(name = "", values=cbPalette)+
  scale_fill_manual(name = "", values=cbPaletteFill)+
  # scale_x_continuous(breaks = seq(2008, cur.yr, by = 4))+
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0))+
  xlab("Year") + ylab("Cone index")+
  labs(title = paste0("Sediment compaction recorded since 2008 as part of the SFGPBM programme"),
       subtitle = "Higher values indicate more compacted sediments.\nRed lines indicate generalised additive model trend")+
  theme(
    legend.position="none",
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12, face = 2),
    axis.title.y = element_text(size = 14),
    strip.text.x = element_text(size = 14),
    strip.text.y = element_text(size = 14),
    strip.text = element_text(face="bold"),
    strip.background = element_rect(color = "black",fill = "grey95", size = 1),
    )
dev.off()

## Statistical comparison ####
## remove Wash data
com_cu_lnc <- droplevels(com_cu %>% filter(.,zone1!="Wash"))
##re-level to compare all with Inside
com_cu_lnc$zone1 <- factor(com_cu_lnc$zone1,levels=c("Inside","Above","Inside2","Below"))

com_cur_mod1 <- lmerTest::lmer(value ~ zone1 + (1|shore) , data = com_cu_lnc,REML=FALSE)
com_cur_mod2 <- lmerTest::lmer(value ~ zone1 + (transect|shore) , data = com_cu_lnc,REML=FALSE)
anova(com_cur_mod3 <- lmerTest::lmer(value ~ zone1 + (1|transect) + (1|shore) , data = com_cu_lnc,REML=FALSE))
anova(com_cur_mod3.1 <- lmerTest::lmer(value ~ zone1 + (shore|transect) , data = com_cu_lnc,REML=TRUE))
anova(com_cur_mod3.2 <- lmerTest::lmer(value ~ zone1*shore + (shore|transect) , data = com_cu_lnc,REML=TRUE))
anova(com_cur_mod4 <- lmerTest::lmer(value ~ zone1*shore + (1|transect), data = com_cu_lnc,REML=FALSE))

anova(com_cur_mod1,com_cur_mod2,com_cur_mod3,com_cur_mod3.1,com_cur_mod4)
AIC(com_cur_mod1,com_cur_mod2,com_cur_mod3,com_cur_mod3.1,com_cur_mod4) %>% 
  as.data.frame() %>%
  arrange(AIC) %>% print()
### Lowest AIC = com_cur_mod3.1

summary(com_cur_mod3.1)
print(com_cur_mod3.1,cor=FALSE)
anova(com_cur_mod3.1)
lmerTest::ls_means(com_cur_mod3.1, test.effs = "Group",pairwise = TRUE)
sjPlot::plot_model(com_cur_mod3.1,show.values=TRUE, show.p=TRUE)
visreg::visreg(com_cur_mod3.1)
performance::check_posterior_predictions(com_cur_mod3.1)
ggplot(data=com_cu, aes(x=zone1,y=value,colour=transect))+
  # geom_point()+facet_wrap(.~shore)+
  geom_jitter(width = 0.2)+
  facet_wrap(.~shore)

## time series ####
# png(file = "output/figs/sed.cone.ts.png",
#     width=12*ppi, height=6*ppi, res=ppi)
# df %>% 
#   filter(.,type=="cone") %>% 
#   mutate(value = as.numeric(value)) %>% 
#   dplyr::select(.,c(year,shore,type,value,zone1)) %>% 
#   ggplot(data=., aes(x=year, y=value, fill=zone1))+
#   geom_boxplot(outlier.shape = NA,
#                aes(group=as.factor(year)), show.legend = FALSE)+
#   geom_jitter(show.legend = FALSE, alpha=0.25)+
#   facet_grid(shore~zone1)+
#   scale_fill_manual(values=cbPaletteFill[c(1:4,8)])+
#   geom_smooth(span=0.9, show.legend = FALSE, col="red",
#               # method = "gam"
#               method = "loess"
#   )+
#   ylab("Cone index")+
#   theme(axis.title.x = element_blank(),
#         axis.text.y = element_text(size = 12),
#         axis.title.y = element_text(size = 14),
#         strip.text.x = element_text(size = 12),
#         strip.text.y = element_text(size = 12),
#         strip.text = element_text(face="bold"))
# dev.off()

# Wave Class #### not interesting!
# df %>% filter(.,type == "WaveClass") -> df_wave
# plot(df_wave$value)

# create summary table for appendices ####

glimpse(cur_dat)

cur_dat %>% 
  filter(.,type=="cone") %>% 
  mutate(value=as.numeric(value)) %>% 
  group_by(.,transect,shore,year) %>% 
  summarise(mean = mean(value, na.rm=TRUE),
            sd = sd(value, na.rm=TRUE)) %>% 
  View(.)

cur_dat %>% 
  filter(.,type=="angle") %>% 
  mutate(value=as.numeric(value)) %>% 
  group_by(.,transect,shore,year) %>% 
  summarise(mean = mean(value, na.rm=TRUE),
            sd = sd(value, na.rm=TRUE)) %>% 
  View(.)

# smry$value <- as.numeric(smry$value)

# smry <- smry %>% 
#   select(.,transect,shore,type,value,zone1) %>%
#   group_by(transect,shore,type,zone1) %>% 
#   summarise(mean=mean(value,na.rm=TRUE),sd=sd(value,na.rm=TRUE)) %>% ungroup()
# 
# smry$ID <- paste0(smry$transect,".",smry$shore)
# 
# smry.cone <- smry %>% 
#   filter(.,type=="cone") %>% ungroup() %>% 
#   rename(., "Mean cone" = mean, "SD cone"=sd)
# smry.ang <- smry %>% 
#   filter(.,type=="angle") %>% 
#   select(.,-(transect:zone1)) %>% ungroup() %>% 
#   rename(., "Mean angle" = mean, "SD angle"=sd)
# 
# oo <- left_join(smry.cone,smry.ang,by="ID")
# write.csv(oo,file="output/forDoc/2021.sed.bulk.csv",row.names = FALSE)

# WIP time series stats WIP ####
## compaction ####
### fix formatting ####
df_tm_com <- df_tm %>% 
  dplyr::filter(type=="cone",
                zone1 != "Wash") %>% 
  dplyr::mutate(value = as.numeric(value))

com_fit <- mgcv::gam(
  # value ~ zone1 + s(year, by = zone1,bs="cr"),
  value ~ s(year, by = zone1,bs="cr"),
  data = df_tm_com,
  method = "REML")
summary(com_fit)  
gratia::draw(com_fit)
# com.pl+  labs(
#     title = "SDSD"
#   )

df_tm_com <- df_tm_com %>%
  mutate(zone_shore = interaction(zone1, shore, drop = TRUE))

com_fit_B <- mgcv::gam(
  value ~ zone1 + shore + s(year, by = zone_shore, bs = "cr"),
  data = df_tm_com,
  method = "REML"
  )

p1 <- draw(com_fit_B,select="s(year):zone_shoreAbove.Upper")
p2 <- draw(com_fit_B,select="s(year):zone_shoreInside.Upper")
p3 <- draw(com_fit_B,select="s(year):zone_shoreInside2.Upper")
p4 <- draw(com_fit_B,select="s(year):zone_shoreBelow.Upper")
p5 <- draw(com_fit_B,select="s(year):zone_shoreAbove.Mid")
p6 <- draw(com_fit_B,select="s(year):zone_shoreInside.Mid")
p7 <- draw(com_fit_B,select="s(year):zone_shoreInside2.Mid")
p8 <- draw(com_fit_B,select="s(year):zone_shoreBelow.Mid")

p1+p2+p3+p4+p5+p6+p7+p8+patchwork::plot_layout(ncol = 4,nrow=2)

# Extract smooth estimates for the zone×shore smooths of year
sm <- gratia::smooth_estimates(com_fit_B) %>%
  # keep only s(year):... smooths
  dplyr::filter(stringr::str_detect(.data$.smooth, "^s\\(year\\):")) %>%
  gratia::add_confint()

# Ensure we have clean zone & shore labels.
# gratia often provides a `by` column (here likely `zone_shore`) OR individual factor columns.
# We'll create robust columns from `.smooth` if needed.
sm <- sm %>%
  mutate(
    # If there is already a zone_shore column (from model matrix), keep it;
    # otherwise derive from `.smooth`
    zone_shore = dplyr::coalesce(.data$zone_shore, sub("^s\\(year\\):", "", .data$.smooth)),
    # Split into zone1 and shore. `interaction(zone1, shore)` defaults to "zone1.shore"
    # Adjust the sep if your interaction used a different separator.
    .group = .data$zone_shore
  ) %>%
  tidyr::separate(.group, into = c("zone1", "shore"), sep = "\\.", remove = FALSE) %>% 
  dplyr::mutate( zone1 = factor(zone1,levels=c("Above","Inside","Inside2","Below")),
                 shore = factor(shore,levels=c("Upper","Mid","Low")))

# Check column names to confirm x/estimate naming varies by gratia version
# names(sm)
# Typical columns: .smooth, x or year, .estimate, .se, .lower_ci, .upper_ci, zone1, shore, zone_shore
# png(file = "figs/sed.ts.mor.pen.gam.trend.png",
#     width=12*ppi, height=6*ppi, res=ppi)
ggplot(sm, aes(x = year, y = .estimate)) +
  geom_vline(xintercept = c(2008:cur.yr),linetype = 2,linewidth = 0.4,col="lightgrey")+
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(aes(ymin = .lower_ci,
                  ymax = .upper_ci,
                  fill = zone1),
              alpha = 0.25, colour = NA) +
  geom_line(aes(colour = zone1), linewidth = 0.9) +
  # Facet by zone and shore for a grid view (readable across management areas and shore levels)
  facet_grid(shore ~ zone1, scales = "free_y") +
  scale_fill_manual(values = cbPaletteFill, guide = "none") +
  # If you want lines to match ribbons exactly, use scale_colour_manual with the same vector
  scale_colour_manual(values = cbPalette, name = "Zone") +
  labs(
    title = "Sediment compaction trends since 2008",
    # subtitle = "GAM smooths of s(Year) by zone × shore\nRibbons show simultaneous confidence intervals",
    x = "Year",
    y = "Smooth effect s(Year)"
  ) +
  ggthemes::theme_few() +
  theme(
    legend.position = "none",
    plot.title = element_text(face=2,size=18),
    plot.subtitle = element_text(face=2,size=12),
    axis.title.y = element_text(face=2),
    axis.title.x=element_blank(),
    axis.text.y = element_text(face=2),
    axis.text.x = element_text(face=2,size = 12),
    strip.text = element_text(face=2,size=14),
    strip.background = element_rect(color = "black",fill = "grey95", size = 1),
    ) ->com.pl.A
# dev.off()

## get derivatives ##
# Get exact smooth labels from the model
terms_year_by <- gratia::smooths(com_fit_B) %>%
  stringr::str_subset("^s\\(year\\):")

# First derivative
deriv1 <- gratia::derivatives(
  com_fit_B,
  term = terms_year_by,
  order = 1,
  interval = "simultaneous",
  n = 200,
  type="central",
  unconditional = TRUE
) %>%
  rename(year = year, deriv = .derivative) %>%
  mutate(zone_shore = sub("^s\\(year\\):", "", .smooth)) %>%
  separate(zone_shore, into = c("zone1", "shore"), sep = "\\.", remove = FALSE) %>% 
  mutate(zone = gsub("^.{0,10}", "", zone1)) %>% 
  mutate(
    zone = factor(zone,levels=c("Above","Inside","Inside2","Below")),
    shore = factor(shore,levels=c("Upper","Mid","Low")),
    )

# png(file = "figs/sed.ts.mor.pen.1stderiv.png",
#     width=12*ppi, height=6*ppi, res=ppi)
deriv1 %>% 
  mutate(change = ifelse(.lower_ci < 0 & .upper_ci > 0, NA_real_, deriv)) %>% 
  mutate(zone1 = ifelse(grepl("Above",.smooth),"Above",
                        ifelse(grepl("Inside2",.smooth),"Inside2",
                               ifelse(grepl("Below",.smooth),"Below","Inside")))) %>% 
  mutate(zone1 = factor(zone1, levels = c("Above","Inside","Inside2","Below"))) %>% 
  ggplot(., aes(x = year, y = deriv)) +
  geom_vline(xintercept = c(2008:cur.yr),linetype = 2,linewidth = 0.4,col="lightgrey")+
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = zone),
              alpha = 0.25, colour = NA) +
  geom_line(aes(colour = zone), linewidth = 0.8,show.legend = FALSE) +
  facet_grid(shore ~ zone, scales = "free_y") +
  scale_fill_manual(values = cbPaletteFill, guide = "none") +
  geom_line(aes(x=year, y= change), col=2,lwd=2)+
  scale_colour_manual(values = cbPalette[1:4]) +
  labs(
    # title = "Sediment compaction",
    subtitle="First derivative of s(Year) by zone × shore\nDerivative > 0 indicates increasing trend; < 0 decreasing",
    y = "ds/dYear"
    ) +
  ggthemes::theme_few() +
  theme(
    plot.title = element_text(face=2,size=18),
    plot.subtitle = element_text(face=2,size=12),
    axis.title.y = element_text(face=2),
    axis.title.x=element_blank(),
    axis.text.y = element_text(face=2),
    axis.text.x = element_text(face=2,size = 12),
    strip.text = element_text(face=2,size=14),
    strip.background = element_rect(color = "black",fill = "grey95", size = 1),
    ) -> com.pl.B
# dev.off()

png(file = "figs/sed.ts.mor.pen_AB.png",
        width=12*ppi, height=10*ppi, res=ppi)
com.pl.A/com.pl.B+plot_layout(guides='collect')+plot_annotation(tag_levels = 'A')
dev.off()

deriv2 <- gratia::derivatives(
  com_fit_B,
  term = terms_year_by,
  order = 2,
  interval = "simultaneous",
  n = 200,
  unconditional = TRUE
) %>%
  rename(year = year, deriv = .derivative) %>%
  mutate(zone_shore = sub("^s\\(year\\):", "", .smooth)) %>%
  separate(zone_shore, into = c("zone1", "shore"), sep = "\\.", remove = FALSE) %>% 
  mutate(zone = gsub("^.{0,10}", "", zone1)) %>% 
  mutate(
    zone = factor(zone,levels=c("Above","Inside","Inside2","Below")),
    shore = factor(shore,levels=c("Upper","Mid","Low")),
  )

# png(file = "figs/sed.ts.mor.pen.2ndderiv.png",
    # width=12*ppi, height=6*ppi, res=ppi)
deriv2 %>% 
  mutate(change = ifelse(.lower_ci < 0 & .upper_ci > 0, NA_real_, deriv)) %>% 
  mutate(zone1 = ifelse(grepl("Above",.smooth),"Above",
                        ifelse(grepl("Inside2",.smooth),"Inside2",
                               ifelse(grepl("Below",.smooth),"Below","Inside")))) %>% 
  mutate(zone1 = factor(zone1, levels = c("Above","Inside","Inside2","Below"))) %>% 
  ggplot(., aes(x = year, y = deriv)) +
  geom_vline(xintercept = c(2008:cur.yr),linetype = 2,linewidth = 0.4,col="lightgrey")+
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = zone),
              alpha = 0.25, colour = NA) +
  geom_line(aes(colour = zone), linewidth = 0.8,show.legend = FALSE) +
  facet_grid(shore ~ zone, scales = "free_y") +
  scale_fill_manual(values = cbPaletteFill, guide = "none") +
  geom_line(aes(x=year, y= change), col=2,lwd=2)+
  scale_colour_manual(values = cbPalette[1:4]) +
  labs(
    # title = "Sediment compaction",
    subtitle="Second derivative of s(Year) by zone × shore\nCurvature: positive = accelerating increase / easing decline; negative = turning down",
    y = "ds/dYear"
  ) +
  ggthemes::theme_few() +
  theme(
    plot.title = element_text(face=2,size=18),
    plot.subtitle = element_text(face=2,size=12),
    axis.title.y = element_text(face=2),
    axis.title.x=element_blank(),
    axis.text.y = element_text(face=2),
    axis.text.x = element_text(face=2,size = 12),
    strip.text = element_text(face=2,size=14),
    strip.background = element_rect(color = "black",fill = "grey95", size = 1),
  )->com.pl.C
# dev.off()

com.pl.B/com.pl.C+plot_layout(guides='collect')+plot_annotation(tag_levels = 'A')

###############
## angle ####
### fix formatting ####
df_tm_ang <- df_tm %>% 
  dplyr::filter(type=="angle",
                zone1 != "Wash") %>% 
  dplyr::mutate(value = as.numeric(value))

ang_fit <- mgcv::gam(
  # value ~ zone1 + s(year, by = zone1,bs="cr"),
  value ~ s(year, by = zone1,bs="cr"),
  data = df_tm_ang,
  method = "REML")
summary(ang_fit)  
gratia::draw(ang_fit)

df_tm_ang <- df_tm_ang %>%
  mutate(zone_shore = interaction(zone1, shore, drop = TRUE))

ang_fit_B <- mgcv::gam(
  value ~ zone1 + shore + s(year, by = zone_shore, bs = "cr"),
  data = df_tm_ang,
  method = "REML"
)

p1 <- draw(ang_fit_B,select="s(year):zone_shoreAbove.Upper")
p2 <- draw(ang_fit_B,select="s(year):zone_shoreInside.Upper")
p3 <- draw(ang_fit_B,select="s(year):zone_shoreInside2.Upper")
p4 <- draw(ang_fit_B,select="s(year):zone_shoreBelow.Upper")
p5 <- draw(ang_fit_B,select="s(year):zone_shoreAbove.Mid")
p6 <- draw(ang_fit_B,select="s(year):zone_shoreInside.Mid")
p7 <- draw(ang_fit_B,select="s(year):zone_shoreInside2.Mid")
p8 <- draw(ang_fit_B,select="s(year):zone_shoreBelow.Mid")

p1+p2+p3+p4+p5+p6+p7+p8+patchwork::plot_layout(ncol = 4,nrow=2)

# Extract smooth estimates for the zone×shore smooths of year
sm <- gratia::smooth_estimates(ang_fit_B) %>%
  # keep only s(year):... smooths
  dplyr::filter(stringr::str_detect(.data$.smooth, "^s\\(year\\):")) %>%
  gratia::add_confint()

# Ensure we have clean zone & shore labels.
# gratia often provides a `by` column (here likely `zone_shore`) OR individual factor columns.
# We'll create robust columns from `.smooth` if needed.
sm <- sm %>%
  mutate(
    # If there is already a zone_shore column (from model matrix), keep it;
    # otherwise derive from `.smooth`
    zone_shore = dplyr::coalesce(.data$zone_shore, sub("^s\\(year\\):", "", .data$.smooth)),
    # Split into zone1 and shore. `interaction(zone1, shore)` defaults to "zone1.shore"
    # Adjust the sep if your interaction used a different separator.
    .group = .data$zone_shore
  ) %>%
  tidyr::separate(.group, into = c("zone1", "shore"), sep = "\\.", remove = FALSE) %>% 
  dplyr::mutate( zone1 = factor(zone1,levels=c("Above","Inside","Inside2","Below")),
                 shore = factor(shore,levels=c("Upper","Mid","Low")))

# Check column names to confirm x/estimate naming varies by gratia version
# names(sm)
# Typical columns: .smooth, x or year, .estimate, .se, .lower_ci, .upper_ci, zone1, shore, zone_shore
# png(file = "figs/sed.ts.mor.ang.gam.trend.png",
    # width=12*ppi, height=6*ppi, res=ppi)
ggplot(sm, aes(x = year, y = .estimate)) +
  geom_vline(xintercept = c(2008:cur.yr),linetype = 2,linewidth = 0.4,col="lightgrey")+
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(aes(ymin = .lower_ci,
                  ymax = .upper_ci,
                  fill = zone1),
              alpha = 0.25, colour = NA) +
  geom_line(aes(colour = zone1), linewidth = 0.9) +
  # Facet by zone and shore for a grid view (readable across management areas and shore levels)
  facet_grid(shore ~ zone1, scales = "free_y") +
  scale_fill_manual(values = cbPaletteFill, guide = "none") +
  # If you want lines to match ribbons exactly, use scale_colour_manual with the same vector
  scale_colour_manual(values = cbPalette, name = "Zone") +
  labs(
    title = "Beach slope trends since 2008",
    # subtitle = "GAM smooths of s(Year) by zone × shore\nRibbons show simultaneous confidence intervals",
    x = "Year",
    y = "Smooth effect s(Year)"
  ) +
  ggthemes::theme_few() +
  theme(
    legend.position = "none",
    plot.title = element_text(face=2,size=18),
    plot.subtitle = element_text(face=2,size=12),
    axis.title.y = element_text(face=2),
    axis.title.x=element_blank(),
    axis.text.y = element_text(face=2),
    axis.text.x = element_text(face=2,size = 12),
    strip.text = element_text(face=2,size=14),
    strip.background = element_rect(color = "black",fill = "grey95", size = 1),
  ) -> ang.pl.A
# dev.off()

## get derivatives ##
# Get exact smooth labels from the model
terms_year_by <- gratia::smooths(ang_fit_B) %>%
  stringr::str_subset("^s\\(year\\):")

# First derivative
deriv1 <- gratia::derivatives(
  ang_fit_B,
  term = terms_year_by,
  order = 1,
  interval = "simultaneous",
  n = 200,
  unconditional = TRUE
  ) %>%
  rename(year = year, deriv = .derivative) %>%
  mutate(zone_shore = sub("^s\\(year\\):", "", .smooth)) %>%
  separate(zone_shore, into = c("zone1", "shore"), sep = "\\.", remove = FALSE) %>% 
  mutate(zone = gsub("^.{0,10}", "", zone1)) %>% 
  mutate(
    zone = factor(zone,levels=c("Above","Inside","Inside2","Below")),
    shore = factor(shore,levels=c("Upper","Mid","Low")),
  )

# png(file = "figs/sed.ts.mor.ang.1stderiv.png",
#     width=12*ppi, height=6*ppi, res=ppi)
ggplot(deriv1, aes(x = year, y = deriv)) +
  geom_vline(xintercept = c(2008:cur.yr),linetype = 2,linewidth = 0.4,col="lightgrey")+
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = zone),
              alpha = 0.25, colour = NA) +
  geom_line(aes(colour = zone), linewidth = 0.8,show.legend = FALSE) +
  facet_grid(shore ~ zone, scales = "free_y") +
  scale_fill_manual(values = cbPaletteFill, guide = "none") +
  scale_colour_manual(values = cbPalette[1:4]) +
  labs(
    # title = "Beach slope",
    # subtitle="First derivative of s(Year) by zone × shore\nDerivative > 0 indicates increasing trend; < 0 decreasing",
    y = "ds/dYear"
  ) +
  ggthemes::theme_few() +
  theme(
    plot.title = element_text(face=2,size=18),
    plot.subtitle = element_text(face=2,size=12),
    axis.title.y = element_text(face=2),
    axis.title.x=element_blank(),
    axis.text.y = element_text(face=2),
    axis.text.x = element_text(face=2,size = 12),
    strip.text = element_text(face=2,size=14),
    strip.background = element_rect(color = "black",fill = "grey95", size = 1),
  )->ang.pl.B
# dev.off()

deriv2 <- gratia::derivatives(
  ang_fit_B,
  term = terms_year_by,
  order = 2,
  interval = "simultaneous",
  n = 200,
  unconditional = TRUE
) %>%
  rename(year = year, deriv = .derivative) %>%
  mutate(zone_shore = sub("^s\\(year\\):", "", .smooth)) %>%
  separate(zone_shore, into = c("zone1", "shore"), sep = "\\.", remove = FALSE) %>% 
  mutate(zone = gsub("^.{0,10}", "", zone1)) %>% 
  mutate(
    zone = factor(zone,levels=c("Above","Inside","Inside2","Below")),
    shore = factor(shore,levels=c("Upper","Mid","Low")),
  )

# png(file = "figs/sed.ts.mor.ang.2ndderiv.png",
#     width=12*ppi, height=6*ppi, res=ppi)
ggplot(deriv2, aes(x = year, y = deriv)) +
  geom_vline(xintercept = c(2008:cur.yr),linetype = 2,linewidth = 0.4,col="lightgrey")+
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = zone),
              alpha = 0.25, colour = NA) +
  geom_line(aes(colour = zone), linewidth = 0.8,show.legend = FALSE) +
  facet_grid(shore ~ zone, scales = "free_y") +
  scale_fill_manual(values = cbPaletteFill, guide = "none") +
  scale_colour_manual(values = cbPalette[1:4]) +
  labs(
    # title = "Sediment compaction",
    # subtitle="Second derivative of s(Year) by zone × shore\nCurvature: positive = accelerating increase / easing decline; negative = turning down",
    y = "ds/dYear"
  ) +
  ggthemes::theme_few() +
  theme(
    plot.title = element_text(face=2,size=18),
    plot.subtitle = element_text(face=2,size=12),
    axis.title.y = element_text(face=2),
    axis.title.x=element_blank(),
    axis.text.y = element_text(face=2),
    axis.text.x = element_text(face=2,size = 12),
    strip.text = element_text(face=2,size=14),
    strip.background = element_rect(color = "black",fill = "grey95", size = 1),
  )
# dev.off()


png(file = "figs/sed.ts.mor.ang_AB.png",
    width=12*ppi, height=10*ppi, res=ppi)
ang.pl.A/ang.pl.B+plot_layout(guides='collect')+plot_annotation(tag_levels = 'A')
dev.off()























############################

# # Extract smooth estimates for the year-by-zone1 smooths and add CIs
# gratia::smooth_estimates(com_fit)  %>% #names()
#   # keep only the s(year):zone1* smooths
#   dplyr::filter(str_detect(.data$.smooth, "^s\\(year\\):")) %>% 
#   gratia::add_confint() %>% #View()
# 
# # Inspect columns: typically `x` (the covariate), `est`, `se`, `lower_ci`, `upper_ci`, `by`
# # Some versions name x as the covariate (e.g., `year`) — if so, adapt aes(x = year)
# ggplot(., aes(x = .data$year, y = .data$.estimate)) +
#   geom_hline(yintercept = 0, lty=2)+
#   geom_ribbon(aes(ymin = .data$.lower_ci,
#                   ymax = .data$.upper_ci,
#                   fill = .data$zone1),
#               alpha = 0.25, colour = NA) +
#   geom_line(aes(colour = .data$zone1), linewidth = 0.8) +
#   facet_wrap(~ zone1, scales = "free_y") +
#   scale_fill_manual(name = "", values=cbPaletteFill)+
#   # scale_fill_brewer(type = "qual", palette = "Set2") +
#   scale_colour_brewer(type = "qual", palette = "Set2") +
#   labs(
#     title="Generalised additive model trends of sediment compaction since 2008",
#     #x = "Year",
#     y = "Smooth effect s(Year)",
#     colour = "Zone",
#     fill   = "Zone"
#     ) +
#   ggthemes::theme_few() +
#   theme(
#     legend.position = "none",
#     axis.text = element_text(face=2),
#     axis.title.y = element_text(face=2),
#     axis.title.x = element_blank(),
#     strip.text = element_text(face=2),
#     )
# 
# #### derivatives ####
# gratia::smooths(com_fit)
# 
# # Explicit vector of smooth terms (from gratia::smooths(com_fit))
# terms_year_by <- c(
#   "s(year):zone1Above",
#   "s(year):zone1Inside",
#   "s(year):zone1Inside2",
#   "s(year):zone1Below"
# )
# 
# # --- First derivative (velocity of change) ---
# 
# # --- First derivative ---
# deriv1 <- gratia::derivatives(
#   com_fit,
#   term = terms_year_by,      # use exact names
#   order = 1,                 # first derivative
#   interval = "simultaneous", # or "confidence" for pointwise CI
#   n = 200,                   # grid density over year
#   unconditional = TRUE       # include parametric uncertainty
# ) %>% 
#   # .x is the covariate value; .derivative is the derivative
#   rename(year = year, deriv = .derivative) %>%
#   # Create a clean zone1 label from the smooth name if not already present
#   mutate(
#     zone1 = if ("zone1" %in% names(.)) .data$zone1 else sub("^s\\(year\\):", "", .data$.smooth),
#     order = 1
#   )
# 
# # --- Second derivative (acceleration / curvature) ---
# deriv2 <- gratia::derivatives(
#   com_fit,
#   term = terms_year_by,      # use exact names
#   order = 2,                 # second derivative
#   interval = "simultaneous", # or "confidence" for pointwise CI
#   n = 200,                   # grid density over year
#   unconditional = TRUE       # include parametric uncertainty
# ) %>% 
#   # .x is the covariate value; .derivative is the derivative
#   rename(year = year, deriv = .derivative) %>%
#   # Create a clean zone1 label from the smooth name if not already present
#   mutate(
#     zone1 = if ("zone1" %in% names(.)) .data$zone1 else sub("^s\\(year\\):", "", .data$.smooth),
#     order = 2
#   )
# 
# # Combine for easier plotting (optional)
# derivs <- bind_rows(deriv1, deriv2)
# 
# 
# ggplot(deriv1, aes(x = year, y = deriv)) +
#   geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
#   geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = zone1),
#               alpha = 0.25, colour = NA) +
#   geom_line(aes(colour = zone1), linewidth = 0.8) +
#   facet_wrap(~ zone1, scales = "free_y") +
#   scale_fill_manual(name = "", values = cbPaletteFill) +
#   scale_colour_brewer(type = "qual", palette = "Set2") +
#   labs(
#     title = "First derivative of s(Year) by zone",
#     subtitle = "Derivative > 0 indicates increasing trend; < 0 decreasing",
#     x = "Year",
#     y = "First derivative ds/dYear"
#   ) +
#   ggthemes::theme_few() +
#   theme(legend.position = "none")
# 
# 
# ggplot(deriv2, aes(x = year, y = deriv)) +
#   geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
#   geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = zone1),
#               alpha = 0.25, colour = NA) +
#   geom_line(aes(colour = zone1), linewidth = 0.8) +
#   facet_wrap(~ zone1, scales = "free_y") +
#   scale_fill_manual(name = "", values = cbPaletteFill) +
#   scale_colour_brewer(type = "qual", palette = "Set2") +
#   labs(
#     title = "Second derivative of s(Year) by zone",
#     subtitle = "Curvature: positive = accelerating increase / easing decline; negative = turning down",
#     x = "Year",
#     y = expression(d^2*s/dYear^2)
#   ) +
#   ggthemes::theme_few() +
#   theme(legend.position = "none")
# 
# #Highlight only significant periods
# ## deriv1
# deriv1_sig <- deriv1 %>%
#   mutate(sig = .lower_ci > 0 | .upper_ci < 0)
# 
# ggplot(deriv1_sig, aes(x = year, y = deriv)) +
#   geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
#   geom_ribbon(aes(ymin = ifelse(sig, .lower_ci, NA_real_),
#                   ymax = ifelse(sig, .upper_ci, NA_real_),
#                   fill = zone1),
#               alpha = 0.35, colour = NA) +
#   geom_line(aes(colour = zone1), linewidth = 0.8) +
#   facet_wrap(~ zone1, scales = "free_y") +
#   scale_fill_manual(values = cbPaletteFill) +
#   scale_colour_brewer(type = "qual", palette = "Set2") +
#   labs(title = "First derivative (significant periods highlighted)",
#        x = "Year", y = "ds/dYear") +
#   ggthemes::theme_few() +
#   theme(legend.position = "none")
# 
# ## deriv2
# deriv2_sig <- deriv2 %>%
#   mutate(sig = .lower_ci > 0 | .upper_ci < 0)
# 
# ggplot(deriv2_sig, aes(x = year, y = deriv)) +
#   geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
#   geom_ribbon(aes(ymin = ifelse(sig, .lower_ci, NA_real_),
#                   ymax = ifelse(sig, .upper_ci, NA_real_),
#                   fill = zone1),
#               alpha = 0.35, colour = NA) +
#   geom_line(aes(colour = zone1), linewidth = 0.8) +
#   facet_wrap(~ zone1, scales = "free_y") +
#   scale_fill_manual(values = cbPaletteFill) +
#   scale_colour_brewer(type = "qual", palette = "Set2") +
#   labs(title = "Second derivative (significant periods highlighted)",
#        x = "Year", y = "ds/dYear") +
#   ggthemes::theme_few() +
#   theme(legend.position = "none")
# 
# ## angle ####
# ### fix formatting ####
# df_tm_ang <- df_tm %>% 
#   dplyr::filter(type=="angle",
#                 zone1 != "Wash") %>% 
#   dplyr::mutate(value = as.numeric(value))
# 
# 
# ang_fit <- mgcv::gam(
#   value ~ zone1 + s(year, by = zone1,bs="cr"),
#   data = df_tm_ang,
#   method = "REML")
# summary(ang_fit)  
# draw(ang_fit)



# tidy up ####
### Unload data
rm(list = ls(pattern = "^ang"))
rm(list = ls(pattern = "^com"))
rm(list = ls(pattern = "^df"))
rm(list = ls(pattern = "^mean_"))
rm(list = ls(pattern = "^cur_"))
rm(list = ls(pattern = "^pl_"))
rm(list = ls(pattern = "^cbP"))
rm(ppi,mn_st_cur, cur.yr, perm,fol,gisfol,projfol,sum_zero)

### Unload packages
detach("package:tidyverse", unload=TRUE)
detach("package:effects", unload=TRUE)
# detach("package:flexplot", unload=TRUE)
detach("package:ggthemes", unload=TRUE)
detach("package:lmerTest", unload=TRUE)
detach("package:sjPlot", unload=TRUE)
detach("package:visreg", unload=TRUE)
