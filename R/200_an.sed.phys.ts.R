# 200_an.sed.phys.ts.R ####
## Time series analysis ####
### import & analyse cone penetrometer and beach angle data

# Set up ####
### load packages ####
ld_pkgs <- c("tidyverse","ggthemes","lmerTest","effects","tictoc",
             "sjPlot","visreg","mgcv","gratia","patchwork","ggnewscale")
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

# WIP time series stats WIP ####

# compaction data ####
## fix formatting ####
df_tm_com <- df %>% 
  dplyr::filter(type=="cone",
                zone1 != "Wash") %>% 
  dplyr::mutate(value = as.numeric(value))

## Run GAM ####
### create interaction term for model ####
df_tm_com <- df_tm_com %>%
  mutate(zone_shore = interaction(zone1, shore, drop = TRUE))

### generate model ####
com_fit_B <- mgcv::gam(
  value ~ zone1 + shore + s(year, by = zone_shore, bs = "cr"),
  data = df_tm_com %>% select(year,zone1,shore,value,zone_shore),
  method = "REML",
)

zone1 <- rep(rep(c("Above","Inside","Inside2","Below"),each=2),250)
shore <- rep(rep(c("Upper","Mid"),4),250)
year <- seq(from=2008, to = 2025, length = 2000)

com_pred <- data.frame(
  year,
  zone1,
  shore
  ) %>% 
  mutate(zone_shore = interaction(zone1, shore, drop = TRUE))

predict(com_fit_B, newdata = com_pred) -> nu
as_tibble(bind_cols(com_pred,as_tibble(nu))) -> com_pred; rm(nu)

# com_fit_tw <- mgcv::gam(
#   value ~ zone1 + shore + s(year, by = zone_shore, bs = "cr"),
#   data = df_tm_com,
#   method = "REML",
#   family = Tweedie()
# )

### initial plots ####
# p1 <- draw(com_fit_B,select="s(year):zone_shoreAbove.Upper")
# p2 <- draw(com_fit_B,select="s(year):zone_shoreInside.Upper")
# p3 <- draw(com_fit_B,select="s(year):zone_shoreInside2.Upper")
# p4 <- draw(com_fit_B,select="s(year):zone_shoreBelow.Upper")
# p5 <- draw(com_fit_B,select="s(year):zone_shoreAbove.Mid")
# p6 <- draw(com_fit_B,select="s(year):zone_shoreInside.Mid")
# p7 <- draw(com_fit_B,select="s(year):zone_shoreInside2.Mid")
# p8 <- draw(com_fit_B,select="s(year):zone_shoreBelow.Mid")
# 
# (com_fit_B_draw <- p1+p2+p3+p4+p5+p6+p7+p8+patchwork::plot_layout(ncol = 4,nrow=2))
# rm(p1,p2,p3,p4,p5,p6,p7,p8)
# 
## generate model estimates ####
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

ilink <- gratia::inv_link(com_fit_B)
sm %>% 
  mutate(meas = ilink(.estimate))

## plot raw values ####

zones <- unique(df$zone1)
bg_cols <- cbPaletteFill[c(1:4,7)]
bg_cols <- setNames(bg_cols[seq_along(zones)], zones)

strip_elems <- lapply(zones, function(z)
  element_rect(fill = bg_cols[[z]], color = "black", linewidth = 1)
  )

# png(file = paste0("figs/sed.ts.mor.pen_vals.png"),
#     width=12*ppi, height=6*ppi, res=ppi)
df %>% 
  filter(.,type=="cone") %>% 
  filter(., zone1 != "Wash") %>%
  droplevels(.) %>% 
  mutate(value = as.numeric(value)) %>% 
  filter(., value >= 0) %>% 
  dplyr::select(.,c(year,shore,type,value,zone1)) %>% 
  ggplot(.,
         aes(y = as.numeric(value), x = year, fill = zone1))+
  geom_vline(xintercept = seq(2008,cur.yr,by=1),linetype=2, colour="lightgrey")+
  geom_boxplot(aes(group=year),outlier.shape = NA,show.legend = FALSE)+
  geom_jitter(width = 0.1, height = 0,alpha=0.3,show.legend = FALSE)+
  # geom_smooth(method = "loess", colour = "red", span = .9)+
  geom_smooth(method = "gam", colour = "red", #span = .9
              show.legend = FALSE
  )+
  # facet_grid(shore~zone1)+
  ggh4x::facet_grid2(
    rows = vars(shore),
    cols = vars(zone1),
    strip = strip_themed(
      background_x = strip_elems)
    )+
  scale_colour_manual(name = "", values=cbPalette)+
  scale_fill_manual(name = "", values=cbPaletteFill[c(1:4,7)])+
  # scale_x_continuous(breaks = seq(2008, cur.yr, by = 4))+
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0))+
  xlab("Year") + ylab("Observed cone index")+
  labs(title = paste0("Sediment compaction recorded since 2008 as part of the SFGPBM programme"),
       #subtitle = "Higher values indicate more compacted sediments.\nRed lines indicate generalised additive model trend"
       )+
  # theme(
  #   legend.position="none",
  #   axis.title.x = element_blank(),
  #   axis.text.y = element_text(size = 12),
  #   axis.text.x = element_text(size = 12, face = 2),
  #   axis.title.y = element_text(size = 14),
  #   strip.text.x = element_text(size = 14),
  #   strip.text.y = element_text(size = 14),
  #   strip.text = element_text(face="bold"),
  #   strip.background = element_rect(color = "black",fill = "grey95", size = 1),
    theme(
      plot.title = element_text(face=2,size=18),
      plot.subtitle = element_text(face=2,size=12),
      plot.caption = element_text(face=2,size=12),
      axis.title.y = element_text(face=2),
      axis.title.x=element_blank(),
      axis.text.y = element_text(face=2),
      axis.text.x = element_text(face=2,size = 12),
      strip.text = element_text(face=2,size=14),
      strip.background = element_rect(color = "black",fill = "grey95", size = 1),
    # )
    
    ) -> pl_com_A
# dev.off()

## plot model smooths ####
# (ggplot(sm, aes(x = year, y = .estimate)) +
#   geom_vline(xintercept = c(2008:cur.yr),linetype = 2,linewidth = 0.4,col="lightgrey")+
#   geom_hline(yintercept = 0, linetype = 2) +
#   geom_ribbon(aes(ymin = .lower_ci,
#                   ymax = .upper_ci,
#                   fill = zone1),
#               alpha = 0.25, colour = NA) +
#   geom_line(aes(colour = zone1), linewidth = 0.9) +
#   # Facet by zone and shore for a grid view (readable across management areas and shore levels)
#   facet_grid(shore ~ zone1, scales = "free_y") +
#   scale_fill_manual(values = cbPaletteFill, guide = "none") +
#   # If you want lines to match ribbons exactly, use scale_colour_manual with the same vector
#   scale_colour_manual(values = cbPalette, name = "Zone") +
#   labs(
#     title = "Sediment compaction trends since 2008",
#     # subtitle = "GAM smooths of s(Year) by zone × shore\nRibbons show simultaneous confidence intervals",
#     x = "Year",
#     y = "Smooth effect s(Year)"
#   ) +
#   ggthemes::theme_few() +
#   theme(
#     legend.position = "none",
#     plot.title = element_text(face=2,size=18),
#     plot.subtitle = element_text(face=2,size=12),
#     axis.title.y = element_text(face=2),
#     axis.title.x=element_blank(),
#     axis.text.y = element_text(face=2),
#     axis.text.x = element_text(face=2,size = 12),
#     strip.text = element_text(face=2,size=14),
#     strip.background = element_rect(color = "black",fill = "grey95", size = 1),
#   ) ->com_pl_A)

## Identify change points ####
### calculate first derivatives ####
# Get exact smooth labels from the model
terms_year_by <- gratia::smooths(com_fit_B) %>%
  stringr::str_subset("^s\\(year\\):")

# Calculate first derivative
com_deriv1 <- gratia::derivatives(
  com_fit_B,
  term = terms_year_by,
  order = 1,
  interval = "simultaneous",
  n = 100,
  type="central",
  unconditional = TRUE
  ) %>%
  rename(year = year, deriv = .derivative) %>%
  mutate(zone_shore = sub("^s\\(year\\):", "", .smooth)) %>%
  separate(zone_shore, into = c("zone1", "shore"), sep = "\\.",
           remove = FALSE) %>% 
  mutate(zone = gsub("^.{0,10}", "", zone1)) %>% 
  mutate(
    zone = factor(zone,levels=c("Above","Inside","Inside2","Below")),
    shore = factor(shore,levels=c("Upper","Mid","Low")),
  ) ->com_deriv1_tmp
  
## bind smooth and derivatives together ####
bind_cols(
  sm %>% rename_with( ~ paste0("sm_", .x)),
  com_deriv1_tmp %>% rename_with( ~ paste0("d1_", .x))
    ) -> com_sm_d1
rm(com_deriv1,com_deriv1_tmp,sm)

names(com_sm_d1)


## calculate changepoints
com_sm_d1 <- com_sm_d1 %>% 
  arrange(sm_zone_shore,sm_year) %>% 
  ## identify direction of change based on derivatives
  mutate(
    sig_pos = d1_.lower_ci >0,
    sig_neg = d1_.upper_ci <0,
    sig_dir = case_when(
      sig_pos ~ "increasing",
      sig_neg ~ "decreasing",
      TRUE ~ "uncertain"
    )
  ) %>% #names()
  ## is this a changepoint?
  mutate(
    direction_prev = lag(sig_dir)) %>% 
  mutate(changepoint = sig_dir != direction_prev &
           sig_dir %in% c("increasing","decreasing")) %>% 
  mutate(d1_zone1 = ifelse(grepl("Above",d1_.smooth),"Above",
                        ifelse(grepl("Inside2",d1_.smooth),"Inside2",
                               ifelse(grepl("Below",d1_.smooth),"Below","Inside")))) %>% 
  mutate(d1_zone1 = factor(d1_zone1, levels = c("Above","Inside","Inside2","Below")))
rm(terms_year_by)

## extract changepoints
com_sm_d1 %>% filter(changepoint == TRUE) %>% 
  arrange(sm_year) %>% slice(-1) -> com_deriv1_ch

## plot ####
# png(file = paste0("figs/sed.ts.mor.pen_changepts.png"),
#     width=12*ppi, height=6*ppi, res=ppi)
com_mean <- mean(as.numeric(df %>% filter(value=="cone")),na.rm=TRUE)
com_sm_d1 %>% #names()
  ggplot(.,
         aes(
           x = sm_year,
           y = sm_.estimate
           # y = sm_.estimate + com_fit_B$coefficients[1]
           )
         ) +
  geom_vline(xintercept = c(2008:cur.yr),linetype = 2,linewidth = 0.4,col="lightgrey")+
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_ribbon(aes(
    ymin = sm_.lower_ci, ymax = sm_.upper_ci,
    # ymin = sm_.lower_ci + com_fit_B$coefficients[1], ymax = sm_.upper_ci+com_fit_B$coefficients[1],
    fill = sm_zone1
    ),
    alpha = 0.25, colour = NA) +
  geom_line(aes(colour = sm_zone1), linewidth = 0.8,show.legend = FALSE) +
  # facet_grid(sm_shore ~ sm_zone1, scales = "free_y") +
  ggh4x::facet_grid2(
    rows = vars(sm_shore),
    cols = vars(sm_zone1),
    strip = strip_themed(
      background_x = strip_elems)
  )+
  scale_fill_manual(values = cbPaletteFill, guide = "none") +
  scale_colour_manual(values = cbPalette[1:4]) +
  ggnewscale::new_scale_fill()+ggnewscale::new_scale_colour()+
  geom_point(data = com_deriv1_ch,
             aes(
               x = sm_year,
               y = sm_.estimate,
               # y = sm_.estimate + com_fit_B$coefficients[1],
               shape = sig_dir,
               fill = sig_dir
             ),
             show.legend = FALSE,
             size = 5,
             colour = 1,
             stroke = 1.2, # control point thickness
             # pch = 21
             )+
  labs(
    #title = "Time series of sediment compaction values recorded since 2008",
    caption = "Symbols indicate change points to increasing (black triangles) and decreasing (white circles) trends
    Lines indicated generalised additive model trends ± 95% CI",
    y = "Model estimate"
  )+
  scale_shape_manual(values = c(21,24))+
  scale_fill_manual(values = c("increasing"="black", "decreasing" = "white")) +
  # scale_colour_manual(values = c("increasing"="blue", "decreasing" = "red")) +
  ggthemes::theme_few() +
  theme(
    plot.title = element_text(face=2,size=18),
    plot.subtitle = element_text(face=2,size=12),
    plot.caption = element_text(face=2,size=12),
    axis.title.y = element_text(face=2),
    axis.title.x=element_blank(),
    axis.text.y = element_text(face=2),
    axis.text.x = element_text(face=2,size = 12),
    strip.text = element_text(face=2,size=14),
    strip.background = element_rect(color = "black",fill = "grey95", size = 1),
    ) -> pl_com_B
# dev.off()

## export patchwork object ####
png(file = "figs/sed.ts.mor.pen_changeptsAB.png",
    width=12*ppi, height=10*ppi, res=ppi)
pl_com_A/pl_com_B + patchwork::plot_annotation(tag_levels = "A") +
  plot_layout(axes = "collect") & 
  theme(plot.tag = element_text(face = 2, size = 18))
dev.off()

## remove items beginning with "com" to avoid cross-contamination of outputs
rm(list = ls(pattern = "^com"))
rm(df_tm_com, pl_com_A,pl_com_B)

# reproduce for angle data ####
# angle data ####
## fix formatting ####
df_tm_ang <- df %>% 
  dplyr::filter(type=="angle",
                zone1 != "Wash") %>% 
  dplyr::mutate(value = as.numeric(value))

## Run GAM ####
### create interaction term for model ####
df_tm_ang <- df_tm_ang %>%
  mutate(zone_shore = interaction(zone1, shore, drop = TRUE))

### generate model ####
com_fit_B <- mgcv::gam(
  value ~ zone1 + shore + s(year, by = zone_shore, bs = "cr"),
  data = df_tm_ang,
  method = "REML",
)

### initial plots ####
# p1 <- draw(com_fit_B,select="s(year):zone_shoreAbove.Upper")
# p2 <- draw(com_fit_B,select="s(year):zone_shoreInside.Upper")
# p3 <- draw(com_fit_B,select="s(year):zone_shoreInside2.Upper")
# p4 <- draw(com_fit_B,select="s(year):zone_shoreBelow.Upper")
# p5 <- draw(com_fit_B,select="s(year):zone_shoreAbove.Mid")
# p6 <- draw(com_fit_B,select="s(year):zone_shoreInside.Mid")
# p7 <- draw(com_fit_B,select="s(year):zone_shoreInside2.Mid")
# p8 <- draw(com_fit_B,select="s(year):zone_shoreBelow.Mid")
# 
# (com_fit_B_draw <- p1+p2+p3+p4+p5+p6+p7+p8+patchwork::plot_layout(ncol = 4,nrow=2))
# rm(p1,p2,p3,p4,p5,p6,p7,p8)

## generate model estimates ####
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

ilink <- gratia::inv_link(com_fit_B)
sm %>% 
  mutate(meas = ilink(.estimate))

## plot raw values ####
df %>% 
  filter(.,type=="angle") %>% 
  filter(., zone1 != "Wash") %>%
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
  # facet_grid(shore~zone1)+
  ggh4x::facet_grid2(
    rows = vars(shore),
    cols = vars(zone1),
    strip = strip_themed(
      background_x = strip_elems)
  )+
  scale_colour_manual(name = "", values=cbPalette)+
  scale_fill_manual(name = "", values=cbPaletteFill)+
  # scale_x_continuous(breaks = seq(2008, cur.yr, by = 4))+
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0))+
  xlab("Year") + ylab("Observed beach slope")+
  labs(title = paste0("Beach slopes recorded since 2008 as part of the SFGPBM programme"),
       # subtitle = "Higher values indicate more compacted sediments.\nRed lines indicate generalised additive model trend"
       )+
  theme(
    plot.title = element_text(face=2,size=18),
    plot.subtitle = element_text(face=2,size=12),
    plot.caption = element_text(face=2,size=12),
    axis.title.y = element_text(face=2),
    axis.title.x=element_blank(),
    legend.position="none",
    axis.text.y = element_text(face=2),
    axis.text.x = element_text(face=2,size = 12),
    strip.text = element_text(face=2,size=14),
    strip.background = element_rect(color = "black",fill = "grey95", size = 1),  # theme(
  #   
  #   legend.position="none",
  #   axis.title.x = element_blank(),
  #   axis.text.y = element_text(size = 12),
  #   axis.text.x = element_text(size = 12, face = 2),
  #   axis.title.y = element_text(size = 14),
  #   strip.text.x = element_text(size = 14),
  #   strip.text.y = element_text(size = 14),
  #   strip.text = element_text(face="bold"),
  #   strip.background = element_rect(color = "black",fill = "grey95", size = 1),
  ) -> pl_ang_A

##########

# ## plot model smooths ####
# (ggplot(sm, aes(x = year, y = .estimate)) +
#    geom_vline(xintercept = c(2008:cur.yr),linetype = 2,linewidth = 0.4,col="lightgrey")+
#    geom_hline(yintercept = 0, linetype = 2) +
#    geom_ribbon(aes(ymin = .lower_ci,
#                    ymax = .upper_ci,
#                    fill = zone1),
#                alpha = 0.25, colour = NA) +
#    geom_line(aes(colour = zone1), linewidth = 0.9) +
#    # Facet by zone and shore for a grid view (readable across management areas and shore levels)
#    facet_grid(shore ~ zone1, scales = "free_y") +
#    scale_fill_manual(values = cbPaletteFill, guide = "none") +
#    # If you want lines to match ribbons exactly, use scale_colour_manual with the same vector
#    scale_colour_manual(values = cbPalette, name = "Zone") +
#    labs(
#      title = "Beach slope trends since 2008",
#      # subtitle = "GAM smooths of s(Year) by zone × shore\nRibbons show simultaneous confidence intervals",
#      x = "Year",
#      y = "Smooth effect s(Year)"
#    ) +
#    ggthemes::theme_few() +
#    theme(
#      legend.position = "none",
#      plot.title = element_text(face=2,size=18),
#      plot.subtitle = element_text(face=2,size=12),
#      axis.title.y = element_text(face=2),
#      axis.title.x=element_blank(),
#      axis.text.y = element_text(face=2),
#      axis.text.x = element_text(face=2,size = 12),
#      strip.text = element_text(face=2,size=14),
#      strip.background = element_rect(color = "black",fill = "grey95", size = 1),
#    ) ->com_pl_A)

## Identify change points ####
### calculate first derivatives ####
# Get exact smooth labels from the model
terms_year_by <- gratia::smooths(com_fit_B) %>%
  stringr::str_subset("^s\\(year\\):")

# Calculate first derivative
ang_deriv1 <- gratia::derivatives(
  com_fit_B,
  term = terms_year_by,
  order = 1,
  interval = "simultaneous",
  n = 100,
  type="central",
  unconditional = TRUE
) %>%
  rename(year = year, deriv = .derivative) %>%
  mutate(zone_shore = sub("^s\\(year\\):", "", .smooth)) %>%
  separate(zone_shore, into = c("zone1", "shore"), sep = "\\.",
           remove = FALSE) %>% 
  mutate(zone = gsub("^.{0,10}", "", zone1)) %>% 
  mutate(
    zone = factor(zone,levels=c("Above","Inside","Inside2","Below")),
    shore = factor(shore,levels=c("Upper","Mid","Low")),
  ) ->ang_deriv1_tmp

## bind smooth and derivatives together ####
bind_cols(
  sm %>% rename_with( ~ paste0("sm_", .x)),
  ang_deriv1_tmp %>% rename_with( ~ paste0("d1_", .x))
) -> ang_sm_d1
rm(ang_deriv1,ang_deriv1_tmp,sm)

names(ang_sm_d1)


## calculate changepoints
ang_sm_d1 <- ang_sm_d1 %>% 
  arrange(sm_zone_shore,sm_year) %>% 
  ## identify direction of change based on derivatives
  mutate(
    sig_pos = d1_.lower_ci >0,
    sig_neg = d1_.upper_ci <0,
    sig_dir = case_when(
      sig_pos ~ "increasing",
      sig_neg ~ "decreasing",
      TRUE ~ "uncertain"
    )
  ) %>% #names()
  ## is this a changepoint?
  mutate(
    direction_prev = lag(sig_dir)) %>% 
  mutate(changepoint = sig_dir != direction_prev &
           sig_dir %in% c("increasing","decreasing")) %>% 
  mutate(d1_zone1 = ifelse(grepl("Above",d1_.smooth),"Above",
                           ifelse(grepl("Inside2",d1_.smooth),"Inside2",
                                  ifelse(grepl("Below",d1_.smooth),"Below","Inside")))) %>% 
  mutate(d1_zone1 = factor(d1_zone1, levels = c("Above","Inside","Inside2","Below")))
rm(terms_year_by)

## extract changepoints
ang_sm_d1 %>% filter(changepoint == TRUE) %>% 
  arrange(sm_year) %>% slice(-1) -> ang_deriv1_ch

## plot ####
png(file = paste0("figs/sed.ts.mor.ang_changepts.png"),
    width=12*ppi, height=6*ppi, res=ppi)
ang_mean <- mean(as.numeric(df %>% filter(value=="angle")),na.rm=TRUE)
ang_sm_d1 %>% #names()
  ggplot(.,
         aes(
           x = sm_year,
           y = sm_.estimate
           # y = sm_.estimate + com_fit_B$coefficients[1]
         )
  ) +
  geom_vline(xintercept = c(2008:cur.yr),linetype = 2,linewidth = 0.4,col="lightgrey")+
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_ribbon(aes(
    ymin = sm_.lower_ci, ymax = sm_.upper_ci,
    # ymin = sm_.lower_ci + com_fit_B$coefficients[1], ymax = sm_.upper_ci+com_fit_B$coefficients[1],
    fill = sm_zone1
  ),
  alpha = 0.25, colour = NA) +
  geom_line(aes(colour = sm_zone1), linewidth = 0.8,show.legend = FALSE) +
  # facet_grid(sm_shore ~ sm_zone1, scales = "free_y") +
  ggh4x::facet_grid2(
    rows = vars(sm_shore),
    cols = vars(sm_zone1),
    strip = strip_themed(
      background_x = strip_elems)
  )+
  scale_fill_manual(values = cbPaletteFill, guide = "none") +
  scale_colour_manual(values = cbPalette[1:4]) +
  ggnewscale::new_scale_fill()+ggnewscale::new_scale_colour()+
  geom_point(data = ang_deriv1_ch,
             aes(
               x = sm_year,
               y = sm_.estimate,
               # y = sm_.estimate + com_fit_B$coefficients[1],
               shape = sig_dir,
               stroke = 1.2, # control point thickness
               fill = sig_dir
             ),
             show.legend = FALSE,
             size = 5,
             colour = 1,
             # pch = 21
  )+
  labs(
    # title = "Time series of sediment compaction values recorded since 2008",
    caption = "Symbols indicate change points to increasing (black triangles) and decreasing (white circles) trends
    Lines indicated generalised additive model trends ± 95% CI",
    y = "Model estimate"
  )+
  scale_shape_manual(values = c(21,24))+
  scale_fill_manual(values = c("increasing"="black", "decreasing" = "white")) +
  # scale_colour_manual(values = c("increasing"="blue", "decreasing" = "red")) +
  ggthemes::theme_few() +
  theme(
    plot.title = element_text(face=2,size=18),
    plot.subtitle = element_text(face=2,size=12),
    plot.caption = element_text(face=2,size=12),
    axis.title.y = element_text(face=2),
    axis.title.x=element_blank(),
    axis.text.y = element_text(face=2),
    axis.text.x = element_text(face=2,size = 12),
    strip.text = element_text(face=2,size=14),
    strip.background = element_rect(color = "black",fill = "grey95", size = 1),
  ) -> pl_ang_B
dev.off()

## export patchwork object ####
png(file = paste0("figs/sed.ts.mor.ang_changeptsAB.png"),
    width=12*ppi, height=10*ppi, res=ppi)
pl_ang_A/pl_ang_B + patchwork::plot_annotation(tag_levels = "A") +
  plot_layout(axes = "collect") & 
  theme(plot.tag = element_text(face = 2, size = 18))
dev.off()



## remove items beginning with "com" to avoid cross-contamination of outputs
rm(list = ls(pattern = "^com"))
rm(list = ls(pattern = "^pl"))
rm(cone.mn, df_tm_com)
rm(list = ls(pattern = "^cb"))
rm(list = ls(pattern = "^ang"))
rm(list = ls(pattern = "^df"))
rm(strip_elems,bg_cols,cur.yr,fol,gisfol,perm,ppi,projfol,
   shore,year,zone1,zones,granstat,ilink,sum_zero)


