## 200_imp.int.sed.R ##

## generate plots from gradistat (G2Sd) outputs

# Set up ####
### load packages ####
ld_pkgs <- c("tidyverse","tictoc","ggthemes","ggh4x")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)

tictoc::tic.clearlog();tic("SET UP")

### load metadata ####
source("R/00_meta_setMeta.R")

# load data
df <- readRDS(file="output/models/grad_all_raw.Rdat")

# Produce plots ####

# Build a strip background list that matches your facet order
zones <- unique(df$sedim$texture$zone1)
bg_cols <- cbPaletteFill[c(1:4,7)]
bg_cols <- setNames(bg_cols[seq_along(zones)], zones)

strip_elems <- lapply(zones, function(z)
  element_rect(fill = bg_cols[[z]], color = "black", linewidth = 1)
)

## Sediment type ####
png(file = "figs/sed.ts.texture.png",
    width=12*ppi, height=10*ppi, res=ppi)
df$sedim$texture %>% as.data.frame(.) %>%
  filter(method == "5cm") %>%
  ggplot(.,
         aes(
           x=year,
           y=texture,
           # colour = (paste0(shore)),
           colour = shore,
           shape = shore,
         )) +
  geom_vline(xintercept = c(min(df$sedim$texture$year-.5):
                              max(df$sedim$texture$year)),
             lty=2,
             colour = "lightgrey"
             )+
  geom_hline(yintercept = c(0.5:length(unique(df$sedim$texture$texture))),
             lty=2,
             colour = "lightgrey"
             )+
  geom_jitter(width = 0.25,
              size = 4, alpha = 0.7
              )+
  scale_x_continuous(breaks=min(df$sedim$texture$year):
                       max(df$sedim$texture$year))+
  scale_shape_manual(values=c(15:18))+
  # facet_wrap(.~zone1,ncol = 2)+
  facet_wrap2(. ~ zone1, ncol = 2,
              strip = strip_themed(background_x = strip_elems)) +
  scale_fill_manual(values=cbPaletteshr) +
  scale_color_manual(values=cbPaletteshr) +
  theme(
    plot.title = element_text(face=2),
    strip.text = element_text(face=2),
    axis.text.x = element_text(face=2,angle=270,vjust=0.5),
    axis.text.y = element_text(face=2),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.direction = "horizontal",
    legend.position = "inside",
    legend.position.inside = c(0.8,0.1),
    legend.title = element_blank(),
    legend.text = element_text(face = 2)
  )
dev.off()

# convert Folk Ward to long then generate plots programatically ####
tictoc::tic("Generate plots programatically")
## convert to tibble
df_fowa <- as_tibble(df$stat$fowa)

## remove non-numeric variable and lengthen
df_fowa %>% 
  select(-sediment) %>% 
  filter(method=="5cm") %>% 
  pivot_longer(cols = -c(samples,year,transect,shore,method,zone1)) %>% 
  ## create name for variable
  mutate(label = case_when(
    name == "mean.fw.um" ~ "Mean (um)",
    name == "sd.fw.um" ~ "SD (um)",
    name == "skewness.fw.um" ~ "Skew (um)",
    name == "kurtosis.fw.um" ~ "Kurtosis (um)",
    name == "mean.fw.phi" ~ "Mean (phi)",
    name == "sd.fw.phi" ~ "SD (phi)",
    name ==  "skewness.fw.phi" ~ "Skew (phi)",
    name == "kurtosis.fw.phi" ~ "Kurtosis (phi)",
    TRUE ~ NA
    )) -> df_fowa_l

for(i in 1:length(unique(df_fowa_l$label))){
  lb <- unique(df_fowa_l$label)[i]
  df_tmp <- df_fowa_l %>% filter(label == lb)
  print(unique(df_tmp$label))
  
  df_tmp %>% 
    ggplot(
      .,
      aes(
        x = year,
        y = value,
        colour = shore,
        shape = shore,
        )
    ) +
    geom_jitter(width = 0.25,
                size = 4, alpha = 0.7
                )+
    labs(
      title = paste0(unique(df_tmp$label)),
      y = (paste0(unique(df_tmp$label))
           ),
      caption = "Lines represent generalised additive model trend (where sufficient data allow)",
      ) +
    scale_shape_manual(values=c(15:18))+
    geom_smooth(inherit.aes = FALSE,
                aes(
                  x = year,
                  # y = texture,
                  y = value,
                  # group = (paste0(shore,zone1)),
                  group = zone1,
                  colour = (paste0(shore)),
                ),
                se=FALSE,
                show.legend = FALSE,
                method = "gam",
    ) +
    scale_x_continuous(breaks=min(df_tmp$year):max(df_tmp$year))+
    # facet_wrap(.~zone1,ncol = 2)+
    facet_wrap2(. ~ zone1, ncol = 2,
                strip = strip_themed(background_x = strip_elems)) +
    scale_fill_manual(values=cbPaletteshr) +
    scale_color_manual(values=cbPaletteshr)+
    theme(
      plot.title = element_text(face=2),
      strip.text = element_text(face=2),
      axis.text.x = element_text(face=2,angle=270,vjust=0.5),
      axis.text.y = element_text(face=2),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face=2),
      legend.direction = "horizontal",
      legend.position = "inside",
      legend.position.inside = c(0.8,0.1),
      legend.title = element_blank(),
      legend.text = element_text(face = 2),
      # strip.background = element_rect(color = "black",fill = "grey95",
      #                                 linewidth = 1),
    ) -> pl
  
  png(file = paste0("figs/sed.ts.",unique(df_tmp$name),".png"),
                    width=12*ppi, height=10*ppi, res=ppi)
  print(pl)
  dev.off()
  
  ## with zero line
  df_tmp %>% 
    ggplot(
      .,
      aes(
        x = year,
        y = value,
        colour = shore,
        shape = shore,
      )
    ) +
    geom_hline(yintercept = 0,lty = 2, colour = "lightgrey")+
    geom_jitter(width = 0.25,
                size = 4, alpha = 0.7
    )+
    labs(
      title = paste0(unique(df_tmp$label)),
      y = (paste0(unique(df_tmp$label))
      ),
      caption = "Lines represent generalised additive model trend (where sufficient data allow)",
    ) +
    scale_shape_manual(values=c(15:18))+
    geom_smooth(inherit.aes = FALSE,
                aes(
                  x = year,
                  # y = texture,
                  y = value,
                  # group = (paste0(shore,zone1)),
                  group = zone1,
                  colour = (paste0(shore)),
                ),
                se=FALSE,
                show.legend = FALSE,
                method = "gam",
    ) +
    scale_x_continuous(breaks=min(df_tmp$year):max(df_tmp$year))+
    # facet_wrap(.~zone1,ncol = 2)+
    facet_wrap2(. ~ zone1, ncol = 2,
                strip = strip_themed(background_x = strip_elems)) +
    scale_fill_manual(values=cbPaletteshr) +
    scale_color_manual(values=cbPaletteshr)+
    theme(
      plot.title = element_text(face=2),
      strip.text = element_text(face=2),
      axis.text.x = element_text(face=2,angle=270,vjust=0.5),
      axis.text.y = element_text(face=2),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face=2),
      legend.direction = "horizontal",
      legend.position = "inside",
      legend.position.inside = c(0.8,0.1),
      legend.title = element_blank(),
      legend.text = element_text(face = 2),
      # strip.background = element_rect(color = "black",fill = "grey95",
      #                                 linewidth = 1),
    ) -> pl2
  
  png(file = paste0("figs/sed.ts.",unique(df_tmp$name),".line.png"),
      width=12*ppi, height=10*ppi, res=ppi)
  print(pl2)
  dev.off()
  
  rm(pl,df_tmp, pl2)
      
}
tictoc::toc(log = TRUE)

# generate table of sediment statistics ####
group_vars <- c("year", "shore", "method", "zone1")

df$index %>%
  ungroup() %>% 
  # optionally drop columns you don't want in the summary
  select(-samples, -transect) %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(
    across(
      # all numeric columns except the grouping vars
      .cols = where(is.numeric) & !any_of(group_vars),
      # compute both mean and sd
      .fns  = list(mean = ~mean(.x, na.rm = TRUE),
                   sd   = ~sd(.x,   na.rm = TRUE)),
      # name like var_mean, var_sd
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  # set NA SDs to 0, keep means untouched
  mutate(across(ends_with("_sd"), ~tidyr::replace_na(.x, 0))) -> stats_out

write.csv(stats_out, file = "output/sed.stats.csv",row.names=FALSE)

df$index %>%
  ungroup() %>% 
  filter(year == cur.yr, method=="5cm") %>% 
  # optionally drop columns you don't want in the summary
  select(-samples, -transect) %>%
  group_by(across(all_of(group_vars))) %>%
  
  ####
  ## sorting calculation not required.  The @`SD parameter in the output
  ## is actually the Sorting value 
  ####
  
  # # calculate sorting value
  # ### Based on Folk & Ward (1957)
  
  # ## First, convert diameters to phi
  # mutate(
  #   D5_phi = -log(D5/1000,base = 2),
  #   D16_phi = -log(D16/1000,base = 2),
  #   D84_phi = -log(D84/1000,base = 2),
  #   D95_phi = -log(D95/1000,base = 2),) %>% 
  # ##Next, calculate sorting
  # mutate(
  #   sort = ((D84_phi-D16_phi)/4)+((D95_phi-D5_phi)/6.6)
  # ) %>%
  summarise(
    across(
      # all numeric columns except the grouping vars
      .cols = where(is.numeric) & !any_of(group_vars),
      # compute both mean and sd
      .fns  = list(mean = ~mean(.x, na.rm = TRUE),
                   sd   = ~sd(.x,   na.rm = TRUE)),
      # name like var_mean, var_sd
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  # set NA SDs to 0, keep means untouched
  mutate(across(ends_with("_sd"),
                ~tidyr::replace_na(.x, 0))) -> stats_out_cur.yr

write.csv(stats_out_cur.yr, file = "output/sed.stats.cur.yr.csv",row.names=FALSE)

df$stat$fowa %>% 
  ungroup() %>% 
  filter(year == cur.yr, method=="5cm") %>% 
  # optionally drop columns you don't want in the summary
  select(-samples, -transect) %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(
    across(
      # all numeric columns except the grouping vars
      .cols = where(is.numeric) & !any_of(group_vars),
      # compute both mean and sd
      .fns  = list(mean = ~mean(.x, na.rm = TRUE),
                   sd   = ~sd(.x,   na.rm = TRUE)),
      # name like var_mean, var_sd
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  # set NA SDs to 0, keep means untouched
  mutate(across(ends_with("_sd"),
                ~tidyr::replace_na(.x, 0))) -> stats_out_bulk_cur.yr

write.csv(stats_out_bulk_cur.yr, file = "output/sed.stats.bulk.cur.yr.csv",row.names=FALSE)
