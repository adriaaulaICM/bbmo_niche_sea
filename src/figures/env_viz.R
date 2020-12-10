library(tidyverse)
library(speedyseq)
library(lomb)
library(patchwork)
library(ggthemes)

source("src/util/backbone_functions.R") 
source("src/util/backbone_params-graphs.R") 
source('src/analysis/helper_prepare_nichepart_data.R')

bl.phy <- readRDS('data/cleaned/blphyloseq.rds')

labels <- readxl::read_xlsx('data/raw/metadata/labels_pretty.xlsx',
                            sheet = 2,
                            col_names = c('Variable', 'Measure'))  %>% 
  mutate(all = str_c(Variable, "~", Measure ))


env.raw <- sample_data(bl.phy) %>% as_tibble(rownames = 'sample')

par.gen <-  c("Day_length", "Temperature", "Salinity", "Secchi",
                "NO2", "NO3", "PO4", "Si", 'Chla_3um',"Chla_total", "BP_FC1.55",
                "Bacteria_joint", "HNA", "LNA", "Synechococcus", "Prochlorococcus_FC",
                "Peuk1", "Peuk2", "PNF_Micro", "PNF2_5um_Micro", "PNF_5um_Micro",
                "Cryptomonas", "Micromonas", "HNF_Micro")

params <- c("Day_length", "Temperature", "Salinity", "Secchi",
    "NO2", "NO3", "PO4", "Si", 'Chla_3um', "Chla_total", 'BP_FC1.55',
    "Bacteria_joint", "Synechococcus", "Prochlorococcus_FC",
    "PNF_Micro", "PNF2_5um_Micro", "PNF_5um_Micro",
    "Cryptomonas", "Micromonas", "HNF_Micro")

lomb.env <- env.raw %>% 
  select(decimal_date, one_of(params)) %>% 
  pivot_longer(names_to = 'param', values_to = 'value', -decimal_date) %>% 
  split(.$param) %>% 
  map(~randlsp( x =.x$value,
                times = .x$decimal_date,
                type = 'period',
                plot = F))

lomb.env.sea <- tibble( parameter = names(lomb.env),
                       pval = map_dbl(lomb.env, ~.x$p.value),
                       peak = map_dbl(lomb.env, ~.x$peak),
                       interval  = map(lomb.env, ~.x$peak.at),
                       int.min = map_dbl(interval, ~.[[2]]),
                       int.max = map_dbl(interval, ~.[[1]])) %>% 
  filter(pval <= 0.01, int.max <= 2)

plot_periodogram <- function(sel, lombobj){
 
  
  map(lombobj[sel], ~tibble( scanned = .x$scanned,
                              power = .x$power)) %>%
    bind_rows(.id = 'sel') %>%
    ggplot(aes(scanned, power)) +
    geom_hline(yintercept = 10, linetype =2 , color = 'grey2') + 
    geom_line(aes(group = sel)) +
    facet_wrap(~sel, scales = 'free_y') + 
    scale_x_continuous(breaks = c(0,1,2,3,4,5), limits = c(0,5)) 
  
}

plot_periodogram(params, lomb.env) + 
  lil.strip

ggsave('results/figures/env_periodograms.pdf', width = 8, height = 6)


# The main plots  ---------------------------------------------------------
env <- env.raw %>% 
  select(year, day_of_year, season, month, one_of(params[1:8])) %>% 
  gather(key = 'parameter', value = 'val',
         -month, -year, -day_of_year, -season,
         factor_key = T) %>%
  # Outlier value the 2008 , at may
  filter(!(year == 2008 & month == '05')) %>% 
  mutate( seasonal = ifelse(parameter %in% lomb.env.sea$parameter, 
                            TRUE,
                            FALSE),
          parameter = factor(parameter, levels = par.gen,
                             labels = labels$all)) 

num.days.mnt <- c(0,31,28,31,30,31,30,31,31,30,31,30)
cumnum <- cumsum(num.days.mnt)

envplot <- env  %>% 
  ggplot( aes(day_of_year, val)) + 
  geom_point(alpha = 0.3)  +
  stat_smooth(aes(group = parameter),
              method = "gam",
              color = 'dodgerblue',
              formula = y ~ s(x, k =12, bs = 'cc'),
              se = F, size = 1.2,
              alpha = 0.7) + 
  facet_wrap(~parameter, ncol = 3, 
             labeller = label_parsed,
             scales = 'free_y') + 
  ylab(NULL) +
  xlab('Day of the year (month label)') +
  scale_x_continuous(breaks = cumnum, labels = str_to_title(month.order)) +
  lab.flip + 
  lil.strip + 
  leg.bottom + 
  theme( strip.text = element_text(size = 8))

ggsave('results/figures/metadata_env.pdf',
       plot = envplot, 
       width = 8, 
       height = 9)

options(scipen=-3)

bio <- env.raw %>% 
  select(year, month, season, day_of_year, one_of(params[c(9:15,18,19)])) %>% 
  gather(key = 'parameter', value = 'val',
         -month, -year, -day_of_year, -season,
         factor_key = T) %>%
  # Outlier value the 2008 , at may
  filter(!(year == 2008 & month == '05')) %>% 
  mutate( seasonal = ifelse(parameter %in% lomb.env.sea$parameter, 
                            TRUE,
                            FALSE),
          parameter = factor(parameter, levels = par.gen,
                             labels = labels$all)) 


bioplot <- bio  %>% 
  ggplot( aes(day_of_year, val)) + 
  geom_point(alpha = 0.3)  +
  stat_smooth(aes(group = parameter),
              method = "gam",
              color = 'dodgerblue',
              formula = y ~ s(x, k =12, bs = 'cc'),
              se = F, size = 1.2,
              alpha = 0.7) + 
  facet_wrap(~parameter, ncol = 3, 
             labeller = label_parsed,
             scales = 'free_y')  + 
  ylab(NULL) +
  xlab('Day of the year (month label)') +
  scale_x_continuous(breaks = cumnum, labels = str_to_title(month.order)) +
  lab.flip + 
  lil.strip + 
  leg.bottom + 
  theme( strip.text = element_text(size = 8))

ggsave('results/figures/metadata_bio.pdf',
       plot = bioplot, 
       width = 8, 
       height = 9)

composite <- envplot + bioplot +
  plot_layout(ncol = 2, guides = 'collect') +
  plot_annotation(tag_levels = 'A') & 
  theme(legend.position = 'bottom')

ggsave('results/figures/metadata_all.pdf',
       plot = composite, 
       width = 14, 
       height = 8)

