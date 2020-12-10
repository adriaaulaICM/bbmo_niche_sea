library(tidyverse)
library(speedyseq)
library(patchwork)
library(ggtext)

# To clean 

source("src/util/backbone_functions.R") 
source("src/util/backbone_params-graphs.R") 
source("src/util/params-graphs.R")
source("src/util/main_functions.R")
source("src/util/timeseries_functions.R")

figpath <- 'results/figures/timeseries/'
dir.create(file.path(figpath), showWarnings = FALSE)


theme_set(theme_classic(base_size = 25))

bl.phy <- readRDS('data/cleaned/blphyloseq.rds') 
lomb.02 <- readRDS('data/analysis/lomball.rds')
lombsea <- readRDS('data/analysis/lombsea.rds')

lomb.sea.02 <- tibble( asv = names(lomb.02),
                       pval = map_dbl(lomb.02, ~.x$p.value),
                       peak = map_dbl(lomb.02, ~.x$peak),
                       interval  = map(lomb.02, ~.x$peak.at),
                       int.min = map_dbl(interval, ~.[[2]]),
                       int.max = map_dbl(interval, ~.[[1]])) %>% 
  mutate( qval = fdrtool::fdrtool(pval, statistic = 'pvalue')$qval) 

asv.sel <- str_c('asv', c(1,9))

geoMean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

bl.phy.relab <- transform_sample_counts(bl.phy, function(x) x / sum(x))
bl.phy.clr <- transform_sample_counts(bl.phy, function(x) x+0.01) %>%
  transform_sample_counts(function(x) log(x/geoMean(x)))

psmelt.clr <- psmelt(bl.phy.relab) %>% as_tibble() %>% 
  filter(OTU %in% asv.sel)

gam.gg <- ggplot(data = psmelt.clr, aes(day_of_year,Abundance)) + 
  geom_jitter(aes(color = OTU), alpha = 0.4) + 
  stat_smooth(aes(x = day_of_year,
                  group = OTU,
                  color = OTU),
              method = "gam",
              formula = y ~ s(x, k =12, bs = 'cc'),
              se = F, size = 2,
              show.legend = F, alpha = 0.7) + 
  facet_wrap(~str_c(str_to_upper(OTU), genus, sep = ', '), scales = 'free_y') + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 2L)) +
  scale_x_continuous(breaks = cumnum,
                     name = 'Month',
                     labels = str_to_title(month.order) %>% str_sub(1,1)) +
  guides(color = FALSE) + 
  ylab('Relative abundance') + 
  lil.strip + 
  leg.bottom 


periodo.gg <- map(lomb.02, ~tibble( scanned = .x$scanned,
                                    power = .x$power)) %>%
  bind_rows(.id = 'asv') %>% 
  filter(asv %in% asv.sel) %>% 
  mutate( desc = ifelse(asv == 'asv1', 'Strong peak with a periodicity of 1 year',
                        'No patterns, random walk') %>% as.factor() %>% fct_inorder()) %>% 
  ggplot(aes(scanned, power)) + 
  geom_line( aes(group = asv)) + 
  ylab('Strength recurrence') + 
  xlab('Periods checked (0.083 ~ 1/12 months)') + 
  facet_wrap(~desc) + lil.strip

gam.gg + periodo.gg + plot_layout(ncol = 1)

ggsave(filename = str_c(figpath, 'example_seasonality.pdf'),
       width = 11, height = 7)
# lomb.sea.02 %>% 
#   mutate( index = str_replace(asv, 'asv', '') %>% as.integer()) %>% 
#   arrange(index) %>% View()