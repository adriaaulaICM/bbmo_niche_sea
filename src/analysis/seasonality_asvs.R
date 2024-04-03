library(tidyverse)
library(phyloseq)
library(lubridate)
library(lomb)

source("src/util/params-graphs.R") 
source("src/util/main_functions.R")


bl.phy <- readRDS('data/cleaned/blphy_filt_7sams.rds') 

physeq.list <- data_transformation(bl.phy)

bl.phy.relab <- physeq.list$relab
bl.phy.raref <- physeq.list$raref
bl.phy.log <- physeq.list$log
bl.phy.clr <- physeq.list$clr
bl.phy.asinh <- transform_sample_counts(bl.phy, function(x) asinh(x))
physeq.list <- NULL


tsdf.0.2 <- bl.phy.asinh %>% 
  psmelt_dplyr() %>% 
  mutate(decimaldat = decimal_date(Date)) 

df.ts <- tsdf.0.2 %>% 
  arrange(Date) %>% 
  select(Date, decimaldat, Abundance, OTU) %>% 
  split(.$OTU, drop = T)

  
lomb.02  <-  df.ts %>% 
  map(~randlsp( x =.x$Abundance,
                times = .x$decimaldat,
                type = 'period',
                plot = F))

lomb.sea.02 <- tibble( asv = names(lomb.02),
                       pval = map_dbl(lomb.02, ~.x$p.value),
                       peak = map_dbl(lomb.02, ~.x$peak),
                       interval  = map(lomb.02, ~.x$peak.at),
                       int.min = map_dbl(interval, ~.[[2]]),
                       int.max = map_dbl(interval, ~.[[1]])) %>% 
  mutate( qval = fdrtool::fdrtool(pval, statistic = 'pvalue')$qval) %>% 
  filter(qval <= 0.01, peak >= 0.154, int.max <= 2)


write_rds(lomb.02, 'data/analysis/lomball.rds')
write_rds(lomb.sea.02, 'data/analysis/lombsea.rds')

results.lomb02 <- lomb.02[lomb.sea.02 %>% pull(asv)]

periodoplots <- map(results.lomb02, ~tibble( scanned = .x$scanned,
                          power = .x$power)) %>% 
  bind_rows(.id = 'asv') %>% 
  split(.$asv) %>% 
  map(~ggplot(.x, aes(scanned, power)) + 
        geom_line(aes(group = asv)) + 
        facet_wrap(~asv) + 
        lil.strip)
