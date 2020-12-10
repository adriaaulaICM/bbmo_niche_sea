library(tidyverse) 
library(speedyseq)
library(mgcv)
library(broom)

source("src/util/backbone_functions.R") 
source("src/util/backbone_params-graphs.R") 

geoMean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Inside this script i have the objects and selection of vars
source('src/analysis/helper_prepare_nichepart_data.R')

genus.interest <-  psmelt(bl.phy.gen.filt) %>% 
  as_tibble() %>% 
  pull(genus) %>% 
  as.character() %>% 
  unique()

geomean_calculation_intragenus <- function(genus.sel){
  
  # Inside each genus calculate the geomean and 
  # divide the values to obtain a CLR measurement
  asv.sel <- tax %>% 
    filter(genus == genus.sel) %>% 
    pull(asv)
  
  bl.phy.genfilt.geomean <- prune_taxa(taxa = asv.sel,
                                       x = bl.phy)
  
  bl.phy.pseudo <- transform_sample_counts(bl.phy.genfilt.geomean,
                                          function(x) x + 0.01)
  bl.phy.check <- transform_sample_counts(bl.phy.pseudo,
                                          function(x) log(x/geoMean(x)))
  
  return(psmelt(bl.phy.check))
  
}

psmelt.clr.intragenus <- genus.interest %>% 
  map(~geomean_calculation_intragenus(.x)) %>% 
  bind_rows() %>% 
  as_tibble()


models.par <- psmelt.clr.intragenus  %>% 
  filter(OTU %in% taxa_names(bl.phy.gen.filt)) %>% 
  filter(Sample_name %in% sample_names(bl.phy.gen.filt)) %>% 
  select(OTU, Abundance, Date, Temperature, Chla_total,
         PO4, NH4, NO2, NO3, Si,
         PNF_Micro, HNF_Micro)  %>% 
  pivot_longer(names_to = 'parameter',
               values_to = 'value',
               cols = -c(OTU, Abundance, Date)) %>% 
  group_by(OTU,parameter) %>% 
  nest() %>% 
  mutate(model = map(data, ~gam(Abundance ~ s(value, k = 2),
                                data = .x, method = 'REML')),
         results = map(model, glance),
         pval.term = map_dbl(model, ~ summary(.)$s.pv) ,
         R.square = map_dbl(model, ~ summary(.)$r.sq)) 

# how many are close to linearlity 
edfs <- models.par %>% 
  inner_join(estimates.all, by = c('OTU' = 'asv', 'parameter' = 'par')) %>% 
  pull(model) %>% 
  map(~summary(.x)) %>% 
  map_dbl(~.x$edf)

table(edfs >= 1.5)

results <- models.par %>% 
  select(OTU, parameter, results, pval.term, R.square) %>% 
  unnest(results)

original.data <- models.par %>% 
  left_join(tax %>% select(asv, genus),
            by = c('OTU' = 'asv')) %>% 
  select(OTU, genus,parameter, data) %>% 
  unnest(data) %>% 
  ungroup() %>% 
  as_tibble()

write_rds(results, 'data/analysis/gams_envresponse_res.rds')
write_rds(original.data, 'data/analysis/gams_originaldata.rds')


# Try it with all the variables at the same time  -------------------------
# models.par <- psmelt.clr.intragenus  %>% 
#   filter(OTU %in% taxa_names(bl.phy.gen.filt)) %>% 
#   filter(Sample_name %in% sample_names(dat)) %>% 
#   select(OTU, Abundance, Date, Temperature, Chla_total,
#          PO4, NH4, NO2, NO3, Si,
#          PNF_Micro, HNF_Micro)  %>% 
#   group_by(OTU)  %>% 
#   nest() %>% 
#   mutate(model = map(data, ~gam(Abundance ~ s(Temperature, k = 2) +
#                                   s(Chla_total, k = 2) +
#                                   s(PO4, k = 2) +
#                                   s(NH4, k = 2) +
#                                   s(NO2, k = 2) +
#                                   s(NO3, k = 2) +
#                                   s(Si, k = 2) +
#                                   s(PNF_Micro, k = 2) +
#                                   s(HNF_Micro, k = 2), 
#                                 data = .x, 
#                                 method = 'REML'))) %>% 
#   mutate(results =map(model, glance),
#          pval.term = map(model, ~ summary(.)$s.pv) ,
#          R.square = map_dbl(model, ~ summary(.)$r.sq)) 
# 
# results.gam <- models.par %>% 
#   select(OTU, results, pval.term, R.square) %>% 
#   unnest(c(results)) %>%
#   unnest(c(pval.term)) %>% 
#   ungroup() %>% 
#   mutate( par = rep(c('Temperature', 'Chla_total', 'PO4',
#                       'NH4', 'NO2', 'NO3', 'Si', 'PNF_Micro',
#                       'HNF_Micro'), 47)) %>% 
#   filter(pval.term <= 0.01) %>% 
#   select(OTU, par)  %>% 
#   rename('asv'= OTU)  %>% 
#   mutate( whole = str_c(asv,par))


# Old ---------------------------------------------------------------------

#plot_gam_response <- function(genus.sel){
#  sm <-  smooths %>% dplyr::filter(genus == genus.sel)
#  org.data <- original.data %>% dplyr::filter(genus == genus.sel)
#  
#  plot <- ggplot(sm, aes(value,est)) + 
#    geom_point(data = org.data, 
#               aes(x =value, y = Abundance, color = OTU),
#               alpha = 0.7) + 
#    geom_ribbon(aes(ymin = est - se,
#                    ymax = est + se, 
#                    group = OTU), alpha = 0.6) + 
#    geom_line(aes(group = OTU, color = OTU) ) + 
#    lil.strip + 
#    facet_wrap(~parameter, scales = 'free')
#  
#  return(plot)
#  
#}
#plots.gams <- genus.interest %>% 
#  map(~plot_gam_response(.x))
#
#names(plots.gams) <- genus.interest
#
#plots.gams$Synechococcus_C
#plots.gams$HIMB59
#plots.gams$Pelagibacter_A
#plots.gams$Pelagibacter
#
#

#smooths <- models.par %>% 
#  filter(pval.term <= 0.001) %>% 
#  mutate(dataforgg = map(model, ~evaluate_smooth(.x, 'value'))) %>% 
#  left_join(tax %>% select(asv, genus), by = c('OTU' = 'asv')) %>% 
#  select(OTU,genus,data, parameter, dataforgg) %>% 
#  unnest(c(dataforgg)) %>% 
#  ungroup() %>% 
#  as_tibble()