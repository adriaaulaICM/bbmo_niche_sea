library(tidyverse)
library(speedyseq)
library(patchwork)
library(corncob)
library(gt)
library(ggthemes)
library(ggtext)
library(ggtree)

source("src/util/backbone_functions.R") 
source("src/util/backbone_params-graphs.R") 
source('src/analysis/helper_prepare_nichepart_data.R')

# Metadata labels 

labels <- readxl::read_xlsx('data/raw/metadata/labels_pretty.xlsx',
                            sheet = 2,
                            col_names = c('Variable', 'Measure'))  %>% 
  mutate(all = str_c(Variable, "~", Measure ))

par.order <-  c("Day_length", "Temperature", "Salinity", "Secchi",
    "NO2", "NO3", "PO4", "Si", "Chla_3um","Chla_total", "BP_FC1.55",
    "Bacteria_joint", "HNA", "LNA", "Synechococcus", "Prochlorococcus_FC",
    "Peuk1", "Peuk2", "PNF_Micro", "PNF2_5um_Micro", "PNF_5um_Micro",
    "Cryptomonas", "Micromonas", "HNF_Micro")

# The ASVs with significant responses 
estimates.all <- readRDS('data/analysis/difftest_par.rds')
estimates.par <-  estimates.all %>% 
  # only the genus with enoughs ASVs to compare
  filter(genus %in% c('AG-337-I02', 'D2472', 'Luminiphilus'))   %>% 
  filter(!par %in% c('Si'))
  # and only the parameters that show both positive and negative trends
  # inside the same genus 
  # filter(par %in% c('Temperature', 'Chla_total', 'NO3', 'PO4', 'NO2'))

# ASV Label creation ----------------------------------------------------------
varnumbering <- psmelt(bl.phy.gen.filt) %>% 
  group_by(genus,OTU) %>% 
  summarise(total = sum(Abundance)) %>% 
  arrange(-total) %>% 
  mutate(varnum = str_c('Var. ',1:n() )) %>% 
  ungroup() %>% 
  select(OTU, varnum)

varpalette <- tibble( varnum = c("Var. 1", "Var. 2", "Var. 3",
                                     "Var. 4", "Var. 5", "Var. 6",
                                     "Var. 7", "Var. 8"),
                          color = c('#50514F', '#F25F5C', '#FFE066',
                                    '#247BA0', '#70C1B3', '#F2D7EE',
                                    '#A5668B', '#0E103D'))

var.df <- varnumbering %>% 
  left_join(varpalette, by = 'varnum') %>% 
  mutate( label = str_c("<b style='color:", color,"'>", str_to_upper(OTU), "</b>"),
          label = ifelse(is.na(color),
                         str_c("<b style='color:", '#696969' ,
                               "'>", str_to_upper(OTU), "</b>"),
                         label),
          color = ifelse(is.na(color),'#696969', color))

asvvarpal <- var.df$color
names(asvvarpal) <- var.df$OTU

asvvarlabels <- var.df$label
names(asvvarlabels) <- var.df$OTU

# Hclust nucdistances ------------------------------------------------------------------

cluster_asvs <- function(genus.asv){
  thesel <- estimates.par  %>% 
    filter(genus %in% genus.asv) %>%
    pull(asv) %>% 
    unique()
  
  hamm.matrix <- hamm %>% 
    filter(asv %in% thesel & asv2 %in% thesel) %>% 
    mutate(asv = factor(asv), 
           asv2 = factor(asv2)) %>% 
    arrange(asv,asv2) %>% 
    select(-genus) %>% 
    pivot_wider(names_from = asv2,
                values_from = hammdist,
                values_fill = list(hammdist = 30))
  
  
  thedat <- hamm.matrix %>% column_to_rownames(var = 'asv') %>% as.matrix()
  thedat <- thedat[ levels(hamm.matrix$asv), levels(hamm.matrix$asv)]
  thedat.dist <- thedat %>% as.dist()
  
  return(hclust(thedat.dist))
  
}

gen_to_cluster <- estimates.par$genus %>% unique() %>% sort()


asv.order.fig <- gen_to_cluster %>% 
  map(~cluster_asvs(.x)) %>% 
  map(~ .x$labels[.x$order]) %>% 
  unlist()


trees <- gen_to_cluster %>% 
  map(~cluster_asvs(.x)) %>% 
  map(~ggtree(.x)) 

filenames.trees <- str_c('results/figures/env_effects/nucdiv_',
                   gen_to_cluster,
                   '.pdf')

map2(filenames.trees, trees, ggsave, width = 2, height =3 )

trees %>% 
  map(~.x + geom_tiplab())

ggsave(filenames.trees[2], trees[[2]], width = 2, height = 8)

# Drawing corncob ---------------------------------------------------------
estimates.par <- estimates.par %>% 
  left_join(varnumbering, by = c('asv' = 'OTU'))

plotres <-  estimates.par %>% 
  mutate( par = factor(par, levels = par.order,
                             labels = labels$all)) %>% 
  arrange(genus, par, -estimate) %>% 
  mutate(asv = factor(asv, levels = asv.order.fig)) %>% 
  mutate(genus = ifelse(genus %in% 'Luminiphilus',
                        'italic(Luminiphilus)',
                        genus)) %>% 
  ggplot( aes(x = estimate,
              y =str_to_upper(asv),
              color = asv)) + 
  geom_vline(xintercept = 0, linetype = 2, color = 'red') + 
  geom_pointrange(aes(xmin = lower, xmax = upper),
                  show.legend = F) + 
  facet_grid(genus~par,
             labeller = label_parsed,
             scales = 'free', space = 'free_y') + 
  theme_minimal(base_size = 15) + 
  # theme(strip.text.y = element_markdown()) + 
  xlab('Coefficient estimate') + 
  ylab('ASV') + 
  scale_x_continuous( guide = guide_axis(n.dodge = 2)) + 
  scale_color_manual(values = asvvarpal) + 
  theme(legend.position = 'bottom') 

# and select what we have to decided to choose the GAMs
ggsave(plot = plotres,
       'results/figures/env_effects/res_corncob.pdf',
       width = 12, height = 10)

selected <- estimates.par %>% 
  select(asv, par) %>% 
  rename('OTU' = asv, 'parameter' = par) 

# Draw GAMs ---------------------------------------------------------------
results.gam <- readRDS('data/analysis/gams_envresponse_res.rds') %>% 
  inner_join(selected, by = c('OTU', 'parameter'))
original.data <- readRDS('data/analysis/gams_originaldata.rds') %>% 
  inner_join(selected, by = c('OTU', 'parameter'))


# dont draw the ones in Pelagibacter, too noise! 
notgrey <- var.df  %>% 
  filter(color != '#696969') %>% 
  pull(OTU)

thenewdatahowpretty <- original.data  %>% 
  mutate( parameter = factor(parameter, levels = par.order,
                             labels = labels$all)) %>% 
  filter(OTU %in% notgrey) %>% 
  left_join(varnumbering, by = 'OTU')  %>% 
  rename('asv' = OTU)

gamsplot <- thenewdatahowpretty %>% 
  mutate(genus = ifelse(genus %in% 'Luminiphilus',
                        'italic(Luminiphilus)',
                        genus)) %>% 
  ggplot(aes(value,Abundance)) + 
  # geom_point(aes(color = asv), alpha = 0.3, show.legend = F) + 
  geom_smooth(method = 'gam',
              aes(color = asv,
                  fill = asv),
              formula = y ~ s(x, k = 1), show.legend = F) +
  facet_grid(genus ~ parameter,
             labeller = label_parsed,
             scales = 'free', switch = 'x') + 
  theme_bw(base_size = 15) + 
  scale_color_manual(values = asvvarpal) + 
  scale_fill_manual(values = asvvarpal) + 
  scale_x_continuous( guide = guide_axis(n.dodge = 2)) + 
  ylab('Centered logarithm ratio ( ASV read count / geoMean)') + 
  xlab('Environmental parameter value') + 
  lil.strip


# Joining the three worlds  -------------------------------------------------

alottotakein <- plotres + 
  gamsplot +
  plot_layout(ncol = 1)


ggsave('results/figures/env_effects/corncob_plus_gams_envresponse_othergroups.pdf',
       alottotakein,
       width = 13,
       height = 16)