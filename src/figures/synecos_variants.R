library(tidyverse)
library(speedyseq)
library(patchwork)


source(("src/util/backbone_functions.R") )
source(("src/util/backbone_params-graphs.R") )
source(("src/util/main_functions.R"))
source(("src/util/timeseries_functions.R"))

theme_set(theme_bw(base_size = 13))
geoMean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

bl.phy <- readRDS(('data/cleaned/blphyloseq.rds') ) %>% 
  transform_sample_counts(., function(x) x+1) %>%
  transform_sample_counts(function(x) log(x/geoMean(x)))

lombsea <- readRDS(('data/analysis/lombsea.rds'))
abunprev <- readRDS(('data/analysis/abund_prev.rds'))

tax <- as(tax_table(bl.phy), 'matrix') %>% as_tibble(rownames = 'asv')
env <-  sample_data(bl.phy)

synes <- prune_taxa(taxa = c('asv1','asv5'), bl.phy) %>% 
  psmelt() %>% 
  mutate(bloom.var = ifelse(Abundance >= 6 & OTU == 'asv5', TRUE, FALSE),
         bloom.count = ifelse(Synechococcus > 49000, TRUE, FALSE))

# plot.variants <- synes %>% 
#   ggplot(aes(Date, Abundance)) + 
#   geom_line() + 
#   geom_point(aes(color = bloom.count), show.legend = F) + 
#   ggrepel::geom_label_repel(data = . %>% filter(bloom.count),
#                             aes(label = as.character(Date))) + 
#   facet_wrap(~OTU, ncol = 1) + 
#   lil.strip + 
#   scale_color_manual(values = c('black','aquamarine4')) + 
#   ylab('Centered log ratio (log10)')

library(viridis)

syneco.counts <- synes %>% 
  select(Date, Synechococcus, Temperature, bloom.count) %>%
  distinct() %>% 
  ggplot(aes(Date, Synechococcus)) + 
  geom_hline(yintercept = 49000, linetype = 2, color = 'grey') + 
  geom_line(color = 'aquamarine4') + 
  geom_point(aes(color = Temperature), size = 3) + 
  ggrepel::geom_label_repel(data = . %>% filter(bloom.count),
                            aes(label = as.character(Date),
                                color = Temperature)) + 
  ylab(expression(italic(Synechococcus)~(cells~ml^{-1}))) + 
  scale_color_viridis_c(option = 'inferno', end = 0.9) 
  


# plot.comparison <- plot.variants + 
#   syneco.counts
# 
# ggsave('results/figures/syneco_comparison.pdf', width = 19, height = 9)


num.days.mnt <- c(0,31,28,31,30,31,30,31,31,30,31,30)
cumnum <- cumsum(num.days.mnt)

month.var <- ggplot(synes, aes(day_of_year, Abundance)) + 
  geom_point(data = . %>% filter(!bloom.count),
             show.legend = F, 
             size = 2,
             alpha = 0.6) + 
  geom_point(data = . %>% filter(bloom.count),
             aes(color = Temperature),
             show.legend = F, 
             size = 2.2) + 
  stat_smooth(aes(x = day_of_year,
                  group = OTU), 
              method = "gam",
              formula = y ~ s(x, k =12, bs = 'cc')) + 
  ggrepel::geom_label_repel(data = . %>% filter(bloom.count),
                            aes(label = as.character(Date),
                                color = Temperature),
                            size = 2.5, show.legend = F ) + 
  facet_wrap(~str_to_upper(OTU), ncol = 2) + 
  lil.strip + 
  scale_x_continuous(breaks = cumnum, labels = str_to_title(month.order),
                     name = 'Month') +
  scale_color_viridis_c(option = 'inferno', end = 0.9) + 
  ylab('Centered log ratio (log10( ASV read count / geoMean))')

plot.month <- month.var + 
  syneco.counts +
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = 'A')
  

ggsave('results/figures/syneco_comparison_month.pdf',
       plot = plot.month, width = 13, height = 9)
  
# infamous poster code 

month.var <- ggplot(synes, aes(day_of_year, Abundance)) + 
  geom_point(data = . %>% filter(!bloom.count),
             show.legend = F, 
             size = 2,
             alpha = 0.6) + 
  geom_point(data = . %>% filter(bloom.count),
             aes(color = Temperature),
             show.legend = F, 
             size = 2.2) + 
  stat_smooth(aes(x = day_of_year,
                  group = OTU), 
              method = "gam",
              formula = y ~ s(x, k =12, bs = 'cc')) + 
  ggrepel::geom_label_repel(data = . %>% filter(bloom.count),
                            aes(label = as.character(Date),
                                color = Temperature),
                            size = 3, show.legend = F ) + 
  facet_wrap(~str_to_upper(OTU), ncol = 1) + 
  lil.strip + 
  scale_x_continuous(breaks = cumnum, labels = str_to_title(month.order),
                     name = 'Month') +
  scale_color_viridis_c(option = 'inferno', end = 0.9) + 
  ggtitle('Niche partitioning between *Synechococcus* ASVs') + 
  theme(plot.title = ggtext::element_markdown()) + 
  ylab('Centered log ratio (log10( ASV read count / geoMean))')

plot.month.pos <- month.var + 
  (syneco.counts +
     theme(strip.text.x = element_blank(),
           strip.background = element_rect(colour="white", fill="white"),
           legend.position=c(.9,.75) 
     ) + 
     ggtitle('Flow cytometry counts',
             subtitle = 'The labelled dates correspond to blooming events')) + 
  plot_layout(ncol = 2,widths = c(0.4, 0.6)) +
  plot_annotation(tag_levels = 'A') 


ggsave('results/figures/syneco_comparison_month_poster.pdf',
       plot = plot.month.pos, width = 19, height = 9)
