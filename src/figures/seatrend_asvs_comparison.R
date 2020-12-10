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

figpath <- 'results/figures/timeseries/niche_part'
dir.create(file.path(figpath), showWarnings = FALSE)

bl.phy <- readRDS('data/cleaned/blphyloseq.rds') 
lomb.02 <- readRDS('data/analysis/lomball.rds')
lombsea <- readRDS('data/analysis/lombsea.rds')


geoMean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

bl.phy.clr <- transform_sample_counts(bl.phy, function(x) x+0.01) %>%
  transform_sample_counts(function(x) log(x/geoMean(x)))

otus99 <- read_tsv('data/cleaned/dada2/04_clustering/bbmo_99correspondence.tsv') %>% 
  rename( 'otu.corr99' = OTU, 'asv' = ASV)  %>% 
  select(asv, otu.corr99)

# we will only work with the seasonal ones
bl.phy <- prune_taxa(lombsea$asv, bl.phy)
hammdists <- read_rds('data/analysis/hammdist.rds')

bl.phy.asinh <- transform_sample_counts(bl.phy, function(x) asinh(x))
psmelt.asinh <- psmelt(bl.phy.asinh) %>% as_tibble()

bl.phy.clr <- prune_taxa(lombsea$asv, bl.phy.clr)
psmelt.clr <- psmelt(bl.phy.clr) %>% as_tibble()

# Number of OTUs presenting more than 1 ASVs
otus99 <- otus99[otus99$asv %in% taxa_names(bl.phy),] %>% 
  group_by(otu.corr99) %>% 
    filter(n() >= 2) %>% 
  ungroup()

# Selection of the genus presenting ASVs with divergence in the distribution
genus.part <- c("AG-337-I02", "TMED189", "MS024-2A", "HIMB59",
                "UBA10364", "Pelagibacter", "Pelagibacter_A","Synechococcus_C",
                "Puniceispirillum", "SAR86A", "Luminiphilus", "Marinisoma")

dataset <- psmelt.clr %>% 
  filter(genus %in% genus.part) %>%
  left_join(otus99, by = c('OTU' = 'asv'))
  

listdat <- dataset %>% 
  split(.$otu.corr99)

num.days.mnt <- c(0,31,28,31,30,31,30,31,31,30,31,30)
cumnum <- cumsum(num.days.mnt)


plots.part <- listdat %>% 
  map(~ ggplot(data = .x, aes(day_of_year,Abundance)) + 
        geom_jitter(aes(color = OTU), alpha = 0.4) + 
        stat_smooth(aes(x = day_of_year,
                        group = OTU,
                        color = OTU),
                    method = "gam",
                    formula = y ~ s(x, k =12, bs = 'cc'),
              se = F, size = 2,
              show.legend = F, alpha = 0.7) + 
        facet_wrap(~str_c(phylum, order, family, genus, otu.corr99, sep = ', ')) + 
        scale_x_continuous(breaks = cumnum,
                           name = 'Month',
                           labels = str_to_title(month.order)) +
        lil.strip + 
        leg.bottom 
  )

filename <- str_c(figpath, '/', names(listdat), '.pdf')

map2(filename, plots.part, ggsave, width = 9, height = 9)


# From looking at the different distributions we will only plot some of them 
# Taking the best ones 
otu.interest <- str_c('otu', c(1,27,30,46, 86,243, 62,90,243))
asvs.interest <- otus99 %>% filter(otu.corr99 %in% otu.interest) %>% pull(asv)


plot_composite <- function(dataset){
  
  distribution <- ggplot(data = dataset, aes(day_of_year, Abundance)) + 
    geom_jitter(aes(color = OTU),
                alpha = 0.6,
                show.legend = F, 
                width = 0.1) + 
    stat_smooth(aes(x = day_of_year,
                    group = OTU,
                    color = OTU),
                method = "gam",
                formula = y ~ s(x, k =12, bs = 'cc'),
                se = F, size = 2,
                show.legend = F, alpha = 0.7) + 
    scale_color_manual(values = asvvarpal, labels = asvvarlabels) + 
    scale_x_continuous(breaks = cumnum,
                       name = 'Month',
                       labels = str_to_title(month.order)) +
    facet_wrap(~str_c(phylum, order, family, genus, otu.corr99, sep = ', ')) + 
    lil.strip + 
    ylab('Centered log ratio (log10( ASV read count / geoMean))') +
    ggtitle('Seasonal trends') + 
    theme(legend.text = element_markdown())
  
  
  hammplot <-   hammdists   %>% 
    filter(asv %in% dataset$OTU & asv2 %in% dataset$OTU  ) %>% 
    ggplot(aes(asv,asv2)) +
    geom_tile(aes(fill = hammdist), show.legend = FALSE) + 
    geom_text(aes(label = hammdist), color = 'white', size = 4) + 
    xlab(NULL) + ylab(NULL) + 
    theme_minimal() +        
    scale_fill_continuous(high = "#132B43", low = "#56B1F7") + 
    scale_x_discrete(expand = c(0.1,0.1), labels = asvvarlabels) + 
    scale_y_discrete(expand = c(0.1,0.1), labels = asvvarlabels) + 
    theme(axis.text.x = element_markdown(),
          axis.text.y = element_markdown()) + 
    ggtitle('Nucleotide divergencies') + 
    coord_equal() 
  
  distribution +
    hammplot +
    plot_layout(widths = c(0.6, 0.4)) +
    plot_annotation(tag_levels = 'A')
  
}

compositedat <- listdat[c('otu1', 'otu30', 'otu46', 'otu243')] %>% 
  bind_rows(.id = 'otu') %>% 
  filter(!(otu %in% 'otu243') | OTU %in% str_c('asv', c(285,337,243, 422))) %>% 
  tax_to_italics()  %>% 
  mutate(wrap_col = str_c(phylum, order, family.ital,
                          genus.ital, str_to_upper(otu), sep = ', ') %>% 
           as.factor() %>% 
           forcats::fct_rev())


varnumbering <- compositedat %>% 
  distinct(OTU, .keep_all = T) %>% 
  group_by(otu) %>% 
  mutate(varnum = str_c('Var. ',1:n() )) %>% 
  ungroup() %>% 
  select(OTU, varnum)

varpalette <- tibble( varnum = c("Var. 1", "Var. 2", "Var. 3",
                                     "Var. 4", "Var. 5", "Var. 6",
                                     "Var. 7", "Var. 8"),
                          color = c('#50514F', '#F25F5C', 'darkgoldenrod1','#70C1B3',
                                    '#F2D7EE', '#A5668B',  '#247BA0', '#0E103D'))


var.df <- varnumbering %>% 
  left_join(varpalette, by = 'varnum') %>% 
  mutate( label = str_c("<b style='color:", color,"'>", OTU, "</b>"))

asvvarpal <- var.df$color
names(asvvarpal) <- var.df$OTU

asvvarlabels <- var.df$label %>% str_replace(pattern = 'asv', replacement = 'ASV')
names(asvvarlabels) <- var.df$OTU 

nichepart.gg <- compositedat %>% 
  ggplot(aes(day_of_year, Abundance)) + 
  geom_jitter(aes(color = OTU),
              alpha = 0.6,
              show.legend = F, 
              width = 0.1) + 
  stat_smooth(aes(x = day_of_year,
                  group = OTU,
                  color = OTU),
              method = "gam",
              formula = y ~ s(x, k =12, bs = 'cc'),
              se = F, size = 2,
              show.legend = F, alpha = 0.7) + 
  scale_color_manual(values = asvvarpal, labels = asvvarlabels) + 
  scale_x_continuous(breaks = cumnum,
                     name = 'Month',
                     labels = str_to_title(month.order)) +
  facet_wrap(~wrap_col, ncol =1, scales = 'free_y' ) + 
  lil.strip + 
  ylab('Centered log ratio (log10( ASV read count / geoMean))') +
  theme(legend.text = element_markdown(),
        strip.text.x = element_markdown())

hammplot <-   hammdists   %>% 
  filter(asv %in% compositedat$OTU & asv2 %in% compositedat$OTU  ) %>% 
  left_join(otus99, by = 'asv') %>% 
  mutate(otu.corr99 = factor(otu.corr99,
                             levels = c('otu30', 'otu243', 'otu1', 'otu46'))) %>% 
  ggplot(aes(asv,asv2)) +
  geom_tile(aes(fill = hammdist), show.legend = FALSE) + 
  geom_text(aes(label = hammdist), color = 'white', size = 4) + 
  xlab(NULL) + ylab(NULL) + 
  theme_minimal() +        
  scale_fill_continuous(high = "#132B43", low = "#56B1F7") + 
  scale_x_discrete(expand = c(0.1,0.1), labels = asvvarlabels) + 
  scale_y_discrete(expand = c(0.1,0.1), labels = asvvarlabels) + 
  facet_wrap(~otu.corr99, ncol =1, scales = 'free')  +
  theme(panel.spacing = unit(2, 'lines'),
    strip.text.x = element_blank(),
    axis.text.x = element_markdown(),
    axis.text.y = element_markdown())

compo <- nichepart.gg +
  hammplot + 
  plot_layout(widths = c(0.7, 0.3)) + 
  plot_annotation(tag_levels = 'A')

ggsave(plot = compo, str_c(figpath, '/nichepart_new.pdf'),
       width = 10, height = 13)

composites <- listdat[otu.interest] %>% 
  map(~plot_composite(.x))

names(composites) <- otu.interest

filenames <- str_c(figpath, '/', 'nichepart_', otu.interest, '.pdf')

map2( filenames, composites, ggsave, width = 11, height = 6)
