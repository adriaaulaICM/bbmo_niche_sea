library(tidyverse)
library(geofacet)
library(speedyseq)
library(lubridate)
library(patchwork)
library(ggtext)

fig.path <- 'results/figures/timeseries/cohesivity/'
dir.create(file.path(fig.path), showWarnings = FALSE)

source("src/util/params-graphs.R") 
source("src/util/main_functions.R")


valtax_to_italics <- function(df){
  df %>% 
    mutate(value = case_when( value %in% c(fams.italy,
                                           genus.italy) ~ str_c('*', value, '*'),
                              TRUE ~  as.character(value)))
  }

gridtax_to_italics <- function(df){
  df %>% 
    mutate(code = case_when( code %in% c(fams.italy,
                                         genus.italy) ~ str_c('*', code, '*'),
                             TRUE ~  as.character(code)),
           name = case_when( name %in% c(fams.italy,
                                         genus.italy) ~ str_c('*', name, '*'),
                             TRUE ~  as.character(name))) 
}


bl.phy <- readRDS('data/cleaned/blphy_filt_7sams.rds') 
bl.phy <- transform_sample_counts(bl.phy, function(x) x / sum(x))
lombsea <- read_rds('data/analysis/lombsea.rds')

grid.alpha <- readxl::read_xlsx('data/raw/cohesivity_grid.xlsx', sheet = 1) %>% 
  gridtax_to_italics()
grid.gamma <- readxl::read_xlsx('data/raw/cohesivity_grid.xlsx', sheet = 2) %>% 
  gridtax_to_italics()
grid.flavo <- readxl::read_xlsx('data/raw/cohesivity_grid.xlsx', sheet = 3) %>% 
  gridtax_to_italics()
  
tax <- as(tax_table(bl.phy), 'matrix') %>% as_tibble(rownames = 'asv')
  
dist.raw <- readRDS('data/analysis/cohesivity_dist.rds') 

distributions <- dist.raw %>% 
  ungroup() %>%
  unnest(cols = values) %>% 
#Some changes in the nomenclature since we have ids repeated 
  mutate( value = case_when(
    rank == 'family' & value == 'D2472' ~ 'D2472_f',
    rank == 'genus' & value == 'D2472' ~ 'D2472_g',
    TRUE ~ value)) %>% 
  valtax_to_italics()


corr <- tax %>% 
  select(class, order, family, genus)  %>% 
  filter(class %in% c('Alphaproteobacteria',
                      'Gammaproteobacteria',
                      'Bacteroidia')) %>% 
  pivot_longer(names_to = 'level', values_to = 'value', cols = -class) %>% 
  mutate( value = case_when(
    level == 'family' & value == 'D2472' ~ 'D2472_f',
    level == 'genus' & value == 'D2472' ~ 'D2472_g',
    TRUE ~ value)) %>% 
  select(-level)

# All!  -------------------------------------------------------------------

# dist.hist <- distributions %>% 
#   mutate(rank = factor(rank, levels = c('class', 'order',
#                                         'family', 'genus')),
#          significance = ifelse(pval < 0.01 & peak >= 10 & int.max <= 2,
#                                'seasonal', 'non seasonal')) %>% 
#   left_join(corr, by = 'value') %>% 
#   mutate(class = ifelse(is.na(class), value, class))
# 
# library(ggridges)
# 
# seas.label <- data.frame(label = c('', 'seasonal'),
#                          y = c('class', 'class'),
#                          x = c(7.3,16.9))
#   
# global_res <- ggplot(dist.hist %>%
#                        mutate(class = factor(class,
#                                              levels = c('Alphaproteobacteria',
#                                                         'Gammaproteobacteria',
#                                                         'Bacteroidia'))),
#                      aes(x = peak,
#                          y = rank,
#                          fill = rank)) + 
#   geom_vline(xintercept = 10, linetype = 2, color = 'tomato3') + 
#   geom_density_ridges(show.legend = F) + 
#   geom_text(data = seas.label,
#             aes(label = label, x = x, y = y),
#             inherit.aes = FALSE,
#             size = 3,
#             nudge_y = -0.3, color = 'tomato3') + 
#   facet_wrap(~class) + 
#   scale_y_discrete(label = str_to_title) + 
#   scale_fill_ptol() + 
#   lil.strip + 
#   xlab('Peak power of recurrence (strength of signal)') + 
#   ylab('Rank level') 
#   
# ggsave(
#   filename = 'results/figures/timeseries/cohesivity/global_v2.pdf',
#   plot = global_res, 
#   width = 8,
#   height = 6)
# 
# global <- dist.hist %>% 
#   ggplot( aes(peak)) + 
#   geom_histogram(aes(fill = rank)) +
#   facet_wrap(~significance, scales = 'free_y', ncol = 1) + 
#   theme_bw() + 
#   lil.strip + 
#   scale_fill_ptol() + 
#   xlab('Peak power of recurrence (strength of signal)') + 
#   ylab('Count of statistic') + 
#   leg.bottom  
# 
# ggsave(
#   filename = 'results/figures/timeseries/cohesivity/global_distr_diffranks.pdf',
#   plot = global, 
#   width = 9,
#   height = 9)
#        
# 
# # Some numbers of that 
# dist.hist %>% 
#   group_by(rank, significance) %>% 
#   summarize(totals = n()) %>% 
#   group_by(rank) %>% 
#   mutate(relab = totals / sum(totals))




without20max <- readRDS('data/analysis/cohesivity_max.rds') %>% 
  ungroup() %>% 
  unnest(cols = values)

normal <- readRDS('data/analysis/cohesivity_normal.rds') %>% 
  ungroup() %>% 
  unnest(cols = values) %>% 
#Some changes in the nomenclature since we have ids repeated 
  mutate( value = case_when(
    rank == 'family' & value == 'D2472' ~ 'D2472_f',
    rank == 'genus' & value == 'D2472' ~ 'D2472_g',
    TRUE ~ value))  %>% 
  valtax_to_italics()


# to save a plot from the lombscargle results 
plot_lombstat <- function(grid.sel, limits){
  
  seldist <- distributions %>% 
    filter(value %in% grid.sel$name)
  
  selnormal <- normal %>% 
    filter(value %in% grid.sel$name)
  
  
  ggplot(seldist, aes(peak)) +
    geom_vline(data = selnormal, aes(xintercept = peak),
               color = 'green', linetype = 2) +
    geom_histogram(aes(fill = pval < 0.01 & peak >= 10 & int.max <= 2)) +
    geom_text(data = selnormal,
              aes(label = str_c('N ASVs = ', nasvs)),
              x=Inf, y = Inf, vjust=1.3, hjust=1.4) + 
    facet_geo(~value, grid = grid.sel, label = 'name') + 
    scale_fill_discrete(name = 'Significance', limits = c(FALSE,TRUE)) + 
    scale_x_continuous(limits = limits) + 
    theme_bw() + 
    lil.strip + 
    xlab('Peak power of recurrence (strength of signal)') + 
    ylab('Count of statistic') + 
    theme(legend.position = c(0.80, 0.9),
          legend.direction = 'horizontal',
          strip.text.x = element_markdown())
  
}

# and the code for the trends 
prepare_dataset_trend <- function(physeq){
  
  tsdf.0.2 <- physeq %>%
    psmelt_dplyr() %>%
    mutate(decimaldat = decimal_date(Date))
  
  df.ts <- tsdf.0.2 %>%
    arrange(Date) %>%
    select(Date, decimaldat,month, Abundance, OTU)
  
  return(df.ts)
}


aggregate_asvs <- function(physeq, keepasvs){
  bl.phy.rem <- prune_taxa(taxa = keepasvs, physeq)
  bl.agg <- tax_glom(bl.phy.rem, taxrank = 'class')
  # bl.agg.asinh <- transform_sample_counts(bl.agg, function(x) asinh(x))
  
  prepare_dataset_trend(bl.agg) %>% 
    return()
}


trend_rank_level <- function(rank.sel, rank.level){
  
  tax.raw <- as(tax_table(bl.phy), 'matrix') %>% as_tibble(rownames = 'asv')
  taxa <- tax.raw %>% 
    pivot_longer(names_to = 'rank', values_to = 'value',  cols = -asv)
  tax.rank.sel <- taxa %>% 
    filter( value == rank.sel, rank == rank.level)
  
  #rank.pos <- tax.rank.sel %>% pull(rank) %>% unique() %>% which(ranks == .)
  
  asv.sel <- tax.rank.sel %>% pull(asv) %>% unique()
  bl.phy.rank <- prune_taxa(taxa = asv.sel, x = bl.phy)
  
  iterations <- 1:iter %>% 
    map(~aggregate_asvs(bl.phy.rank,
                        sample(taxa_names(bl.phy.rank),
                               size = ntaxa(bl.phy.rank) * 0.8 ) 
    )) %>% 
    bind_rows(.id = 'iteration') %>% 
    as_tibble()
  
  return(iterations)
}

iter <- 10

trends.df <- dist.raw %>% 
  select(rank, value)  %>% 
  mutate( values = map2(.x = value, .y = rank,
                        ~trend_rank_level(.x, .y)))  %>% 
  ungroup()  %>%  
  unnest(c(values)) %>% 
#Some changes in the nomenclature since we have ids repeated 
  mutate( value = case_when(
    rank == 'family' & value == 'D2472' ~ 'D2472_f',
    rank == 'genus' & value == 'D2472' ~ 'D2472_g',
    TRUE ~ value))  %>% 
  valtax_to_italics()


plot_trend <- function(grid.sel){
  
  seltrends <- trends.df %>% 
    filter(value %in% grid.sel$name)
  
  ggplot(seltrends, aes(x = month, y = Abundance)) + 
    geom_boxplot(outlier.colour = 'transparent') + 
    geom_smooth(aes(x = month,
                    group = iteration,
                    color = iteration),
                method = "lm",
                formula = y ~ poly(x, 3),
                se = F, size = 0.7,
                show.legend = F, alpha = 0.5) +
    facet_geo(~value, scales = 'free_y', grid = grid.sel, label = 'name') + 
    scale_x_discrete(name = 'Month',
                    labels = str_to_title(month.order) %>% str_sub(1,1)) + 
    scale_y_continuous(labels = scales::percent,
                       name = 'Relative abundance') +
    lil.strip + 
    theme(strip.text.x = element_markdown())
  
  
}



# Alphaproteobacteria -----------------------------------------------------

lombalpha <- plot_lombstat(grid.sel = grid.alpha, 
              limits =  c(0,37))

ggsave(str_c(fig.path, 'lombstat_alphaproteobacteria_new', '.pdf'),
       plot = lombalpha,
       width = 9, height = 8)

trendalpha <- plot_trend(grid.sel = grid.alpha)

ggsave(str_c(fig.path, 'trend_alphaproteobacteria_new', '.pdf'),
       plot = trendalpha,
       width = 9, height = 8)

# Gammaproteobacteria 
lombgamma <- plot_lombstat(grid.sel = grid.gamma, 
              limits =  c(0,30))

ggsave(str_c(fig.path, 'lombstat_gammaproteobacteria_new', '.pdf'),
       plot = lombgamma,
       width = 9, height = 8)

trendgamma <- plot_trend(grid.sel = grid.gamma)

ggsave(str_c(fig.path, 'trend_gammaproteobacteria_new', '.pdf'),
       plot = trendgamma,
       width = 11, height = 8)


# Flavobacteria  
lombflavo <- plot_lombstat(grid.sel = grid.flavo, 
              limits =  c(0,25))

ggsave(str_c(fig.path, 'lombstat_flavo_new', '.pdf'),
       plot = lombflavo,
       width = 9, height = 8)

trendflavo<- plot_trend(grid.sel = grid.flavo)

ggsave(str_c(fig.path, 'trend_flavo_new', '.pdf'),
       plot = trendflavo,
       width = 9, height = 8)


