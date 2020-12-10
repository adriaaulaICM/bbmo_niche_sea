library(tidyverse)
library(phyloseq)
library(pheatmap)
library(paletteer)

source("src/util/backbone_functions.R") 
source("src/util/backbone_params-graphs.R") 
source('src/util/params-graphs.R')
source("src/util/main_functions.R")

figpath <- 'results/figures/cluster'
dir.create(file.path(figpath), showWarnings = FALSE)


bl.phy <- readRDS('data/cleaned/blphyloseq.rds') 
lomb.02 <- readRDS('data/analysis/lomball.rds')
lomb.sea <- readRDS('data/analysis/lombsea.rds') 

geoMean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

bl.phy.clr <- transform_sample_counts(bl.phy, function(x) x+0.01) %>%
  transform_sample_counts(function(x) log(x/geoMean(x)))

bl.phy <- prune_taxa(lomb.sea$asv, bl.phy)
bl.phy.clr <- prune_taxa(lomb.sea$asv, bl.phy.clr)

bl.phy.relab <- transform_sample_counts(bl.phy, function(x) x / sum(x))
bl.phy.asinh <- transform_sample_counts(bl.phy, function(x) asinh(x))


# Clusters of seasonal ASVs -----------------------------------------------
library(cluster)
#Takes quite a while I have hid this part of the code since it is unnecesary now
#pamfun = function(x, k) list(cluster = pam(x, k, cluster.only = TRUE))
#
#gss = clusGap(t(otu_table(bl.phy.asinh)), 
#              FUN = pamfun, K.max = 8, B = 250,
#              verbose = FALSE)
#
#plot_gap = function(x) {
#  gstab = data.frame(x$Tab, k = seq_len(nrow(x$Tab)))
#  ggplot(gstab, aes(k, gap)) + geom_line() +
#    geom_errorbar(aes(ymax = gap + SE.sim,
#                      ymin = gap - SE.sim), width=0.1) +
#    geom_point(size = 3, col =  "red")
#}
#
#plot_gap(gss)
#
#
#ggsave( 'results/figures/cluster/heatmap_clusgap.pdf',
#        height = 10, width = 10)
#
#k2 = maxSE(gss$Tab[, "gap"], gss$Tab[, "SE.sim"],
#           method = "Tibs2001SEmax")
#
## The number of clusters taken into account are the following: 
#k2
k2 = 3


# Visualization -----------------------------------------------------------

env <- sample_data(bl.phy) %>% data.frame() %>% 
  select(season)  %>% 
  rownames_to_column() %>% 
  mutate(season = factor(season, levels = season.order)) %>% 
  column_to_rownames( )

# we will put the family 
tax <- as(tax_table(bl.phy), 'matrix') %>% as_tibble(rownames = 'asv')  %>% 
  mutate(family = fct_lump(family, n = 8)) %>% 
  mutate(genus = fct_lump(genus, n = 7)) %>% 
  select(asv, family)  %>% 
  # mutate( family = ifelse(family %in% fams.italy,
  #                         bquote(italic(.(family))),
  #                         family)) %>% 
  column_to_rownames(var = 'asv') 


family.pal <- ptol_pal()(9)
names(family.pal) <- tax$family %>% levels()

# and the genus 

# gen.pal <- paletteer_c("scico::tokyo", n = 8)
# names(gen.pal) <- tax$genus %>% levels()

temp.colors <- list(season = sea.col, 
                    family = family.pal)
                    # genus = gen.pal)


matrix.asv <- t(otu_table(bl.phy.clr))

rownames(matrix.asv)

pheatmap(matrix.asv, 
         annotation_col = env,
         annotation_row = tax,
         labels_row = rownames(matrix.asv) %>% str_to_upper(),
         cutree_rows = k2,
         annotation_colors = temp.colors, 
         cellheight = 9, cellwidth = 9,
         cluster_rows = T, cluster_cols = F,
         filename = str_c(figpath, '/heatmap_clr.pdf'))



pheatmap(t(otu_table(subset_taxa(bl.phy.clr, family == "Flavobacteriaceae"))),
         annotation_col = env,
         annotation_row = tax,
         # cutree_rows = k2,
         annotation_colors = temp.colors, 
         cellheight = 7, cellwidth = 7,
         cluster_rows = T, cluster_cols = F,
         filename = str_c(figpath, '/heatmap_flavos_clr.pdf'))

 
# A shorter version for the poster 

library(mgcv)

sel.top.ab <- lomb.sea %>% 
  mutate(index =  str_replace(asv, 'asv', '') %>% as.integer()) %>% 
  arrange(index) %>% 
  head(n = 25) %>% 
  pull(asv)

bl.phy.clr <- prune_taxa(sel.top.ab, bl.phy.clr)

megafile.agg <- psmelt(bl.phy.clr)

gam.gene.model <- megafile.agg %>% 
  group_by(OTU) %>% 
  nest()  %>% 
  mutate( model = map(data,
                      ~gam(Abundance ~ s(day_of_year, k =12, bs = 'cc'), data = .)),
          predict = map(model, ~predict(.))) %>% 
  unnest(c(data,predict)) 

pred.matrix <- gam.gene.model %>% 
  select(OTU, Sample, predict) %>% 
  pivot_wider(names_from = 'OTU', values_from = 'predict') %>% 
  arrange(Sample) %>% 
  column_to_rownames('Sample') %>% 
  as.matrix() %>% 
  .[1:12,]

pheatmap((pred.matrix), 
         annotation_row = env,
         annotation_col = tax,
         labels_col = colnames((pred.matrix)) %>% str_to_upper(),
         show_rownames = FALSE,
         annotation_colors = temp.colors, 
         cellheight = 10, cellwidth = 10,
         cluster_rows = F, cluster_cols = T,
         filename = str_c(figpath, '/heatmap_clr_poster.pdf'))

