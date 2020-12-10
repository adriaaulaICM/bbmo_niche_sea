library(tidyverse)
library(ggthemes)
library(speedyseq)
library(ggforce)
library(patchwork)
library(scales)

bl.phy <- readRDS('data/cleaned/blphyloseq.rds')
abunprev.crt <- read_rds('data/analysis/abund_prev.rds')
abuntotal <- abunprev.crt %>% 
  group_by(behavior) %>% 
  summarize(total = n()) %>% 
  mutate(prevalence = c(0.25, 0.1, 0.75,1)) %>% 
  mutate( label = str_c(behavior, total, sep = ':\n' ))
  
tax <- as(tax_table(bl.phy), 'matrix') %>% 
  as_tibble(rownames = 'asv') 

abund.class.df <- abunprev.crt %>% 
  left_join(tax, by = 'asv') %>% 
  filter(!is.na(class)) %>% 
  select(prevalence, mean.reads, behavior, class, crt)  %>% 
  mutate(behavior.num = case_when(behavior == 'Broad' ~ 1,
                              behavior == 'Intermediate' ~ 0.75,
                              behavior == 'Narrow' ~ 0.10, 
                              behavior == 'CRT' ~ 0.25,
                              TRUE ~ 0),
         class = ifelse(class %in%  c("Alphaproteobacteria", "Cyanobacteriia",
                                      "Gammaproteobacteria", "Bacteroidia",
                                      "Verrucomicrobiae", "Acidimicrobiia"),
                        class,
                        'Other class') %>% 
           as.factor())

# Plotting ----------------------------------------------------------------
# prev.pal <- c('gray16', 'firebrick', 'coral2', 'darkslateblue')
# names(prev.pal) <- abunprev.crt$behavior %>% unique()

abundplot <- ggplot(abund.class.df, aes(prevalence, mean.reads)) +
  geom_mark_hull(data = abuntotal %>% filter(behavior == 'CRT'),
                 aes(label = behavior, description = 81), 
                 linetype = 2, 
                 expand = unit(3, "mm"),
                 show.legend = T) +
  geom_vline(data = abuntotal %>% filter(behavior != 'CRT'),
             aes(xintercept = prevalence), linetype = 2) + 
  geom_label(data = abuntotal %>% filter(behavior != 'CRT'),
             aes(x = prevalence,
                 y = 0.0000002,
                 label = label),
             inherit.aes = F) + 
  geom_point(aes(color = class, shape = crt), alpha = 0.8) + 
  coord_flip() +
  scale_y_log10(labels = scales::percent) + 
  scale_x_continuous(labels = scales::percent,
                     limits = c(0,1.1),
                     breaks = c(0,0.1,0.3,0.6,0.75,0.9,1)) + 
  guides(color = FALSE, shape = FALSE) + 
  theme_minimal() + 
  ylab(' Mean relative abundance (%, log10 scaled)') + 
  xlab('Occurrence') + 
  theme(text = element_text(size = 13),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.background = element_blank()) 

# Distribution of classes  ------------------------------------------------

solpal <- solarized_pal()(7)[c(2,1,4,6,5,3,7)]
names(solpal) <-  abund.class.df$class %>% levels
solpal[6] <- '#FED766'

class.distro <- abund.class.df %>% 
  ggplot(aes(y = behavior.num, fill = class)) +
  geom_hline(data = abuntotal %>% filter(behavior != 'CRT'),
             aes(yintercept = prevalence),
             linetype = 2, show.legend = F) + 
  geom_bar(position = 'fill') + 
  scale_x_continuous(labels = scales::percent) + 
  scale_y_continuous(labels = scales::percent, limits = c(0,1.1)) + 
  theme_minimal() + 
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.background = element_blank()) + 
  xlab(' Number of ASVs (%)') 

psmelt.phy <- psmelt(bl.phy)

gen.sel <- c('HIMB59', 'AG-337-I02', 'Pelagibacter',
             'Pelagibacter_A','Luminiphilus', 'Synechococcus_C',
             'Prochlorococcus_A', 'Glaciecola', 'HIMB11', 'Amylibacter',
             'Nereida', 'SAR86A', 'HIMB59_g', 'UBA7446')

fam.sel <- c('D2472', 'HIMB59', 'HIMB59_f', 'Flavobacteriaceae')

order.sel <- c('Flavobacteriales', 'Opitutales',
               'Verrucomicrobiales', 'Rhodobacterales')
class.sel <- c('Acidimicrobiia')

phy.class.ord <- tax  %>% 
  mutate_if(is.factor, as.character) %>% 
  filter(genus %in% gen.sel | family %in% fam.sel |
           order %in% order.sel | class %in% class.sel)

#-------------
sil <- readRDS('data/cleaned/dada2/03_taxonomy/bbmo_silva/bbmo_silva_tax_assignation.rds')%>% as_tibble(rownames = 'seq')
thechosenones <- psmelt.phy %>% filter(is.na(phylum)) %>% pull(seq) %>% unique()
silunkows <- sil %>% 
  filter(seq %in% thechosenones)

seqtotals <- psmelt.phy %>%
  select(seq,Abundance) %>%
  group_by(seq) %>%
  summarise(total = sum(Abundance)) %>% 
  mutate( seq = as.character(seq))

silunkows %>% 
  left_join(seqtotals) %>% 
  arrange(-total) %>% 
  View()

#---------------
selphy <- phy.class.ord$phylum %>% unique()
selclass <- phy.class.ord$class %>% unique() %>% na.omit()
selord <- phy.class.ord$order %>% unique() %>% na.omit()

fam <- tax  %>% 
  mutate_if(is.factor, as.character) %>% 
  filter(genus %in% gen.sel | family %in% fam.sel) %>% 
  pull(family) %>% 
  unique()

strdata <- psmelt.phy %>% 
  mutate_if(is.factor, as.character) %>% 
  group_by(phylum, class, order, family, genus) %>% 
  # filter(!is.na(class)) %>%
  # filter(!is.na(phylum)) %>%
  summarise(value = sum(Abundance)) %>% 
  ungroup() %>% 
  mutate( value = value / sum(value)) %>% 
  arrange(phylum) %>% 
  mutate(phylum = fct_inorder(phylum) %>%
           fct_other(keep = selphy, other_level = 'Other phylum'),
         class = fct_inorder(class) %>%
           fct_other(keep = selclass, other_level = 'Other class'),
         order = fct_inorder(order) %>%
           fct_other(keep = selord , other_level = 'Other order'),
         family = ifelse(family == 'HIMB59', 'HIMB59_f', family) %>% 
           fct_inorder() %>%
           fct_other(keep = c(fam, fam.sel), other_level = 'Other family'),
         genus = ifelse(genus == 'HIMB59', 'HIMB59_g', genus) %>% 
           fct_inorder() %>%
           fct_other(keep = gen.sel, other_level = 'Other genus')) %>% 
  select(-phylum)


paralelldf <- gather_set_data(strdata, 1:4) %>% 
  mutate(x = factor(x, levels = c('class', 'order',
                                  'family', 'genus')))  %>% 
  group_by(y) %>% 
  mutate(total = (sum(value) * 100) %>% round(digits = 1)) %>% 
  mutate(ynew = str_c(y, '; ', total, '%')) %>% 
  filter(!is.na(genus))  %>% 
  ungroup()   


prettylabs <- paralelldf %>% 
  select(y, ynew) %>% 
  distinct() %>% 
  arrange(y)

paralo <- paralelldf %>%
  mutate( y = factor(y, levels = levels(y), labels = prettylabs$ynew))

  
parallel.plot <- ggplot(paralo , aes(x, id = id, split = y, value = value)) +
  geom_parallel_sets(aes(fill = class),
                     show.legend = F, alpha = 0.6) +
  geom_parallel_sets_labels(angle = 360, size = 3)  + 
  scale_x_discrete(labels = str_to_title) + 
  theme_minimal() + 
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        axis.line.y = element_blank()) + 
  ylab('Total relative abundance (%)') + 
  xlab('Taxonomic rank')

mainplot <-  (parallel.plot / (abundplot + class.distro)) & 
  plot_layout(guides = 'collect') &
  plot_annotation(tag_levels = 'A') & 
  scale_fill_manual(values = solpal, name = 'Class') &
  scale_color_manual(values = solpal) & 
  theme(legend.position = 'bottom')


ggsave(mainplot, filename = 'results/figures/abun_preval_mainplot.pdf',
       width = 9, height = 10)
