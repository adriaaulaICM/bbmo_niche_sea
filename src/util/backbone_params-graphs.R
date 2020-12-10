
#### General #####

library(ggthemes)
library(tidyverse)

#Black and white theme and the text size bigger
### If the work is a presentation, its best to give a base size of...around 13
theme_set(theme_bw(base_size = 11))

# Lets also set the seed for random processes
set.seed(42)


#### Palettes phylogeny (BACTERIA) #####

superpalette <- c("yellow1", #B_Flavobacteriia
                  "#ffce00", #B_Sphingobacteriia
                  "#ff9a00",  #Other_Bacteroidetes
                  "limegreen", #Cya_Synechococcus
                  "springgreen4", #Other_Cyanobacteria
                  "lightcyan2", #P_Alp_Rhodobacterales
                  "#78aaff", #P_Alp_Rhodospirillales
                  "steelblue1", #P_Alp_Rickettsiales
                  "#4188ff", #"P_Alp_SAR11_clade
                  "steelblue4", #P_Alp_Sphingomonadales
                  "royalblue4", #Alp_Proteobacteria
                  "#efbbff", #P_Gam_Alteromonadales
                  "#d896ff", #P_Gam_Cellvibrionales
                  "#b76bce", #P_Gam_Oceanospirillales
                  "#be29ec", #P_Gam_Pseudomonadales
                  "#660066",  #Gam_Proteobacteria
                  "#330033", #Other_Proteobacteria
                  "#f7cfcf", #O_Actinobacteria
                  "indianred1", #O_Firmicutes
                  "indianred", #O_Planctomycetes
                  "#a03e3e", #"O_Verrucomicrobia
                  "#492a2a") #"Other_Bacteria


tax.order <-c("B_Flavobacteriia",
              "B_Sphingobacteriia",
              "Other_Bacteroidetes",
              "Cya_Synechococcus",
              "Other_Cyanobacteria",
              "P_Alp_Rhodobacterales",
              "P_Alp_Rhodospirillales",
              "P_Alp_Rickettsiales",
              "P_Alp_SAR11_clade",
              "P_Alp_Sphingomonadales",
              "Alp_Proteobacteria",
              "P_Gam_Alteromonadales",
              "P_Gam_Cellvibrionales",
              "P_Gam_Oceanospirillales",
              "P_Gam_Pseudomonadales",
              "Gam_Proteobacteria",
              "Other_Proteobacteria",
              "O_Actinobacteria",
              "O_Firmicutes",
              "O_Planctomycetes",
              "O_Verrucomicrobia",
              "Other_Bacteria")

names(superpalette) <- tax.order


bac.colScale <- scale_color_manual(name = "16S taxonomy",
                                   values = superpalette,
                                   limits = names(superpalette))
bac.fillScale <- scale_fill_manual(name = "16S taxonomy",
                                   values = superpalette,
                                   limits = names(superpalette))


### Palettes EUKARYA #####

superpalette.euk <- c("olivedrab3", #Alveolata/Dinoflagellata
                      "springgreen4", #Alveolata/MALV-I
                      "yellow",  #Alveolata/MALV-II
                      "#ffce00", #Other MALVs
                      "#ff9a00", #Other Alveolates
                      "#ff4d00", #Amoebozoa
                      "#78aaff" , #Archaeplastida/Chlorodendrophyceae
                      "steelblue1" , #Archaeplastida/Mamiellophyceae
                      "steelblue4" , #Archaeplastida/Prasinophyceae
                      "royalblue4" , #Other Archaeplastida
                      "#b76bce", #Hacrobia/Cryptomonadales
                      "#efbbff", #Hacrobia/Picozoa
                      "#be29ec", #Other Hacrobia
                      "#660066", #Opisthokonta
                      "indianred", #Rhizaria
                      "#a03e3e",  #Stramenopiles/Diatomea
                      "#492a2a", #Stramenopiles/Labyrinthulomycetes
                      "indianred1", #Other Stramenopiles
                      "#666463"  #Other Eukaryotes 
)

tax.order.euk <- c(   
  "Alveolata/Dinoflagellata",
  "Alveolata/MALV-I",
  "Alveolata/MALV-II",
  "Other MALVs",
  "Other Alveolates",
  "Amoebozoa",
  "Archaeplastida/Chlorodendrophyceae",
  "Archaeplastida/Mamiellophyceae",
  "Archaeplastida/Prasinophyceae",
  "Other Archaeplastida",
  "Hacrobia/Cryptomonadales",
  "Hacrobia/Picozoa",
  "Other Hacrobia",
  "Opisthokonta",
  "Rhizaria",
  "Stramenopiles/Diatomea",
  "Stramenopiles/Labyrinthulomycetes",
  "Other Stramenopiles",
  "Other Eukaryotes")
  
names(superpalette.euk) <- tax.order.euk

euk.colScale <- scale_color_manual(name = "18S taxonomy",
                                   values = superpalette.euk,
                                   limits = names(superpalette.euk))
euk.fillScale <- scale_fill_manual(name = "18S taxonomy",
                                   values = superpalette.euk,
                                   limits = names(superpalette.euk))

###### VARIOUS AESTHETICS #####

##### GGPlot ####

# Label related 
leg.bottom <- theme(legend.direction="horizontal",legend.position="bottom")
lab.flip <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
labelsx.diag <- theme(axis.text.x = element_text(angle = -45,vjust = 0.5))

# Strip problems
lil.strip <- theme(strip.background = element_blank(),
                   strip.text.x =element_text(margin = margin(.05, 0, .1, 0, "cm")))


#### $$ TIME $$ #####

## ~~~ Season~~~

season.order <- c("winter","spring","summer","autumn")

sea.shape <- c(15,16,17,25)
sea.col <- c("#63bfb8","#52b216","#ecee21","#e3a40c")

names(sea.shape) <- season.order
names(sea.col) <- season.order

sea.shapeScale <- scale_shape_manual(name = "Season",
                                     values = sea.shape,
                                     limits = season.order)
sea.colScale <- scale_color_manual(name = "Season",
                                   values = sea.col,
                                   limits = season.order)
sea.fillScale <- scale_fill_manual(name = "Season",
                                   values = sea.col,
                                   limits = season.order)

## ~~~Month~~~~

month.num <- c("01", "02",  "03", "04",
              "05", "06", "07", "08",
              "09", "10", "11", "12")
               
month.order <- c("jan","feb","mar","apr","may","jun",
                 "jul","aug","sep","oct","nov","dec")

mon.colScale <- scale_colour_manual(name = "Month",
                                    values=c(hc_pal()(10)[-5],
                                             "#d11919","#B37400","#4C4F8B"),
                                    limits=month.num,
                                    labels = str_to_title(month.order))

mon.fillScale <- scale_fill_manual(name = "Month",
                                    values=c(hc_pal()(10)[-5],
                                             "#d11919","#B37400","#4C4F8B"),
                                    limits=month.num,
                                    labels = str_to_title(month.order))

scale_x_month <- scale_x_discrete(name = 'Month',
                                  labels = str_to_title(month.order))
## ~~~Decade plots~~~ 

darkzone <- data.frame(x = c(2003:2013) + 0.75,
                       x2 = c(2004:2014) + 0.25)

thedark <- geom_rect(data=darkzone,
                     fill="grey",alpha=0.4,
                     aes(xmin=x,
                         xmax=x2),
                     ymin=-Inf,ymax=Inf, inherit.aes = FALSE)

theme_timeseries <- list(
  thedark ,
  scale_x_continuous(breaks = c(2005, 2007, 2009, 2011, 2013)) ,
  xlab("Year") ,
  ylab("Relative abundance"),
  lil.strip)

#### Day of the year 

num.days.mnt <- c(0,31,28,31,30,31,30,31,31,30,31,30)
cumnum <- cumsum(num.days.mnt)

scale_dayyear <- scale_x_continuous(breaks = cumnum,
                                 labels = str_to_title(month.order)) 
  
scale_dayyear_shrt <- scale_x_continuous(breaks = cumnum,
                                 labels = str_to_title(str_sub(month.order,1,1))) 

#### ~~~~~Location~~~~~~~ ####

# Implemented for gradients analysis, around 20180606
location.order <-  c('Girona','Barcelona',
                     'Tarragona','Palma')

colors.location <- c("green4", "#F8766D", "#C71CFF", "dodgerblue")

names(colors.location) <- location.order

loc.colScale <- scale_color_manual(name = "Location",
                                   values = colors.location,
                                   limits = location.order)

#### Fraction ####

# Implemented for gradients analysis, around 20180606
frac.linetype <- list(scale_linetype_manual(name = 'Fraction',
                                       values = c(1,5),
                                       breaks = c('pico','nano'),
                                       limits = c('pico', 'nano')),
  guides(linetype = F))


