# Aggregation of most of the functions used in the time series analysis

######## DATA PREP #########

reorder_chrono <- function(df){
  require(lubridate)
  # Reorders the time columns, adding them into dat, 
  # and creating a decimal date tab
  
  mutate(df, year= str_c("20",year),
         month=str_to_title(month)) %>% 
    unite(dat, year,month,day,sep = "-",remove = F) %>% 
    mutate(dat =lubridate::ymd(dat), 
           decimaldat=decimal_date(dat)) %>% 
    arrange(dat)
  
}

ts_conversor <- function(df, s=c(2004, 1), f=12){
  # Not the pretiest function in the world..
  # I used it for creating the time series object in a situation 
  # when everything is in order (no missing data).
  
  ts(df$Abundance,start=s,frequency = f) 
  
}

######### EDA time series ############

## Some fast functions to visuazlie the data 

plot_timeseries <- function(asvs, psmelt = psmelt.relab,
                            title = NULL, color = "Phylum",
                            wrap = NULL, percent = TRUE, scale = "fixed"){
  # A simple plotting of a timeseries. A vector of the asvs for plotting is given
  # and the results are plotted through the psmelt.
  
  require(lubridate)
  
  basic <-   psmelt %>% 
    filter(OTU %in% asvs) %>% 
    mutate(decimaldat = decimal_date(Date) ) %>% 
    ggplot(aes(decimaldat, Abundance)) + 
    geom_line(aes_string(group = "OTU", color = color)) +
    ggtitle(title) +
    theme_timeseries  
  
  if (!is.null(wrap)) {
    basic = basic +
      facet_wrap(wrap, scale = scale)
  }
  
  if (percent) {
    basic = basic + 
      ylab('Relative abundance') +
      scale_y_continuous(labels = scales::percent)
  }
  
  if (color == "phylogroup") {
    basic = basic +
      phy.colScale
  }
  
  basic
}

# plot_ts <- function(asvs, physeq, wrap){
#   psmelt_dplyr(physeq) %>% 
#     
#   
#   }

check_periodogram <- function(ls,wrap="phylogroup"){
  
  # From a periodogram output, it plots the resulting periodograms in a nice way. 
  # It should be improved (quite a lot)
  
  ls %>% 
    map(~left_join(.x,bl.psmelt.raw,by=c("asv"="OTU"))) %>% 
    map(~ggplot(.) +
          geom_line(aes(whole.freq,whole.spec)) +
          facet_wrap(wrap) +
          theme(strip.background = element_blank(),
                strip.text.x = element_text(margin = margin(.05, 0,
                                                            .1, 0, "cm")),
                strip.text = element_text(size = 7),
                axis.text.y = element_text(size = 7)) +
          ggtitle(paste0("Periodogram of values for phylogroup level"))
    )
  
}

timedecay_plot <- function(physeq, distance="bray", lineareq = T, is.diss = T){
  
  ## We calculate the distance 
  dist.df <- as.data.frame(as.matrix(distance(physeq, method=distance)))
  
  ## Add some flavours to the mixxxxx
  dist.df$rwname <- rownames(dist.df)
  dist.df$rwseason <- get_variable(physeq, "season")
  dist.df$julian <- julian(get_variable(physeq, "Date"))
  
  # Put it in a nice way
  dist.df.tidied = dist.df %>%
    gather(clname,value,-rwname,-rwseason, -julian)
  
  match.julian <- data.frame(rwname =dist.df$rwname,
                             cljulian= dist.df$julian)
  
  dist.df.tidied <- merge(dist.df.tidied,match.julian,
                          by.x='clname',by.y='rwname')
  
  julian.df <- dist.df.tidied[,c(6,4,5)]
  
  # Calculate the day difference, and convert it to month
  # Any value above 0 is converted to 0 month!
  julian.df$diff_julian <- floor((julian.df$cljulian -
                                    julian.df$julian)/30)
  
  # And calculate from each the final values
  julian.df <- julian.df %>%
    filter(diff_julian > 0) %>% 
    arrange(diff_julian) %>% 
    mutate(x=diff_julian,
           y  = value)
  
  if(is.diss){
    julian.df$value = 1 - julian.df$value
    julian.df$y = 1 - julian.df$y
  }
  
  
  # Add the error bars
  diss.sum <- julian.df %>%
    group_by(diff_julian) %>% 
    dplyr::summarize(mean = mean(value),
                     ci = qt(0.95,
                             df=length(value)-1)*
                       sd(value)/sqrt(length(value)))
  
  
  # Linear regression equation?
  lm_eqn <- function(df){
    m <- lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x), 
                     list(a = format(coef(m)[[1]], digits = 2), 
                          b = format(coef(m)[[2]], digits = 2), 
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
  }
  
  label.y = case_when( 
    distance == 'bray' ~ 'Bray Curtis similarity',
    distance == 'euclidean' ~ 'Euclidean distance',
    distance == 'unifrac' ~ 'Unifrac similarity',
    distance == 'wunifrac' ~ 'Weigthed Unifrac similarity',
    TRUE ~ "Other distance"
  )
    
    
  # PLot everything nice
  dissplot <-   ggplot(julian.df, aes(x=diff_julian,y=value)) +
    geom_point(shape=1,size=1, color='grey',alpha=0.3) + 
    geom_linerange(data=diss.sum,
                   aes(ymin=mean-ci, ymax=mean+ci, x=diff_julian),
                   colour='gray64',
                   inherit.aes = FALSE) +
    geom_point(data=diss.sum,
               aes(x=diff_julian,y=mean),
               shape=1, size=3) +
    geom_line(data=diss.sum,
              aes(x=diff_julian,y=mean, group=1),
              color=ifelse(distance == 'bray',
                           tableau_color_pal('Tableau 20')(1),
                           'darkorange2')) +
    geom_smooth(method='lm',size=0.5) + 
    scale_x_continuous(breaks=c(1,12,24,36,48,60,72,84,96,108,120),
                       labels = c(0:10)) + 
    ylab(label.y) +
    xlab('Number of years between samples')
  
  if(lineareq){
    dissplot <- dissplot + 
      annotate(geom = "text",
               x=24 , y= mean(julian.df$y) / 1.8,
               label = lm_eqn(julian.df), size=3, parse = TRUE) 
  }
  
  
  #Save it!!
  return( dissplot + coord_cartesian( ylim = c(min(julian.df$y),
                                               max(julian.df$y))))
}

########## STATISTICAL CHECKING ########

check_seasonality <- function(ls, pval.t=0.01){
  
  # From a time series list, it check for each ts the seasonality
  # with a fisher.g.test of significance.
  # Additionally, the function outputs the maximum frequency from the periodogram
  # And corrects by FDR.
  
  ls %>% 
    map(~ts_conversor(.x)) %>% 
    map(~data.frame(per=GeneCycle::periodogram(.)$freq[
      which.max(GeneCycle::periodogram(.)$spec)],
      whole=GeneCycle::periodogram(.),
      pval=GeneCycle::fisher.g.test(.)))  %>% 
    bind_rows(.id = "asv") %>% 
    mutate(qval=fdrtool::fdrtool(pval, statistic="pvalue")$qval) %>% 
    filter(qval < pval.t)
}

check_seasonality_lomb <- function(ls, pval.t, reps){
  
  
  }

plot_monthts <- function(physeq, asvs){
 
  psmelt(physeq) %>% 
    filter(OTU %in% asvs) %>% 
    ggplot(aes(month, Abundance)) + 
    geom_jitter(alpha = 0.6) + 
    stat_smooth(aes(x = month,
                    group = OTU ), # continuous x-axis
              method = "lm",
              formula = y ~ poly(x, 3), se = F, size = 2, show.legend = F) + 
    facet_wrap(~OTU)
  
  }
