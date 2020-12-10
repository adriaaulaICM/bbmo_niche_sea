# Functions related to the gradients project

## Creator: Adria Auladell


# phyloseq ----------------------------------------------------------------

tax_df <- function(physeq){
  
  as( tax_table(physeq), 'matrix') %>% 
    data.frame() %>% 
    return() 
  
}
  

######## dbRDA ##########

#From a selection of variables made through PERMANOVA (function adonis) or 
#other means, with distance and ordination an dbRDA is performed finding the explanatory
#variables of the distribution.


select_variables_adonis <- function(physeq, pval = 0.01){
  
  # Selects from an adonis test the significant variables and returns them 
  otu.adonis <- as(otu_table(physeq), "matrix") %>%
    t() %>%
    data.frame() 
  
  envdata <- data.frame(sample_data(physeq)) %>% 
    # This is a MONSTER AND it has to be changed
    select(-Fraction)
  
  ad <- adonis(otu.adonis ~ . , data = envdata, permutations = 999)
  print(ad)
  
  bestEnvVariables <- data.frame(ad$aov.tab) %>% 
    rownames_to_column(var = 'varinterest') %>% 
    filter(Pr..F. <= pval) %>% 
    .$varinterest
  
  return(bestEnvVariables)
  
}

distance_and_ordination <- function(physeq, bestEnvVariables){
  #Calculate the distances
  cohort.distance <- phyloseq::distance(physeq,
                                        method = "bray",
                                        type = "samples")
  
  #And ordinate by the function, extracting the biplot
  cohort.distance.ord <- do.call("ordinate",
                                 list(physeq,
                                      "CAP",
                                      cohort.distance,
                                      as.formula(paste0("~",
                                                        str_c(bestEnvVariables,
                                                              collapse = "+"))),
                                      na.action = na.exclude))
  
  
  names(cohort.distance.ord)
  
  print(anova(cohort.distance.ord, by = "terms"))
  
  print(summary(eigenvals(cohort.distance.ord)))
  
  vectors <- data.frame(cohort.distance.ord$CCA$biplot)
  # The find the length (magnitude) of the vectors in the first two dimensions
  vectors$length1and2 <- sqrt(vectors$CAP1^2 + vectors$CAP2^2)
  
  #Better names for the variables
  vector.data <- data.frame(cohort.distance.ord$CCA$biplot) %>% 
    rownames_to_column(var = 'rowname')
  
  return( list( ordination = cohort.distance.ord, 
                vectordata = vector.data))
}

plot_dbrda <- function(physeq, ordination, vector.data, color, shape){
  
  circleFun <- function(center = c(0,0), radius = 1, npoints = 100){
    r = radius
    tt <- seq(0, 2*pi, length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }
  
  plot <- plot_ordination(physeq,
                          ordination,
                          "samples", 
                          color = color, 
                          shape = shape) + 
    geom_point( size = 2.5) +
    geom_path(aes(x = x, y = y),
              data = circleFun(c(0,0), 1, 100),
              color = "lightgrey",
              alpha = 0.7,
              inherit.aes = FALSE) + 
    geom_segment(data = vector.data,
                 aes(x = 0, xend = CAP1,
                     y = 0, yend = CAP2),
                 arrow = arrow(length = unit(0.25, "cm")),
                 inherit.aes = FALSE,
                 color = 'grey') + 
    geom_text(data = vector.data, 
              aes(x = CAP1*1.2, y = CAP2*1.2,
                  label = str_replace(rowname,"_", " ")), 
              size = 3.5,
              inherit.aes = FALSE) +
    coord_fixed() + 
    theme_bw() +
    theme(plot.margin = unit(c(0,0,0,0),"mm"),
          legend.text = element_text(size = 12)) 
  
  plot
  
}

# This does all of it. Wonderful eh?
dbrda_general_pipeline <- function(physeq, envdata, color, shape){
  
  best.environmental.var <- select_variables_adonis(physeq)
  ordination_vector <- distance_and_ordination(physeq, best.environmental.var)
  
  plot_dbrda(physeq, ordination = ordination_vector$ordination, 
             vector.data =  ordination_vector$vectordata, 
             color = color,
             shape = shape)
  
  }




###### OTHER #######

write_fasta <- function(physeq, outfile, totals = FALSE){
  # This cool function writes easily any fasta
  # from a phyloseq object.  
  # one column should be named fasta, with the sequences 
  
  require('glue')
  require('readr')
  
  tax <- data.frame(as(tax_table(physeq), 'matrix'),
                    taxsum = taxa_sums(physeq)) %>% 
    select(-asv) %>% 
    rownames_to_column(var = 'asv') %>% 
    arrange(-taxsum)
  
  # The rownames will always be the asv name! 
  # If the id is the seq you will have quite a redundant fasta...
  
  if (totals) {
    
    fasta <- (tax) %>% glue_data(" >{asv};size={taxsum} 
                        {fasta}")
    write_lines(fasta, outfile)   
    
  }else{
    fasta <- (tax) %>% glue_data(" >{asv} 
                        {fasta}")
    write_lines(fasta, outfile)   
  }

}

summary_column <- function(df){
  df %>% 
  summarize(min = min(Abundance), 
            Q1 = quantile(Abundance)[[2]],
            mean = mean(Abundance),
            median = median(Abundance),  
            Q3 = quantile(Abundance)[[4]],
            max = max(Abundance)) 
}

psmelt_dplyr = function(physeq) {
  #Implementation of the psmelt from phyloseq with dplyr
  # It is indeed faster, without further complications or differences. 
  # (well, in fact, its way more prone to give errors if the output is not well established, 
  # it doesn't check anything)

  sd = data.frame(sample_data(physeq)) %>% 
    rownames_to_column("Sample")
  TT = data.frame(as(tax_table(physeq),'matrix')) %>%
    rownames_to_column("OTU")
  if(taxa_are_rows(physeq)){
    otu.table = data.frame(as(otu_table(physeq),"matrix"),
                           check.names = FALSE) %>%
      rownames_to_column("OTU")
  } else {
    otu.table = data.frame(t(as(otu_table(physeq),"matrix")),
                           check.names = FALSE) %>%
      rownames_to_column("OTU")
  }
  
  all <- otu.table %>%
    left_join(TT, by = "OTU") %>%
    gather_("Sample", "Abundance", setdiff(colnames(otu.table), "OTU")) %>%
    left_join(sd, by = "Sample") %>% 
    select(Sample, Abundance, OTU, everything())
  
  return(all)
}

## Rarefaction curve, ggplot style
ggrare <- function(physeq, step = 10,
                   label = NULL, color = NULL,
                   plot = TRUE, parallel = FALSE, se = TRUE) {
  ## Args:
  ## - physeq: phyloseq class object, from which abundance data are extracted
  ## - step: Step size for sample size in rarefaction curves
  ## - label: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - color: (Optional). Default ‘NULL’. Character string. The name of the
  ##          variable to map to colors in the plot. This can be a sample
  ##          variable (among the set returned by
  ##          ‘sample_variables(physeq)’ ) or taxonomic rank (among the set
  ##          returned by ‘rank_names(physeq)’).
  ##
  ##          Finally, The color scheme is chosen automatically by
  ##          ‘link{ggplot}’, but it can be modified afterward with an
  ##          additional layer using ‘scale_color_manual’.
  ## - color: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - plot:  Logical, should the graphic be plotted.
  ## - parallel: should rarefaction be parallelized (using parallel framework)
  ## - se:    Default TRUE. Logical. Should standard errors be computed.
  ## require vegan
  x <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) { x <- t(x) }

  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)

  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)

  ## Get sample data
  if (!is.null(sample_data(physeq, FALSE))) {
    sdf <- as(sample_data(physeq), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }

  ## Add, any custom-supplied plot-mapped variables
  if( length(color) > 1 ){
    data$color <- color
    names(data)[names(data)=="color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  if( length(label) > 1 ){
    labels$label <- label
    names(labels)[names(labels)=="label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }

  p <- ggplot(data = data, aes_string(x = "Size", y = ".S", group = "Sample", color = color))
  p <- p + labs(x = "Sample Size", y = "Species Richness")
  if (!is.null(label)) {
    p <- p + geom_text(data = labels, aes_string(x = "x", y = "y", label = label, color = color),
                       size = 4, hjust = 0)
  }
  p <- p + geom_line()
  if (se) { ## add standard error if available
    p <- p + geom_ribbon(aes_string(ymin = ".S - .se", ymax = ".S + .se", color = NULL, fill = color), alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}

data_transformation <- function(physeq){
  require(tidyverse)
  require(phyloseq)
  
  
  
  # Relative abundance of OTUs
  phy.relab = transform_sample_counts(physeq, function(x) x / sum(x))
  
  # Log10 abundance of OTUs
  phy.log = transform_sample_counts(physeq, function(x) log10(x+1))
  
  # Rarefied dataset
  phy.raref = rarefy_even_depth(physeq,rngseed = 42)
  
  # Log centered ratio (Aitchison)
  geoMean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  
  phy.clr <- transform_sample_counts(physeq, function(x) x+1) %>%
    transform_sample_counts(function(x) log(x/geoMean(x)))
  
  print('if this is not close to zero or zero, something is wrong! Beware')
  print(sum(otu_table(phy.clr)[1,]))
  
  print("Starting a psmelt of a large dataset, this could take a while!")
  
  return(list( relab = phy.relab, 
               log = phy.log,
               raref = phy.raref,
               clr = phy.clr))
}

save_mult_plot <- function(plots, vec_files, common, endfile, width, height){
  
  # Saves a multiple plot generated with map_ procedures 
  # Both plots and vec_files should have the same length
  
  filenames <- str_c(common,
                     vec_files ,
                     endfile)
  
  map2(filenames, plots, ggsave, width = width , height = height )
  
}

correct_samnames <- function(df){
  
  df$sample <-  str_replace_all(df$sample,
                                   c( "BL100413" = "BL100414",
                                      "BL110704" = "BL110705",
                                      "BL120518" = "BL120511",
                                      "BL131204" = "BL131215"))
  return(df)
}

