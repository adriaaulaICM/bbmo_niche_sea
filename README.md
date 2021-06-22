# Scripts for *Seasonal niche differentiation among closely related marine bacteria*

Author: Adrià Auladell

Inside `src` there are all the scripts with the following structure: 

```
src/
    data/ # processing of raw files with DADA2 + phyloseq
    analysis/ # all the analysis generating new data or statistics
    figures/  # the visualizations, mainly through ggplot
    utils/ # scripts called from the abovementioned scripts
```

Inside the `data/cleaned` folder there is a phyloseq R object with the sequence data, abundance table and sample metadata. 
To execute most of the analysis, create an R project in the main folder and execute the desired scripts. Most of the `figures` scripts 
require to have executed scripts in the `analysis` folder for the statistics to be saved in the correspondent folders. 

You will need some minimal packages for running the scripts. For most analysis: 

- [tidyverse](https://www.tidyverse.org/): data wrangling, visualization and general purpose package.
- [phyloseq](https://joey711.github.io/phyloseq/): specific microbiome experiment functions. 
    - [speedyseq](https://github.com/mikemc/speedyseq): faster processes.
- [propr](https://github.com/tpq/propr): proportionality calculation. 
- [lomb](https://cran.r-project.org/web/packages/lomb/index.html): seasonality testing.
- [DECIPHER](http://www2.decipher.codes/index.html): nucelotide calculations.

Regretabbly, there are also tons of other packages that I used for small processes. Sorry for that. Package-addict here :^/

Some of there are: `janitor, geofacet, patchwork, gt, ggthemes, ggtext, ggraph...`. 

If there is an issue dealing with the code or something else either post an issue in the project or write an email to aauladell@icm.csic.es :)

