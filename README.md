# Scripts for **Seasonal niche differentiation between evolutionary closely related marine bacteria**

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



If there is an issue dealing with the code or something else either post an issue in the project or write an email to aauladell@icm.csic.es :)

