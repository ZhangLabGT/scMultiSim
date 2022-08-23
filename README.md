# scMultiSim
scMultiSim is an in-silico simulator that generates multi-modality data of single cells, including gene-expression data and chromatin accessibility data, guided by gene regulatory networks. In particular, we generate both unspliced counts and spliced counts as measures of gene-expression, as well as the true RNA velocity.  Technical noise can be added to generate observed RNA-seq counts (which can be read counts or UMI counts) from true transcript counts.   Technical noise considers capture efficiency, amplification bias, sequencing depth and batch effects.  

Refer to the [scMultiSim vignette](https://github.com/squiresmf/scMultiSim/blob/master/vignettes/scMultiSimTutorial.Rmd) for examples of using scMultiSim to simulate datasets and for a description of tool functionalities.

### Install from Github
This package can be installed with R package devtools.  The following installs scMultiSim and other necessary packages:
```{r, message=F, warning=F, eval=T}
library("devtools")
devtools::install_github("squiresmf/scMultiSim")
library("scMultiSim")
list.of.packages <- c("reshape", "ape", "phytools", "repr", "KernelKnn")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
```

### References

Xiuwei Zhang &ast;, Chenling Xu &ast;, Nir Yosef. **Simulating multiple faceted variability in Single Cell RNA sequencing**. _Nature Communications_, 10:2611, 2019. (https://www.nature.com/articles/s41467-019-10500-w).
