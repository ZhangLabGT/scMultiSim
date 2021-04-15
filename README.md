# SymSim
SymSim (Synthetic model of multiple variability factors for Simulation) is an R package for simulation of single cell RNA-Seq data. 

### Install from Github
This package can be installed with R package devtools. First, pull the package with git clone to your working directory. Make sure that you have installed the packages listed in the DESCRIPTION file.

The required Bioconductor packages can be installed as follows in R:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("IRanges")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biobase")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("S4Vectors")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")
```
To install SymSim, run:
```{r}
library("devtools")
devtools::install_github("YosefLab/SymSim")
```

### Generating datasets with SymSim

Refer to the [SymSim vignette](https://github.com/YosefLab/SymSim/blob/master/vignettes/SymSimTutorial.Rmd) for examples of using SymSim to simulate datasets.

[New feature] SymSim can now generate co-expressed gene modules. A new argument to the function `SimulateTrueCounts()` which is `gene_module_prop` represent the proportion of genes the user wants to be in gene modules. The smallest module size is 10 genes. For example, if the total number of genes in modules is 100, and the total number of populations is 5, then the 100 genes are randomy assigned to 5 gene modules. A gene module often correspond to a population, meaning that that set of genes tend to be highly expressed in a given population, but there is no guarantee on this. In the output of the function `SimulateTrueCounts()` there is an element `is_in_module` which gives the module ID for each gene (0 means the gene is not in any module).

### References

Xiuwei Zhang &ast;, Chenling Xu &ast;, Nir Yosef. **Simulating multiple faceted variability in Single Cell RNA sequencing**. _Nature Communications_, 10:2611, 2019. (https://www.nature.com/articles/s41467-019-10500-w).
