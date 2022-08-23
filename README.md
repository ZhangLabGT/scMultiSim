# scMultiSim

scMultiSim is an in silico simulator that generates multi-modality data of single-cells, including gene expression, chromatin accessibility, RNA velocity, and spatial location of cells. It takes a cell differential tree and a gene regulatory network (GRN) as input, and simulates spliced and unspliced counts while accounting for the relationships between modalities. The output single cell gene expression data is determined by three factors: cell-cell interactions, within-cell GRNs and chromatin accessibility. Users can tune the effect of each factor on the output data and set various parameters for the underlying model. Furthermore, the GRN can be set in a time-varying mode where the network's structure changes temporally to reflect the dynamic nature of biological networks. We also provide options to simulate technical variations such as batch effects. scMultiSim can be used to benchmark challenging computational tasks on single-cell multi-omics data, including the inference of GRNs, estimation of RNA velocity, integration of single-cell datasets from multiple batches and modalities, and analysis of cell-cell interaction using the cell spatial location data.

Refer to the [vignette](https://zhanglabgt.github.io/scMultiSim/vignettes/sim_new.nb.html) for examples of using scMultiSim to simulate datasets and for a description of tool functionalities.

![Overview](https://github.com/ZhangLabGT/scMultiSim/raw/img/img/scMultisim.png)

## Results

The following figure briefly shows results from the same cell differential tree:

1. Connected scATAC-seq and scRNA-seq, in continuous or discrete mode. Visualized by t-SNE.
2. GRN correlation heatmap, where genes regulated by the same regulator have similar correlations with others.
3. Unspliced counts and RNA velocity ground truth visualized by t-SNE.
4. Spatial cell locations and cell-cell interaction ground truth.
5. Discrete cell population with added batch effects.

![Results](https://github.com/ZhangLabGT/scMultiSim/raw/img/img/results.png)

