**Table of contents**

* [Results](#results)
* [Installation](#installation)
* [Basic Workflow](#basic-workflow)
* [Options](#options)
  + [Options: General](#options-general)
  + [Options: Genes](#options-genes)
  + [Options: Cells](#options-cells)
  + [Options: CIF](#options-cif)
  + [Options: Simulation - ATAC](#options-simulation---atac)
  + [Options: Simulation - RNA](#options-simulation---rna)
  + [Options: Simulation - Spatial Cell-Cell Interaction](#options-simulation---spatial-cell-cell-interaction)
* [Technical Noise and Batch Effects](#technical-noise-and-batch-effects)
* [Output](#output)
* [FAQ](#faq)
* [Contact](#contact)


# scMultiSim

scMultiSim is an in silico simulator that generates multi-modality data of single-cells, including gene expression, chromatin accessibility, RNA velocity, and spatial location of cells. It takes a cell differential tree and a gene regulatory network (GRN) as input, and simulates spliced and unspliced counts while accounting for the relationships between modalities. The output single cell gene expression data is determined by three factors: cell-cell interactions, within-cell GRNs and chromatin accessibility. Users can tune the effect of each factor on the output data and set various parameters for the underlying model. Furthermore, the GRN can be set in a time-varying mode where the network's structure changes temporally to reflect the dynamic nature of biological networks. We also provide options to simulate technical variations such as batch effects. scMultiSim can be used to benchmark challenging computational tasks on single-cell multi-omics data, including the inference of GRNs, estimation of RNA velocity, integration of single-cell datasets from multiple batches and modalities, and analysis of cell-cell interaction using the cell spatial location data.

Please refer to the following vignettes for examples and tutorials.
- [scMultisim Basics](https://zhanglabgt.github.io/scMultiSim/doc/basics.html)
- [Simulating Spatial Cell-Cell Interactions](https://zhanglabgt.github.io/scMultiSim/doc/spatialCCI.html)

Alternatively, you can also browse the documentation and vignettes after installing the package.

![Overview](https://github.com/ZhangLabGT/scMultiSim/raw/img/img/scMultisim.png)

#### Please Cite

```
Li, H., Zhang, Z., Squires, M., Chen, X., & Zhang, X. (2023). scMultiSim: simulation of multi-modality single cell data guided by cell-cell interactions and gene regulatory networks. BioRxiv, 2022.10.15.512320. https://doi.org/10.1101/2022.10.15.512320
```

## Results

The following figure briefly shows results from the same cell differential tree:

1. Connected scATAC-seq and scRNA-seq, in continuous or discrete mode. Visualized by t-SNE.
2. GRN correlation heatmap, where genes regulated by the same regulator have similar correlations with others.
3. Unspliced counts and RNA velocity ground truth visualized by t-SNE.
4. Spatial cell locations and cell-cell interaction ground truth.
5. Discrete cell population with added batch effects.

![Results](https://github.com/ZhangLabGT/scMultiSim/raw/img/img/results.png)

## Installation

Please install from this repo and specify the `main` branch:

```R
devtools::install_github("ZhangLabGT/scMultiSim@main")
```

## Basic Workflow

scMultiSim provides several example GRN and trajectory data:
- GRN: `GRN_params_100` and `GRN_params_1139`
- Tree: `Phyla5()`, `Phyla3()` and `Phyla1()`

A typical workflow consists these three main steps:

1. Simulate the true count:
   ```R
   results <- sim_true_counts(list2(
       GRN = GRN_params_100,
       tree = Phyla5(),
       num.cells = 1000,
       # optional options
       num.cif = 50,
       discrete.cif = F,
       cif.sigma = 1,
       do.velocity = T,
       # ... other options
   ))
   ```
2. Add technical noise to the dataset:
   ```R
   add_expr_noise(results)
   ```
3. Add batch effects:
   ```R
   divide_batches(results, nbatch = 4)
   ```

## Options

scMultiSim requires users to provide the following options:

- [`GRN`](#grn): The Gene Regulatory Network.
  If the GRN effect is unneeded, Provide `NA`.
- [`tree`](#tree): The cell differential tree.
  If the cell trajectory is unimportant,
  simply use `Phyla1()` to generate a linear tree.

Typically, you may also want to adjust the following options
to control the basic cell population:

- [`num.cells`](#numcells): Specify the number of cells.
- [`unregulated.gene.ratio`](#unregulatedgeneratio) or
  [`num.genes`](#numgenes): Control the total number of genes.
- [`discrete.cif`](#discretecif): Whether generating discrete or
  continuous cell population.
- [`diff.cif.fraction`](#diffciffraction):
  Control the contribution of the trajectory/cluster specified by the tree.
- [`cif.sigma`](#cifcenter-cifsigma):
  Control the variation of cells along the trajectory.

Here, we provide a quick guide on adjusting the effect of each biological factors:

- `diff.cif.fraction` controls the balance between trajectory and GRN.
  | Higher Value | Lower Value |
  | ------------------------------------------- | --------------------------------------------------- |
  | **Clear** trajectories / cluster boundaries | **Vague** trajectories / cluster boundaries |
  | Trajectory determined by the tree | Trajectory determined by other factors, e.g., GRN |
- `cif.sigma` control the noise along trajectory.
  | Higher Value | Lower Value |
  | ------------------------------------------- | --------------------------------------------------- |
  | **Vague** trajectories / cluster boundaries | **Clear** trajectories / cluster boundaries |
  | Cells stick to the trajectory | More randomness in the trajectory |
- `intrisic.noise` controls the randomness added during sampling mRNA counts,
  therefore also affecting the trajectory shape.
  | Higher Value | Lower Value |
  | ------------------------------------------- | --------------------------------------------------- |
  | **Clear** trajectories / cluster boundaries | **Vague** trajectories / cluster boundaries |
  | Expression contains more noise | Expression directly computed from kinectic parameters |
- `atac.effect` controls the contribution of the chtomatin accessibility.
  | Higher Value | Lower Value |
  | ------------------------------------------- | --------------------------------------------------- |
  | Higher correlation between ATAC and RNA-seq | Lower correlation between ATAC and RNA-seq |

### Options: General

#### rand.seed

> integer (default: `0`)

scMultiSim should produce the same result if all other parameters are the same.

#### threads

> integer (default: `1`)

Use multithreading only when generating the CIF matrix.
It will not speed up the simulation a lot, thus not recommended.

### Options: Genes

#### GRN

> A data frame with 3 columns as below.
> Supply `NA` to disable the GRN effect. (required)

| Column | Value                                      |
| ------ | ------------------------------------------ |
| 1      | target gene ID: `integer or character`;    |
| 2      | regulator gene ID: `integer or character`; |
| 3      | effect: `number`.                          |

If `num.genes` presents, the gene IDs should not exceed this number.
The gene IDs should start from 1 and should not ship any intermidiate numbers.

Two sample datasets `GRN_params_100` and `GRN_params_1000` from
[Dibaeinia, P., &amp; Sinha, S. (2020)](https://doi.org/10.1016/j.cels.2020.08.003) are provided for testing and inspection.

#### num.genes

> integer (default: `NULL`)

If a GRN is supplied, override the total number of genes.
It should be larger than the largest gene ID in the GRN.
Otherwise, the number of genes will be determined by `N_genes * (1 + r_u)`,
where `r_u` is `unregulated.gene.ratio`.

If GRN is disabled,
this option specifies the total number of genes.

#### unregulated.gene.ratio

> number > 0 (default: `0.1`)

Ratio of unreulated to regulated genes.
When a GRN is supplied with `N` genes,
scMultiSim will simulate `N * r_u` extra (unregulated) genes.

#### giv.mean, giv.sd, giv.prob

> (default: `0, 1, 0.3`)

The parameters used to sample the GIV matrix.
With probability `giv.prob`, the value is sampled from N(`giv.mean`, `giv.sd`).
Otherwise the value is 0.

#### dynamic.GRN

> list (default: `NULL`)

Enables dynamic (cell-specific GRN).
Run `scmultisim_help("dynamic.GRN")` to see more explaination.

#### hge.prop, hge.mean, hge.sd

> (default: `0, 5, 1`)

Treat some random genes as highly-expressed (house-keeping) genes.
A proportion of `hge.prop` genes will have expression scaled by a
multiplier sampled from N(`hge.mean`, `hge.sd`).

#### hge.range

> integer (default: `1`)

When selecting highly-expressed genes, only choose genes with ID > `hge.range`.

#### hge.max.var

> number (default: `500`)

When selecting highly-expressed genes, only choose genes
with variation < `hge.max.var`.

### Options: Cells

#### num.cells

> integer (default: `1000`)

The number of cells to be simulated.

#### tree

> phylo (default: `Phyla5()`)

The cell differential tree,
which will be used to generate cell trajectories (if `discrete.cif = T`)
or clusters (if `discrete.cif = F`).
In discrete population mode, only the tree tips will be used.
Three demo trees, `Phyla5()`, `Phyla3()` and `Phyla1()`, are provided.

#### discrete.cif

> logical (default: `FALSE`)

Whether the cell population is discrete (continuous otherwise).

#### discrete.min.pop.size, discrete.min.pop.index

> integer, integer (default: `70, 1`)

In discrete population mode, specify one cluster to have the
smallest cell population.
The cluster will contain `discrete.min.pop.size` cells.
`discrete.min.pop.index` should be a valid cluster index (tree tip number).

#### discrete.pop.size

> integer vector (default: `NA`); e.g. `c(200, 250, 300)`

Manually specify the size of each cluster.

### Options: CIF

#### num.cifs

> integer (default: `50`)

Total number of differential and non-differential CIFs,
which can be viewed as latent representation of cells.

#### diff.cif.fraction

> number (default: `0.9`)

Fraction of differential CIFs.
Differential CIFs encode the cell type information,
while non-differential CIFs are randomly sampled for each cell.

#### cif.center, cif.sigma

> (default: `1, 0.1`)

The distribution used to sample CIF values.

#### use.impulse

> logical (default: `FALSE`)

In continuous population mode, when sampling CIFs along the tree,
use the impulse model rather than the default gaussian random walk.

### Options: Simulation - ATAC

#### atac.effect

> number ∈ [0, 1] (default: `0.5`)

The influence of chromatin accessability data on gene expression.

#### region.distrib

> vector of length 3, should sum to 1 (default: `c(0.1, 0.5, 0.4)`)

The probability that a gene is regulated by 0, 1, 2
consecutive regions, respectively.

#### atac.p_zero

> number ∈ [0, 1] (default: `0.8`)

The proportion of zeros we see in the simulated scATAC-seq data.

#### riv.mean, riv.sd, riv.prob

> (default: `0, 1, 0.3`)

The parameters used to sample the RIV (Region Identity Vectors).
With probability `riv.prob`, the value is sampled from N(`riv.mean`, `riv.sd`).
Otherwise the value is 0.

### Options: Simulation - RNA

#### do.velocity

> logical (default: `FALSE`)

When set to `TRUE`,
simulate using the full kinetic model and generate RNA velocity data.
Otherwise, the Beta-Poission model will be used.

#### beta

> number (default: `0.4`)

The splicing rate of each gene in the kinetic model.

#### d

> number (default: `1`)

The degradation rate of each gene in the kinetic model.

#### num.cycles

> number (default: `3`)

The number of cycles run before sampling the gene expression of a cell.

#### cycle.len

> number (default: `1`)

In velocity mode, a multiplier for the cell cycle length.
It is multiplied by the expected time to
transition from k_on to k_off and back to form the the length of a cycle.

### Options: Simulation - Spatial Cell-Cell Interaction

#### cci

Enables cell-cell interaction. See `scmultisim_help("cci")` for details.

## Technical Noise and Batch Effects

### `add_expr_noise`

Options:

- `protocol: "umi"/"nonUMI"`: Whether simulate the UMI protocol.
- `alpha_mean`, `alpha_sd`: Mean and deviation of rate of subsampling of transcripts during capture step.
- `alpha_gene_mean`, `alpha_gene_sd`: `alpha` parameters, but gene-wise.
- `depth_mean`, `depth_sd`: Mean and deviation of sequencing depth.
- `gene_len`: A vector with lengths of all genes.
- `atac.obs.prob`: For each integer count of a particular region for a particular cell, the probability the count will be observed.
- `atac.sd.frac`: The fraction of ATAC-seq data value used as the standard deviation of added normally distrubted noise.
- `randseed`: random seed.

### `divide_batches`

Options:

- `nbatch`: Number of batches.
- `effect`: The batch effect size.

## Output

scMultiSim returns an environment with the following fields:

- `counts`: Gene-by-cell scRNA-seq counts.
- `atac_counts`: Region-by-cell scATAC-seq counts.
- `region_to_gene`: Region-by-gene 0-1 marix indicating the corresponding relationship between chtomatin regions and genes.
- `atacseq_data`: The "clean" scATAC-seq counts without added intrinsic noise.
- `cell_meta`: A dataframe containing cell type labels and pseudotime information.
- `cif`: The CIF used during the simulation.
- `giv`: The GIV used during the simulation.
- `kinetic_params`: The kinetic parameters used during the simulation.
- `.grn`: The GRN used during the simulation.
  - `.grn$regulators`: The list of TFs used by all gene-by-TF matrices.
  - `.grn$geff`: Gene-by-TF matrix representing the GRN used during the simulation.
- `.n`: Other metadata, e.g. `.n$cells` is the number of cells.

If `do.velocity` is enabled, it has these additional fields:

- `unspliced_counts`: Gene-by-cell unspliced RNA counts.
- `velocity`: Gene-by-cell RNA velocity ground truth.
- `cell_time`: The pseudotime at which the cell counts were generated.

If dynamic GRN is enabled, it has these additional fields:

- `cell_specific_grn`: A list of length `n_cells`. Each element is a gene-by-TF matrix, indicating the cell's GRN.

If cell-cell interaction is enabled, it has these additional fields:

- `grid`: The grid object used during the simulation.
  - `grid$get_neighbours(i)`: Get the neighbour cells of cell `i`.
- `cci_locs`: A dataframe containing the X and Y coordinates of each cell.
- `cci_cell_type_param`: A dataframe containing the CCI network ground truth: all ligand-receptor pairs between each pair of cell types.
- `cci_cell_types`: For continuous cell population, the sub-divided cell types along the trajectory used when simulating CCI.

If it is a debug session (`debug = TRUE`), a `sim` field is available,
which is an environment contains all internal states and data structures.

### Reference

Hechen Li, Ziqi Zhang, Michael Squires, Xi Chen, and Xiuwei Zhang. 2023. “scMultiSim: Simulation of Multi-Modality Single Cell Data Guided by Cell-Cell Interactions and Gene Regulatory Networks.” bioRxiv.

## FAQ

### Running Speed

Simulations should finish in a reasonable time in most cases. On a machine with an i7-12700K CPU and 64GB RAM, using 1000 cells, 100 genes and 50 CIFs, the simulation took under 1 mimute to generate both scRNA-seq and scATAC-seq data. If also generating unspliced and spliced counts, or enabling cell-cell interactions, the running time is longer (~3 minutes when RNA velocity is enabled, and 30 minutes for 500 cells with spatial cell-cell interaction enabled).

## Contact

GitHub issues are welcomed.
It is also possible to send email to the main author
`Hechen Li (hli691 at gatech.edu)`.
