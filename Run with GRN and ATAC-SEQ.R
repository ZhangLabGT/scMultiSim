library("devtools")
load_all("scMultiSim")

# load the file target_gene_GRN_params
# target_gene_GRN_params is a matrix where:
  #    - column 1 is the target gene ID,
  #    - column 2 is the gene ID which acts as a transcription factor for the target (regulated) gene
  #    - column 3 is the effect of the column 2 gene ID on the column 1 gene ID
load("scMultiSim/data/100_gene_GRN.RData")
# load("scMultiSim/data/1139_gene_GRN.RData")

# if TRUE, saves a variety of metrics in the script directory
metrics = TRUE

results <- SimulateTrueCounts(num_cells = 200, unregulated_to_regulated_gene_ratio = 0.1, num_evfs = 500, diffEVF_fraction = 0.9, Sigma = 0.1,
                              atac_effect = 0.5, beta = 0.4, d = 1, num_cycles = 2, cycle_length = 1.0, intrinsic_noise = 1, phyla = Phyla5(), randseed = 0, do_velocity = TRUE, metrics = TRUE)

counts <- results$counts
unspliced_counts <- results$unspliced_counts
velocity <- results$velocity
cell_time <- results$cell_time
atacseq_data <- results$atacseq_data
num_genes <- results$num_genes
