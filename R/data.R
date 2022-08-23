#' distribution of kinetic parameters learned from the Zeisel UMI cortex datasets
#' @name match_params
#' @docType data
#' @usage data(param_realdata.zeisel.imputed)
#' @format a data frame.
#' @keywords datasets
#' @examples 
#' data(param_realdata.zeisel.imputed)
"match_params"


#' a pool of gene lengths to sample from
#' @name gene_len_pool
#' @docType data
#' @usage data(gene_len_pool)
#' @format a vector.
#' @keywords datasets
#' @examples 
#' data(gene_len_pool)
"gene_len_pool"


#' from transcript length to number of fragments (for the nonUMI protocol)
#' @name len2nfrag
#' @docType data
#' @usage data(len2nfrag)
#' @format a vector.
#' @keywords datasets
#' @examples 
#' data(len2nfrag)
"len2nfrag"


#' this is the density function of log(x+1), where x is the non-zero values for ATAC-SEQ data
#' @name dens_nonzero
#' @docType data
#' @usage data(dens_nonzero)
#' @format a vector.
#' @keywords datasets
#' @examples
#' data(dens_nonzero)
"dens_nonzero"


#' 100_gene_GRN is a matrix of GRN params consisting of 100 genes where: #    - column 1 is the target gene ID, #    - column 2 is the gene ID which acts as a transcription factor for the target (regulated) gene #    - column 3 is the effect of the column 2 gene ID on the column 1 gene ID
#' @name GRN_params_100
#' @docType data
#' @usage data(GRN_params_100)
#' @format a vector.
#' @keywords datasets
#' @examples
#' data(GRN_params_100)
"GRN_params_100"


#' GRN_params_1139 is a matrix of GRN params consisting of 1139 genes where: #    - column 1 is the target gene ID, #    - column 2 is the gene ID which acts as a transcription factor for the target (regulated) gene #    - column 3 is the effect of the column 2 gene ID on the column 1 gene ID
#' @name GRN_params_1139
#' @docType data
#' @usage data(GRN_params_1139)
#' @format a vector.
#' @keywords datasets
#' @examples
#' data(GRN_params_1139)
"GRN_params_1139"

