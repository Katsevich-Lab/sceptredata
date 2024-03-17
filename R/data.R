#' Example high-MOI data
#'
#' The example high-MOI CRISPRi data are taken from the paper "A genome-wide framework for mapping gene regulation via cellular genetic screens" by Gasperini et al., 2019. The authors used a CRISPRi-based assay to target 5,000+ putative enhancers in a population of K562 cells. The authors additionally targeted 200+ gene transcription start sites (TSSs) and designed a library of 100 non-targeting gRNAs to serve as negative controls. The genes, gRNAs, and cells are downsampled for inclusion in the example data.
#'
#' The example data are stored in a list containing the following components:
#' - `response_matrix`: the gene-by-cell expression matrix
#' - `grna_matrix`: the gRNA-by-cell expression matrix
#' - `extra_covariates`: a data frame containing a single column, `batch`, specifying the batch in which each cell was sequenced (`batch_1` or `batch_2`)
#' - `gene_names`: the human-readable name of each gene
#'
#' @usage data(highmoi_example_data)
"highmoi_example_data"

#' gRNA target data frame
#'
#'
#' The gRNA target data frame corresponding to the example high MOI data. The data frame contains the columns `grna_id` (ID of an individual gRNA), `grna_target` (genomic target of the gRNA), `chr` (target chromosome), `start` (start coordinate of target location), and `end` (end coordinate of target location).
#' @usage data(grna_target_data_frame_highmoi)
"grna_target_data_frame_highmoi"

#' Example low-MOI data
#'
#' The example low MOI data are taken from the paper "Characterizing the molecular regulation of inhibitory immune checkpoints with multimodal single-cell screens" by Papalexi et al., 2021. The authors used CRISPRko to target 26 genes in 20,729 THP-1 cells and measured the effects of these perturbations via single-cell sequencing. The authors designed three or four gRNAs to target each gene as well as a library of nine non-targeting gRNAs. The genes were downsampled for inclusion in the example data.
#'
#' The example data are stored in a list containing four components:
#' - `response_matrix`: the gene-by-cell expression matrix
#' - `grna_matrix`: the gRNA-by-cell expression matrix
#' - `grna_target_data_frame`:  a data frame containing the columns `grna_id` (ID of an individual gRNA) and `grna_target` (genomic target of the gRNA)
#' - `extra_covariates`: a data frame with a single column, `bio_rep`, specifying the biological replicate of each cell (`rep_1`, `rep_2`, or `rep_3`)
#'
#' @usage data(lowmoi_example_data)
"lowmoi_example_data"
