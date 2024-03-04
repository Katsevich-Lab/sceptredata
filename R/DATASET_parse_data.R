setwd("~/parse_data/")
library(sceptre)

sceptre_object <- import_data_from_parse(
  gene_mat_fp = "gene_mat.mtx",
  grna_mat_fp = "grna_mat.mtx",
  all_genes_fp = "all_genes.csv",
  all_grnas_fp = "all_grnas.csv",
  moi = "low",
  grna_target_data_frame = data.frame(grna_id = c("guide_A", "guide_B", "guide_C"),
                                      grna_target = c(NA, NA, NA))
)

# take a subset of genes and cells
response_matrix <- get_response_matrix(sceptre_object)
grna_matrix <- get_grna_matrix(sceptre_object)
gene_ids <- rownames(response_matrix)[1:50]
n_cells <- 5000
response_matrix_sub <- response_matrix[gene_ids, seq(1, n_cells)]
grna_matrix_sub <- grna_matrix[,seq(1, n_cells)]

write_ondisc_backed_sceptre_object()
