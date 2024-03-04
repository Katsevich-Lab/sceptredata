library(sceptre)
library(Matrix)

process_matrix_and_ids <- function(mat_fp, id_fp, gene = TRUE) {
  n_cells <- 5000
  feature_df <- read.csv(file = id_fp, col.names = c("gene_id", "gene_name", "genome"))
  matrix_metadata <- sceptre:::get_mtx_metadata(mat_fp, c("n_cells", "n_features", "n_nonzero"))
  dt <- data.table::fread(file = mat_fp, skip = matrix_metadata$n_to_skip,
                          col.names = c("cell_idx", "feature_idx", "x"),
                          colClasses = c("integer", "integer", "double"),
                          showProgress = FALSE, nThread = 2L)
  dt <- dt |> dplyr::filter(cell_idx %in% seq(1, n_cells))

  if (gene) {
    mt_gene_idxs <- grep(pattern = "^MT-", x = feature_df$gene_name)
    rows_to_keep <- c(sample(x = seq(1L, nrow(feature_df)), size = 245), sample(mt_gene_idxs, 5)) |> sort()
    feature_df <- feature_df[rows_to_keep,]
    dt <- dt |> dplyr::filter(feature_idx %in% rows_to_keep)
    old_to_new_feature_map <- data.frame(feature_idx = rows_to_keep,
                                         new_feature_idx = seq(1L, length(rows_to_keep)))
    dt <- dplyr::left_join(x = dt, y = old_to_new_feature_map, by = "feature_idx") |>
      dplyr::select(cell_idx, feature_idx = new_feature_idx, x)
  }

  n_features <- nrow(feature_df)
  data.table::setorderv(dt, cols = c("cell_idx", "feature_idx"))
  dt_with_header <- rbind(data.frame(cell_idx = n_cells,
                                     feature_idx = nrow(feature_df),
                                     x = nrow(dt)), dt)
  comments <- c("%%MatrixMarket matrix coordinate integer general",
                paste0("%Rows=cells (", n_cells,"), Cols=genes (", n_features,")"))
  ex_dir <- "~/research_code/sceptredata/inst/extdata/parse_example/"
  file_path_gene_mtx <- paste0(ex_dir, if (gene) "gene_mat.mtx" else "grna_mat.mtx")
  file_path_all_genes <- paste0(ex_dir, if (gene) "all_genes.csv" else "all_grnas.csv")

  # write to directory
  writeLines(comments, con = file_path_gene_mtx)
  write.table(x = dt, file = file_path_gene_mtx, append = TRUE, row.names = FALSE, col.names = FALSE)
  rownames(feature_df) <- NULL
  write.csv(x = feature_df, file = file_path_all_genes, row.names = FALSE)
}

# call function on gene and grna modalities
set.seed(6)
mat_fp <- "~/parse_data/gene_mat.mtx"
id_fp <- "~/parse_data/all_genes.csv"
process_matrix_and_ids(mat_fp, id_fp, TRUE)
mat_fp <- "~/parse_data/grna_mat.mtx"
id_fp <- "~/parse_data/all_grnas.csv"
process_matrix_and_ids(mat_fp, id_fp, FALSE)
