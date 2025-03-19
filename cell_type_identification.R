library(rjson)
library(Seurat)
library(tidyverse)


cell_type_scoring <- function(seurat_object, markers, assay = NULL, layer = NULL) {
  markers_positive <- markers$positive_markers
  markers_negative <- markers$negative_markers

  # In the Seurat data, consider only marker genes
  markers_to_use <- intersect(rownames(seurat_object), unique(unlist(markers)))
  seurat_object_subset <- subset(seurat_object, features = markers_to_use)
  markers_expression_matrix <- as.matrix(GetAssayData(seurat_object_subset, layer = layer, assay = assay))

  # z-scale
  z_score_matrix <- t(scale(t(markers_expression_matrix)))
  z_score_matrix <- na.omit(z_score_matrix)

  # Consider only markers found in Seurat data
  for (cell_type in names(markers_positive)) {
    markers_positive[[cell_type]] <- intersect(markers_positive[[cell_type]], rownames(z_score_matrix))
    markers_negative[[cell_type]] <- intersect(markers_negative[[cell_type]], rownames(z_score_matrix))
  }

  # marker sensitivity
  marker_occurrences <- sort(table(unlist(markers_positive)), decreasing = T)
  marker_sensitivity <- data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_occurrences),
                                                                             to = c(0, 1),
                                                                             from = c(length(markers_positive), 1)),
                                  gene = names(marker_occurrences),
                                  stringsAsFactors = FALSE)
  rownames(marker_sensitivity) <- marker_sensitivity$gene

  # multiple by marker sensitivity
  for (gene in rownames(z_score_matrix)) {
    z_score_matrix[gene,] <- z_score_matrix[gene,] * marker_sensitivity[gene, 'score_marker_sensitivity']
  }

  # Combine scores for each cell type
  marker_to_cell_type_matrix <- matrix(0,
                                       nrow = nrow(z_score_matrix),
                                       ncol = length(markers_positive))
  rownames(marker_to_cell_type_matrix) <- rownames(z_score_matrix)
  colnames(marker_to_cell_type_matrix) <- names(markers_positive)

  for (cell_type in names(markers_positive)) {
    marker_to_cell_type_matrix[markers_positive[[cell_type]], cell_type] <- 1 / sqrt(max(length(markers_positive[[cell_type]]), 1))
    if (length(markers_negative[[cell_type]]) > 0) {
      marker_to_cell_type_matrix[markers_negative[[cell_type]], cell_type] <- -1 / sqrt(max(length(markers_negative[[cell_type]]), 1))
    }
  }

  cell_type_scores <- t(t(z_score_matrix) %*% marker_to_cell_type_matrix)
  cell_type_scores
}


cluster_cell_type_clasification <- function(seurat_object, cell_type_scores) {
  cluster_scores <- do.call("rbind", lapply(unique(seurat_object@meta.data$seurat_clusters), function(cluster_id) {
    cluster_cell_indices <- which(seurat_object@meta.data$seurat_clusters == cluster_id)
    cluster_cell_type_scores <- sort(rowSums(cell_type_scores[, cluster_cell_indices]), decreasing = TRUE)
    data.frame(seurat_clusters = cluster_id,
               cell_type = names(cluster_cell_type_scores),
               scores = cluster_cell_type_scores,
               ncells = sum(seurat_object@meta.data$seurat_clusters == cluster_id))
  }))

  clusters_cell_type <- cluster_scores |>
    group_by(seurat_clusters) |>
    slice_max(order_by = scores, n = 1, with_ties = FALSE) |>
    mutate(cell_type = ifelse(scores < ncells / 4, "Unknown", cell_type)) |>
    select(seurat_clusters, cell_type)

  clusters_cell_type
}
