lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
library(rjson)

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

path_prefix <- "/cellfile/datapublic/jkoubele/ercc1/data/"
input_folder <- file.path(path_prefix, "clustering_output")
output_folder <- file.path(path_prefix, "cell_type_annotation")

seurat_object <- readRDS(file.path(input_folder, "seurat_object.rds"))

mouse_markers <- fromJSON(file='/cellfile/datapublic/jkoubele/sc-type-refactored/cell_type_markers_mouse.json')


markers_kidney <- mouse_markers$Kidney
markers_to_use <- intersect(rownames(seurat_object), unique(unlist(markers_kidney)))

markers_positive <- mouse_markers$Kidney$positive_markers
markers_negative <- mouse_markers$Kidney$negative_markers

seurat_object_subset <- subset(seurat_object, features = markers_to_use)

seurat_object_subset <- subset(seurat_object_subset, cells = colnames(seurat_object_subset)[1:1000])

markers_expression_matrix <- as.matrix(GetAssayData(seurat_object_subset, layer = 'data', assay = 'SCT'))

markers_positive_uppercase <- list()
for(name in names(markers_positive)){
  markers_positive_uppercase[[name]] <- toupper(markers_positive[[name]])
}


## Re-implementation of sctype_score:

# marker sensitivity
marker_occurrences = sort(table(unlist(markers_positive)), decreasing = T); 
marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_occurrences),
                                                                           to = c(0,1),
                                                                           from = c(length(markers_positive),1)),
                                gene = names(marker_occurrences), 
                                stringsAsFactors = FALSE)
rownames(marker_sensitivity) <- marker_sensitivity$gene

# z-scale
z_score_matrix <- t(scale(t(markers_expression_matrix)))
z_score_matrix <- na.omit(z_score_matrix)

# subselect genes only found in data
for(cell_type in names(markers_positive)){
  print(cell_type)
  print(length(markers_positive[[cell_type]]))
  markers_positive[[cell_type]] <- intersect(markers_positive[[cell_type]], rownames(z_score_matrix))
  markers_negative[[cell_type]] <- intersect(markers_negative[[cell_type]], rownames(z_score_matrix))
  print(length(markers_positive[[cell_type]]))
}

# multiple by marker sensitivity
for(gene in rownames(z_score_matrix)){
  z_score_matrix[gene,] <- z_score_matrix[gene,] * marker_sensitivity[gene, 'score_marker_sensitivity']
}

# Combine scores for each cell type
marker_to_cell_type_matrix <- matrix(0, nrow = nrow(z_score_matrix), ncol = length(markers_positive))
rownames(marker_to_cell_type_matrix) <- rownames(z_score_matrix)
colnames(marker_to_cell_type_matrix) <- names(markers_positive)
for(cell_type in names(markers_positive)){
  marker_to_cell_type_matrix[markers_positive[[cell_type]], cell_type] <- 1 / sqrt(max(length(markers_positive[[cell_type]]), 1))
  if(length(markers_negative[[cell_type]]) > 0){
    marker_to_cell_type_matrix[markers_negative[[cell_type]], cell_type] <- marker_to_cell_type_matrix[markers_negative[[cell_type]], cell_type] - 
      1 / sqrt(max(length(markers_negative[[cell_type]]), 1))
  }
}

es <- t(t(z_score_matrix) %*% marker_to_cell_type_matrix)
# Combine scores for each cell type
es_original <- do.call(rbind, lapply(names(markers_positive), function(cell_type) { 
  vapply(1:ncol(z_score_matrix), function(j) {
    
    # Get positive and negative markers
    pos_genes <- markers_positive[[cell_type]]
    neg_genes <- markers_negative[[cell_type]]
    
    # Compute scores only if genes are available
    sum_t1 <- if (length(pos_genes) > 0) sum(z_score_matrix[pos_genes, j]) / sqrt(length(pos_genes)) else 0
    sum_t2 <- if (length(neg_genes) > 0) sum(z_score_matrix[neg_genes, j] * -1) / sqrt(length(neg_genes)) else 0
    
    sum_t1 + sum_t2
    
  }, numeric(1))  # Ensures output is a numeric vector
}))

# # Assign row names (cell types)
# rownames(es) <- names(markers_positive)
# 
# 
# 
# # Assign row names (cell types)
# rownames(es) <- names(markers_positive)
