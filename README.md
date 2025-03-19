# sc-type-refactored

Re-implementation of the [sc-type](https://github.com/IanevskiAleksandr/sc-type) tool for cell type annotation. 
When using the original implementation,
I encountered few minor issues, which are addressed here:
* For mouse data, I used [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) package
to map the marker genes to corresponding mouse orthologs. Only uniquely mapped genes are kept (one-to-many mapping genes are discarded).
The mapping is done in the script [prepare_marker_files.R](./prepare_marker_files.R).
  * The marker genes are stored as JSON files ([human](./cell_type_markers_human.json) and [mouse](./cell_type_markers_mouse.json)).
  * The cell type names were slightly renamed, by replacing special characters which caused issues when e.g. using the cell type as name of file.
  * Duplicate cell types were removed (as this was affecting the marker sensitivity calculation).
* Fixed issue when genes with zero expression caused NA values in the normalization step and subsequently propagated NAs to the cell type scores.

## Example usage

```r
library(Seurat)
library(rjson)
library(tidyverse)

# Source the code in this repository
source("https://raw.githubusercontent.com/jkoubele/sc-type-refactored/main/cell_type_identification.R")

# Prepare your Seurat object
seurat_object <- readRDS(file = "path_to_your_seurat_object.rds")

# Load the gene markers
mouse_markers <- fromJSON(file="https://raw.githubusercontent.com/jkoubele/sc-type-refactored/main/cell_type_markers_mouse.json")
human_markers <- fromJSON(file="https://raw.githubusercontent.com/jkoubele/sc-type-refactored/main/cell_type_markers_human.json")

# You can display available tissues:
print(names(mouse_markers))
# > [1] "Immune system" "Pancreas"      "Liver"         "Eye"           "Kidney"        "Brain"        
# > [7] "Lung"          "Adrenal"       "Heart"         "Intestine"     "Muscle"        "Placenta"     
# > [13] "Spleen"        "Stomach"       "Thymus"        "Hippocampus" 

# Load the markers for specific tissue
tissue <- 'Kidney'
markers <- mouse_markers[[tissue]]

# At first, we compute cell type scores for each cell
cell_type_scores <- cell_type_scoring(seurat_object, kidney_markers)

# After that, we aggregate the scores by Seurat clusters, and predict cell type for each cluster
clusters_cell_type <- cluster_cell_type_clasification(seurat_object, cell_type_scores)

# We may store the results to metadata of the Seurat object 
seurat_object@meta.data <- seurat_object@meta.data |> left_join(clusters_cell_type , by = "seurat_clusters")

```

## Copyright and license
This repository is re-implementation of a code provided in [sc-type](https://github.com/IanevskiAleksandr/sc-type) GitHub repository, originally licensed under GPLv3.
When using, please also refer to the published [article](https://www.nature.com/articles/s41467-022-28803-w).
