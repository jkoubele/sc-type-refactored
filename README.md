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
* Fixed issue when genes with zero expression caused NA values in the normalization step.

## Usage



## Copyright and license
This repository is re-implementation of a code provided in [sc-type](https://github.com/IanevskiAleksandr/sc-type) GitHub repository, originally licensed under GPLv3.
When using, please also refer to the published [article](https://www.nature.com/articles/s41467-022-28803-w).
