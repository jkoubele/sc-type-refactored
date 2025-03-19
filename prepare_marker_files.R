library(openxlsx)
library(rjson)
library(HGNChelper)
library(biomaRt)


repo_folder <- if (interactive()) {
  dirname(rstudioapi::getActiveDocumentContext()$path)  # Works in RStudio
} else {
  dirname(sys.frames()[[1]]$ofile)  # Works when sourced
}

sanitize_names <- function(cell_type_names) {
  cell_type_names |>
    trimws() |>
    gsub("\u00A0", " ", x = _) |>
    gsub("\\s+", " ", x = _) |>
    gsub("[()]", "", x = _) |>
    gsub("[ /]", "_", x = _) |>
    gsub("\\+", "_plus", x = _) |>
    gsub("_+", "_", x = _) |>
    gsub("α", "Alpha", x = _) |>
    gsub("β", "Beta", x = _) |>
    gsub("γδ", "Gamma_delta", x = _)
}


source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")

markers_excel_path <- file.path(repo_folder, 'ScTypeDB_full.xlsx')
cell_markers_df = read.xlsx(markers_excel_path)
markers_original <- list()

for (tissue in unique(cell_markers_df$tissueType)) {
  gs_list <- gene_sets_prepare(markers_excel_path, tissue)
  names(gs_list$gs_positive) <- sanitize_names(names(gs_list$gs_positive))
  names(gs_list$gs_negative) <- sanitize_names(names(gs_list$gs_negative))


  # Drop multiple occurrences of the same cell type
  gs_list$gs_positive <- gs_list$gs_positive[!duplicated(names(gs_list$gs_positive))]
  gs_list$gs_negative <- gs_list$gs_negative[!duplicated(names(gs_list$gs_negative))]

  markers_original[[tissue]] <- list(positive_markers = gs_list$gs_positive,
                                     negative_markers = gs_list$gs_negative)
}


biomart_human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
biomart_mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

gene_name_mapping <- getLDS(
  attributes = c("hgnc_symbol"),
  filters = "hgnc_symbol",
  values = unique(unlist(markers_original)),
  mart = biomart_human,
  attributesL = c("mgi_symbol"),
  martL = biomart_mouse,
  uniqueRows = TRUE
)

gene_name_mapping_unique <- gene_name_mapping[!gene_name_mapping$HGNC.symbol %in% gene_name_mapping$HGNC.symbol[duplicated(gene_name_mapping$HGNC.symbol)],]
rownames(gene_name_mapping_unique) <- gene_name_mapping_unique$HGNC.symbol

markers_mouse <- list()
for (tissue in names(markers_original)) {
  positive_markers_mouse <- list()
  negative_markers_mouse <- list()

  for (cell_type in names(markers_original[[tissue]]$positive_markers)) {
    human_markers <- markers_original[[tissue]]$positive_markers[[cell_type]]
    positive_markers_mouse[[cell_type]] <- unique(na.omit(gene_name_mapping_unique[human_markers, "MGI.symbol"]))
  }
  for (cell_type in names(markers_original[[tissue]]$negative_markers)) {
    human_markers <- markers_original[[tissue]]$negative_markers[[cell_type]]
    negative_markers_mouse[[cell_type]] <- unique(na.omit(gene_name_mapping_unique[human_markers, "MGI.symbol"]))
  }
  markers_mouse[[tissue]] <- list(positive_markers = positive_markers_mouse,
                                  negative_markers = negative_markers_mouse)
}

write(toJSON(markers_original), file = file.path(repo_folder, "cell_type_markers_human.json"))
write(toJSON(markers_mouse), file = file.path(repo_folder, "cell_type_markers_mouse.json"))
