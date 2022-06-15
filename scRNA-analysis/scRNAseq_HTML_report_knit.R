#' # Gill Lab `r meta_data[i,1]` report

#' ## Standard SCT normalization performed

#' ## QC summary plot
#' ---
#+ fig.width=16, fig.height=9
qc
#' ## UMAP cell types L1
#' ---
#+ fig.width=16, fig.height=9
UMAP_celltype_l1
#' ## UMAP cell types L2
#' ---
#+ fig.width=16, fig.height=9
UMAP_celltype_l2
#' ## UMAP doublets detection*
#' ---
#+ fig.width=16, fig.height=9
UMAP_doublets
#' ## A UMAP of clustering
#' ---
#+ fig.width=16, fig.height=9
UMAP_clusters
#' ## A QC Gene counts summary plot
#' ---
#+ fig.width=16, fig.height=9
Heatmap_cluster_Genes
#' ## Top 5 highest expressed genes per cluster
#' ---
kbl(cbind(top5))%>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "250px")

#' ## Notes
#'
#' 1. Notes
#' 1. Notes