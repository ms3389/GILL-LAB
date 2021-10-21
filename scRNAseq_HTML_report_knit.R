#' # Gill Lab `r meta_data[i,1]` report

#' ## Standard SCT normalization performed

#' ## Sample Info
#' ---
#' #### Total cells post filtering: `r total_cells`
#' ## Sample Metadata Information
#' ---
kbl(cbind(sample_meta))%>%
  kable_styling()
#' ## QC summary plot
#' ---
#+ fig.width=16, fig.height=9
qc
#' ## UMAP cell types L1
#' ---
#+ fig.width=16, fig.height=9
UMAP_celltype_l1
#' ## Cell Type L1 count
#' ---
kbl(cbind(cellcount_l1))%>%
  kable_styling()
#' ## UMAP cell types L2
#' ---
#+ fig.width=16, fig.height=9
UMAP_celltype_l2
#' ## Cell Type L2 count
#' ---
kbl(cbind(cellcount_l2))%>%
  kable_styling()
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