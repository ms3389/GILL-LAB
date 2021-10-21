library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(garnett)
library(org.Hs.eg.db)
library(DESeq2)
library(DoubletFinder)
library(stringr)
library(data.table)
# install.packages("reticulate")
# library(reticulate)
# conda_create("CellPhoneDB", python_version = "3.7")
# # conda_install("CellPhoneDB", "CellPhoneDB")
# 
# reticulate::pip_install("CellPhoneDB")
# 
# 
# 
# 
# use_condaenv("CellPhoneDB")
# 
# os <- import("os")
# os$getcwd()
# 
# 
# py_run_string("cellphonedb")
# 
# 
# py_install("cellphonedb")
# 
# 
# shell("dir", wait = TRUE, intern = TRUE)
fig_names = "AML_AUG2021_NOIMPUTE"

setwd("C:/Users/Max/Dropbox/Work Saar/AML_August_2021")
nCoV.integrated = readRDS(paste(fig_names, '_', "nCoV.integrated.rds",sep = ''))

Idents(nCoV.integrated) = 'orig.ident'

id = unique(nCoV.integrated@meta.data$orig.ident)[1]

nCoV.integrated@meta.data$BTM_Others <- plyr::mapvalues(
  x = nCoV.integrated@meta.data$predicted.celltype.l2, 
  from = c('HSPC', "CD8 TEM",   "CD4 TCM",   "CD4 Naive", "CD14 Mono", "CD16 Mono",      "NK",      "NK Proliferating",   "CD8 Naive",    "CD4 CTL",   "CD4 TEM",   "CD8 TCM",   "CD4 Proliferating",   "CD8 Proliferating", "Treg"),
  to = c(  'HSPC', "T cells",   "T cells",   "T cells",   "Monocytes", "Monocytes",      "NK",      "NK"              ,   "T cells"  ,    "T cells",   "T cells",   "T cells",   "T cells"          ,   "T cells",           "T cells")
)

unique(nCoV.integrated@meta.data$BTM_Others)
unique(Sample@meta.data$BTM_Others)



Filenames = list()

for(id in unique(nCoV.integrated@meta.data$orig.ident)) {
  
  Sample = subset(nCoV.integrated, ident = id)
  
  Idents(Sample) = 'BTM_Others'
  
  Sample = subset(Sample, idents = c("NK","T cells","Monocytes","HSPC"))
  
  Sample = subset(Sample, ident = id)
  Filename = paste(unique(Sample@meta.data$Patient), unique(Sample@meta.data$Day),sep = '_'  )
  Filenames = c(Filenames, Filename)


 colz = colnames(Sample@assays$SCT@counts)

 colnames(Sample@assays$SCT@counts) = gsub("-", "_", colnames(Sample@assays$SCT@counts) )

 colz2 = colnames(Sample@assays$SCT@counts)

 # olnames(Sample@assays$RNA@counts))-2)
 # colnames(Sample@assays$RNA@scale.data) = str_sub(colnames(Sample@assays$RNA@scale.data), 1, str_lecolnames(Sample@assays$RNA@counts) = str_sub(colnames(Sample@assays$RNA@counts), 1, str_length(cngth(colnames(Sample@assays$RNA@scale.data))-2)


 # take raw data and normalise it

 # colnames(count_norm) = str_sub(colnames(count_norm), 1, str_length(colnames(count_norm))-2)

 # count_norm <- apply(count_raw, 2, function(x) (x/sum(x))*10000)
 # a = as.data.frame(Sample@assays$SCT@counts)
 # rownames(a)
 fwrite(as.data.frame(Sample@assays$SCT@counts), paste(Filename,".txt",sep = ''), sep="\t", quote=F, row.names=TRUE)

 # generating meta file
 meta_data <- cbind(colnames(Sample@assays$SCT@counts), Sample@meta.data[,'BTM_Others', drop=F])

 meta_data = arrange(meta_data, BTM_Others)
 #####  cluster is the user's specific cluster column
 write.table(meta_data, paste(Filename,"_meta.txt",sep = ''), sep="\t", quote=F, row.names=F)


   file_ = Filenames[[1]]

}

saveRDS(Filenames, "Filenames.rds")
Filenames = readRDS("Filenames.rds")


for(file_ in Filenames) {
  

system("cmd.exe", input = paste('cd "C:\\Users\\Max\\Dropbox\\Work Saar\\AML_August_2021" && conda activate cpdb && cellphonedb method statistical_analysis "',file_,'_meta.txt" "',file_,'.txt" --counts-data=hgnc_symbol --threads=6 --output-path="',file_,'"',sep = ''))


  }
# system("cmd.exe", input = paste('cd "C:\\Users\\Max\\Dropbox\\Work Saar\\AML_August_2021" && conda activate C:\\Users\\Max\\AppData\\Local\\r-miniconda\\envs\\CellPhoneDB && cellphonedb method statistical_analysis "Patient 13_d14_meta.txt" "Patient 13_d14.txt" --counts-data=hgnc_symbol --threads=6 --output-path=',Filename,sep = ''))

og_wd = getwd()
library(ggplot2)




for(file_ in Filenames) {

setwd(file_)

dot_plot = function(selected_rows = NULL,
                    selected_columns = NULL,
                    filename = 'plot.pdf',
                    width = 16,
                    height = 100,
                    means_path = 'means.txt',
                    pvalues_path = 'pvalues.txt',
                    means_separator = '\t',
                    pvalues_separator = '\t',
                    output_extension = '.pdf'
){
  
  all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
  all_means = read.table(means_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)
  
  intr_pairs = all_pval$interacting_pair
  all_pval = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]
  
  if(is.null(selected_rows)){
    selected_rows = intr_pairs
  }
  
  if(is.null(selected_columns)){
    selected_columns = colnames(all_pval)
  }
  
  sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]
  
  df_names = expand.grid(selected_rows, selected_columns)
  pval = unlist(sel_pval)
  pval[pval==0] = 0.0009
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  pr[pr==0] = 1
  plot.data = cbind(plot.data,log2(pr))
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
  
  my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
  
  ggplot(plot.data,aes(x=clusters,y=pair)) +
    geom_point(aes(size=-log10(pvalue),color=mean)) +
    scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text=element_text(size=14, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_blank(),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
  
  if (output_extension == '.pdf') {
    ggsave(filename, width = width, height = height, device = cairo_pdf, limitsize=F)
  }
  else {
    ggsave(filename, width = width, height = height, limitsize=F)
  }
}

dot_plot()

setwd(og_wd)


}

for(file_ in Filenames) {
  
  
  system("cmd.exe", input = paste('cd "C:\\Users\\Max\\Dropbox\\Work Saar\\AML_August_2021" && conda activate cpdb && cellphonedb plot heatmap_plot "',file_,'_meta.txt" --output-path="',file_,'"',sep = ''))
  
  
}












library(pheatmap)
heatmaps_plot = function(meta_file, pvalues_file, count_filename, log_filename, count_network_filename, interaction_count_filename, count_network_separator, interaction_count_separator, show_rownames = T, show_colnames = T,
                         scale="none", cluster_cols = T,border_color='white', cluster_rows = T, fontsize_row=11,
                         fontsize_col = 11, main = '',treeheight_row=0, family='Arial', treeheight_col = 0,
                         col1 = "dodgerblue4", col2 = 'peachpuff', col3 = 'deeppink4', meta_sep='\t', pvalues_sep='\t', pvalue=0.05){
  #######   Network
  
  meta = read.csv(meta_file, comment.char = '', sep=meta_sep)
  
  all_intr = read.table(pvalues_file, header=T, stringsAsFactors = F, sep=pvalues_sep, comment.char = '', check.names = F)
  intr_pairs = all_intr$interacting_pair
  all_intr = all_intr[,-c(1:11)]
  
  
  split_sep = '\\|'
  join_sep = '|'
  
  pairs1_all = unique(meta[,2])
  
  pairs1 = c()
  for (i in 1:length(pairs1_all))
    for (j in 1:length(pairs1_all))
      pairs1 = c(pairs1,paste(pairs1_all[i],pairs1_all[j],sep=join_sep))
  
  all_count = matrix(ncol=3)
  colnames(all_count) = c('SOURCE','TARGET','count')
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]
    
    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    
    if(p1!=p2)
      count1 = length(unique(n1))+length(unique(n2))
    else
      count1 = length(unique(n1))
    
    new_count = c(p1,p2,count1)
    names(new_count) = c('SOURCE','TARGET','count')
    all_count = rbind(all_count, new_count)
  }
  
  all_count = all_count[-1,]
  write.table(all_count, count_network_filename, sep=count_network_separator, quote=F, row.names = F)
  
  #######   count interactions
  
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]
    
    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    if(p1!=p2)
      count1 = c(count1,length(unique(n1))+length(unique(n2)))
    else
      count1 = c(count1,length(unique(n1)))
    
  }
  
  if (any(count1)>0)
  {
    count_matrix = matrix(count1, nrow=length(unique(meta[,2])), ncol=length(unique(meta[,2])))
    rownames(count_matrix)= unique(meta[,2])
    colnames(count_matrix)= unique(meta[,2])
    
    all_sum = rowSums(count_matrix)
    all_sum = cbind(names(all_sum), all_sum)
    write.table(all_sum, file=interaction_count_filename, quote=F, sep=count_network_separator, row.names=F)
    
    col.heatmap <- colorRampPalette(c(col1,col2,col3 ))( 1000 )
    
    pheatmap(count_matrix, show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
             border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
             main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = count_filename)
    
    pheatmap(log(count_matrix+1), show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
             border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
             main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = log_filename)
  } else {
    stop("There are no significant results using p-value of: ", pvalue, call.=FALSE)
  }
}


heatmaps_plot()
############################### CELLPHONEDB #####################
# Sample <- RunPCA(Sample, features = VariableFeatures(object = Sample))
# Sample <- RunUMAP(Sample, dims = 1:10)
# 
# Sample <- FindNeighbors(Sample, dims = 1:10)
# Sample <- FindClusters(Sample, resolution = 0.2)
# 
# 
# 
# 
# 
# 
# DimPlot(object = Sample, reduction = "umap", group.by = "ident", pt.size = 1, do.return = TRUE)
# FeaturePlot(Sample, reduction = "umap" ,features = c("huCD19scFV","CD3D","CD19","CD79A","CD14","SERPINA1")) 
# 
# Sample@meta.data$CellType <- plyr::mapvalues(
#   x = Sample@meta.data$orig.ident, 
#   from = c("Sample1", "Sample2", "Sample4", "Sample5","Sample6", "Sample7", "Sample14", "Sample15", "Sample17", "Sample18"), 
#   to = c("Day10", "Day10", "Day0", "Day7", "Day7", "Day10", "Day0", "Day7", "Day0", "Day7")
# )
# 
# DimPlot(object = Sample, reduction = "umap", group.by = "CellType", pt.size = 1, do.return = TRUE)

###############################
###############################
###############################
Filename = 'counts_AML1'
setwd('E:/iCloud Drive/Work Saar/11_12_2020_New_Data_AML_GvH/AML')


Sample = readRDS(file = 'ALL_GENES_SP_AML_SCT_NOIMPUTE_nCoV.integrated.rds')

Idents(Sample) = 'Patient'

Sample = subset(Sample, idents = 'AML_1')

Idents(Sample) = 'cluster_ext_type'

Sample = subset(Sample, idents = 'Unknown', invert = TRUE)

Idents(Sample) = 'DoubletFinder'

Sample = subset(Sample, idents = 'Singlet')


saveRDS(Sample, file = 'aml1_Singlets_noUknown.rds')

Sample = readRDS('aml1_Singlets_noUknown.rds')

colz = colnames(Sample@assays$SCT@counts)

colnames(Sample@assays$SCT@counts) = gsub("-", "_", colnames(Sample@assays$SCT@counts) )

colz2 = colnames(Sample@assays$SCT@counts)

# olnames(Sample@assays$RNA@counts))-2)
# colnames(Sample@assays$RNA@scale.data) = str_sub(colnames(Sample@assays$RNA@scale.data), 1, str_lecolnames(Sample@assays$RNA@counts) = str_sub(colnames(Sample@assays$RNA@counts), 1, str_length(cngth(colnames(Sample@assays$RNA@scale.data))-2)


# take raw data and normalise it

# colnames(count_norm) = str_sub(colnames(count_norm), 1, str_length(colnames(count_norm))-2)

# count_norm <- apply(count_raw, 2, function(x) (x/sum(x))*10000)
write.table(Sample@assays$SCT@counts, paste(Filename,".txt",sep = ''), sep="\t", quote=F)

# generating meta file
meta_data <- cbind(colnames(Sample@assays$SCT@counts), Sample@meta.data[,'cluster_ext_type', drop=F])  

meta_data = arrange(meta_data, cluster_ext_type)
#####  cluster is the user's specific cluster column
write.table(meta_data, paste(Filename,"_meta.txt",sep = ''), sep="\t", quote=F, row.names=F)