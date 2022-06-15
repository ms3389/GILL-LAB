install.packages('pheatmap')

n_clusters = 8

setwd("E:/iCloud Drive/Work Saar/11_12_2020_New_Data_AML_GvH/AML2")

cluster_names = c('B_cells','CD34','CD4_T_cells','CD8_T_cells','Dendritic_cells','Monocytes','NK_cells','T_cells')

file.create("dotplot.sh")
u = 1
for (i in 1:(n_clusters)) {
  pvalues = read.csv('cellphoneDB_out/pvalues.txt', sep = '\t')
  
  pvalues_base = pvalues[, 1:11]
  
  pvalues_clusters = pvalues[, 12:length(colnames(pvalues))]
  
  pval_coln = colnames(pvalues_clusters)
  
  # cluster_name = strtrim(pval_coln[u+i-1],floor(nchar(pval_coln[u+i-2])/2))
  cluster_name = cluster_names[i]
  
  dataframe_out = cbind(pvalues_base, pvalues_clusters[,u:(u+n_clusters-1)])
  
  write.table(dataframe_out, paste('cellphoneDB_out/','pvalues_',cluster_name,'.txt', sep = ''), sep = '\t')
  
  
  dataframe_condition = dataframe_out
  
  dataframe_condition$condition = rowSums(dataframe_out[ , c(12:19)])
  
  dataframe_condition = dataframe_condition[!(dataframe_condition$condition==8),]
  
  write.table(dataframe_condition$interacting_pair, paste('cellphoneDB_out/rows_in_',cluster_name,'.txt',sep = ''), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  
  
  
  means = read.csv('cellphoneDB_out/means.txt', sep = '\t')
  
  means_base = means[, 1:11]
  
  means_clusters = means[, 12:length(colnames(means))]
  
  means_coln = colnames(means_clusters)
  
  # cluster_name = strtrim(means_coln[u+i-1],floor(nchar(means_coln[u+i-2])/2))
  
  cluster_name = cluster_names[i]
  
  dataframe_out2 = cbind(means_base, means_clusters[,u:(u+n_clusters-1)])
  
  write.table(dataframe_out2, paste('cellphoneDB_out/','means_',cluster_name,'.txt', sep = ''), sep = '\t') 
  
  u = u + n_clusters
  
  write(paste('cellphonedb plot dot_plot --means-path cellphoneDB_out/means_',cluster_name,'.txt --rows ', 'cellphoneDB_out/rows_in_',cluster_name,'.txt', ' --pvalues-path cellphoneDB_out/pvalues_',cluster_name,'.txt --output-name ',cluster_name,'_dot_plot.pdf' , sep = ''), "dotplot.sh",append = TRUE)
  
  cat(i)
}



