





library("devtools")
devtools::install_github("broadinstitute/infercnv")



BiocManager::install("infercnv")
install.packages('rjags')
library(cowplot)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(dplyr)
library(stringr)
library(rjags)
library(infercnv)
library(data.table)

#C:\Users\Max\

setwd("C:/Users/Max/Dropbox/Work Saar/AML_August_2021")

setwd('E:\\Dropbox\\Work Saar\\11_12_2020_New_Data_AML_GvH\\AML2')
setwd('C:\\Users\\Max\\Dropbox\\Work Saar\\11_12_2020_New_Data_AML_GvH\\AML2')

# seurat_obj = readRDS('SP_AML_SCT_NOIMPUTE_nCoV.integrated.rds')
# 
# Idents(seurat_obj) = 'Patient'
# 
# seurat_obj = subset(seurat_obj, idents = c('AML_2'))
# 
# Idents(seurat_obj) = 'cluster_ext_type'
# 
# seurat_obj = subset(seurat_obj, idents = c('T cells','CD8 T cells','CD4 T cells', 'CD34+'))
# 
# Patient = 'AML2'

Patient = ("Patient_10")
seurat_obj = readRDS('nCoV.integrated.Patient 10.rds')


FeaturePlot(seurat_obj, features = c('CD3D','CD34','CD8','CD4','CD19'))

a = as.data.frame(rownames(seurat_obj))

Idents(seurat_obj) = 'BTM_Others'

seurat_obj = subset(x = seurat_obj, subset = CD34 > 0 | CD3D > 0)

seurat_obj_T <- subset(seurat_obj, subset = CD3D > 0)

## Get cell names
cellNames <- rownames(seurat_obj_T@meta.data)
rm(seurat_obj_T)
gc()

## Mutate a column in original SO metadata
seurat_obj$barcode <- rownames(seurat_obj@meta.data)
seurat_obj@meta.data <- seurat_obj@meta.data %>% mutate(inferCNV_Groups = ifelse((seurat_obj$barcode %in% cellNames), "T cells",  "CD34+"))


# seurat_obj <- subset(seurat_obj, cells = sample(colnames(seurat_obj), size = 100))



# seurat_obj@meta.data$Proper_Patient <- plyr::mapvalues(
#   x = seurat_obj@meta.data$Patient,
#   from = c("Patient 3", "Patient 19"),
#   to = c("Patient 3", "Healthy")
# )




# counts_matrix = as.matrix(seurat_obj@raw.data[,seurat_obj@cell.names])

# counts_matrix = as.data.frame(seurat_obj@assays$SCT@counts)

# write.table(rep,round(counts_matrix, digits=3), file='sc.10x.counts.matrix', quote=F, sep="\t")        

# use more palatable column names (cell identifiers)            
# cell.names <- sapply(seq_along(colnames(counts_matrix)), function(i) paste0("cell_", i), USE.NAMES = F)      
# colnames(counts_matrix) = colnames(seurat_obj@assays$RNA@counts)    
# rownames(counts_matrix) = rownames(seurat_obj@assays$RNA@counts)    

# write the output table 
# write.table(round(counts_matrix, digits=3), file='sc.10x.counts.matrix', quote=F, sep="\t")        
fwrite(round(as.data.frame(seurat_obj@assays$SCT@counts), digits=3), paste(Patient, 'sc.10x.counts.matrix.txt',sep = '.'), quote=FALSE, sep = "\t", row.names = TRUE)
##################################################################################################################################################################

seurat_obj_cells = WhichCells(seurat_obj)
seurat_obj_annot = seurat_obj@meta.data$inferCNV_Groups

seurat_obj_df <- data.frame(seurat_obj_cells,seurat_obj_annot)

# write.table(seurat_obj_df, "cellAnnotations.txt", sep="\t")
fwrite(seurat_obj_df, paste(Patient, "cellAnnotations.txt",sep = '.'), quote=FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

##################################################################################################################################################################
seurat_obj_cells = WhichCells(seurat_obj)
seurat_obj_annot = seurat_obj@meta.data$cluster_ext_type

seurat_obj_df <- data.frame(seurat_obj_cells,seurat_obj_annot)

# write.table(seurat_obj_df, "cellAnnotations.txt", sep="\t")
fwrite(seurat_obj_df, paste(Patient, "cellAnnotations.txt",sep = '.'), quote=FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
##################################################################################################################################################################


# fwrite(round(counts_matrix[0:100,0:100], digits=3), 'sc.10x.counts.mat.txt', quote=FALSE, sep = "\t", row.names = TRUE)
hg38_pos = read.table("hg38_pos.txt",sep="\t", header=FALSE)
# 
hg38_pos$V2 = strsplit(as.character(hg38_pos$V2), "|", fixed=TRUE)

hg38_pos$V2 = strsplit(as.character(hg38_pos$V2), "_", fixed=TRUE)


hg38_pos$V2 = str_extract(hg38_pos$V2, regex("\\w\\w\\w(\\d{1,2}|[XY])"))

# 
# hg38_pos2 = hg38_pos[1:4]
# 
fwrite(hg38_pos, "hg38_pos_chr.txt", quote=FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste(Patient, "sc.10x.counts.matrix.txt",sep = '.'),
                                    annotations_file=paste(Patient, "cellAnnotations.txt",sep = '.'),
                                    delim="\t",
                                    gene_order_file="hg38_pos_chr.txt",
                                    ref_group_names=c('T cells'),
                                    chr_exclude = c("chrM"))
rm(seurat_obj)
gc()
# fact = unique(as.character(infercnv_obj@gene_order$chr))

infercnv_obj@gene_order$chr = factor(infercnv_obj@gene_order$chr, levels=c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'))


# install.packages("igraph")
# 
# gc() 
# 
# library(reticulate)
# py_config()
# reticulate::conda_install("r-reticulate", "python-igraph", channel = "vtraag")
# 
# 
# 



cnvdir = paste("InferCNV_Out_",Patient,sep = '')

infercnv_obj = infercnv::run(infercnv_obj,
                             no_plot = T,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir=cnvdir,  # dir is auto-created for storing outputs
                             cluster_by_groups=FALSE,   # cluster
                             denoise=T,
                             HMM=T,
                             analysis_mode='subclusters',
                             no_prelim_plot=TRUE,
                             plot_steps=FALSE,
                             num_threads=4)


copy_infercnv_obj = infercnv_obj
new_gene_order = data.frame()
for (chr_name in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")) {
  new_gene_order = rbind(new_gene_order, infercnv_obj@gene_order[which(infercnv_obj@gene_order[["chr"]] == chr_name) , , drop=FALSE])
}
names(new_gene_order) <- c("chr", "start", "stop")
copy_infercnv_obj@gene_order = new_gene_order
copy_infercnv_obj@expr.data = infercnv_obj@expr.data[rownames(new_gene_order), , drop=FALSE]


plot_cnv(copy_infercnv_obj,
         cluster_by_groups=FALSE,
         cluster_references=TRUE,
         out_dir=cnvdir,
         output_filename=paste(Patient,"_rev_chr_infercnv.png",sep='')
)







seurat_obj = add_to_seurat(seurat_obj = seurat_obj,infercnv_output_path = cnvdir)

colnames(seurat_obj@meta.data)

DimPlot(seurat_obj, group.by = "has_cnv_chr11",pt.size = 1)
DimPlot(seurat_obj, group.by = "has_cnv_chr1",pt.size = 1)
DimPlot(seurat_obj, group.by = "top_loss_4",pt.size = 1)


DimPlot(seurat_obj, group.by = "cluster_ext_type",pt.size = 1)

##########################################################################################################################################
##########################################################################################################################################

library(devtools)
install_github("navinlabcode/copykat")
library(copykat)



library(cowplot)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(dplyr)
library(stringr)
library(infercnv)
library(data.table)

#C:\Users\Max\
setwd('E:\\Dropbox\\Work Saar\\11_12_2020_New_Data_AML_GvH\\AML2')
# setwd('C:\\Users\\Max\\Dropbox\\Work Saar\\11_12_2020_New_Data_AML_GvH\\AML2')

seurat_obj = readRDS('SP_AML_SCT_NOIMPUTE_nCoV.integrated.rds')

Idents(seurat_obj) = 'Patient'

seurat_obj = subset(seurat_obj, idents = c('AML_2'))

Idents(seurat_obj) = 'cluster_ext_type'

seurat_obj = subset(seurat_obj, idents = c('T cells','CD8 T cells','CD4 T cells', 'CD34+'))

seurat_obj_norm_names = colnames(subset(seurat_obj, idents = c('T cells','CD8 T cells','CD4 T cells')))






exp.rawdata <- as.matrix(seurat_obj@assays$RNA@counts)

rm(seurat_obj)
rm(copykat.test)
gc()

save.image(file = 'pre_copykat_data.rds')
load(file = 'pre_copykat_data.rds')

copykat.test <- copykat(rawmat=exp.rawdata, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, sam.name="test", distance="euclidean", norm.cell.names=seurat_obj_norm_names, n.cores=1)

save.image(file = 'copykat_aml2.rds')
load(file = 'copykat_aml2.rds')

pred.test <- data.frame(copykat.test$prediction)
CNA.test <- data.frame(copykat.test$CNAmat)

meta = seurat_obj@meta.data

jointdataset <- merge(meta, pred.test, by = 0, all.x= TRUE)

rownames(jointdataset) = jointdataset$Row.names

seurat_obj@meta.data = jointdataset

DimPlot(object = seurat_obj, group.by = 'cluster_ext_type', label = TRUE ,pt.size = 1.5, raster=FALSE, repel = TRUE)
DimPlot(object = seurat_obj, group.by = 'copykat.pred', label = TRUE ,pt.size = 1.5, raster=FALSE, repel = TRUE)


my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

chr <- as.numeric(CNA.test$chrom) %% 2+1
rbPal1 <- colorRampPalette(c('black','grey'))
CHR <- rbPal1(2)[as.numeric(chr)]
chr1 <- cbind(CHR,CHR)

rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
com.preN <- pred.test$copykat.pred
pred <- rbPal5(2)[as.numeric(factor(com.preN))]

cells <- rbind(pred,pred)
col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))

heatmap.3(t(CNA.test[,4:ncol(CNA.test)]),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
          ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
          notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
          keysize=1, density.info="none", trace="none",
          cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
          symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))

legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=0.6, bty="n")





tumor.cells <- pred.test$cell.names[which(pred.test$copykat.pred=="aneuploid")]
tumor.mat <- CNA.test[, which(colnames(CNA.test) %in% tumor.cells)]
hcc <- hclust(parallelDist::parDist(t(tumor.mat),threads =4, method = "euclidean"), method = "ward.D2")
hc.umap <- cutree(hcc,2)

rbPal6 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4])
subpop <- rbPal6(2)[as.numeric(factor(hc.umap))]
cells <- rbind(subpop,subpop)

heatmap.3(t(tumor.mat),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
          ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
          notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
          keysize=1, density.info="none", trace="none",
          cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
          symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))

legend("topright", c("c1","c2"), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4], cex=0.9, bty='n')
