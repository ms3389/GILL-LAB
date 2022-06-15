install.packages("immunarch") 
library(immunarch)
library(plyr)
library(readr)
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
# library(garnett)
library(org.Hs.eg.db)
library(DESeq2)
library(DoubletFinder)
library(stringr)

# data(immdata) 
# repOverlap(immdata_10x$data) %>% vis() 
# 
# repDiversity(immdata$data) %>% vis(.by = "Status", .meta = immdata$meta) 
path1 = 'C:\\Users\\Max\\Dropbox\\Work Saar\\AML_August_2021\\data_vdj\\Patient_2'
path1 = 'C:\\Users\\Max\\Dropbox\\Work Saar\\11_12_2020_New_Data_AML_GvH\\GvH_TCR_Immunarch2\\New folder'
setwd(path1)
immdata_10x <- repLoad(path1)

immdata_10x

immdata_10x$data = immdata_10x$data[c(2,1,3)]
immdata_10x$meta = immdata_10x$meta[c(2,1,3),]


repExplore(immdata_10x$data, "lens") %>% vis()

repClonality(immdata_10x$data, "homeo") %>% vis()


rO = repOverlap(immdata_10x$data) 


vis(rO) 

repDiversity(immdata_10x$data) %>% vis(.by = "Sample", .meta = immdata_10x$meta) 

exp_cnt <- repExplore(immdata_10x$data, .method = "count")

vis(exp_cnt)

tc1 <- trackClonotypes(immdata_10x$data, list('Patient_2_d14', 5), .col = "aa")

vis(tc1)

tc2 <- trackClonotypes(immdata_10x$data, list('Patient_2_d-21 to -7', 5), .col = "aa")

vis(tc2, .order = c(2, 1, 3))

pr.nt <- pubRep(immdata_10x$data, "aa", .verbose = F)

pr.nt[is.na(pr.nt)] <- 0

sample_counts = pr.nt[,3:ncol(pr.nt)]

pr.nt$total_count = rowSums(sample_counts[])



setwd('E:\\Dropbox\\Work Saar\\11_12_2020_New_Data_AML_GvH\\GvH_TCR_Immunarch2')

# pr.nt = arrange(pr.nt, desc(total_count))

pr.nt = pr.nt[with(pr.nt, order(-total_count)), ]

pr.nt.10 = pr.nt[1:10,]

top10_clones = as.list(pr.nt.10[['CDR3.aa']])

i = 1

write.csv(pr.nt, 'pr.nt.csv')


nCoV.integrated.Patient = nCoV.integrated.Patients$`Patient 2`
# GvH = readRDS("E:\\Dropbox\\Work Saar\\11_12_2020_New_Data_AML_GvH\\GvH\\SP_GvH_SCT_NOIMPUTE_nCoV.integrated.rds")

seurat_meta = nCoV.integrated.Patient@meta.data

myfiles = list.files(path=path1, pattern="*.csv", full.names=TRUE)

cells = ldply(myfiles, read_csv)

cells = cells[c('barcode','cdr3')]

cells$clonotype = paste0("Other", cells$clonotype)

cells_list = as.list(cells[['barcode']])

aaa = as.list(tc1$CDR3.aa)

cells_list_seurat = as.list(colnames(nCoV.integrated.Patient))

for (clone in aaa) {
  
  sequences = as.character(unlist(strsplit(clone, ";", fixed = TRUE)))
  
  clonotype = clone
  
    
  for (cell_ in cells_list) {
      
      if (length(intersect(cells[cells$barcode == cell_,'cdr3'],sequences))>0) {
      
        cells$clonotype[cells$barcode == cell_] <- clonotype
      
      }
    
  }
}

# write.csv(cells, 'cells.csv')

nCoV.integrated.Patient$predicted.celltype.l1

cells_unique = distinct(cells, barcode, .keep_all = TRUE)

rownames(cells_unique) = cells_unique$barcode

cells_unique_clonotype = cells_unique

cells_unique_clonotype$barcode = NULL
cells_unique_clonotype$cdr3 = NULL
table(cells_unique_clonotype$clonotype)


seurat_meta_clono = merge(seurat_meta, cells_unique_clonotype, by = 0, all.x= TRUE)

rownames(seurat_meta_clono) = seurat_meta_clono$Row.names
seurat_meta_clono$Row.names = NULL
# seurat_meta_clono_copy = seurat_meta_clono
# seurat_meta_clono = seurat_meta_clono_copy
# write.csv(seurat_meta_clono, 'seurat.csv')


seurat_meta_clono = seurat_meta_clono %>% replace(is.na(.), 'N/A')

seurat_meta_clono_actual = seurat_meta_clono[, c('clonotype','predicted.celltype.l1')]
seurat_meta_clono_actual$predicted.celltype.l1 = NULL

unique(seurat_meta_clono_actual$clonotype)

GvH = AddMetaData(nCoV.integrated.Patient, seurat_meta_clono_actual, col.name = 'Clonotype_forward')
###################################################################################################
###################################################################################################
###################################################################################################
seurat_meta = GvH@meta.data

myfiles = list.files(path=path1, pattern="*.csv", full.names=TRUE)

cells = ldply(myfiles, read_csv)

cells = cells[c('barcode','cdr3')]

cells$clonotype = paste0("Other", cells$clonotype)

cells_list = as.list(cells[['barcode']])

aaa2 = as.list(tc2$CDR3.aa)

cells_list_seurat = as.list(colnames(GvH))

for (clone in aaa2) {
  
  sequences = as.character(unlist(strsplit(clone, ";", fixed = TRUE)))
  
  clonotype = clone
  
  
  for (cell_ in cells_list) {
    
    if (length(intersect(cells[cells$barcode == cell_,'cdr3'],sequences))>0) {
      
      cells$clonotype[cells$barcode == cell_] <- clonotype
      
    }
    
  }
}

# write.csv(cells, 'cells.csv')



cells_unique = distinct(cells, barcode, .keep_all = TRUE)

rownames(cells_unique) = cells_unique$barcode

cells_unique_clonotype = cells_unique

cells_unique_clonotype$barcode = NULL
cells_unique_clonotype$cdr3 = NULL

seurat_meta_clono = merge(seurat_meta, cells_unique_clonotype, by = 0, all.x= TRUE)

rownames(seurat_meta_clono) = seurat_meta_clono$Row.names
seurat_meta_clono$Row.names = NULL
# seurat_meta_clono_copy = seurat_meta_clono
# seurat_meta_clono = seurat_meta_clono_copy
# write.csv(seurat_meta_clono, 'seurat.csv')


seurat_meta_clono = seurat_meta_clono %>% replace(is.na(.), 'N/A')

seurat_meta_clono_actual = seurat_meta_clono[, c('clonotype','cluster_ext_type')]
seurat_meta_clono_actual$cluster_ext_type = NULL



GvH = AddMetaData(GvH, seurat_meta_clono_actual, col.name = 'Clonotype_backward')




###################################################################################################
###################################################################################################
###################################################################################################
# GvH@meta.data = seurat_meta_clono

extra = c("N/A", "Other")

order_aaa = as.character(aaa)

order_aaa2 = as.character(aaa2)

order_new = append(extra, order_aaa)

order_new2 = append(extra, order_aaa2)

colors_new = c("grey50", "grey40","#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3")

Idents(GvH) = 'cluster_ext_type'

GvH_CD4 = subset(GvH, idents = c('CD4 T cells'))
GvH_CD8 = subset(GvH, idents = c('CD8 T cells'))



dpi = 200
png(file='CLONOTYPE_GvHCD4_START5_3.png', width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
print(DimPlot(object = GvH_CD4, group.by = 'Clonotype_forward', split.by = 'Timepoint',pt.size = 1.5, order = rev(order_new), cols = colors_new))
dev.off()

dpi = 200
png(file='CLONOTYPE_GvHCD4_END5_3.png', width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
print(DimPlot(object = GvH_CD4, group.by = 'Clonotype_backward', split.by = 'Timepoint',pt.size = 1.5, order = rev(order_new2), cols = colors_new))
dev.off()

dpi = 200
png(file='CLONOTYPE_GvHCD8_START5_3.png', width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
print(DimPlot(object = GvH_CD8, group.by = 'Clonotype_forward', split.by = 'Timepoint',pt.size = 1.5, order = rev(order_new), cols = colors_new))
dev.off()

dpi = 200
png(file='CLONOTYPE_GvHCD8_END5_3.png', width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
print(DimPlot(object = GvH_CD8, group.by = 'Clonotype_backward', split.by = 'Timepoint',pt.size = 1.5, order = rev(order_new2), cols = colors_new))
dev.off()


Idents(GvH_CD8) = 'Clonotype_forward'


GvH_CD8_CVV = subset(GvH_CD8, idents = c("CVVDDAASKLTF;CASSAGTVYEQYF"))

Idents(GvH_CD8_CVV) = 'Timepoint'


mm_CVV = FindMarkers(GvH_CD8_CVV, ident.1 = c('GVHD TP_2'), ident.2 = c('GVHD TP_1'))

write.csv(mm_CVV, file = "mm_CVV.csv")




colors = c("#F8766D", "#DB8E00", "#AEA200", "#64B200", "#00BD5C", "#00C1A7", "#00BADE", "#00A6FF", "#B385FF", "#EF67EB", "grey50")
order = c("N/A", "Other", "CVVSDVSYSGGGADGLTF;CASSLTPVGDEQYF", "CATCLYGNNRLAF;CASSFGPYTEAFF", "CAGHLNNNARLMF;CASSFHGEKLFF", "CVVDDAASKLTF;CASSAGTVYEQYF", "CANTRFWTPNTGELFF", "CAISGTGETEKLFF", "CAESETGANNLFF;CAISGTGETEKLFF", "CASSQYEGGPYNEQFF", "CASSLGWGHTSSYEQYF")  



rm(clonotype)
unique(GvH@meta.data$Clonotype)

