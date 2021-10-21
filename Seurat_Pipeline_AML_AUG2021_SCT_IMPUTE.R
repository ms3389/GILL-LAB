library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
# library(garnett)
library(org.Hs.eg.db)
library(DESeq2)
library(DoubletFinder)
library(cowplot)
library(kableExtra)
# setwd("C:/Users/Max/Dropbox/Work Saar/AML_August_2021")

checkpoint_files <- list.files(pattern = "\\.checkpoint$")

starting_wd = getwd()
setwd("data/")

checkpoint_files <- list.files(pattern = "\\.checkpoint$")

unlink(checkpoint_files)
setwd(starting_wd)
options(future.globals.maxSize = 100000 * 1024^2)


fig_names = "AML_AUG2021_NOIMPUTE"


if (file.exists(paste(fig_names, '_', "nCoV.integrated.RData",sep = ''))) {
  
  sink(paste('Loading checkpoint integrated samples','.checkpoint',sep = ''))
  cat('Loading checkpoint integrated samples')
  sink()
  
  
  load(paste(fig_names, '_', "nCoV.integrated.RData",sep = '')) 
  
} else {
  
  if (file.exists(paste(fig_names, '_', "nCoV.list.RData",sep = ''))) {
    
    sink(paste('Loading checkpoint processed samples','.checkpoint',sep = ''))
    cat('Loading checkpoint processed samples')
    sink()
    
    
    load(paste(fig_names, '_', "nCoV.list.RData",sep = ''))
    
  } else {
    
    #### Start here for file list only ####
    
    
    Imputation_TRUE_FALSE = TRUE
    SCTransform_TRUE_FALSE = TRUE
    Marker_genes = c("CD3D","CD19","MS4A1","CD79A","CD14","SERPINA1")
    split_list = c("Day","Patient")
    use_reference = TRUE
    cell_type_ident = TRUE
    remove_genes_HLA_IG = TRUE
    remove_genes_MT = FALSE
    remove_genes_RIBO = FALSE
    reference_name = '104860_filtered_feature_bc_matrix.h5'
    doublet_detection = TRUE
    subset_size_fig = 20000
    subset_TRUE_FALSE = FALSE
    subset_size = 2000
    SAVE_TRUE_FALSE = TRUE
    SCTransform_integration_features = 3000
    nfeatures_for_find_var_feat = 3000
    nFeature_RNA_low = 100
    nFeature_RNA_high = 30000
    min_features = 200
    min_cells = 3
    percent_mito_filter = 30
    nCount_RNA_filter = 5
    # setwd("E:/Dropbox/Work Saar/Georgia_Collaboration")
    starting_wd = getwd()
    
    
    
    setwd("data/")
    h5_files <- list.files(pattern = "\\.h5$")
    rds_files <- list.files(pattern = "\\.rds$")
    
    rds_files = rds_files[!rds_files %in% c('reference_SCT.rds')]
    
    
    data_folders <- list.dirs(recursive = FALSE)
    samples_list = c(data_folders, h5_files, rds_files)
    meta_data_file = list.files(pattern = "\\.csv$")
    if (!identical(meta_data_file, character(0))) {
      meta_data = read.csv(meta_data_file, stringsAsFactors = F)
      colnames(meta_data)[1] = 'ID'
      meta_data_list = colnames(meta_data)
      meta_data$ID = as.character(meta_data$ID)
    }
    
    if (cell_type_ident == TRUE) {
      # classifier <- readRDS("hsPBMC_20191017.RDS")
      
      cell_type_reference = readRDS("reference_SCT.rds")
      
      # cell_type_reference@assays$RNA = cell_type_reference@assays$SCT
      # 
      # cell_type_reference <- SCTransform(cell_type_reference, verbose = FALSE)
      
    }
    
    
    
    nCoV.list = list()
    feat_uni = as.character()
    options(future.globals.maxSize = 100000 * 1024^2)
    i = 0
    for(sample_s in samples_list){
      i = i + 1
      sink(paste(starting_wd, '/Processing sample ', i, ' out of ', as.character(length(samples_list)),'.checkpoint',sep = ''))
      cat('Processing sample', i, ' out of ', as.character(length(samples_list)))
      sink()
      
      
      
      
      # }
      if(endsWith(sample_s, ".h5")) {
        sample.tmp = Read10X_h5(sample_s)
      } else if(endsWith(sample_s, ".rds")) {
        sample_s_seurat = readRDS(sample_s)
        sample.tmp = sample_s_seurat@assays$RNA@counts
        sample_s_meta = sample_s_seurat@meta.data
        meta_data_list = colnames(sample_s_seurat@meta.data)
        remove(sample_s_seurat)
        gc()
      } else  {
        sample.tmp = Read10X(sample_s)
        
      }
      if (subset_TRUE_FALSE == TRUE) {
        sample.tmp = sample.tmp[,sample(colnames(sample.tmp),subset_size)]
      }
      
      if(Imputation_TRUE_FALSE == TRUE) {
        sample.tmp.foralra=CreateSeuratObject(counts = sample.tmp, min.cells = 0, min.features = 0)
        sample.tmp.foralra = RunALRA(sample.tmp.foralra)
        sample.tmp = sample.tmp.foralra@assays$alra@data
        sample.tmp = round(sample.tmp*5)
      }
      
      #### Cell Type identification ##### OLD
      # if (cell_type_ident == TRUE) {
      #   
      #   sample.tmp.seurat.reference <- CreateSeuratObject(counts = sample.tmp, min.cells = min_cells, min.features = min_features)
      #   
      #   pbmc3k <- SCTransform(pbmc3k, verbose = FALSE)
      #   
      #   
      #   
      #   
      #   gene_names = as.data.frame(row.names(sample.tmp))
      #   colnames(gene_names) = c('gene_short_name')
      #   rownames(gene_names) = row.names(sample.tmp)
      #   fd <- new("AnnotatedDataFrame", data = gene_names)
      #   sample.tmp_cds <- newCellDataSet(as.matrix(sample.tmp), featureData = fd)
      #   
      #   sample.tmp_cds <- estimateSizeFactors(sample.tmp_cds)
      #   
      #   sample.tmp_cds <- classify_cells(sample.tmp_cds, classifier,
      #                                    db = org.Hs.eg.db,
      #                                    cluster_extend = TRUE,
      #                                    cds_gene_id_type = "SYMBOL")
      #   
      #   meta_data_for_seurat = pData(sample.tmp_cds)
      #   
      #   
      # }
      allGenes <- rownames(sample.tmp)
      if (remove_genes_HLA_IG == TRUE) {
        ### Remove HLA IGK IGH IGL genes from analysis #####
        RemoveGenes <- allGenes[grepl("^[Ii][Gg][Hh]", allGenes) & !grepl("^IGHMBP", allGenes) | grepl("^[Ii][Gg][Kk]", allGenes)
                                | grepl("^[Ii][Gg][Ll]", allGenes) | grepl("^[Hh][Ll][Aa]", allGenes)]
        sample.tmp <- sample.tmp[allGenes[!allGenes %in% RemoveGenes],]
      }
      allGenes <- rownames(sample.tmp)
      if (remove_genes_MT == TRUE) {
        ### Remove Mitochodrial genes from analysis #####
        RemoveGenes <- allGenes[grepl("^MT-", allGenes)]
        sample.tmp <- sample.tmp[allGenes[!allGenes %in% RemoveGenes],]
      }
      allGenes <- rownames(sample.tmp)
      if (remove_genes_RIBO == TRUE) {
        ### Remove ribosomal genes from analysis #####
        RemoveGenes <- allGenes[grepl("^[Rr][Pp][Ss]", allGenes) | grepl("^[Rr][Pp][Ll]", allGenes)]
        sample.tmp <- sample.tmp[allGenes[!allGenes %in% RemoveGenes],]
      }
      
      
      if(endsWith(sample_s, ".rds")) {
        meta_data_for_seurat = sample_s_meta
      } 
      
      
      if (endsWith(sample_s, ".rds")) {
        sample.tmp.seurat <- CreateSeuratObject(counts = sample.tmp, min.cells = min_cells, min.features = min_features,
                                                project = sample_s, meta.data = meta_data_for_seurat)
      } else {
        sample.tmp.seurat <- CreateSeuratObject(counts = sample.tmp, min.cells = min_cells, min.features = min_features,
                                                project = meta_data[i,1])
      }
      sample.tmp.seurat[['percent.mito']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^MT-")
      sample.tmp.seurat[["percent.rb"]] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^RPS|RPL")
      sample.tmp.seurat[["percent.hg"]] <- PercentageFeatureSet(sample.tmp.seurat, "^HB[^(P)]")
      
      
      if (!identical(meta_data_file, character(0))) {
        u = 0
        for(meta in meta_data_list) {
          u = u + 1
          sample.tmp.seurat@meta.data[, meta] <- meta_data[i,u]
        }
        remove(u)
        
      }
      
      sample.tmp.seurat <- subset(x = sample.tmp.seurat, subset = nFeature_RNA > nFeature_RNA_low & nFeature_RNA < nFeature_RNA_high 
                                  & nCount_RNA > nCount_RNA_filter & percent.mito < percent_mito_filter)
      
      if (SCTransform_TRUE_FALSE == TRUE) {
        cat('Performing SCTransform on sample', i, 'out of', as.character(length(samples_list)))
        sample.tmp.seurat <- SCTransform(sample.tmp.seurat, verbose = FALSE) 
        sample.tmp.seurat <- RunPCA(sample.tmp.seurat)
        sample.tmp.seurat <- RunUMAP(sample.tmp.seurat, dims = 1:10)
        feat_uni = union(feat_uni,sample.tmp.seurat@assays$SCT@var.features)
      } else {
        cat('Performing normalization and scaling on sample', i, 'out of', as.character(length(samples_list)))
        sample.tmp.seurat <- NormalizeData(sample.tmp.seurat, verbose = FALSE)
        sample.tmp.seurat <- FindVariableFeatures(sample.tmp.seurat, selection.method = "vst", nfeatures = nfeatures_for_find_var_feat,verbose = FALSE)
        sample.tmp.seurat <- ScaleData(sample.tmp.seurat)
        sample.tmp.seurat <- RunPCA(sample.tmp.seurat)
        sample.tmp.seurat <- RunUMAP(sample.tmp.seurat, dims = 1:10)
      }
      
      
      if (cell_type_ident == TRUE) {
        
        anchors <- FindTransferAnchors(
          reference = cell_type_reference,
          query = sample.tmp.seurat,
          normalization.method = "SCT",
          reference.reduction = "spca",
          dims = 1:50
        )
        
        sample.tmp.seurat <- MapQuery(
          anchorset = anchors,
          query = sample.tmp.seurat,
          reference = cell_type_reference,
          refdata = list(
            celltype.l1 = "celltype.l1",
            celltype.l2 = "celltype.l2",
            predicted_ADT = "ADT"
          ),
          reference.reduction = "spca", 
          reduction.model = "wnn.umap"
        )
        
        
      }
      
      ##### DoubletsFinder #####
      
      if (doublet_detection == TRUE) {
        
        if (SCTransform_TRUE_FALSE == TRUE) {
          sweep.res.list <- paramSweep_v3(sample.tmp.seurat, PCs = 1:10, sct = TRUE)
        } else {
          sweep.res.list <- paramSweep_v3(sample.tmp.seurat, PCs = 1:10, sct = FALSE)
        }
        
        
        
        sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
        bcmvn <- find.pK(sweep.stats)
        
        
        ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
        # homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- sample.tmp.seurat@meta.data$ClusteringResults
        TenxGenomics_Doublet_Estimator = length(colnames(sample.tmp.seurat))*0.00078/100
        nExp_poi <- round(TenxGenomics_Doublet_Estimator*length(colnames(sample.tmp.seurat)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
        # nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
        
        ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
        if (SCTransform_TRUE_FALSE == TRUE) {
          sample.tmp.seurat <- doubletFinder_v3(sample.tmp.seurat, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
        } else {
          sample.tmp.seurat <- doubletFinder_v3(sample.tmp.seurat, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
        }
        names(sample.tmp.seurat@meta.data)[names(sample.tmp.seurat@meta.data)==
                                             colnames(sample.tmp.seurat@meta.data)[startsWith(colnames(sample.tmp.seurat@meta.data)
                                                                                              , 'DF.classifications')]] <- "DoubletFinder"
        
        ## add processed sample to list
      }
      
      nCoV.list[sample_s] = sample.tmp.seurat
      
      # rm(sample.tmp.seurat)
      # rm(sample.tmp)
      
      
      #### Add report stuff here  ##################################################################
      # sample.tmp.seurat = nCoV.list[[1]]
      
      
      
      a = VlnPlot(sample.tmp.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.rb", "percent.hg"), ncol = 5, group.by = "orig.ident", pt.size = 0)
      b = FeatureScatter(sample.tmp.seurat, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
      
      qc = ggdraw() +
        draw_plot(a, x = 0, y = 1/2, width = 1, height = 1/2) +
        draw_plot(b, x = 0, y = 0, width = 1, height = 1/2)
      
      sample_meta = meta_data[i,]
      
      total_cells = length(colnames(sample.tmp.seurat))
      
      cellcount_l1 = table(sample.tmp.seurat@meta.data$predicted.celltype.l1)
      cellcount_l2 = table(sample.tmp.seurat@meta.data$predicted.celltype.l2)
      
      
      sample.tmp.seurat <- FindNeighbors(sample.tmp.seurat, dims = 1:10)
      sample.tmp.seurat <- FindClusters(sample.tmp.seurat, set.seed(length(samples_list)))
      
      
      UMAP_celltype_l1 = DimPlot(object = sample.tmp.seurat, group.by = 'predicted.celltype.l1', label = TRUE ,pt.size = 1.5, raster=FALSE, repel = TRUE)
      UMAP_celltype_l2 = DimPlot(object = sample.tmp.seurat, group.by = 'predicted.celltype.l2', label = TRUE ,pt.size = 1.5, raster=FALSE, repel = TRUE)
      
      
      if (doublet_detection == TRUE) {
        UMAP_doublets = DimPlot(object = sample.tmp.seurat, group.by = 'DoubletFinder', label = TRUE ,pt.size = 1.5, raster=FALSE, repel = TRUE)
      } else {
        UMAP_doublets = DimPlot(object = sample.tmp.seurat, group.by = 'orig.ident', label = TRUE ,pt.size = 1.5, raster=FALSE, repel = TRUE)
      }
      
      UMAP_clusters = DimPlot(object = sample.tmp.seurat, group.by = 'seurat_clusters', label = TRUE ,pt.size = 1.5, raster=FALSE, repel = TRUE)
      
      sample.tmp.seurat.markers <- FindAllMarkers(sample.tmp.seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
      
      top5 <- sample.tmp.seurat.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
      feat_hm <- top5$gene
      
      sample.tmp.seurat <- ScaleData(sample.tmp.seurat, verbose = FALSE, assay = "SCT", features = feat_hm)
      
      Heatmap_cluster_Genes = DoHeatmap(sample.tmp.seurat, features = feat_hm, assay = "SCT", slot = 'scale.data') + NoLegend()
      
      current = getwd()
      setwd(starting_wd)
      
      files = c(paste(meta_data[i,1],'_report.R',sep = ''),paste(meta_data[i,1],'_report.md',sep = ''),paste(meta_data[i,1],'_report.html',sep = ''))
      
      for (fn in files) {
        if (file.exists(fn)) {
          #Delete file if it exists
          file.remove(fn)
        }
      }
      
      test_knit = file.copy('scRNAseq_HTML_report_knit.R',paste(meta_data[i,1],'_report.R',sep = ''))
      
      if (test_knit == TRUE) {
        knitr::spin(paste(meta_data[i,1],'_report.R',sep = ''))
      }
      setwd(current)
      
      ###################################################
      
      
      
      
      
      
      # rm(aml2_)
      gc()
      
    }
    
    
    gc()
    
    
    sink(paste('Saving processed individual samples','.checkpoint',sep = ''))
    cat('Saving processed individual samples')
    sink()
    
    setwd(starting_wd)
    
    if (SAVE_TRUE_FALSE == TRUE) {
      save.image(paste(fig_names, '_', "nCoV.list.RData",sep = ''))
      load(paste(fig_names, '_', "nCoV.list.RData",sep = ''))
      
      saveRDS(nCoV.list, paste(fig_names, '_nCoV.list.rds', sep = ''))
      nCoV.list = readRDS(paste(fig_names, '_nCoV.list.rds', sep = ''))
    }
    
    
    
    
    
    #### End here for file list only ####
    
    
    
  }  # START OF ANALYSIS IF LIST FILE FOUND ######
  #####################################################
  #####################################################
  #####################################################
  #####################################################
  #####################################################
  
  
  setwd(starting_wd)
  
  # rm(sample.tmp.foralra,sample.tmp.seurat,sample.tmp,fuckme)
  # gc()
  rm(i)
  sink(paste('Starting Integration','.checkpoint',sep = ''))
  cat('Starting Integration')
  sink()
  
  if(SCTransform_TRUE_FALSE == TRUE) {
    # for (i in names(nCoV.list)) {
    #   nCoV.list[[i]] <- SCTransform(nCoV.list[[i]], verbose = FALSE)
    # }
    
    nCoV.features <- SelectIntegrationFeatures(object.list = nCoV.list, nfeatures = SCTransform_integration_features)
    nCoV.list <- PrepSCTIntegration(object.list = nCoV.list, anchor.features = nCoV.features, 
                                    verbose = TRUE)
    if (use_reference == TRUE) {
      reference_dataset <- which(names(nCoV.list) == reference_name)
      
      nCoV.anchors <- FindIntegrationAnchors(object.list = nCoV.list, normalization.method = "SCT", 
                                             anchor.features = nCoV.features, verbose = TRUE, reference = reference_dataset)
    } else {
      nCoV.anchors <- FindIntegrationAnchors(object.list = nCoV.list, normalization.method = "SCT", 
                                             anchor.features = nCoV.features, verbose = TRUE)
      
    }
    nCoV.integrated <- IntegrateData(anchorset = nCoV.anchors, normalization.method = "SCT", 
                                     verbose = TRUE)
  } else {
    
    nCoV <- FindIntegrationAnchors(object.list = nCoV.list, dims = 1:50)
    nCoV.integrated <- IntegrateData(anchorset = nCoV, dims = 1:50,features.to.integrate = rownames(nCoV))
    
    
  }
  
  
  
  nCoV.integrated@meta.data$BTM_Others <- plyr::mapvalues(
    x = nCoV.integrated@meta.data$predicted.celltype.l1, 
    from = c("B",       "CD4 T",      "CD8 T",       "other T", "Mono",      "other",   "DC",      "NK"),
    to = c(  "B cells", "T cells",    "T cells",     "T cells", "Monocytes", "Others",  "Others",  "Others")
  )
  
  sink(paste('Starting standard analysis on integrated raw RNA data','.checkpoint',sep = ''))
  cat('Starting standard analysis on integrated raw RNA data')
  sink()
  
  
  all.genes.COVID <- rownames(nCoV.integrated)
  
  feat = intersect(rownames(nCoV.integrated),feat_uni)
  
  
  if (SCTransform_TRUE_FALSE == TRUE) {
    DefaultAssay(nCoV.integrated) <- "SCT"
    nCoV.integrated <- RunPCA(nCoV.integrated, verbose = TRUE, features = feat, assay = "SCT")
    
  } else {
    DefaultAssay(nCoV.integrated) <- "RNA"
    nCoV.integrated <- ScaleData(nCoV.integrated, features = all.genes.COVID)
    nCoV.integrated <- FindVariableFeatures(nCoV.integrated, selection.method = "vst")
    nCoV.integrated <- RunPCA(nCoV.integrated, features = VariableFeatures(object = nCoV.integrated))
  }
  
  
  nCoV.integrated <- RunUMAP(nCoV.integrated, dims = 1:10, assay = "SCT")
  
  
  nCoV.integrated <- FindNeighbors(nCoV.integrated, dims = 1:10)
  nCoV.integrated <- FindClusters(nCoV.integrated, resolution = 0.01, set.seed(length(samples_list)))
  nCoV.integrated <- FindClusters(nCoV.integrated, resolution = 0.05, set.seed(length(samples_list)))
  nCoV.integrated <- FindClusters(nCoV.integrated, resolution = 0.1, set.seed(length(samples_list)))
  nCoV.integrated <- FindClusters(nCoV.integrated, resolution = 0.2, set.seed(length(samples_list)))
  nCoV.integrated <- FindClusters(nCoV.integrated, resolution = 0.4, set.seed(length(samples_list)))
  nCoV.integrated <- FindClusters(nCoV.integrated, resolution = 0.8, set.seed(length(samples_list)))
  nCoV.integrated <- FindClusters(nCoV.integrated, resolution = 1.2, set.seed(length(samples_list)))
  nCoV.integrated <- FindClusters(nCoV.integrated, resolution = 1.8, set.seed(length(samples_list)))
  nCoV.integrated <- FindClusters(nCoV.integrated, resolution = 2.4, set.seed(length(samples_list)))
  nCoV.integrated <- FindClusters(nCoV.integrated, resolution = 3.0, set.seed(length(samples_list)))
  # meta_data_list = c('predicted.celltype.l2', 'DoubletFinder')
  nCoV.integrated@meta.data$cell_type_l1_confidence_bins <- cut(round(nCoV.integrated@meta.data$predicted.celltype.l1.score,2), breaks=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1), labels=c("0-10","11-20","21-30","31-40","41-50","51-60","61-70","71-80","81-90","91-100"))
  
  # DimPlot(aml2_, reduction = "umap", combine = TRUE, group.by = 'cell_type_confidence_bins', pt.size = 1.5, label = TRUE, repel = TRUE)
  
  nCoV.integrated@meta.data$cell_type_l2_confidence_bins <- cut(round(nCoV.integrated@meta.data$predicted.celltype.l2.score,2), breaks=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1), labels=c("0-10","11-20","21-30","31-40","41-50","51-60","61-70","71-80","81-90","91-100"))
  
  # DimPlot(aml2_, reduction = "umap", combine = TRUE, group.by = 'cell_type_confidence_bins', pt.size = 1.5, label = TRUE, repel = TRUE)
  
  
  
  if (SAVE_TRUE_FALSE == TRUE) {
    
    sink(paste('Saving integrated dataset','.checkpoint',sep = ''))
    cat('Saving integrated dataset')
    sink()
    save.image(paste(fig_names, '_', "nCoV.integrated.RData",sep = ''))
    load(paste(fig_names, '_', "nCoV.integrated.RData",sep = ''))
    
    saveRDS(nCoV.integrated, paste(fig_names, '_nCoV.integrated.rds', sep = ''))
    nCoV.integrated = readRDS(paste(fig_names, '_nCoV.integrated.rds', sep = ''))
  }
  
} # START OF ANALYSIS IF INTEGRATED FILE FOUND ######
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################

# unique(nCoV.integrated@meta.data$dataset)

sink(paste('Making plots','.checkpoint',sep = ''))
cat('Making plots')
sink()


a = VlnPlot(nCoV.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.rb", "percent.hg"), ncol = 6, group.by = "orig.ident", pt.size = 0)
b = FeatureScatter(nCoV.integrated, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)

c = ggdraw() +
  draw_plot(a, x = 0, y = 1/2, width = 1, height = 1/2) +
  draw_plot(b, x = 0, y = 0, width = 1, height = 1/2)


dpi = 200
png(file=paste(fig_names, '_QC','.png', sep = ''), width = dpi*16, height = dpi*12, units = "px",res = dpi,type='cairo')
print(c)
dev.off()

# rm(aml2_)
# gc()
# nCoV.integrated = readRDS(file = 'SP_Georgia_SCT_NOIMPUTE_nCoV.integrated.rds')

# nCoV.integrated@meta.data$Condition <- plyr::mapvalues(
#   x = nCoV.integrated@meta.data$Sample_ID,
#   from = c("T1", "T2","T3","T4","T5","T6","T7","T8"),
#   to = c("Steady State","Mobilized","Steady State","Mobilized","Steady State","Mobilized","Steady State","Mobilized")
# )
# 
# nCoV.integrated@meta.data$Patient <- plyr::mapvalues(
#   x = nCoV.integrated@meta.data$Sample_ID,
#   from = c("T1", "T2","T3","T4","T5","T6","T7","T8"),
#   to = c("Donor 1","Donor 1","Donor 2","Donor 2","Donor 3","Donor 3","Donor 4","Donor 4")
# )

meta_data_list_plots = colnames(nCoV.integrated@meta.data)
# list_unique <- lapply(aml2_@meta.data, unique)
count_unique <- rapply(nCoV.integrated@meta.data, function(x) length(unique(x))) < 40
meta_data_list_plots = meta_data_list_plots[count_unique]

for(meta in meta_data_list_plots) {
  Idents(nCoV.integrated) <- meta
  dpi = 200
  png(file=paste(fig_names, '_',paste(meta,'.png', sep = ''), sep = ''), width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
  print(DimPlot(object = nCoV.integrated, group.by = "ident", label = TRUE ,pt.size = 1.5, raster=FALSE, repel = TRUE))
  dev.off()
}

for(split_ in split_list) {
  for(meta in meta_data_list_plots) {
    Idents(nCoV.integrated) <- meta
    dpi = 200
    png(file=paste('Split_',split_,'_', fig_names, '_',paste(meta,'.png', sep = ''), sep = ''), width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
    print(DimPlot(object = nCoV.integrated, group.by = "ident",split.by = split_, label = TRUE ,pt.size = 1.5, raster=FALSE, repel = TRUE))
    dev.off()
  }
}







dpi = 200
png(file=paste(fig_names, '_',paste('Marker_Genes','.png', sep = ''), sep = ''), width = dpi*32, height = dpi*16, units = "px",res = dpi,type='cairo')
print(FeaturePlot(nCoV.integrated, reduction = "umap", features = Marker_genes))
dev.off()

nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, assay = "SCT")

FindClusters_res_list = c('SCT_snn_res.0.01', 'SCT_snn_res.0.05', 'SCT_snn_res.0.1', 'SCT_snn_res.0.2', 'SCT_snn_res.0.4', 'SCT_snn_res.0.8', 'SCT_snn_res.1.2', 'SCT_snn_res.1.8', 'SCT_snn_res.2.4', 'SCT_snn_res.3')

Idents(nCoV.integrated) = 'DoubletFinder'

nCoV.integrated.Singlets = subset(nCoV.integrated, idents = c("Singlet"))

nCoV.integrated.Singlets.subset = subset(nCoV.integrated.Singlets, cells = sample(colnames(nCoV.integrated.Singlets),subset_size_fig))

rownames(nCoV.integrated.Singlets.subset)


for (res in FindClusters_res_list) {
  if (length(unique(as.character(nCoV.integrated.Singlets@meta.data[, res]))) < 17 & length(unique(as.character(nCoV.integrated.Singlets@meta.data[, res]))) > 6) {
    Idents(nCoV.integrated.Singlets) <- res
    DefaultAssay(nCoV.integrated.Singlets) <- "SCT"
    Idents(nCoV.integrated.Singlets.subset) <- res
    DefaultAssay(nCoV.integrated.Singlets.subset) <- "SCT"
    if (file.exists(paste(fig_names, '_', res, "_AllClusterMarkers.csv", sep = ''))) {
      
      nCoV.integrated.Singlets.markers2 = read.csv2(paste(fig_names, '_', res, "_AllClusterMarkers.csv", sep = ''),sep = ',', row.names = 1)
      
    } else {
      
      nCoV.integrated.Singlets.markers <- FindAllMarkers(nCoV.integrated.Singlets, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
      
      write.csv(nCoV.integrated.Singlets.markers, file = paste(fig_names, '_', res, "_AllClusterMarkers.csv", sep = ''))
      
    }
    nCoV.integrated.Singlets.markers = read.csv2(paste(fig_names, '_', res, "_AllClusterMarkers.csv", sep = ''),sep = ',', row.names = 1)
    top5 <- nCoV.integrated.Singlets.markers %>% group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC)
    feat_hm <- top5$gene
    
    nCoV.integrated.Singlets.subset <- ScaleData(nCoV.integrated.Singlets.subset, verbose = FALSE, assay = "SCT", features = feat_hm)
    
    dpi = 200
    png(file=paste(fig_names, paste('_cluster_heatmap_',res,'.png', sep = ''), sep = ''), width = dpi*16, height = dpi*12, units = "px",res = dpi,type='cairo')
    print(DoHeatmap(nCoV.integrated.Singlets.subset, features = feat_hm, assay = "SCT", slot = 'scale.data') + NoLegend())
    dev.off()
    
  } else if (length(unique(as.character(nCoV.integrated.Singlets@meta.data[, res]))) < 30 & length(unique(as.character(nCoV.integrated.Singlets@meta.data[, res]))) >= 17) {
    Idents(nCoV.integrated.Singlets) <- res
    DefaultAssay(nCoV.integrated.Singlets) <- "SCT"
    Idents(nCoV.integrated.Singlets.subset) <- res
    DefaultAssay(nCoV.integrated.Singlets.subset) <- "SCT"
    if (file.exists(paste(fig_names, '_', res, "_AllClusterMarkers.csv", sep = ''))) {
      
      nCoV.integrated.Singlets.markers = read.csv2(paste(fig_names, '_', res, "_AllClusterMarkers.csv", sep = ''),sep = ',', row.names = 1)
      
    } else {
      
      nCoV.integrated.Singlets.markers <- FindAllMarkers(nCoV.integrated.Singlets, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
      
      write.csv(nCoV.integrated.Singlets.markers, file = paste(fig_names, '_', res, "_AllClusterMarkers.csv", sep = ''))
      
    }
    
    
    nCoV.integrated.Singlets.markers = read.csv2(paste(fig_names, '_', res, "_AllClusterMarkers.csv", sep = ''),sep = ',', row.names = 1)
    top5 <- nCoV.integrated.Singlets.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)
    feat_hm <- top5$gene
    
    nCoV.integrated.Singlets.subset <- ScaleData(nCoV.integrated.Singlets.subset, verbose = FALSE, assay = "SCT", features = feat_hm)
    
    dpi = 200
    png(file=paste(fig_names, paste('_cluster_heatmap_',res,'.png', sep = ''), sep = ''), width = dpi*16, height = dpi*12, units = "px",res = dpi,type='cairo')
    print(DoHeatmap(nCoV.integrated.Singlets.subset, features = feat_hm, assay = "SCT", slot = 'scale.data') + NoLegend())
    dev.off()
    
  }
}







sink(paste('Saving final workspace','.checkpoint',sep = ''))
cat('Saving final workspace')
sink()


save.image(paste(fig_names, '_', "nCoV.integrated.Final_Workspace.RData",sep = ''))

#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html
### CELLCHAT ###


nCoV.integrated = readRDS('AML_AUG2021_NOIMPUTE_nCoV.integrated.rds')

rm(nCoV.integrated.Patient)
rm(nCoV.integrated.Patient.10)
gc()


nCoV.integrated@meta.data$BTM_Others <- plyr::mapvalues(
  x = nCoV.integrated@meta.data$predicted.celltype.l1, 
  from = c("B",       "CD4 T",      "CD8 T",       "other T", "Mono",      "other",   "DC",      "NK"),
  to = c(  "B cells", "T cells",    "T cells",     "T cells", "Monocytes", "Others",  "Others",  "Others")
)

# unique(nCoV.integrated@meta.data$Patient)

unique(nCoV.integrated@meta.data$Day)

nCoV.integrated@meta.data$Day2 <- plyr::mapvalues(
  x = nCoV.integrated@meta.data$Day, 
  from = c("d-21 to-7",   "d-1",        "d14",        "d-2 to -7",  "d28",        "d7",         "d-21 to -7"),
  to =   c("d-21 to -7",  "d-1",        "d14",        "d-2 to -7",  "d28",        "d7",         "d-21 to -7")
)

unique(nCoV.integrated@meta.data$Day2)
nCoV.integrated@meta.data$Day = nCoV.integrated@meta.data$Day2

nCoV.integrated@meta.data$Day = factor(nCoV.integrated@meta.data$Day, levels = c("d-21 to -7", "d-2 to -7", "d-1", "d7", "d14", "d28"))

Idents(nCoV.integrated) = "Patient"

nCoV.integrated.Patients = list()


# aml2_$Condition =  factor(aml2_$Condition, levels = c("Untreated CD45n", "UTD CD45n", "CARMA CD45n", "Untreated CD45p", "UTD CD45p", "CARMA CD45p"))


for(Patient in unique(nCoV.integrated@meta.data$Patient)) {
  
  nCoV.integrated.Patient = subset(nCoV.integrated, idents = Patient, invert = FALSE)
  
  # saveRDS(nCoV.integrated.Patient, paste('nCoV.integrated.',Patient,'.rds',sep = ''))
  
  nCoV.integrated.Patients[Patient] = nCoV.integrated.Patient
  
}

meta_data_list_plots = colnames(nCoV.integrated@meta.data)

count_unique <- rapply(nCoV.integrated@meta.data, function(x) length(unique(x))) < 40
# count_unique2 <- rapply(nCoV.integrated@meta.data, function(x) length(unique(x))) > 1

# count_unique*count_unique2
# 
# meta_data_list_plots = meta_data_list_plots[count_unique*count_unique2]


meta_data_list_plots = meta_data_list_plots[count_unique]

meta_data_list_plots = meta_data_list_plots[-8]
split_list = c("Day")


for(name in names(nCoV.integrated.Patients)) {
  
  nCoV.integrated.Patient = nCoV.integrated.Patients[[name]]
  
for(split_ in split_list) {
  for(meta in meta_data_list_plots) {
    Idents(nCoV.integrated.Patient) <- meta
    dpi = 200
    png(file=paste(name,'_Split_',split_,'_', fig_names, '_',paste(meta,'.png', sep = ''), sep = ''), width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
    print(DimPlot(object = nCoV.integrated.Patient, group.by = "ident",split.by = split_, label = TRUE ,pt.size = 1.5, raster=FALSE, repel = TRUE))
    dev.off()
  }
}

}

for(name in names(nCoV.integrated.Patients)) {
  
  nCoV.integrated.Patient = nCoV.integrated.Patients[[name]]
  
  
  length(as.character(unique(nCoV.integrated.Patient@meta.data$Day)))
  
  Idents(nCoV.integrated.Patient) = "BTM_Others"
  
  for(celltype in as.character(unique(Idents(nCoV.integrated.Patient)))) {
  
    
    nCoV.integrated.Patient.celltype = subset(nCoV.integrated.Patient, idents = celltype, invert = FALSE)
    
    Idents(nCoV.integrated.Patient.celltype) = "Day"
    
    Days = as.character(unique(Idents(nCoV.integrated.Patient.celltype)))
   for(i in 1:(length(as.character(unique(nCoV.integrated.Patient@meta.data$Day)))-1)) {
    
     
     markers = FindMarkers(nCoV.integrated.Patient.celltype, ident.1 = Days[i+1], ident.2 = Days[i])
     write.csv(markers, file = paste(fig_names, '_',name,'_',celltype,'_',Days[i+1],'_vs_ctrl_',Days[i], "_Markers.csv", sep = ''))
     
    
    
    
  }
  
  }
  
  for(split_ in split_list) {
    for(meta in meta_data_list_plots) {
      Idents(nCoV.integrated.Patient) <- meta
      dpi = 200
      png(file=paste(name,'_Split_',split_,'_', fig_names, '_',paste(meta,'.png', sep = ''), sep = ''), width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
      print(DimPlot(object = nCoV.integrated.Patient, group.by = "ident",split.by = split_, label = TRUE ,pt.size = 1.5, raster=FALSE, repel = TRUE))
      dev.off()
    }
  }
  
}
split_list = c("Day")
marker_genes = c('MKI67','IL2','TNF','IFNG','CSF2','IL6','IL10','TGFB1','TGFB2','TGFB3','IL2RA','TNFRSF1A','IFNGR1','CSF2RA','IL6R','IL10RA','IL10RB','TGFBR1','TGFBR2')

for(name in names(nCoV.integrated.Patients)) {
  
  nCoV.integrated.Patient = nCoV.integrated.Patients[[name]]
  
  nCoV.integrated.Patient@meta.data$Day = factor(nCoV.integrated.Patient@meta.data$Day, levels = c("d-21 to -7", "d-2 to -7", "d-1", "d7", "d14", "d28"))
  

  
  for(split_ in split_list) {
    for(gene in marker_genes) {
      Idents(nCoV.integrated.Patient) <- 'predicted.celltype.l2'
      dpi = 200
      png(file=paste(name,'_Split_',split_,'_', fig_names, '_',paste(gene,'_Imputed.png', sep = ''), sep = ''), width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
      print(FeaturePlot(object = nCoV.integrated.Patient, features = gene, split.by = split_, label = TRUE ,pt.size = 1.5, raster=FALSE, repel = TRUE, order = TRUE))
      dev.off()
    }
  }
  
}
genes = as.data.frame(rownames(nCoV.integrated))

marker_genes = c('MKI67','IL2','TNF','IFNG','CSF2','IL6','IL10','TGFB1','TGFB2','TGFB3','IL2RA','TNFRSF1A','IFNGR1','CSF2RA','IL6R','IL10RA','IL10RB','TGFBR1','TGFBR2')



library(SeuratWrappers)

split_ = c("Day")

FeaturePlot(object = nCoV.integrated.Patient, features = marker_genes, label = TRUE ,pt.size = 1.5, raster=FALSE, repel = TRUE, order = TRUE)

nCoV.integrated.Patient.10 = readRDS('nCoV.integrated.Patient 10.rds')
nCoV.integrated.Patient.5 = nCoV.integrated.Patients[['Patient 5']]

nCoV.integrated.Patient.10 = RunALRA(nCoV.integrated.Patient.10, assay = 'SCT')
# 
nCoV.integrated.Patient = nCoV.integrated.Patient.10

name = 'Patient 10'

rm(nCoV.integrated.Patients)
rm(nCoV.integrated.Patient)
gc()
# 
# nCoV.integrated.Patient
# 
# 
# meta_data_list_plots
# 
# 
# Idents(nCoV.integrated) = "orig.ident"
# 
# nCoV.integrated.104848 = subset(nCoV.integrated, idents = '104848')
# 
# nCoV.integrated.104849 = subset(nCoV.integrated, idents = '104849')
# 
# nCoV.integrated.104850 = subset(nCoV.integrated, idents = '104850')
# 
# 
# Idents(nCoV.integrated.104848) = "predicted.celltype.l1"
# Idents(nCoV.integrated.104849) = "predicted.celltype.l1"
# Idents(nCoV.integrated.104850) = "predicted.celltype.l1"
# 
# 
# 
# 
# nCoV.integrated.104848 = subset(nCoV.integrated.104848, idents = 'other', invert = TRUE)
# 
# nCoV.integrated.104849 = subset(nCoV.integrated.104849, idents = 'other', invert = TRUE)
# 
# nCoV.integrated.104850 = subset(nCoV.integrated.104850, idents = 'other', invert = TRUE)
# 
# saveRDS(nCoV.integrated.104848, 'nCoV.integrated.104848.rds')
# saveRDS(nCoV.integrated.104849, 'nCoV.integrated.104849.rds')
# saveRDS(nCoV.integrated.104850, 'nCoV.integrated.104850.rds')
# 
# nCoV.integrated.104848 = readRDS('nCoV.integrated.104848.rds')
# nCoV.integrated.104849 = readRDS('nCoV.integrated.104849.rds')
# nCoV.integrated.104850 = readRDS('nCoV.integrated.104850.rds')
# 
# 
# Idents(nCoV.integrated.104848) = "predicted.celltype.l1"
# 
# nCoV.integrated.104848 = subset(nCoV.integrated.104848, idents = 'other', invert = TRUE)
# 
# 
# 
# install.packages('NMF')
# devtools::install_github("jokergoo/circlize")
# devtools::install_github("jokergoo/ComplexHeatmap")
# 
# devtools::install_github("sqjin/CellChat")
# 
# library(CellChat)
# library(patchwork)
# options(stringsAsFactors = FALSE)
# 
# 
# CellChat.104848 = createCellChat(nCoV.integrated.104850, group.by = "predicted.celltype.l1")
# cellchat = CellChat.104848
# 
# CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
# showDatabaseCategory(CellChatDB)
# 
# # use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# # use all CellChatDB for cell-cell communication analysis
# # CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# 
# # set the used database in the object
# cellchat@DB <- CellChatDB.use
# 
# 
# 
# 
# 
# # subset the expression data of signaling genes for saving computation cost
# cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
# future::plan("multiprocess", workers = 4) # do parallel
# #> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
# #> in future (>= 1.13.0) when running R from RStudio, because it is
# #> considered unstable. Because of this, plan("multicore") will fall
# #> back to plan("sequential"), and plan("multiprocess") will fall back to
# #> plan("multisession") - not plan("multicore") as in the past. For more details,
# #> how to control forked processing or not, and how to silence this warning in
# #> future R sessions, see ?future::supportsMulticore
# cellchat <- identifyOverExpressedGenes(cellchat)
# cellchat <- identifyOverExpressedInteractions(cellchat)
# # project gene expression data onto PPI network (optional)
# cellchat <- projectData(cellchat, PPI.human)
# 
# # help(subsetData)
# 
# cellchat <- computeCommunProb(cellchat)
# # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
# cellchat <- filterCommunication(cellchat, min.cells = 5)
# 
# cellchat <- computeCommunProbPathway(cellchat)
# 
# cellchat <- aggregateNet(cellchat)
# 
# groupSize <- as.numeric(table(cellchat@idents))
# par(mfrow = c(1,2), xpd=TRUE)
# netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",arrow.size = 1,arrow.width = 1,alpha.edge = 1)
# netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",arrow.size = 1,arrow.width = 1,alpha.edge = 1)
# 
# mat <- cellchat@net$weight
# par(mfrow = c(3,4), xpd=TRUE)
# for (i in 1:nrow(mat)) {
#   mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#   mat2[i, ] <- mat[i, ]
#   netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
# }
# 
# pathways.show <- c("MIF")
# # Hierarchy plot
# # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells
# vertex.receiver = seq(2,3) # a numeric vector.
# netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# 
# 
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# 
# # Chord diagram
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# #> Note: The first link end is drawn out of sector 'Inflam. FIB'.
# 
# # Heatmap
# par(mfrow=c(1,1))
# netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
# #> Do heatmap based on a single object
# 
# pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
# LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# # Hierarchy plot
# vertex.receiver = seq(1,4) # a numeric vector
# netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# 
# 
# netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
# 
# 
# # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# netVisual_bubble(cellchat, sources.use = 1, targets.use = c(5:11), remove.isolate = FALSE)
# #> Comparing communications on a single object
# 
# 
# 
# # Access all the signaling pathways showing significant communications
# pathways.show.all <- cellchat@netP$pathways
# # check the order of cell identity to set suitable vertex.receiver
# levels(cellchat@idents)
# vertex.receiver = seq(1,4)
# for (i in 1:length(pathways.show.all)) {
#   # Visualize communication network associated with both signaling pathway and individual L-R pairs
#   netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
#   # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
#   gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
#   ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
# }
# 
# 
# plotGeneExpression(cellchat, signaling = "MIF")
# 
# 
# # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# # show all the interactions sending from Inflam.FIB
# netVisual_chord_gene(cellchat, sources.use = 1, targets.use = c(1:11), lab.cex = 0.5,legend.pos.y = 30)
# 
# netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 10)
# 
# 
# 
# # Compute the network centrality scores
# cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
# netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
# 
# 
# # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
# ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
# ht1 + ht2
# 
# ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = cellchat@netP$pathways)




# cat('Saving processed integrated data')
# 
# Idents(nCoV.integrated) = 'BTM_Others'
# 
# if (SAVE_TRUE_FALSE == TRUE) {
# 
#     nCoV.integrated.T = subset(nCoV.integrated, idents = c("T cells"))
#     nCoV.integrated.B = subset(nCoV.integrated, idents = c("B cells"))
#     nCoV.integrated.M = subset(nCoV.integrated, idents = c("Monocytes"))
# 
# }
# 
# if (SAVE_TRUE_FALSE == TRUE) {
#   nCoV.integrated.T = subset(nCoV.integrated, idents = c("T cells"))
#   nCoV.integrated.B = subset(nCoV.integrated, idents = c("B cells"))
#   nCoV.integrated.M = subset(nCoV.integrated, idents = c("Monocytes"))
# }
# 
# saveRDS(nCoV.integrated.T, paste(fig_names, '_nCoV.integrated.T.rds', sep = ''))
# saveRDS(nCoV.integrated.B, paste(fig_names, '_nCoV.integrated.B.rds', sep = ''))
# saveRDS(nCoV.integrated.M, paste(fig_names, '_nCoV.integrated.M.rds', sep = ''))
# 
# if (file.exists(paste(fig_names, '_nCoV.integrated.T.rds', sep = ''))) {
#   nCoV.integrated.T = readRDS(paste(fig_names, '_nCoV.integrated.T.rds', sep = ''))
#   Idents(nCoV.integrated.T) = Marker_MetaData
#   Markers.T = FindMarkers(nCoV.integrated.T, ident.1 = Marker_Condition, ident.2 = Marker_Control)
#   write.csv(Markers.T, file = paste(fig_names, '_nCoV.Markers.T.csv', sep = ''))
# }
# if (file.exists(paste(fig_names, '_nCoV.integrated.B.rds', sep = ''))) {
#   nCoV.integrated.B = readRDS(paste(fig_names, '_nCoV.integrated.B.rds', sep = ''))
#   Idents(nCoV.integrated.B) = Marker_MetaData
#   Markers.B = FindMarkers(nCoV.integrated.B, ident.1 = Marker_Condition, ident.2 = Marker_Control)
#   write.csv(Markers.B, file = paste(fig_names, '_nCoV.Markers.B.csv', sep = ''))
# }
# if (file.exists(paste(fig_names, '_nCoV.integrated.M.rds', sep = ''))) {
#   nCoV.integrated.M = readRDS(paste(fig_names, '_nCoV.integrated.M.rds', sep = ''))
#   Idents(nCoV.integrated.M) = Marker_MetaData
#   Markers.M = FindMarkers(nCoV.integrated.M, ident.1 = Marker_Condition, ident.2 = Marker_Control)
#   write.csv(Markers.M, file = paste(fig_names, '_nCoV.Markers.M.csv', sep = ''))
# }




# Marker_Control = 'RLN'
# Marker_Condition = 'HL'
# 
# 
# 
# # fig_names = "SP_CLL_SCT_2000_features_remove_genes"
# 
# nCoV.integrated = readRDS(paste(fig_names, '_nCoV.integrated.rds', sep = ''))
# 
# cat('Starting ALRA on integrated data')
# 
# 
# nCoV.integrated = RunALRA(nCoV.integrated)
# 
# cat('Starting standard analysis on integrated ALRA data')
# all.genes.COVID <- rownames(nCoV.integrated)
# 
# DefaultAssay(nCoV.integrated) <- "alra"
# nCoV.integrated <- ScaleData(nCoV.integrated, features = all.genes.COVID, assay = "alra")
# nCoV.integrated <- FindVariableFeatures(nCoV.integrated, selection.method = "vst", assay = "alra")
# nCoV.integrated <- RunPCA(nCoV.integrated, features = VariableFeatures(object = nCoV.integrated), assay = "alra")
# nCoV.integrated <- RunUMAP(nCoV.integrated, dims = 1:10, assay = "alra",reduction.name = "umap_alra")
# 
# 
# for(meta in meta_data_list_plots) {
#   Idents(nCoV.integrated) <- meta
#   dpi = 200
#   png(file=paste(fig_names, paste('Dimplot_ALRA_',meta,'.png', sep = ''), sep = ''), width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
#   print(DimPlot(object = nCoV.integrated, reduction = "umap_alra", group.by = "ident", label = TRUE ,pt.size = 1.5))
#   dev.off()
# }
# paste('Saving processed integrated ALRA data')
# 
# saveRDS(nCoV.integrated, paste(fig_names, '_nCoV.integrated.imputed.rds', sep = ''))
# 
# Idents(nCoV.integrated) = 'BTM_Others'
# 
# nCoV.integrated.imputed.T = subset(nCoV.integrated, idents = c("T cells"))
# nCoV.integrated.imputed.B = subset(nCoV.integrated, idents = c("B cells"))
# nCoV.integrated.imputed.M = subset(nCoV.integrated, idents = c("Monocytes"))
# 
# saveRDS(nCoV.integrated.imputed.T, paste(fig_names, '_nCoV.integrated.imputed.T.rds', sep = ''))
# saveRDS(nCoV.integrated.imputed.B, paste(fig_names, '_nCoV.integrated.imputed.B.rds', sep = ''))
# saveRDS(nCoV.integrated.imputed.M, paste(fig_names, '_nCoV.integrated.imputed.M.rds', sep = ''))
# 
# 
# SP_CLL_SCT_2000_features_remove_genes_nCoV.integrated <- readRDS("SP_CLL_SCT_2000_features_remove_genes_nCoV.integrated.rds")
# nCoV.integrated.CD3D <- readRDS("E:/iCloud Drive/Work Saar/Puneeth_HL/nCoV.integrated.CD3D.rds")
# 
# DimPlot(nCoV.integrated.CD3D, group.by = "cluster_ext_type", split.by = "Type", label = TRUE)
# DimPlot(nCoV.integrated.CD3D, group.by = "RNA_snn_res.0.8", split.by = "Type", label = TRUE)
# 









# nCoV.integrated.CD3D.markers <- FindAllMarkers(nCoV.integrated.CD3D, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# 
# Markers_by_Cluster <- nCoV.integrated.CD3D.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
# 
# top10 <- nCoV.integrated.CD3D.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# feat <- top10$gene
# DoHeatmap(nCoV.integrated.CD3D, features = feat) + NoLegend()
# 
# Idents(nCoV.integrated) = 'BTM_Others'
# 
# nCoV.integrated.8.15 <- FindMarkers(nCoV.integrated.CD3D, min.pct = 0.25, logfc.threshold = 0.25, ident.1 = c(8,15), ident.2 = c(1,2,3,4,5,6,7,9,11,12,13))
# 
# 
# write.csv(nCoV.integrated.8.15, file = paste('nCoV.integrated.8.15.csv', sep = ''))
# 
# write.csv(top10, file = paste('top10.csv', sep = ''))
# 
# write.csv(nCoV.integrated.CD3D.markers, file = paste('nCoV.integrated.CD3D.markers.csv', sep = ''))

# SP_CLL_SCT_2000_features_remove_genes_nCoV.integrated.B <- readRDS("E:/iCloud Drive/Work Saar/CLL/SP_CLL_SCT_2000_features_remove_genes_nCoV.integrated.B.rds")
# 
# meta_data_list_plots
# 
# for(meta in meta_data_list_plots) {
#   Idents(nCoV.integrated) <- meta
#   dpi = 200
#   png(file=paste(fig_names, paste('Dimplot_SCT_',meta,'.png', sep = ''), sep = ''), width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
#   print(DimPlot(object = nCoV.integrated, reduction = "umap_SCT", group.by = "ident", label = TRUE ,pt.size = 1.5))
#   dev.off()
# }

