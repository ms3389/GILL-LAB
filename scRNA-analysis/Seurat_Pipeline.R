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
library(SeuratWrappers)
library(ggrepel)
library(ggvenn)
library(topGO)
library(stringr)
fig_names = "AML_Batch2_3_CAR123"
setwd("C:/Users/Max/Dropbox/Work Saar/AML_Batch_2_3/")
# setwd("C:/Users/Max/Dropbox/Work Saar/AML_August_2021")
# remotes::install_github('satijalab/seurat-wrappers')
# install.packages("R.utils")
# setwd("C:/Users/Max/Dropbox/Work Saar/AML_August_2021")
# setwd("C:/Users/Max/Dropbox/Work Saar/Pipeline Testing")

# setwd("C:/Users/Max/Dropbox/Work Saar/Kim_CD38_Project/PipelineRun2")


data_directory = 'data'
# reference_directory = NA
reference_directory = '/home/mshestov/R_big_memory/CellType_References/'
# reference_directory = 'C:/Users/Max/Dropbox/Work Saar/Seurat_References/'



checkpoint_files <- list.files(pattern = "\\.checkpoint$")
unlink(checkpoint_files)

starting_wd = getwd()
setwd(data_directory)

checkpoint_files <- list.files(pattern = "\\.checkpoint$")

unlink(checkpoint_files)
setwd(starting_wd)
options(future.globals.maxSize = 100000 * 1024^2)



sink(paste('Packages loaded - checking for checkpoints','.checkpoint',sep = ''), split = TRUE)
closeAllConnections()


if (file.exists(paste(fig_names, '_', "nCoV.integrated.RData",sep = ''))) {
  
  sink(paste('Loading checkpoint integrated samples','.checkpoint',sep = ''))
  cat('Loading checkpoint integrated samples')
  closeAllConnections()
  
  
  
  load(paste(fig_names, '_', "nCoV.integrated.RData",sep = '')) 
  
} else {
  
  if (file.exists(paste(fig_names, '_', "nCoV.list.RData",sep = ''))) {
    
    sink(paste('Loading checkpoint processed samples','.checkpoint',sep = ''), split = TRUE)
    cat('Loading checkpoint processed samples')
    closeAllConnections()
    
    
    
    load(paste(fig_names, '_', "nCoV.list.RData",sep = ''))
    
  } else {
    
    #### Start here for file list only ####
    
    
    Imputation_TRUE_FALSE = FALSE
    SCTransform_TRUE_FALSE = TRUE
    Marker_genes = c("CD3D","CD19","MS4A1","CD79A","CD14","SERPINA1")
    split_list = c("Day","Patient")
    factor_1 = NA #c("Healthy_BM","d-21 to -7", "d-1", "d7", "d14", "d28")
    factor_2 = NA
    use_reference = TRUE
    report_TRUE_FALSE = TRUE
    cell_type_ident_PBMC = TRUE
    cell_type_ident_BM = TRUE
    cell_type_ident_BM_AML = TRUE
    cell_type_ident_TabulaS_PBMC = TRUE
    cell_type_ident_TabulaS_BM = TRUE
    remove_genes_HLA_IG = FALSE
    remove_genes_MT = TRUE
    remove_genes_RIBO = TRUE
    reference_name = '42-GEX-NDBM2_filtered_feature_bc_matrix.h5'
    doublet_detection = TRUE
    subset_TRUE_FALSE = FALSE
    subset_size = 30000
    subset_size_fig = 20000
    subset_celltype_TRUE_FALSE = FALSE
    subset_celltype_ident = "predicted.celltype.l1"
    subset_celltype = 'CD8 T'
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
    
    
    
    setwd(data_directory)
    h5_files <- list.files(pattern = "\\.h5$")
    rds_files <- list.files(pattern = "\\.rds$")
    
    rds_files = rds_files[!rds_files %in% c('reference_SCT_PBMC.rds','reference_SCT_BM.rds','reference_SCT_BM_AML.rds')]
    
    
    data_folders <- list.dirs(recursive = FALSE)
    samples_list = c(data_folders, h5_files, rds_files)
    meta_data_file = list.files(pattern = "\\meta_data.csv$")
    if (!identical(meta_data_file, character(0))) {
      meta_data = read.csv(meta_data_file, stringsAsFactors = F)
      colnames(meta_data)[1] = 'ID'
      meta_data_list = colnames(meta_data)
      meta_data$ID = as.character(meta_data$ID)
    }
    
    split_list = colnames(meta_data)[2:length(colnames(meta_data))]
    
    sink(paste('Loading references','.checkpoint',sep = ''), split = TRUE)
    cat('Loading references')
    
    if (cell_type_ident_PBMC == TRUE) {
      # classifier <- readRDS("hsPBMC_20191017.RDS")
      if (is.na(reference_directory)) {
        
          cell_type_reference_PBMC = readRDS("reference_SCT_PBMC.rds")
          
      } else {
        
        cell_type_reference_PBMC = readRDS(paste(reference_directory,"reference_SCT_PBMC.rds", sep = ''))
        
      }
      # cell_type_reference@assays$RNA = cell_type_reference@assays$SCT
      # 
      # cell_type_reference <- SCTransform(cell_type_reference, verbose = FALSE)
      
    }
    if (cell_type_ident_BM == TRUE) {
      if (is.na(reference_directory)) {
        
          reference_SCT_BM = readRDS('reference_SCT_BM.rds')
          
      } else {
        
        reference_SCT_BM = readRDS(paste(reference_directory,'reference_SCT_BM.rds', sep = ''))
        
      }
      

      
    }    
    if (cell_type_ident_BM_AML == TRUE) {

      if (is.na(reference_directory)) {
          
        reference_SCT_BM_AML = readRDS('reference_SCT_BM_AML.rds')
        
      } else {
        
        reference_SCT_BM_AML = readRDS(paste(reference_directory,'reference_SCT_BM_AML.rds', sep = ''))
        
      }

      
    }
    
    if (cell_type_ident_TabulaS_PBMC == TRUE) {
      
      if (is.na(reference_directory)) {
        
        reference_SCT_TabulaS_PBMC = readRDS('reference_SCT_TabulaS_PBMC.rds')
        
      } else {
        
        reference_SCT_TabulaS_PBMC = readRDS(paste(reference_directory,'reference_SCT_TabulaS_PBMC.rds', sep = ''))
        
      }
      
      
    }
    
    if (cell_type_ident_TabulaS_BM == TRUE) {
      
      if (is.na(reference_directory)) {
        
        reference_SCT_TabulaS_BM = readRDS('reference_SCT_TabulaS_BM.rds')
        
      } else {
        
        reference_SCT_TabulaS_BM = readRDS(paste(reference_directory,'reference_SCT_TabulaS_BM.rds', sep = ''))
        
      }
      
      
    }
    
    
    nCoV.list = list()
    feat_uni = as.character()
    options(future.globals.maxSize = 100000 * 1024^2)
    i = 0
    for(sample_s in samples_list){
      
      
      
      data_wd = getwd()
      setwd(starting_wd)

      if (file.exists(paste(fig_names,'_' ,sample_s, '_sample.rds', sep = ''))) {
        
        
        
        i = i + 1
        sink(paste(starting_wd, '/Loading sample ', i,' (',sample_s ,') out of ', as.character(length(samples_list)),'.checkpoint',sep = ''), split = TRUE)
        cat('Loading sample', i, ' out of ', as.character(length(samples_list)))
        closeAllConnections()
        
        
        
        sample.tmp.seurat = readRDS(paste(fig_names,'_' ,sample_s, '_sample.rds', sep = ''))
        
        setwd(data_wd)
        
        nCoV.list[sample_s] = sample.tmp.seurat
        
      } else {
      
      setwd(data_wd)
        
      i = i + 1
      sink(paste(starting_wd, '/Processing sample ', i,' (',sample_s ,') out of ', as.character(length(samples_list)),'.checkpoint',sep = ''), split = TRUE)
      cat('Processing sample', i, ' out of ', as.character(length(samples_list)))
      closeAllConnections()
      
      
      
      
      
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
        if (ncol(sample.tmp) > subset_size) {
        sample.tmp = sample.tmp[,sample(colnames(sample.tmp),subset_size)]
        }
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
      
      saveRDS(sample.tmp.seurat,file = 'SCTransform_crash.rds')
      sample.tmp.seurat = readRDS(file = 'SCTransform_crash.rds')
      
      if (SCTransform_TRUE_FALSE == TRUE) {
        cat('Performing SCTransform on sample', i, 'out of', as.character(length(samples_list)))
        sample.tmp.seurat <- SCTransform(sample.tmp.seurat, verbose = TRUE) 
        sample.tmp.seurat <- RunPCA(sample.tmp.seurat)
        sample.tmp.seurat <- RunUMAP(sample.tmp.seurat, dims = 1:20)
        feat_uni = union(feat_uni,sample.tmp.seurat@assays$SCT@var.features)
      } else {
        cat('Performing normalization and scaling on sample', i, 'out of', as.character(length(samples_list)))
        sample.tmp.seurat <- NormalizeData(sample.tmp.seurat, verbose = TRUE)
        sample.tmp.seurat <- FindVariableFeatures(sample.tmp.seurat, selection.method = "vst", nfeatures = nfeatures_for_find_var_feat,verbose = FALSE)
        sample.tmp.seurat <- ScaleData(sample.tmp.seurat)
        sample.tmp.seurat <- RunPCA(sample.tmp.seurat)
        sample.tmp.seurat <- RunUMAP(sample.tmp.seurat, dims = 1:20)
      }
      
      
      if (cell_type_ident_PBMC == TRUE) {
        
        anchors <- FindTransferAnchors(
          reference = cell_type_reference_PBMC,
          query = sample.tmp.seurat,
          normalization.method = "SCT",
          reference.reduction = "spca",
          dims = 1:50
        )
        
        sample.tmp.seurat <- MapQuery(
          anchorset = anchors,
          query = sample.tmp.seurat,
          reference = cell_type_reference_PBMC,
          refdata = list(
            celltype.l1 = "celltype.l1",
            celltype.l2 = "celltype.l2",
            predicted_ADT = "ADT"
          ),
          reference.reduction = "spca", 
          reduction.model = "wnn.umap"
        )
        
      }
      if (cell_type_ident_BM == TRUE) {
        
        anchors <- FindTransferAnchors(
          reference = reference_SCT_BM,
          query = sample.tmp.seurat,
          normalization.method = "SCT",
          reference.reduction = "pca",
          dims = 1:50
        )
        
        sample.tmp.seurat <- MapQuery(
          anchorset = anchors,
          query = sample.tmp.seurat,
          reference = reference_SCT_BM,
          refdata = list(
            Healthy = "Prediction_Healthy",
            HCA = "Prediction_HCA",
            BM3 = "Prediction_BM3",
            CellType_BM = "ct"
          ),
          reference.reduction = "pca"
        )
        
        
      }
      if (cell_type_ident_BM_AML == TRUE) {
        
        anchors <- FindTransferAnchors(
          reference = reference_SCT_BM_AML,
          query = sample.tmp.seurat,
          normalization.method = "SCT",
          reference.reduction = "pca",
          dims = 1:50
        )
        
        sample.tmp.seurat <- MapQuery(
          anchorset = anchors,
          query = sample.tmp.seurat,
          reference = reference_SCT_BM_AML,
          refdata = list(
            Ind = "Prediction_Ind",
            CellType_BM_AML = "ct"
          ),
          reference.reduction = "pca"
        )
        
        
      }
      if (cell_type_ident_TabulaS_PBMC == TRUE) {
        
        anchors <- FindTransferAnchors(
          reference = reference_SCT_TabulaS_PBMC,
          query = sample.tmp.seurat,
          normalization.method = "SCT",
          reference.reduction = "pca",
          dims = 1:50
        )
        
        sample.tmp.seurat <- MapQuery(
          anchorset = anchors,
          query = sample.tmp.seurat,
          reference = reference_SCT_TabulaS_PBMC,
          refdata = list(
            TabulaS_Blood_cell_ontology_class = "cell_ontology_class",
            TabulaS_Blood_free_annotation = "free_annotation"
          ),
          reference.reduction = "pca"
        )
        
        
      }
      if (cell_type_ident_TabulaS_BM == TRUE) {
        
        anchors <- FindTransferAnchors(
          reference = reference_SCT_TabulaS_BM,
          query = sample.tmp.seurat,
          normalization.method = "SCT",
          reference.reduction = "pca",
          dims = 1:50
        )
        
        sample.tmp.seurat <- MapQuery(
          anchorset = anchors,
          query = sample.tmp.seurat,
          reference = reference_SCT_TabulaS_BM,
          refdata = list(
            TabulaS_BM_cell_ontology_class = "cell_ontology_class",
            TabulaS_BM_free_annotation = "free_annotation"
          ),
          reference.reduction = "pca"
        )
        
        
      }
      
      sample.tmp.seurat@meta.data$BTM_Others <- plyr::mapvalues(
        x = sample.tmp.seurat@meta.data$predicted.celltype.l1, 
        from = c("B",       "CD4 T",      "CD8 T",       "other T", "Mono",      "other",   "DC",      "NK"),
        to = c(  "B cells", "T cells",    "T cells",     "T cells", "Monocytes", "Others",  "Others",  "Others")
      )
      
      
      if (subset_celltype_TRUE_FALSE == TRUE){
        
        Idents(sample.tmp.seurat) = subset_celltype_ident
        sample.tmp.seurat = subset(sample.tmp.seurat, idents = subset_celltype)
        
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
      data_wd = getwd()
      setwd(starting_wd)
      saveRDS(sample.tmp.seurat, paste(fig_names,'_' ,sample_s, '_sample.rds', sep = ''))
      setwd(data_wd)
      
      
      
      
      nCoV.list[sample_s] = sample.tmp.seurat
      
      }
      # rm(sample.tmp.seurat)
      # rm(sample.tmp)
      
      
      #### Add report stuff here  ##################################################################
      # sample.tmp.seurat = nCoV.list[[1]]
      if (report_TRUE_FALSE == TRUE) {
      
      
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
      
      setwd(starting_wd)
      
      if (file.exists(paste(fig_names, '_', sample_s, "_AllClusterMarkers.csv", sep = ''))) {
        
        sample.tmp.seurat.markers = read.csv2(paste(fig_names, '_', sample_s, "_AllClusterMarkers.csv", sep = ''),sep = ',', row.names = 1)
        
      } else {
        
        sample.tmp.seurat.markers <- FindAllMarkers(sample.tmp.seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
        
        write.csv(sample.tmp.seurat.markers, file = paste(fig_names, '_', sample_s, "_AllClusterMarkers.csv", sep = ''))
      
      }
      
      setwd(data_directory)
      
      # sample.tmp.seurat.markers <- FindAllMarkers(sample.tmp.seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
      
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
      
      }
      
      ###################################################
      
      
      
      
      
      
      # rm(aml2_)
      gc()
      
    } ### Sample list loop closed ###
    
    
    gc()
    
    
    sink(paste('Saving processed individual samples','.checkpoint',sep = ''), split = TRUE)
    cat('Saving processed individual samples')
    
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
  rm(reference_SCT_BM)
  rm(reference_SCT_BM_AML)
  rm(cell_type_reference_PBMC)
  rm(sample.tmp)
  rm(sample.tmp.seurat)
  rm(reference_SCT_TabulaS_PBMC)
  rm(reference_SCT_TabulaS_BM)
  
  
  gc()
  
  setwd(starting_wd)
  
  # rm(sample.tmp.foralra,sample.tmp.seurat,sample.tmp,fuckme)
  # gc()
  rm(i)
  sink(paste('Starting Integration','.checkpoint',sep = ''), split = TRUE)
  cat('Starting Integration')
  
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
  
  # if (!is.na(factor_1) == TRUE) {
  # try(nCoV.integrated@meta.data[,split_list[1]] = factor(nCoV.integrated@meta.data[,split_list[1]], levels = factor_1))
  # }
  # if (!is.na(factor_2) == TRUE) {
  # try(nCoV.integrated@meta.data[,split_list[2]] = factor(nCoV.integrated@meta.data[,split_list[2]], levels = factor_2))
  # }
  
  
  # nCoV.integrated@meta.data
  
  # nCoV.integrated@meta.data$BTM_Others <- plyr::mapvalues(
  #   x = nCoV.integrated@meta.data$predicted.celltype.l1, 
  #   from = c("B",       "CD4 T",      "CD8 T",       "other T", "Mono",      "other",   "DC",      "NK"),
  #   to = c(  "B cells", "T cells",    "T cells",     "T cells", "Monocytes", "Others",  "Others",  "Others")
  # )
  
  sink(paste('Starting standard analysis on integrated raw RNA data','.checkpoint',sep = ''), split = TRUE)
  cat('Starting standard analysis on integrated raw RNA data')

  
  
  all.genes.COVID <- rownames(nCoV.integrated)
  
  feat = intersect(rownames(nCoV.integrated),feat_uni)
  
  
  if (SCTransform_TRUE_FALSE == TRUE) {
    DefaultAssay(nCoV.integrated) <- "SCT"
    nCoV.integrated <- RunPCA(nCoV.integrated, verbose = TRUE, features = all.genes.COVID, assay = "SCT")
    
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
    
    sink(paste('Saving integrated dataset','.checkpoint',sep = ''), split = TRUE)
    cat('Saving integrated dataset')
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

# Idents(nCoV.integrated) = 'BTM_Others'
# nCoV.integrated.Tcells = subset(nCoV.integrated, idents = c("T cells"))
# 
# saveRDS(nCoV.integrated.Tcells, paste(fig_names, '_nCoV.integrated.Tcells.rds', sep = ''))
# 
# nCoV.integrated = readRDS(paste(fig_names, '_nCoV.integrated.Tcells.rds', sep = ''))





sink(paste('Making plots','.checkpoint',sep = ''), split = TRUE)
cat('Making plots')

################################################################## QC PLOT ##################################################################
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
################################################################## BASIC PLOTS ##################################################################

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






################################################################## MARKER GENES PLOTS ##################################################################

dpi = 200
png(file=paste(fig_names, '_',paste('Marker_Genes','.png', sep = ''), sep = ''), width = dpi*32, height = dpi*16, units = "px",res = dpi,type='cairo')
print(FeaturePlot(nCoV.integrated, reduction = "umap", features = Marker_genes))
dev.off()


################################################################## HEATMAP CLUSTERS PLOTS ##################################################################



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

quit()
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

nCoV.integrated = readRDS(paste(fig_names, '_nCoV.integrated.rds', sep = ''))

nCoV.integrated@meta.data$Day = factor(nCoV.integrated@meta.data$Day, levels = c("Healthy_BM","d-21 to -7", "d-1", "d7", "d14", "d28"))

nCoV.integrated@meta.data$Day_New <- plyr::mapvalues(
  x = nCoV.integrated@meta.data$Day, 
  from = c("Healthy_BM","d-21 to -7", "d-1", "d7", "d14", "d28"),
  to = c("Healthy_BM","Pre Chemo", "Pre Infusion", "Post Infusion", "Post Infusion", "Post Infusion")
)


nCoV.integrated@meta.data$Day_New = factor(nCoV.integrated@meta.data$Day_New, levels = c("Healthy_BM","Pre Chemo", "Pre Infusion", "Post Infusion"))

unique(nCoV.integrated@meta.data$Day_New)
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
###### GO ANALYSIS ###############################################################################################################
# BiocManager::install("topGO")
topGOterms = function( fg.genes = NULL,
                       bg.genes = NULL,
                       organism = "Human", 
                       ontology.use = "BP",
                       stats.use = "fisher",
                       algorithm.use = "weight01",
                       topnodes.print=20,
                       num.char=100){
  
  if (is.null(fg.genes) | is.null(bg.genes)){
    stop("Error : Both gene lists are empty")
  }
  
  require(topGO)
  if (organism == "Mouse"){
    mapping.use = "org.Mm.eg.db"
    library(org.Mm.eg.db)
  } else if (organism == "Human"){
    mapping.use = "org.Hs.eg.db"
    library(org.Hs.eg.db)
  } else {
    stop("Error : Organisms other than mouse not supported currently")
  }
  
  n = length(bg.genes)
  geneList = integer(n)
  names(geneList) = bg.genes
  geneList[intersect(names(geneList), fg.genes)]=1
  print(paste0("Total ", length(geneList), " genes. ", sum(geneList), " genes in the foreground"))
  geneList = factor(geneList)
  
  if (ontology.use %in% c("BP", "CC", "MF")){
    print(paste0("Using Ontology : ", ontology.use))
  } else {
    stop("Error: Ontology not available. Should be one of BP, CC or MF")
  }
  # Make GO object
  GOdata <- new("topGOdata",
                description = "GOanalysis",
                ontology = ontology.use,
                allGenes = geneList,
                annot = annFUN.org,
                mapping = mapping.use,
                ID = "SYMBOL",
                nodeSize = 10)
  print(paste0("Using the ", stats.use, " statistic with the ", algorithm.use, " algorithm"))
  res.result <- runTest(GOdata, statistic = stats.use, algorithm = algorithm.use)
  to.return = list()
  to.return$GOdata = GOdata
  to.return$res.table <- GenTable(GOdata, pval = res.result, topNodes = topnodes.print, numChar = num.char)
  return(to.return)
}
############################# GO FIGURE FUNCTION#######################################################################################

GO_plot = function(GO_results_df, file_name = '', title = '')
  
{
  graphics.off()
  GO_results_df$pval = as.numeric(GO_results_df$pval)
  GO_results_df = filter(GO_results_df, pval < 0.05)
  data1 <- GO_results_df
  
  colnames(GO_results_df)[1] = 'Cluster ID'
  
  S1 <- ggplot(data1, aes(x = id, y = Term , color = -log(pval), size = log(Significant))) + 
    geom_point() + xlab("Dataset") + #scale_x_continuous(breaks = 1:12) +
    ylab(NULL) + 
    labs(colour = "-Log10(P.adjust)", size = "Log10(Number of Genes)") +
    theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          axis.text=element_text(size=10,colour ="black"),
          axis.title=element_text(size=14,face="bold",colour ="black"),
          legend.text=element_text(size=12),
          legend.title = element_text(colour="black", size=14, face="bold"),
          axis.line = element_line(colour = "black", size =1),
          panel.grid.major = element_line(colour = "lightgray", size = 0.5), 
          # panel.grid.minor = element_line(colour = "lightgray", size = 0.5),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, size = 12)) + ggtitle(title)
  
  # S1 = S1+scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab", limit = c(min(data1$pval), max(data1$pval)), n.breaks = 5)
  S1+scale_size(range = c(0, length(unique(data1$id))))
  S1
  
  if(str_sub(file_name, -4, -1) != '.png') {
  
    file_name = paste(file_name, '.png', sep = '')
    
  }
  
  dpi = 200
  png(file=file_name, width = dpi*12, height = dpi*12, units = "px",res = dpi,type='cairo')
  print(S1)
  graphics.off()
  # try({
  #   
  #   dev.off()
  #   dev.off()
  #   dev.off()}, silent = T)
  
  
  
}


GO_plot(GO_results_df, file_name = 'whatever', title = 'Whatever Test')

######################################### VOLCANO_PLOT FUNCTION ######################################################################################

markers_list_of_lists = readRDS('AML_AUG2021_CAR123_Patient_predicted.celltype.l2_Day_Healthy_TRUE.rds')


DGE_table = markers_list_of_lists[['Patient 2']][['HSPC']][["d14_vs_ctrl_d-1"]][,2:7]

volcano_plot(DGE_table)


volcano_plot = function(DGE_table, gene_list = NA, Control_name = "Control", Condition_name = 'Condition', Comparison_name = '', FileName = "Volcano"){
  
  results = DGE_table
  
  results$p_val_adj = as.numeric(results$p_val_adj)
  results$avg_log2FC = as.numeric(results$avg_log2FC)
  
  if (is.na(gene_list)) {
    
    results = mutate(results, Genes=ifelse(results$p_val_adj<0.05, "Significant", "Non-Significant"))
    res = results[order(results$p_val_adj, decreasing = FALSE), ]
    list = res$Gene[1:5]
    
  } else{
    
    results = mutate(results, Genes=ifelse(results$Gene %in% gene_list, "List", "Others"))
    list = gene_list
    
  }
  
  
  
  results4 = results[order(results$Genes, decreasing = FALSE), ]
  
  results4$p_val_adj = as.numeric(results4$p_val_adj)
  results4$avg_log2FC = as.numeric(results4$avg_log2FC)
  
  pvals = results4$p_val_adj
  pvals[pvals == 0] = 2
  
  pvals[pvals == 2] = runif(sum(pvals == 2), 0, min(pvals))
  
  
  results4$p_val_adj = pvals
  
  
  p = ggplot(results4, aes(avg_log2FC, -log10(as.numeric(p_val_adj))),group=Genes) +
    geom_jitter(aes(col=Genes,size=Genes,shape = Genes)) + ylim(-10, round(-log10(min(results4$p_val_adj))/10)*10) + xlim(-(round(max(abs(DGE_table$avg_log2FC)))+2),(round(max(abs(DGE_table$avg_log2FC)))+2)) +
    scale_size_manual(values=c(0.5,2)) + scale_shape_manual(values=c(19,19))+
    scale_color_manual(values=c("gray30","red"))+ xlab('Log2 Fold Change') + ylab('Log10 Adjusted Pvalue') + ggtitle(paste(Control_name, "              ",Comparison_name,"                          ", Condition_name)) +theme(plot.title = element_text(hjust = 0.5))
  
  theme_set(theme_gray(base_size = 18)) 
  # dev.off()
  
  p = p+geom_vline(xintercept = 2, linetype="dotted", 
                   color = "black", size=1)+
    geom_vline(xintercept = -2, linetype="dotted", 
               color = "black", size=1)+
    geom_hline(yintercept=20, linetype="dotted", color = "black", size=1)+
    # geom_text_repel(color = 'black', size = 5,data=filter(results, Gene %in% results2$V1 & results$Adjusted.P.Value<1e-50),
    geom_text_repel(color = 'black', size = 5,data=filter(results4, Gene %in% list ),#& results$adj.P.Val<1e-80),
                    direction ='both',xlim=c(-10,1000), ylim=c(20,100), aes(label=list),min.segment.length = 0,
                    point.padding = 0.25, box.padding = 2.00,force=5) #,segment.color = "grey50"
  
  
  print(p)
  
  dpi = 100
  png(file=paste(FileName,'.png', sep = ''), width = dpi*16, height = dpi*9, units = "px",res = dpi,type='cairo')
  print(p)
  dev.off()
  
}

###############################################################################################################################

#https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html
### CELLCHAT ###


# nCoV.integrated = readRDS('AML_AUG2021_IMPUTE_nCoV.integrated_CAR123.rds')
# 
# 
# 
# nCoV.integrated = readRDS('AML_AUG2021_IMPUTE_nCoV.integrated.rds')
# 
# 
# 
# rm(nCoV.integrated.Patient)
# rm(nCoV.integrated.Patients)
# gc()


# nCoV.integrated@meta.data$BTM_Others <- plyr::mapvalues(
#   x = nCoV.integrated@meta.data$predicted.celltype.l1, 
#   from = c("B",       "CD4 T",      "CD8 T",       "other T", "Mono",      "other",   "DC",      "NK"),
#   to = c(  "B cells", "T cells",    "T cells",     "T cells", "Monocytes", "Others",  "Others",  "Others")
# )

# nCoV.integrated@meta.data$BTM_Others_CAR = nCoV.integrated@meta.data$BTM_Others
# 
# poscells <- WhichCells(nCoV.integrated, expression = scvfCAR123 > 0)
# 
# nCoV.integrated@meta.data[poscells, 'BTM_Others_CAR'] = 'scvfCAR123 T cells'
# 
# 
# nCoV.integrated@meta.data$BTM_Others_CD38 = nCoV.integrated@meta.data$BTM_Others
# 
# 
# 
# poscells <- WhichCells(nCoV.integrated, expression = CD38 > 0)
# 
# nCoV.integrated@meta.data[poscells, 'BTM_Others_CD38'] = 'CD38+ cells'
# 
# 
# Idents(nCoV.integrated) = 'BTM_Others_CD38'
# 
# 
# DEG_CD38 = FindMarkers(nCoV.integrated, ident.1 = 'CD38+ cells')
# 
# 
# table(nCoV.integrated@meta.data$BTM_Others_CAR)
# 
# poscells <- WhichCells(nCoV.integrated, expression = scvfCAR123 > 0)
# 
# # so$FOXP3_exp<- ifelse(colnames(so) %in% poscells, "Pos", "Neg")
# 
# 
# # unique(nCoV.integrated@meta.data$Patient)
# 
# unique(nCoV.integrated@meta.data$Day)
# 
# nCoV.integrated@meta.data$Day2 <- plyr::mapvalues(
#   x = nCoV.integrated@meta.data$Day, 
#   from = c("d-21 to-7",   "d-1",        "d14",        "d-2 to -7",  "d28",        "d7",         "d-21 to -7"),
#   to =   c("d-21 to -7",  "d-1",        "d14",        "d-21 to -7",  "d28",        "d7",         "d-21 to -7")
# )

# unique(nCoV.integrated@meta.data$Day2)
# nCoV.integrated@meta.data$Day = nCoV.integrated@meta.data$Day2

nCoV.integrated@meta.data$Day = factor(nCoV.integrated@meta.data$Day, levels = c("Healthy_BM","d-21 to -7", "d-1", "d7", "d14", "d28"))

# unique(nCoV.integrated@meta.data$Day)

# FeaturePlot(object = nCoV.integrated, features = 'CD38', split.by = 'Patient', label = TRUE ,pt.size = 1.5, raster=FALSE, repel = TRUE, order = TRUE)




# aml2_$Condition =  factor(aml2_$Condition, levels = c("Untreated CD45n", "UTD CD45n", "CARMA CD45n", "Untreated CD45p", "UTD CD45p", "CARMA CD45p"))
###############
### Remove ribosomal genes from analysis #####
allGenes <- rownames(nCoV.integrated)
RemoveGenes <- allGenes[grepl("^[Rr][Pp][Ss]", allGenes) | grepl("^[Rr][Pp][Ll]", allGenes)]
nCoV.integrated <- nCoV.integrated[allGenes[!allGenes %in% RemoveGenes],]
allGenes2 <- rownames(nCoV.integrated)

Idents(nCoV.integrated) = "predicted.celltype.l2"

table(nCoV.integrated$predicted.celltype.l2)

nCoV.integrated = subset(nCoV.integrated, idents = 'Doublet', invert = TRUE)
table(nCoV.integrated$predicted.celltype.l2)
levels(nCoV.integrated)

rm(nCoV.integrated.Patients)
rm(nCoV.integrated.Patient)
gc()
Idents(nCoV.integrated) = "Patient"

nCoV.integrated.Patient.5 = subset(nCoV.integrated, idents = c("Patient 5"), invert = FALSE)

Idents(nCoV.integrated.Patient.5) = "BTM_Others_CAR"

DEG_CAR123 = FindMarkers(nCoV.integrated.Patient.5, ident.1 = 'scvfCAR123 T cells', ident.2 = 'T cells')


#####################################
markers = DEG_CAR123
markers1 = cbind(rownames(markers), markers)
colnames(markers1)[1] = 'Gene'
rownames(markers1) = c(1:nrow(markers1))
markers1 = subset(markers1, p_val_adj < 0.05)
markers1 = cbind(c(1:nrow(markers1)), markers1)
markers1$`c(1:nrow(markers1))` = sub('\\d+','_', markers1$`c(1:nrow(markers1))`)
colnames(markers1)[1] = 'Original out p_adj < 0.05 >>'


markers2 = markers1
markers2$pct.1 = NULL
markers2$pct.2 = NULL
markers2$p_val = NULL
markers2 = markers2[order(-markers2$avg_log2FC), ]

colnames(markers2)[1] = 'Sorted logFC asc >>'


markers3 = markers2
markers3 = markers3[order(markers3$avg_log2FC), ]
colnames(markers3)[1] = 'Sorted logFC desc >>'

markers = cbind(markers1,markers2,markers3)

write.csv(markers, file = 'DEG_CAR123_T_cells_vs_T_cells_Ctrl.csv')




#####################################

nCoV.integrated.Patients = list()
gc()

Idents(nCoV.integrated) = "Patient"
for(Patient in unique(nCoV.integrated@meta.data$Patient)) {
  
  nCoV.integrated.Patient = subset(nCoV.integrated, idents = c(Patient,"Healthy_BM"), invert = FALSE)
  
  # saveRDS(nCoV.integrated.Patient, paste('nCoV.integrated.',Patient,'.rds',sep = ''))
  
  nCoV.integrated.Patients[Patient] = nCoV.integrated.Patient
  
}

nCoV.integrated.Patients['Healthy_BM'] = NULL
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
nCoV.integrated@meta.data$pANN_0.25_0.09_5 = NULL
meta_data_list_plots = colnames(nCoV.integrated@meta.data)

count_unique <- rapply(nCoV.integrated@meta.data, function(x) length(unique(x))) < 40
# count_unique2 <- rapply(nCoV.integrated@meta.data, function(x) length(unique(x))) > 1

# count_unique*count_unique2
# 
# meta_data_list_plots = meta_data_list_plots[count_unique*count_unique2]


meta_data_list_plots = meta_data_list_plots[count_unique]


# meta_data_list_plots = meta_data_list_plots[-8]

unique(nCoV.integrated@meta.data$Day)
# nCoV.integrated@meta.data$Day = factor(nCoV.integrated@meta.data$Day, levels = c("Healthy_BM","d-21 to -7", "d-1", "d7", "d14", "d28"))
nCoV.integrated@meta.data$Day = factor(nCoV.integrated@meta.data$Day, levels = c("Normal BM","d-21 to -7", "d-1", "d14", "28d", "Infusion product"))


split_list = c("Day")

Idents(nCoV.integrated) = "Patient"

for(Patient in unique(nCoV.integrated@meta.data$Patient)) {

  nCoV.integrated.Patient = subset(nCoV.integrated, idents = c(Patient, "Normal BM"), invert = FALSE)

for(split_ in split_list) {
  for(meta in meta_data_list_plots) {
    Idents(nCoV.integrated.Patient) <- meta
    dpi = 200
    png(file=paste(Patient,'_Split_',split_,'_', fig_names, '_',paste(meta,'.png', sep = ''), sep = ''), width = dpi*24, height = dpi*8, units = "px",res = dpi,type='cairo')
    print(DimPlot(object = nCoV.integrated.Patient, group.by = "ident",split.by = split_, label = TRUE ,pt.size = 1.5, raster=FALSE, repel = TRUE))
    dev.off()
  }
}

}

unique(nCoV.integrated@meta.data$pANN_0.25_0.09_5)

for(Patient in unique(nCoV.integrated@meta.data$Patient)) {
  
  nCoV.integrated.Patient = subset(nCoV.integrated, idents = c(Patient, "Normal BM"), invert = FALSE)
  
  for(split_ in split_list) {
      Idents(nCoV.integrated.Patient) <- 'predicted.TabulaS_BM_free_annotation'
      dpi = 200
      png(file=paste(Patient,'scvfCAR123_Split_',split_,'_', fig_names, '.png', sep = ''), width = dpi*24, height = dpi*8, units = "px",res = dpi,type='cairo')
      print(FeaturePlot(object = nCoV.integrated.Patient, features = "scvfCAR123",split.by = split_, label = TRUE, pt.size = 1.5, raster=FALSE, repel = TRUE, order = TRUE)
)
      dev.off()
    
  }
  
}
unique(nCoV.integrated$Patient)

nCoV.integrated.Patient = subset(nCoV.integrated, idents = c('Patient 1'), invert = FALSE)


Idents(nCoV.integrated.Patient) <- 'predicted.TabulaS_BM_free_annotation'
FeaturePlot(object = nCoV.integrated.Patient, features = "scvfCAR123",split.by = split_, label = TRUE ,pt.size = 1.5, raster=FALSE, repel = TRUE, order = TRUE)
# meta_data_list_plots = c('predicted.ct_BM','predicted.ct_BM_AML')
#########################################################################################################################################################

################################################################## FIND MARKERS AND GO BETWEEN TIMEPOINTS ##################################################################
Idents(nCoV.integrated) = 'predicted.celltype.l2'

unique(Idents(nCoV.integrated))

nCoV.integrated = subset(nCoV.integrated, idents = c('HSPC'))
# Idents(nCoV.integrated) = CellType_level
# 
# nCoV.integrated@meta.data$BTM_Others <- plyr::mapvalues(
#   x = nCoV.integrated@meta.data$predicted.celltype.l1, 
#   from = c("B",       "CD4 T",      "CD8 T",       "other T", "Mono",      "other",   "DC",      "NK"),
#   to = c(  "B cells", "T cells",    "T cells",     "T cells", "Monocytes", "Others",  "Others",  "Others")
# )
# 
# unique(nCoV.integrated@meta.data$Day)
# 
# setwd("GO_test")
# 
# nCoV.integrated = subset(nCoV.integrated, idents = c('Mono'))


analysis_folder_name = "GO_Analysis_CellType_L2_Updated_Volcano_2"



dir.create(analysis_folder_name, showWarnings = TRUE, recursive = FALSE)


setwd(analysis_folder_name)


Healthy = TRUE
Healthy_Sample = "Healthy_BM"
First_Split = "Patient"
Timepoints_meta = "Day_New"
CellType_level = "predicted.celltype.l2"



if (file.exists(paste(fig_names, '_',First_Split,'_',CellType_level,'_',Timepoints_meta,'_Healthy_',Healthy,'.rds', sep = ''))) {

  markers_list_of_lists = readRDS(paste(fig_names, '_',First_Split,'_',CellType_level,'_',Timepoints_meta,'_Healthy_',Healthy,'.rds', sep = '')) 
  
} else {

markers_list_of_lists = list()
markers_list_of_lists_raw = list()
in_list = unique(nCoV.integrated@meta.data[, "Patient"])

if (Healthy == TRUE){
  
  in_list = in_list[!(in_list %in% Healthy_Sample)]
  
}

for(First_Split_ in in_list) {
  
  Idents(nCoV.integrated) = First_Split
  
  
  if (Healthy == TRUE) {
    nCoV.integrated.Dataset = subset(nCoV.integrated, idents = c(First_Split_, Healthy_Sample), invert = FALSE)
  } else {
    nCoV.integrated.Dataset = subset(nCoV.integrated, idents = c(First_Split_), invert = FALSE)
    
  }
  
  ##length(as.character(unique(nCoV.integrated.Patient@meta.data$Day)))
  
  Timepoints = as.character(unique(nCoV.integrated.Dataset@meta.data[, Timepoints_meta]))
  
  Timepoints_No_Control = Timepoints[!Timepoints %in% Healthy_Sample]
  
  Idents(nCoV.integrated.Dataset) = CellType_level
  
  
  for(celltype in as.character(unique(Idents(nCoV.integrated.Dataset)))) { # loop through all celltypes
    
    
    nCoV.integrated.Dataset.celltype = subset(nCoV.integrated.Dataset, idents = celltype, invert = FALSE) 
    
    Idents(nCoV.integrated.Dataset.celltype) = Timepoints_meta
    # 
    # Idents(nCoV.integrated.Dataset) = Timepoints_meta
    
    
    
    
    for(i in 1:(length(Timepoints_No_Control)-1)) { # loop through all timepoints
      
      
      
      if (length(unique(levels(Idents(nCoV.integrated.Dataset.celltype)))) == length(Timepoints)) { # check for cells in each time point
      
      # WhichCells(nCoV.integrated.Dataset.celltype, idents = Timepoints[i+1])
      # WhichCells(nCoV.integrated.Dataset.celltype, idents = Timepoints[i])
      
      
      
      if(length(WhichCells(nCoV.integrated.Dataset.celltype, idents = Timepoints[i+1]))>=10 && length(WhichCells(nCoV.integrated.Dataset.celltype, idents = Timepoints[i])) >= 10)
        # check for cells for more than 10 cells in each timepoint
      {
        
      markers = FindMarkers(nCoV.integrated.Dataset.celltype, ident.1 = Timepoints[i+1], ident.2 = Timepoints[i],logfc.threshold = 0)
      markers_list_of_lists_raw[[First_Split_]][[celltype]][[paste(Timepoints[i+1],'_vs_ctrl_',Timepoints[i], sep = '')]] = markers
      
      try({
      markers1 = cbind(rownames(markers), markers)
      colnames(markers1)[1] = 'Gene'
      rownames(markers1) = c(1:nrow(markers1))
      markers1 = subset(markers1, p_val_adj < 0.05)
      markers1 = cbind(c(1:nrow(markers1)), markers1)
      markers1$`c(1:nrow(markers1))` = sub('\\d+','_', markers1$`c(1:nrow(markers1))`)
      colnames(markers1)[1] = 'Original out p_adj < 0.05 >>'
      
      
      markers2 = markers1
      markers2$pct.1 = NULL
      markers2$pct.2 = NULL
      markers2$p_val = NULL
      markers2 = markers2[order(-markers2$avg_log2FC), ]
      
      colnames(markers2)[1] = 'Sorted logFC asc >>'
      
      
      markers3 = markers2
      markers3 = markers3[order(markers3$avg_log2FC), ]
      colnames(markers3)[1] = 'Sorted logFC desc >>'
      
      markers = cbind(markers1,markers2,markers3)
      }, silent = FALSE)
      markers_list_of_lists[[First_Split_]][[celltype]][[paste(Timepoints[i+1],'_vs_ctrl_',Timepoints[i], sep = '')]] = markers
      
      write.csv(markers, file = paste(fig_names, '_',First_Split_,'_',CellType_level,'_',celltype,'_',Timepoints[i+1],'_vs_ctrl_',Timepoints[i], '_Markers.csv', sep = ''))
      
      } else {
        
        markers = data.frame()
        markers_list_of_lists[[First_Split_]][[celltype]][[paste(Timepoints[i+1],'_vs_ctrl_',Timepoints[i], sep = '')]] = markers
        
      }

        
        
      }
      
      if (Healthy == TRUE) {
      
      for(Timepoint in Timepoints_No_Control) { # loop through all timepoints
        

        if (length(unique(levels(Idents(nCoV.integrated.Dataset.celltype)))) == length(Timepoints)) { # check for cells in each time point
          
          # WhichCells(nCoV.integrated.Dataset.celltype, idents = Timepoints[i+1])
          # WhichCells(nCoV.integrated.Dataset.celltype, idents = Timepoints[i])
          
          
          
          if(length(WhichCells(nCoV.integrated.Dataset.celltype, idents = Timepoint))>=10 && length(WhichCells(nCoV.integrated.Dataset.celltype, idents = Healthy_Sample)) >= 10)
            # check for cells for more than 10 cells in each timepoint
          {
            
            markers = FindMarkers(nCoV.integrated.Dataset.celltype, ident.1 = Timepoint, ident.2 = Healthy_Sample,logfc.threshold = 0)
            markers_list_of_lists_raw[[First_Split_]][[celltype]][[paste('z_',Timepoint,'_vs_ctrl_',Healthy_Sample, sep = '')]] = markers
            
            try({
              markers1 = cbind(rownames(markers), markers)
              colnames(markers1)[1] = 'Gene'
              rownames(markers1) = c(1:nrow(markers1))
              markers1 = subset(markers1, p_val_adj < 0.05)
              markers1 = cbind(c(1:nrow(markers1)), markers1)
              markers1$`c(1:nrow(markers1))` = sub('\\d+','_', markers1$`c(1:nrow(markers1))`)
              colnames(markers1)[1] = 'Original out p_adj < 0.05 >>'
              
              
              markers2 = markers1
              markers2$pct.1 = NULL
              markers2$pct.2 = NULL
              markers2$p_val = NULL
              markers2 = markers2[order(-markers2$avg_log2FC), ]
              
              colnames(markers2)[1] = 'Sorted logFC asc >>'
              
              
              markers3 = markers2
              markers3 = markers3[order(markers3$avg_log2FC), ]
              colnames(markers3)[1] = 'Sorted logFC desc >>'
              
              markers = cbind(markers1,markers2,markers3)
            }, silent = FALSE)
            markers_list_of_lists[[First_Split_]][[celltype]][[paste('z_',Timepoint,'_vs_ctrl_',Healthy_Sample, sep = '')]] = markers
            
            write.csv(markers, file = paste(fig_names, '_',First_Split_,'_',CellType_level,'_',celltype,'_',Timepoint,'_vs_healthy_ctrl_',Healthy_Sample, '_Markers.csv', sep = ''))
            
          } else {
            
            markers = data.frame()
            markers_list_of_lists[[First_Split_]][[celltype]][[paste(Timepoint,'_vs_ctrl_',Healthy_Sample, sep = '')]] = markers
            
          }
          
          
          
        }
      }  # /// loop through all timepoints ///
      }
      
      
      
      
    }
  }
  
  
  sink(paste('Sample_',First_Split_,'_Markers_Complete','.checkpoint',sep = ''), split = TRUE)
  
  
}

saveRDS(markers_list_of_lists, file = paste(fig_names, '_',First_Split,'_',CellType_level,'_',Timepoints_meta,'_Healthy_',Healthy,'.rds', sep = ''))
saveRDS(markers_list_of_lists_raw, file = paste(fig_names, '_',First_Split,'_',CellType_level,'_',Timepoints_meta,'_Healthy_',Healthy,'_RAW.rds', sep = ''))

}


markers_list_of_lists = readRDS(paste(fig_names, '_',First_Split,'_',CellType_level,'_',Timepoints_meta,'_Healthy_',Healthy,'.rds', sep = ''))
markers_list_of_lists_raw = readRDS(paste(fig_names, '_',First_Split,'_',CellType_level,'_',Timepoints_meta,'_Healthy_',Healthy,'_RAW.rds', sep = ''))

############### Plot volcanos ##################
starting_wd_volcano = getwd()

for (names_1 in names(markers_list_of_lists_raw)) {
  setwd(starting_wd_volcano)
  
  if (!file.exists(names_1)) {
    
    dir.create(names_1, showWarnings = TRUE, recursive = FALSE)
    
  }
  
  setwd(names_1)
  
  starting_wd_names_1 = getwd()
  
  for (names_2 in names(markers_list_of_lists_raw[[names_1]])) {
    
    if (!file.exists(names_2)) {
      
      dir.create(names_2, showWarnings = TRUE, recursive = FALSE)
      
    }
    
    setwd(names_2)
    
    for (names_3 in names(markers_list_of_lists_raw[[names_1]][[names_2]])){
      
      DGE_df = markers_list_of_lists_raw[[names_1]][[names_2]][[names_3]]
    
    DGE_df$Gene = rownames(DGE_df)
    
    volcano_plot(DGE_df, FileName = names_3, Comparison_name = names_3)
    }
    setwd(starting_wd_names_1)
    
  }
  
  
  
}

setwd(starting_wd_volcano)

###############
# DGE_df = markers_list_of_lists[['Patient 2']][['HSPC']][["d14_vs_ctrl_d-1"]][,2:7]






# Go_list_of_lists_UP_test = Go_list_of_lists_UP
# 
# aaa = Go_list_of_lists_UP_test$`Patient 2`$`CD4 T`
# 
# 
# aaa = aaa[order(names(aaa))]





Go_list_of_lists_UP = list()
Go_list_of_lists_DOWN = list()



for(First_Split_ in names(markers_list_of_lists)) {
  
  if (file.exists(paste(fig_names, '_GO_list_UP_',First_Split_,'_',CellType_level,'_',Timepoints_meta,'_Healthy_',Healthy,'.rds', sep = '')) && file.exists(paste(fig_names, '_GO_list_DOWN_',First_Split_,'_',CellType_level,'_',Timepoints_meta,'_Healthy_',Healthy,'.rds', sep = ''))) {
    
    Go_list_of_lists_UP[[First_Split_]] = readRDS(paste(fig_names, '_GO_list_UP_',First_Split_,'_',CellType_level,'_',Timepoints_meta,'_Healthy_',Healthy,'.rds', sep = '')) 
    Go_list_of_lists_DOWN[[First_Split_]] = readRDS(paste(fig_names, '_GO_list_DOWN_',First_Split_,'_',CellType_level,'_',Timepoints_meta,'_Healthy_',Healthy,'.rds', sep = '')) 
    
    
    
    
  } else {
  
  
  for(celltype in names(markers_list_of_lists[[First_Split_]])) {
    
      for(timepoint in names(markers_list_of_lists[[First_Split_]][[celltype]])) {
        
        if (length(colnames(markers_list_of_lists[[First_Split_]][[celltype]][[timepoint]])) > 0) {
        
        markers = markers_list_of_lists[[First_Split_]][[celltype]][[timepoint]][,1:4]
        markers_UP = dplyr::filter(markers, avg_log2FC > 0.25)$Gene
        markers_DOWN = dplyr::filter(markers, avg_log2FC < -0.25)$Gene
        
        if (length(markers_UP) != 0) {
          
          GOterms_UP = topGOterms(fg.genes = markers_UP, bg.genes = rownames(nCoV.integrated), organism = "Human")
          
          
        }
        
        if (length(markers_DOWN) != 0) {
          
          GOterms_DOWN = topGOterms(fg.genes = markers_DOWN, bg.genes = rownames(nCoV.integrated), organism = "Human")
          
          
        }
        
        
        
        Go_list_of_lists_UP[[First_Split_]][[celltype]][[timepoint]] = GOterms_UP$res.table  %>% top_n(n = 15, wt = pval)
        
        Go_list_of_lists_DOWN[[First_Split_]][[celltype]][[timepoint]] = GOterms_DOWN$res.table  %>% top_n(n = 15, wt = pval)
        
        
        } else {
          
          Go_list_of_lists_UP[[First_Split_]][[celltype]][[timepoint]] = data.frame()
          
          Go_list_of_lists_DOWN[[First_Split_]][[celltype]][[timepoint]] = data.frame()
          
        }
        
      }
    
    
    
  }
  
  sink(paste('Sample_',First_Split_,'_GO_enrichment_Complete','.checkpoint',sep = ''), split = TRUE)
  
  checkpoint = Go_list_of_lists_UP[[First_Split_]]
  
  saveRDS(checkpoint, file = paste(fig_names, '_GO_list_UP_',First_Split_,'_',CellType_level,'_',Timepoints_meta,'_Healthy_',Healthy,'.rds', sep = ''))
  
  checkpoint2 = Go_list_of_lists_DOWN[[First_Split_]]
  
  saveRDS(checkpoint2, file = paste(fig_names, '_GO_list_DOWN_',First_Split_,'_',CellType_level,'_',Timepoints_meta,'_Healthy_',Healthy,'.rds', sep = ''))
  }
  
}


saveRDS(Go_list_of_lists_UP, file = paste(fig_names, '_GO_list_UP_',First_Split,'_',CellType_level,'_',Timepoints_meta,'_Healthy_',Healthy,'.rds', sep = ''))
saveRDS(Go_list_of_lists_DOWN, file = paste(fig_names, '_GO_list_DOWN_',First_Split,'_',CellType_level,'_',Timepoints_meta,'_Healthy_',Healthy,'.rds', sep = ''))


Go_list_of_lists_UP = readRDS(file = paste(fig_names, '_GO_list_UP_',First_Split,'_',CellType_level,'_',Timepoints_meta,'_Healthy_',Healthy,'.rds', sep = ''))
Go_list_of_lists_DOWN = readRDS(file = paste(fig_names, '_GO_list_DOWN_',First_Split,'_',CellType_level,'_',Timepoints_meta,'_Healthy_',Healthy,'.rds', sep = ''))

# rm(GO_results_df)


# GO_results_df = bind_rows(Go_list_of_lists_UP$`Patient 2`$`CD4 T`, .id = "id")
# GO_results_df = bind_rows(Go_list_of_lists_DOWN[["Cancer"]], .id = "id")
# 
# 
# GO_results_df$pval = as.numeric(GO_results_df$pval)
# GO_results_df = filter(GO_results_df, pval < 0.05)

for(First_Split_ in names(Go_list_of_lists_UP)) {
  
  for(celltype in names(Go_list_of_lists_UP[[First_Split_]])) {
    
    file_name = paste('GO_Plot_',First_Split_, '_',celltype)
    
    GO_results_df_UP = bind_rows(Go_list_of_lists_UP[[First_Split_]][[celltype]], .id = "id")
    GO_results_df_DOWN = bind_rows(Go_list_of_lists_DOWN[[First_Split_]][[celltype]], .id = "id")
    try({
    GO_plot(GO_results_df_UP, file_name = paste(file_name,'_UP', sep = ''), title = paste(First_Split_,' ',celltype,' Upregulated', sep = ''))
    dev.off()
    },silent = T)
    try({
    GO_plot(GO_results_df_DOWN, file_name = paste(file_name,'_DOWN', sep = ''), title = paste(First_Split_,' ',celltype,' Downregulated', sep = ''))
    dev.off()
    },silent = T)
    
    
  }
}

dev.off()

Go_list_of_lists_UP_flipped = list()
Go_list_of_lists_DOWN_flipped = list()

for(First_Split_ in names(Go_list_of_lists_UP)) {
  
  for(celltype in names(Go_list_of_lists_UP[[First_Split_]])) {
    
    for(timepoint in names(Go_list_of_lists_UP[[First_Split_]][[celltype]])) {
      
      up = Go_list_of_lists_UP[[First_Split_]][[celltype]][[timepoint]]
      down = Go_list_of_lists_DOWN[[First_Split_]][[celltype]][[timepoint]]
      
      # timepoint = gsub('d14', 'post_inf',timepoint)
      # timepoint = gsub('d28', 'post_inf',timepoint)
      # timepoint = gsub('d7', 'post_inf',timepoint)
      # timepoint = gsub('d-21 to -7', 'pre_chem',timepoint)
      # timepoint = gsub('d-1', 'pre_inf',timepoint)
      
    
      Go_list_of_lists_UP_flipped[[timepoint]][[celltype]][[First_Split_]] = up
      Go_list_of_lists_DOWN_flipped[[timepoint]][[celltype]][[First_Split_]] = down
    
    }
    
  }
  
}

# setwd('GO_Analysis_CellType_L1')

for(First_Split_ in names(Go_list_of_lists_UP_flipped)) {
  
  for(celltype in names(Go_list_of_lists_UP_flipped[[First_Split_]])) {
    
    file_name = paste('GO_Plot_Timepoint_',First_Split_, '_',celltype)
    
    GO_results_df_UP = bind_rows(Go_list_of_lists_UP_flipped[[First_Split_]][[celltype]], .id = "id")
    GO_results_df_DOWN = bind_rows(Go_list_of_lists_DOWN_flipped[[First_Split_]][[celltype]], .id = "id")
    
    try({
    GO_plot(GO_results_df_UP, file_name = paste(file_name,'_UP', sep = ''), title = paste(First_Split_,' ',celltype,' Upregulated', sep = ''))
    },silent = F)
    try({
    GO_plot(GO_results_df_DOWN, file_name = paste(file_name,'_DOWN', sep = ''), title = paste(First_Split_,' ',celltype,' Downregulated', sep = ''))
    },silent = F)
    
    
  }
}

closeAllConnections()


gsub('d14', 'post_inf',timepoint)


GO_results = list()

clusters = as.character(sort(as.numeric(unique(object_AXOSpatial_seurat@meta.data$SCT_snn_res.0.8))-1))


clusters = as.character(unique(nCoV.integrated@meta.data$CD38_PDCD1))

Idents(nCoV.integrated) = "CD38_PDCD1"


de_markers = FindAllMarkers(nCoV.integrated, only.pos = TRUE)



de_markers$cluster = as.character(de_markers$cluster)





for (cluster_ in clusters) {
  
  de_markers_cluster = dplyr::filter(de_markers, cluster == cluster_)$gene
  
  if (length(de_markers_cluster) != 0) {
    
    GOterms = topGOterms(fg.genes = de_markers_cluster, bg.genes = rownames(nCoV.integrated), organism = "Human")
    
    GO_results[[cluster_]] = GOterms$res.table  %>% top_n(n = 10, wt = pval)
  }
  
}
################################################################## FIND MARKERS BETWEEN CONDITIONS ##################################################################


for(Patient in unique(nCoV.integrated@meta.data$Patient)) {
  Idents(nCoV.integrated) = "Patient"
  
  nCoV.integrated.Patient = subset(nCoV.integrated, idents = c(Patient,"Healthy_BM"), invert = FALSE)
  
  
  ##length(as.character(unique(nCoV.integrated.Patient@meta.data$Day)))
  
  Idents(nCoV.integrated.Patient) = "BTM_Others"

  for(celltype in as.character(unique(Idents(nCoV.integrated.Patient)))) {


    nCoV.integrated.Patient.celltype = subset(nCoV.integrated.Patient, idents = celltype, invert = FALSE)

    Idents(nCoV.integrated.Patient.celltype) = "Day"

    Days = as.character(unique(Idents(nCoV.integrated.Patient.celltype)))
   for(i in 1:(length(as.character(unique(nCoV.integrated.Patient@meta.data$Day)))-1)) {


     markers = FindMarkers(nCoV.integrated.Patient.celltype, ident.1 = Days[i+1], ident.2 = Days[i])
     
     markers1 = cbind(rownames(markers), markers)
     colnames(markers1)[1] = 'Gene'
     rownames(markers1) = c(1:nrow(markers1))
     markers1 = subset(markers1, p_val_adj < 0.05)
     markers1 = cbind(c(1:nrow(markers1)), markers1)
     markers1$`c(1:nrow(markers1))` = sub('\\d+','_', markers1$`c(1:nrow(markers1))`)
     colnames(markers1)[1] = 'Original out p_adj < 0.05 >>'
     
     
     markers2 = markers1
     markers2$pct.1 = NULL
     markers2$pct.2 = NULL
     markers2$p_val = NULL
     markers2 = markers2[order(-markers2$avg_log2FC), ]
     
     colnames(markers2)[1] = 'Sorted logFC asc >>'
     
     
     markers3 = markers2
     markers3 = markers3[order(markers3$avg_log2FC), ]
     colnames(markers3)[1] = 'Sorted logFC desc >>'
     
     markers = cbind(markers1,markers2,markers3)
     
     write.csv(markers, file = paste(fig_names, '_',name,'_',celltype,'_',Days[i+1],'_vs_ctrl_',Days[i], "_Markers.csv", sep = ''))




   }
  }
  
  Idents(nCoV.integrated.Patient) = 'predicted.celltype.l2'
  
  for(celltype in as.character(unique(Idents(nCoV.integrated.Patient)))) {
    
    
    nCoV.integrated.Patient.celltype = subset(nCoV.integrated.Patient, idents = celltype, invert = FALSE)
    
    Idents(nCoV.integrated.Patient.celltype) = "Day"
    
    Days = as.character(unique(Idents(nCoV.integrated.Patient.celltype)))
    for(i in 1:(length(as.character(unique(nCoV.integrated.Patient@meta.data$Day)))-1)) {
      
      try({
      markers = FindMarkers(nCoV.integrated.Patient.celltype, ident.1 = Days[i+1], ident.2 = Days[i],)
      
      markers1 = cbind(rownames(markers), markers)
      colnames(markers1)[1] = 'Gene'
      rownames(markers1) = c(1:nrow(markers1))
      markers1 = subset(markers1, p_val_adj < 0.05)
      markers1 = cbind(c(1:nrow(markers1)), markers1)
      markers1$`c(1:nrow(markers1))` = sub('\\d+','_', markers1$`c(1:nrow(markers1))`)
      colnames(markers1)[1] = 'Original out p_adj < 0.05 >>'
      
      
      markers2 = markers1
      markers2$pct.1 = NULL
      markers2$pct.2 = NULL
      markers2$p_val = NULL
      markers2 = markers2[order(-markers2$avg_log2FC), ]
      
      colnames(markers2)[1] = 'Sorted logFC asc >>'
      
      
      markers3 = markers2
      markers3 = markers3[order(markers3$avg_log2FC), ]
      colnames(markers3)[1] = 'Sorted logFC desc >>'
      
      markers = cbind(markers1,markers2,markers3)
      
      
      write.csv(markers, file = paste(fig_names, '_',name,'_',celltype,'_',Days[i+1],'_vs_ctrl_',Days[i], "_Markers.csv", sep = ''))
      },silent = T)
      
      
      
    }
  }
  
  
  
  
}

################################################################## FIND MARKERS WITHIN CONDITIONS ##################################################################
nCoV.integrated@meta.data$CD38_exp = nCoV.integrated@meta.data$predicted.celltype.l1

poscells <- WhichCells(nCoV.integrated, expression = CD38 > 0)

nCoV.integrated@meta.data[poscells, 'CD38_exp'] = 'CD38+'

negcells <- WhichCells(nCoV.integrated, expression = CD38 == 0)

nCoV.integrated@meta.data[negcells, 'CD38_exp'] = 'CD38-'
###################################################################################################################################


nCoV.integrated@meta.data$scvfCAR123_exp = nCoV.integrated@meta.data$BTM_Others

poscells <- WhichCells(nCoV.integrated, expression = scvfCAR123 > 0)

nCoV.integrated@meta.data[poscells, 'scvfCAR123_exp'] = 'CAR T Cells'

# negcells <- WhichCells(nCoV.integrated, expression = scvfCAR123 == 0)
# 
# nCoV.integrated@meta.data[negcells, 'CD38_exp'] = 'CD38-'

nCoV.integrated

First_Split = "Dataset_Type"#"Patient"
Second_Split = "Dataset"#"Day"

unique(nCoV.integrated@meta.data[, First_Split])
unique(nCoV.integrated@meta.data[, Second_Split])
unique(nCoV.integrated@meta.data[, celltype_level])


Idents(nCoV.integrated) = Second_Split

nCoV.integrated = subset(nCoV.integrated, idents = c('d7','d14','d28'))

Within_Ident = "CD38_exp"
Condition = 'CD38+'
Control = 'CD38-'
celltype_level = "predicted.celltype.l1"

markers_list_of_lists = list()


for(First_Split_ in unique(nCoV.integrated@meta.data[, First_Split])) {
  
  Idents(nCoV.integrated) = First_Split
  
  nCoV.integrated.Dataset = subset(nCoV.integrated, idents = c(First_Split_), invert = FALSE)
  
  
  ##length(as.character(unique(nCoV.integrated.Patient@meta.data$Day)))

  
  Idents(nCoV.integrated.Dataset) = celltype_level
  
  for(celltype in as.character(unique(Idents(nCoV.integrated.Dataset)))) {
    
    
    nCoV.integrated.Dataset.celltype = subset(nCoV.integrated.Dataset, idents = celltype, invert = FALSE)
    
    Idents(nCoV.integrated.Dataset.celltype) = Second_Split
    
    # Days = as.character(unique(Idents(nCoV.integrated.Dataset.celltype)))
    
    for(Second_Split_ in (unique(nCoV.integrated.Dataset.celltype@meta.data[, Second_Split]))) {
      
      # Idents(nCoV.integrated.Dataset) = Second_Split
      
      
      nCoV.integrated.Dataset.celltype.2nd_split = subset(nCoV.integrated.Dataset.celltype, idents = Second_Split_, invert = FALSE)
      
      # markers = data.frame()
      # 
      # markers_list_of_lists[[First_Split_]][[celltype]][[Second_Split_]] = markers
      
      
       try({
         
        Idents(nCoV.integrated.Dataset.celltype.2nd_split) = Within_Ident
         
        
        markers = FindMarkers(nCoV.integrated.Dataset.celltype.2nd_split, ident.1 = Condition, ident.2 = Control)
        
        markers1 = cbind(rownames(markers), markers)
        colnames(markers1)[1] = 'Gene'
        rownames(markers1) = c(1:nrow(markers1))
        markers1 = subset(markers1, p_val_adj < 0.05)
        markers1 = cbind(c(1:nrow(markers1)), markers1)
        markers1$`c(1:nrow(markers1))` = sub('\\d+','_', markers1$`c(1:nrow(markers1))`)
        colnames(markers1)[1] = 'Original out p_adj < 0.05 >>'
        
        
        markers2 = markers1
        markers2$pct.1 = NULL
        markers2$pct.2 = NULL
        markers2$p_val = NULL
        markers2 = markers2[order(-markers2$avg_log2FC), ]
        
        colnames(markers2)[1] = 'Sorted logFC asc >>'
        
        
        markers3 = markers2
        markers3 = markers3[order(markers3$avg_log2FC), ]
        colnames(markers3)[1] = 'Sorted logFC desc >>'
        
        markers = cbind(markers1,markers2,markers3)
        
        markers_list_of_lists[[First_Split_]][[celltype]][[Second_Split_]] = markers
        
        
        write.csv(markers, file = paste(fig_names, '_',First_Split_,'_',Second_Split_,'_',celltype,'_',Within_Ident,'_',Condition,'_vs_ctrl_',Control, "_Markers.csv", sep = ''))
        
      },silent = F)
      
      
      
      
      
    }
  }
  
  
  
  
}



Go_list_of_lists_UP = list()
Go_list_of_lists_DOWN = list()



for(First_Split_ in names(markers_list_of_lists)) {
  
  for(celltype in names(markers_list_of_lists[[First_Split_]])) {
    
    for(Second_Split_ in names(markers_list_of_lists[[First_Split_]][[celltype]])) {
      
  markers = markers_list_of_lists[[First_Split_]][[celltype]][[Second_Split_]][,1:4]
  markers_UP = dplyr::filter(markers, avg_log2FC > 0.25)$Gene
  markers_DOWN = dplyr::filter(markers, avg_log2FC < -0.25)$Gene
  
  if (length(markers_UP) != 0) {
    
    GOterms = topGOterms(fg.genes = markers_UP, bg.genes = rownames(nCoV.integrated), organism = "Human")
    
    Go_list_of_lists_UP[[First_Split_]][[celltype]][[Second_Split_]] = GOterms$res.table  %>% top_n(n = 15, wt = pval)
    
  }
  
  if (length(markers_DOWN) != 0) {
    
    GOterms = topGOterms(fg.genes = markers_DOWN, bg.genes = rownames(nCoV.integrated), organism = "Human")
    
    Go_list_of_lists_DOWN[[First_Split_]][[celltype]][[Second_Split_]] = GOterms$res.table  %>% top_n(n = 15, wt = pval)
    
  }
    
    }
    
  }
  
}

for(First_Split_ in names(Go_list_of_lists_UP)) {
  
  for(celltype in names(Go_list_of_lists_UP[[First_Split_]])) {
    
    file_name = paste('GO_Plot_',First_Split_, '_',celltype,'_',Within_Ident,'_',Condition,'_vs_Control_',Control, sep = '')
    
    GO_results_df_UP = bind_rows(Go_list_of_lists_UP[[First_Split_]][[celltype]], .id = "id")
    GO_results_df_DOWN = bind_rows(Go_list_of_lists_DOWN[[First_Split_]][[celltype]], .id = "id")
    try({
      GO_plot(GO_results_df_UP, file_name = paste(file_name,'_UP', sep = ''), title = paste(First_Split_,' ',celltype,' ', Condition,' vs Control ',Control,' Upregulated', sep = ''))
      dev.off()
    },silent = T)
    try({
      GO_plot(GO_results_df_DOWN, file_name = paste(file_name,'_DOWN', sep = ''), title = paste(First_Split_,' ',celltype,' ', Condition,' vs Control ',Control,' Downregulated', sep = ''))
      dev.off()
    },silent = T)
    try({
    dev.off()
    dev.off()
    }, silent = T)
  }
}


Within_Ident = "CD38_exp"
Condition = 'CD38+'
Control = 'CD38-'


Go_list_of_lists_UP_collapsed = list()


for(First_Split_ in names(Go_list_of_lists_UP)) {
  
  for(celltype in names(Go_list_of_lists_UP[[First_Split_]])) {
    
    for(Second_Split_ in names(Go_list_of_lists_UP[[First_Split_]][[celltype]])) {

      Go_list_of_lists_UP_collapsed[[paste(First_Split_ ,celltype, Second_Split_,sep = ' ')]] = Go_list_of_lists_UP[[First_Split_]][[celltype]][[Second_Split_]]
      
    }
  }
}
GO_results_df = bind_rows(Go_list_of_lists_UP_collapsed, .id = "id")


GO_plot(GO_results_df, file_name = 'scvfCAR123 T cells vs Ctrl T cells AML', title = 'scvfCAR123 T cells vs Ctrl T cells AML')


Go_list_of_lists_DOWN_collapsed = list()


for(First_Split_ in names(Go_list_of_lists_DOWN)) {
  
  for(celltype in names(Go_list_of_lists_DOWN[[First_Split_]])) {
    
    for(Second_Split_ in names(Go_list_of_lists_DOWN[[First_Split_]][[celltype]])) {
      
      Go_list_of_lists_DOWN_collapsed[[paste(First_Split_ ,celltype, Second_Split_,sep = ' ')]] = Go_list_of_lists_DOWN[[First_Split_]][[celltype]][[Second_Split_]]
      
    }
  }
}
GO_results_df = bind_rows(Go_list_of_lists_DOWN_collapsed, .id = "id")


GO_plot(GO_results_df, file_name = 'scvfCAR123 T cells vs Ctrl T cells AML DOWN', title = 'scvfCAR123 T cells vs Ctrl T cells AML DOWN')


rm(GO_results_df)



GO_results_df = bind_rows(Go_list_of_lists_UP, .id = "id")



GO_results_df = bind_rows(Go_list_of_lists_DOWN[["Cancer"]], .id = "id")


GO_results_df$pval = as.numeric(GO_results_df$pval)
GO_results_df = filter(GO_results_df, pval < 0.05)





GO_results = list()

clusters = as.character(sort(as.numeric(unique(object_AXOSpatial_seurat@meta.data$SCT_snn_res.0.8))-1))


clusters = as.character(unique(nCoV.integrated@meta.data$CD38_PDCD1))

Idents(nCoV.integrated) = "CD38_PDCD1"


de_markers = FindAllMarkers(nCoV.integrated, only.pos = TRUE)



de_markers$cluster = as.character(de_markers$cluster)





for (cluster_ in clusters) {
  
  de_markers_cluster = dplyr::filter(de_markers, cluster == cluster_)$gene
  
  if (length(de_markers_cluster) != 0) {
    
    GOterms = topGOterms(fg.genes = de_markers_cluster, bg.genes = rownames(nCoV.integrated), organism = "Human")
    
    GO_results[[cluster_]] = GOterms$res.table  %>% top_n(n = 10, wt = pval)
  }
  
}



GO_results_df = bind_rows(GO_results, .id = "id")
GO_results_df$pval = as.numeric(GO_results_df$pval)
GO_results_df = filter(GO_results_df, pval < 0.05)

# GO_results_df$id = as.numeric(GO_results_df$id)
# 
# 
# GO_results_df = GO_results_df[order(GO_results_df$id), ]

############################# GO FIGURE #######################################################################################
# GO_results_df$pval = as.numeric(GO_results_df$pval)
# GO_plot = function(DGE_table, file_name = '', gene_list = NA, Control_name = "Control", Condition_name = 'Condition', Comparison_name = '', FileName = "Volcano"){
# }

GO_results_df$pval = as.numeric(GO_results_df$pval)
GO_results_df = filter(GO_results_df, pval < 0.05)
data1 <- GO_results_df

colnames(GO_results_df)[1] = 'Cluster ID'

S1 <- ggplot(data1, aes(x = id, y = Term , color = -log(pval), size = log(Significant))) + 
  geom_point() + xlab("Dataset") + #scale_x_continuous(breaks = 1:12) +
  ylab(NULL) + 
  labs(colour = "-Log10(P.adjust)", size = "Log10(Number of Genes)") +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.text=element_text(size=10,colour ="black"),
        axis.title=element_text(size=14,face="bold",colour ="black"),
        legend.text=element_text(size=12),
        legend.title = element_text(colour="black", size=14, face="bold"),
        axis.line = element_line(colour = "black", size =1),
        panel.grid.major = element_line(colour = "lightgray", size = 0.5), 
        # panel.grid.minor = element_line(colour = "lightgray", size = 0.5),
        panel.background = element_blank())

# S1 = S1+scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab", limit = c(min(data1$pval), max(data1$pval)), n.breaks = 5)
S1+scale_size(range = c(0, length(unique(data1$id))))
S1

############################# GO FIGURE FUNCTION#######################################################################################

GO_plot = function(GO_results_df, file_name = '', gene_list = NA, Control_name = "Control", Condition_name = 'Condition', Comparison_name = '', FileName = "Volcano")
  
  {

data1 <- GO_results_df

colnames(GO_results_df)[1] = 'Cluster ID'

S1 <- ggplot(data1, aes(x = id, y = Term , color = -log(pval), size = log(Significant))) + 
  geom_point() + xlab("Dataset") + #scale_x_continuous(breaks = 1:12) +
  ylab(NULL) + 
  labs(colour = "-Log10(P.adjust)", size = "Log10(Number of Genes)") +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.text=element_text(size=10,colour ="black"),
        axis.title=element_text(size=14,face="bold",colour ="black"),
        legend.text=element_text(size=12),
        legend.title = element_text(colour="black", size=14, face="bold"),
        axis.line = element_line(colour = "black", size =1),
        panel.grid.major = element_line(colour = "lightgray", size = 0.5), 
        # panel.grid.minor = element_line(colour = "lightgray", size = 0.5),
        panel.background = element_blank())

# S1 = S1+scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab", limit = c(min(data1$pval), max(data1$pval)), n.breaks = 5)
S1+scale_size(range = c(0, length(unique(data1$id))))
S1

}

############################# SPLIT FIGURES ###############################################################################################








split_list = c("Day")

for(Patient in unique(nCoV.integrated@meta.data$Patient)) {
  Idents(nCoV.integrated) = "Patient"
  
  nCoV.integrated.Patient = subset(nCoV.integrated, idents = c(Patient,"Healthy_BM"), invert = FALSE)
  
  
  
  for(split_ in split_list) {
    for(meta in meta_data_list_plots) {
      Idents(nCoV.integrated.Patient) <- meta
      dpi = 200
      png(file=paste(Patient,'_Split_',split_,'_', fig_names, '_',paste(meta,'.png', sep = ''), sep = ''), width = dpi*24, height = dpi*8, units = "px",res = dpi,type='cairo')
      print(DimPlot(object = nCoV.integrated.Patient, group.by = "ident",split.by = split_, label = TRUE ,pt.size = 1.5, raster=FALSE, repel = TRUE))
      dev.off()
    }
  }
  
}

################################################################## UMAP - DOTPLOTS = VIOLIN PLOTS ##################################################################


split_list = c("Day")
marker_genes = c('scvfCAR123','MKI67','IL2','TNF','IFNG','CSF2','IL6','IL10','TGFB1','TGFB2','TGFB3','IL2RA','TNFRSF1A','IFNGR1','CSF2RA','IL6R','IL10RA','IL10RB','TGFBR1','TGFBR2','SLC2A5','SLC2A8','ALDOB','ALDOC','KHK','AKR1B1','SORD')
marker_genes = c('IL3','KITLG','CXCR4','CXCL12')

marker_genes = c('FAS','FASLG')

FeaturePlot(object = nCoV.integrated, features = marker_genes)

DotPlot(object = nCoV.integrated, features = marker_genes, split.by = "Day_New", cols = 'Blues')



levels_l1 = unique(nCoV.integrated$predicted.celltype.l1)
levels_l2 = unique(nCoV.integrated$predicted.celltype.l2)


unique(nCoV.integrated$predicted.ct_BM_AML)


First_Split = c("Patient")
split_list = c("Day_New")

violin_plot = FALSE
dot_plot = TRUE

name = fig_names

Idents(nCoV.integrated) = First_Split

for(first_split_ in unique(nCoV.integrated@meta.data[, First_Split])) {
  
  
  nCoV.integrated.Patient = subset(nCoV.integrated, idents = first_split_)
  

  
  for(split_ in split_list) {
    for(gene in marker_genes) {
      Idents(nCoV.integrated.Patient) <- 'predicted.celltype.l2'
      dpi = 200
      png(file=paste(name,'_',first_split_ ,'_Split_',split_,'_', fig_names, '_',paste(gene,'.png', sep = ''), sep = ''), width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
      print(FeaturePlot(object = nCoV.integrated.Patient, features = gene, split.by = split_, label = TRUE ,pt.size = 1.5, raster=FALSE, repel = TRUE, order = TRUE))
      dev.off()
      
      if (violin_plot == TRUE) {
      
      png(file=paste(name,'_',first_split_ ,'_Violin_l1_Split_',split_,'_', fig_names, '_',paste(gene,'.png', sep = ''), sep = ''), width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
      print(VlnPlot(object = nCoV.integrated.Patient, features = gene, split.by = split_, group.by = 'predicted.celltype.l1',log = TRUE, adjust = 2))
      dev.off()
      png(file=paste(name,'_',first_split_ ,'_Violin_l2_Split_',split_,'_', fig_names, '_',paste(gene,'.png', sep = ''), sep = ''), width = dpi*32, height = dpi*9, units = "px",res = dpi,type='cairo')
      print(VlnPlot(object = nCoV.integrated.Patient, features = gene, split.by = split_, group.by = 'predicted.celltype.l2',log = TRUE, adjust = 2))
      dev.off()
    }
    }
    
    ##################################
    
    
    Idents(nCoV.integrated.Patient) <- 'Day'
    
    
    if (dot_plot == TRUE) {
    
    for (Day in unique(Idents(nCoV.integrated.Patient))) {
      
      nCoV.integrated.Patient.Day = subset(nCoV.integrated.Patient, idents = Day)
      
      # setdiff(levels_l2,as.character(unique(Idents(nCoV.integrated.Patient.Day))))
      
      missing_idents = setdiff(levels_l1,as.character(unique(nCoV.integrated.Patient.Day$predicted.celltype.l1)))
      
      for (missing_ident in missing_idents) {
        
        nCoV.integrated.Patient.Day$predicted.celltype.l1[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l1[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l1[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l1[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l1[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l1[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l1[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l1[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l1[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l1[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l1[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l1[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l1[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l1[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l1[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
      }
      
      
      # setdiff(levels_l2,as.character(unique(Idents(nCoV.integrated.Patient.Day))))
    
      Idents(nCoV.integrated.Patient.Day) = 'predicted.celltype.l1'
      levels(nCoV.integrated.Patient.Day) = levels_l1
      
      
      
      png(file=paste(name,'_',first_split_ ,'_',Day,'_Dotplot_l1_Split_', fig_names, '_',paste('Marker_Genes','.png', sep = ''), sep = ''), width = dpi*16, height = dpi*9, units = "px",res = dpi,type='cairo')
      print(DotPlot(object = nCoV.integrated.Patient.Day, features = marker_genes) + theme(axis.text.x = element_text(angle = 90)))
      dev.off()

    }
    
    for (Day in unique(Idents(nCoV.integrated.Patient))) {
      
      nCoV.integrated.Patient.Day = subset(nCoV.integrated.Patient, idents = Day)
      
      # setdiff(levels_l2,as.character(unique(nCoV.integrated.Patient.Day$predicted.celltype.l2)))
      
      missing_idents = setdiff(levels_l2,as.character(unique(nCoV.integrated.Patient.Day$predicted.celltype.l2)))
      
      for (missing_ident in missing_idents) {
      
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
        nCoV.integrated.Patient.Day$predicted.celltype.l2[sample(1:nrow(nCoV.integrated.Patient.Day), 1)] = missing_ident
      }
      
      Idents(nCoV.integrated.Patient.Day) = 'predicted.celltype.l2'
      levels(nCoV.integrated.Patient.Day) = levels_l2
      

      # setdiff(levels_l2,as.character(unique(nCoV.integrated.Patient.Day$predicted.celltype.l2)))
      

      
      
      png(file=paste(name,'_',first_split_ ,'_',Day,'_Dotplot_l2_Split_', fig_names, '_',paste('Marker_Genes','.png', sep = ''), sep = ''), width = dpi*16, height = dpi*9, units = "px",res = dpi,type='cairo')
      print(DotPlot(object = nCoV.integrated.Patient.Day, features = marker_genes) + theme(axis.text.x = element_text(angle = 90)))
      dev.off()
    }
    
    unique(nCoV.integrated.Patient$Day)
    }
  }
  
}
try(dev.off(), silent = TRUE)
try(dev.off(), silent = TRUE)
try(dev.off(), silent = TRUE)


####################################################################################################################################



############################################################# FIND MARKERS POS/ NEG gene expression #######################################################################


nCoV.integrated@meta.data$scvfCAR123_exp = nCoV.integrated@meta.data$BTM_Others

poscells <- WhichCells(nCoV.integrated, expression = scvfCAR123 > 0)

nCoV.integrated@meta.data[poscells, 'scvfCAR123_exp'] = 'CAR T Cells'



nCoV.integrated@meta.data$CD38_exp = nCoV.integrated@meta.data$BTM_Others

poscells <- WhichCells(nCoV.integrated, expression = scvfCAR123 > 0)

nCoV.integrated@meta.data[poscells, 'scvfCAR123_exp'] = 'CAR T Cells'

negcells <- WhichCells(nCoV.integrated, expression = scvfCAR123 == 0)

nCoV.integrated@meta.data[negcells, 'CD38_exp'] = 'CD38-'


###########################


Idents(nCoV.integrated) = "Dataset"

fig_names = 'Kim_CD38'



if (subset_celltype_TRUE_FALSE == FALSE) {
  
  subset_celltype = 'All_Cells'
  
}


Idents(nCoV.integrated) = "Dataset_Type"

as.character(unique(Idents(nCoV.integrated)))

nCoV.integrated.subset = subset(nCoV.integrated, idents = dataset, invert = FALSE)


Idents(nCoV.integrated.subset) = "Dataset"


for (type in as.character(unique(Idents(nCoV.integrated)))) {
  
  nCoV.integrated.subset = subset(nCoV.integrated, idents = type, invert = FALSE)
  gene = 'CD38'
  Idents(nCoV.integrated.subset) = 'predicted.celltype.l1'
  dpi = 200
  png(file=paste(fig_names,'_',paste(unlist(as.character(unique(Idents(nCoV.integrated.subset)))), collapse = '_'),'_feature_plot_',type,'_',subset_celltype,'_', gene,'.png',sep = ''), width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
  print(FeaturePlot(object = nCoV.integrated.subset, features = gene, split.by = "Dataset", label = TRUE ,pt.size = 1.5, raster=FALSE, repel = TRUE, order = TRUE))
  dev.off()
  
  
  Idents(nCoV.integrated.subset) = "Dataset"
  
  markers_cd38_list = list()
  
  markers_cd38_list_sig_genes_UP = list()
  
  markers_cd38_list_sig_genes_DOWN = list()
  
  for(dataset in as.character(unique(Idents(nCoV.integrated.subset)))) {
    
    if (file.exists(paste(fig_names, '_',dataset,'_',subset_celltype,'_CD38p_vs_ctrl_CD38n_Markers.csv', sep = ''))) {
      
      nCoV.integrated.dataset = subset(nCoV.integrated.subset, idents = dataset, invert = FALSE)
      
      
      
      markers = read.csv2(paste(fig_names, '_',dataset,'_',subset_celltype,'_CD38p_vs_ctrl_CD38n_Markers.csv', sep = ''), sep = ',', row.names = 1)
      
      markers_volcano = read.csv2(paste(fig_names, '_',dataset,'_',subset_celltype,'_CD38p_vs_ctrl_CD38n_Markers_Volcano.csv', sep = ''), sep = ',', row.names = 1)
      
      
      
      markers_cd38_list[[dataset]] = markers
      
      markers_cd38_list_sig_genes_UP[[dataset]] = markers$Gene[which(markers$avg_log2FC>0)]
      
      
      markers_cd38_list_sig_genes_DOWN[[dataset]] = markers$Gene[which(markers$avg_log2FC<0)]
      
      
    } else {
      
      nCoV.integrated.dataset = subset(nCoV.integrated.subset, idents = dataset, invert = FALSE)
      
      Idents(nCoV.integrated.dataset) = "CD38_exp"
      
      
      markers = FindMarkers(nCoV.integrated.dataset, ident.1 = 'CD38+', ident.2 = 'CD38-')
      
      markers_volcano = markers 
      
      markers_volcano$Gene = rownames(markers_volcano)
      
      write.csv(markers_volcano, file = paste(fig_names, '_',dataset,'_',subset_celltype,'_CD38p_vs_ctrl_CD38n_Markers_Volcano.csv', sep = ''))
      
      markers = subset(markers, p_val_adj < 0.05)
      
      
      markers_cd38_list_sig_genes_UP[[dataset]] = row.names(markers)[which(markers$avg_log2FC>0)]
      
      
      markers_cd38_list_sig_genes_DOWN[[dataset]] = row.names(markers)[which(markers$avg_log2FC<0)]
      
      markers1 = cbind(rownames(markers), markers)
      colnames(markers1)[1] = 'Gene'
      rownames(markers1) = c(1:nrow(markers1))
      markers1 = cbind(c(1:nrow(markers1)), markers1)
      markers1$`c(1:nrow(markers1))` = sub('\\d+','_', markers1$`c(1:nrow(markers1))`)
      colnames(markers1)[1] = 'Original out p_adj < 0.05 >>'
      
      
      markers2 = markers1
      markers2$pct.1 = NULL
      markers2$pct.2 = NULL
      markers2$p_val = NULL
      markers2 = markers2[order(-markers2$avg_log2FC), ]
      
      colnames(markers2)[1] = 'Sorted logFC asc >>'
      
      
      markers3 = markers2
      markers3 = markers3[order(markers3$avg_log2FC), ]
      colnames(markers3)[1] = 'Sorted logFC desc >>'
      
      markers = cbind(markers1,markers2,markers3)
      
      write.csv(markers, file = paste(fig_names, '_',dataset,'_',subset_celltype,'_CD38p_vs_ctrl_CD38n_Markers.csv', sep = ''))
      
      markers_cd38_list[[dataset]] = markers
      
      
    }
      
    volcano_plot(markers_volcano, FileName = paste(fig_names, '_',dataset,'_',subset_celltype,'_CD38p_vs_ctrl_CD38n_Volcano', sep = ''), Control_name = 'CD38-', Condition_name = 'CD38+')
    
      # row.names(markers1)[which(markers1$avg_log2FC<0)]
      
      # 
      # png(file=paste(name,'_',Day,'_Dotplot_l2_Split_', fig_names, '_',paste('Marker_Genes','.png', sep = ''), sep = ''), width = dpi*16, height = dpi*9, units = "px",res = dpi,type='cairo')
      # print(DotPlot(object = nCoV.integrated.Patient.Day, features = marker_genes) + theme(axis.text.x = element_text(angle = 90)))
      # dev.off()
      
      
      # rm(nCoV.integrated.dataset)
      # gc()
      
      
      
      
      
      
      
    
    
  }
  
  if (length(markers_cd38_list_sig_genes_UP)<5) {
    
    # ggvenn(markers_cd38_list_sig_genes_UP,show_elements=FALSE,stroke_color="Black", stroke_linetype="solid")
    dpi = 200
    png(file=paste(fig_names,'_',paste(unlist(as.character(unique(Idents(nCoV.integrated.subset)))), collapse = '_'),'_',subset_celltype,'_Venn_Diagram_UP','_CD38p_vs_ctrl_CD38n_Markers_UP.png',sep = ''), width = dpi*16, height = dpi*9, units = "px",res = dpi,type='cairo')
    print(ggvenn(markers_cd38_list_sig_genes_UP,show_elements=FALSE,stroke_color="Black", stroke_linetype="solid"))
    dev.off()
    
    png(file=paste(fig_names,'_',paste(unlist(as.character(unique(Idents(nCoV.integrated.subset)))), collapse = '_'),'_',subset_celltype,'_Venn_Diagram_DOWN','_CD38p_vs_ctrl_CD38n_Markers_Down.png',sep = ''), width = dpi*16, height = dpi*9, units = "px",res = dpi,type='cairo')
    print(ggvenn(markers_cd38_list_sig_genes_DOWN,show_elements=FALSE,stroke_color="Black", stroke_linetype="solid"))
    dev.off()
    
    
    
  }
}








install.packages('ggvenn')
library("ggvenn")

# use list as input 
H <-list('Bus'=c(6,7,3),'Truck'=c(4,3,9),
         'Cycle'=c(10,3,2,8),'Car'=c(7,5,4,3))




# create customised venn diagram
ggvenn(markers_cd38_list_sig_genes_UP[1:5],show_elements=FALSE,stroke_color="Red",
       stroke_linetype="solid")


results  <- read.csv(file="AML_vs_Normal_tcellsF.csv", header=TRUE, sep=",", fileEncoding="UTF-8-BOM")

######################################################### VOLCANO PLOT FUNCTION ###########################################################################

volcano_plot = function(DGE_table, file_name = '', gene_list = NA, Control_name = "Control", Condition_name = 'Condition', Comparison_name = '', FileName = "Volcano"){
  
  results = DGE_table
  
  results$p_val_adj = as.numeric(results$p_val_adj)
  results$avg_log2FC = as.numeric(results$avg_log2FC)
  
  if (is.na(gene_list)) {
    
    results = mutate(results, Genes=ifelse(results$p_val_adj<0.05, "Significant", "Non-Significant"))
    res = results[order(results$p_val_adj, decreasing = FALSE), ]
    list = res$Gene[1:5]
    
  } else{
    
    results = mutate(results, Genes=ifelse(results$Gene %in% gene_list, "List", "Others"))
    list = gene_list
    
  }
  
  
  
  results4 = results[order(results$Genes, decreasing = FALSE), ]
  
  results4$p_val_adj = as.numeric(results4$p_val_adj)
  results4$avg_log2FC = as.numeric(results4$avg_log2FC)
  
  pvals = results4$p_val_adj
  pvals[pvals == 0] = 2
  
  pvals[pvals == 2] = runif(sum(pvals == 2), 0, min(pvals))
  
  
  results4$p_val_adj = pvals
  
  
  p = ggplot(results4, aes(avg_log2FC, -log10(as.numeric(p_val_adj))),group=Genes) +
    geom_jitter(aes(col=Genes,size=Genes,shape = Genes)) + ylim(-10, round(-log10(min(results4$p_val_adj))/10)*10) + xlim(-12,12) +
    scale_size_manual(values=c(0.5,2)) + scale_shape_manual(values=c(19,19))+
    scale_color_manual(values=c("gray30","red"))+ xlab('Log2 Fold Change') + ylab('Log10 Adjusted Pvalue') + ggtitle(paste(Control_name, "              ",Comparison_name,"                          ", Condition_name)) +theme(plot.title = element_text(hjust = 0.5))
  
  theme_set(theme_gray(base_size = 18)) 
  # dev.off()
  
  p = p+geom_vline(xintercept = 2, linetype="dotted", 
               color = "black", size=1)+
    geom_vline(xintercept = -2, linetype="dotted", 
               color = "black", size=1)+
    geom_hline(yintercept=20, linetype="dotted", color = "black", size=1)+
    # geom_text_repel(color = 'black', size = 5,data=filter(results, Gene %in% results2$V1 & results$Adjusted.P.Value<1e-50),
    geom_text_repel(color = 'black', size = 5,data=filter(results4, Gene %in% list ),#& results$adj.P.Val<1e-80),
                    direction ='both',xlim=c(-10,1000), ylim=c(20,100), aes(label=list),min.segment.length = 0,
                    point.padding = 0.25, box.padding = 2.00,force=5) #,segment.color = "grey50"
  
  
  print(p)
  
  dpi = 100
  png(file=paste(FileName,'.png', sep = ''), width = dpi*16, height = dpi*9, units = "px",res = dpi,type='cairo')
  print(p)
  dev.off()
  
}
####################################################################################################################################

volcano_plot(markers_volcano)

#list = c("ITK","BTK","SAMHD1")
# list = c("ARG1","NOS2","NOX2","TGFB1","TGFB2","TGFB3","IL10","PGE2","IDO1","IDO2","IL4R","IRF8","CEBPB","S100A8","S100A9","RORC","DDIT3","OLR1","RB1")
list = c('PD1', 'TIGIT', 'CTLA4', 'VISTA', 'TIM3', 'LAG3', 'CD226', 'ICOS','TNFRSF4')

# results2  <- read.csv(file="immuno_Sup.txt", header=FALSE, sep=",", fileEncoding="UTF-8-BOM")

# list = results2$V1

results = markers



#gene_list = NA

if (is.na(gene_list)) {
  
  results = mutate(results, Genes=ifelse(results$p_val_adj<0.05, "Significant", "Non-Significant"))
  res = results[order(results$p_val_adj, decreasing = FALSE), ]
  list = res$Gene[1:10]
  
  
} else{
  
  results = mutate(results, Genes=ifelse(results$Gene %in% gene_list, "List", "Others"))
  list = gene_list
}


results4 = results[order(results$Genes, decreasing = TRUE), ]

results4$p_val_adj = as.numeric(results4$p_val_adj)
results4$avg_log2FC = as.numeric(results4$avg_log2FC)

pvals = results4$p_val_adj
pvals[pvals == 0] = 2

pvals[pvals == 2] = runif(sum(pvals == 2), 0, min(pvals))


results4$p_val_adj = pvals


p = ggplot(results4, aes(avg_log2FC, -log10(as.numeric(p_val_adj))),group=Genes) +
  geom_jitter(aes(col=Genes,size=Genes,shape = Genes)) + ylim(-10, round(-log10(min(results4$p_val_adj))/10)*10) + xlim(-12,12) +
  scale_size_manual(values=c(2,0.5)) + scale_shape_manual(values=c(19,19))+
  scale_color_manual(values=c("red","gray30"))+ xlab('Log2 Fold Change') + ylab('Log10 Adjusted Pvalue')
theme_set(theme_gray(base_size = 18)) 
# dev.off()

p+geom_vline(xintercept = 2, linetype="dotted", 
             color = "black", size=1)+
  geom_vline(xintercept = -2, linetype="dotted", 
             color = "black", size=1)+
  geom_hline(yintercept=20, linetype="dotted", color = "black", size=1)+
  # geom_text_repel(color = 'black', size = 5,data=filter(results, Gene %in% results2$V1 & results$Adjusted.P.Value<1e-50),
  geom_text_repel(color = 'black', size = 5,data=filter(results4, Gene %in% list ),#& results$adj.P.Val<1e-80),
                  direction ='both',xlim=c(-10,1000), ylim=c(20,100), aes(label=list),min.segment.length = 0,
                  point.padding = 0.25, box.padding = 2.00, label.padding = 0.25,force=5) + #,segment.color = "grey50"
  ggtitle(paste(Control_name, "              ",Comparison_name,"                          ", Condition_name)) +theme(plot.title = element_text(hjust = 0.5))



####################################################################################################################################
#### MAP geneset score UMAP ##############################

for (marker_gene in y_chr_genes_vector_filtered) {
dpi = 200
png(file=paste(fig_names, '_',paste(marker_gene,'.png', sep = ''), sep = ''), width = dpi*32, height = dpi*16, units = "px",res = dpi,type='cairo')
print(FeaturePlot(seurat_obj, reduction = "umap", features = marker_gene, label = TRUE, repel = TRUE,order = TRUE, pt.size = 2))
dev.off()

}


y_chr_genes = read.csv2("Y_chr_coding_genes.csv", sep = ',')
y_chr_genes_vector = y_chr_genes$symbol

y_chr_genes_vector_filtered = Reduce(intersect,list(y_chr_genes_vector,rownames(seurat_obj)))




y_chr_genes_vector_filtered = y_chr_genes_vector(y_chr_genes_vector %in% rownames(seurat_obj))

ident_1 = "Patient 10"

Idents(nCoV.integrated) = "Patient"

seurat_obj = subset(nCoV.integrated, idents = ident_1)

FeaturePlot(seurat_obj, features = c('KDM5D','ZFY','TTTY13','TMSB4Y',"UTY"), split.by = ident_split)

DefaultAssay(seurat_obj) = 'RNA'

ident_map = 'predicted.CellType_BM_AML'

ident_split = 'Day_New'

seurat_obj$predicted.CellType_BM_AML

seurat_obj = AddModuleScore(seurat_obj, features = list(y_chr_genes_vector_filtered), search = FALSE, nbin = 1,name = 'Y_chr__')


seurat_obj = AddModuleScore(seurat_obj, features = list(ae_enriched), search = FALSE, nbin = 1,name = 'Test___')


FeaturePlot(seurat_obj, split.by = ident_split, group.by = 'Y_chr_1',)


Idents(seurat_obj) = 'predicted.CellType_BM_AML'


FeaturePlot(seurat_obj, features = 'Y_chr_1', label = TRUE, repel = TRUE,order = TRUE, pt.size = 2) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 8, name = "RdBu")))


seurat_obj$Test_1_Log = log(seurat_obj$Test_1)

FeaturePlot(seurat_obj, features = 'Test___1', label = TRUE, repel = TRUE,order = TRUE, pt.size = 2) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 8, name = "RdBu")))

Idents(seurat_obj) = 'predicted.CellType_BM_AML'
library(tidyverse)
library(RColorBrewer)
# library(SeuratData)

ae_enriched <- FindMarkers(seurat_obj, ident.1 = c("Aberrant erythroid"), verbose = TRUE) %>%
  arrange(-avg_log2FC) %>%
  rownames_to_column(var = "gene") %>%
  pull(gene) %>% 
  .[1:50]





FeaturePlot(seurat_obj, features = 'Y_chr__1', label = TRUE, repel = TRUE, order = TRUE, pt.size = 2) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 20, name = "RdBu")))


Idents(seurat_obj) <- 'predicted.CellType_BM_AML'
dpi = 200
png(file=paste(ident_1,'_Split_',ident_split,'_',paste(ident_map,'.png', sep = ''), sep = ''), width = dpi*24, height = dpi*8, units = "px",res = dpi,type='cairo')
print(DimPlot(object = nCoV.integrated.Patient, group.by = ident_map,split.by = ident_split, label = TRUE ,pt.size = 1.5, raster=FALSE, repel = TRUE))
dev.off()



#################################################### HEATMAP DUAL GENES +/- ################################################################################



nCoV.integrated@meta.data$CD38_PDCD1 = nCoV.integrated@meta.data$BTM_Others

poscells <- WhichCells(nCoV.integrated, expression = sct_CD38 > 0)

nCoV.integrated@meta.data[poscells, 'CD38_exp'] = 'CD38+'


nCoV.integrated.CD38p.PD1p = colnames(subset(nCoV.integrated, subset = CD38 > 0 & PDCD1 > 0))
nCoV.integrated.CD38n.PD1p = colnames(subset(nCoV.integrated, subset = CD38 == 0 & PDCD1 > 0))
nCoV.integrated.CD38p.PD1n = colnames(subset(nCoV.integrated, subset = CD38 > 0 & PDCD1 == 0))
nCoV.integrated.CD38n.PD1n = colnames(subset(nCoV.integrated, subset = CD38 == 0 & PDCD1 == 0))

nCoV.integrated@meta.data$CD38_PD1 = nCoV.integrated@meta.data$BTM_Others

CD38p_PD1p_cells <- WhichCells(nCoV.integrated, expression = sct_CD38 > 0 & sct_PD1 > 0)
CD38p_PD1n_cells <- WhichCells(nCoV.integrated, expression = sct_CD38 > 0 & sct_PD1 == 0)
CD38n_PD1p_cells <- WhichCells(nCoV.integrated, expression = sct_CD38 == 0 & sct_PD1 > 0)
CD38n_PD1n_cells <- WhichCells(nCoV.integrated, expression = sct_CD38 == 0 & sct_PD1 == 0)






nCoV.integrated@meta.data[nCoV.integrated.CD38p.PD1p, 'CD38_PDCD1'] = 'CD38+ PD1+'
nCoV.integrated@meta.data[nCoV.integrated.CD38n.PD1p, 'CD38_PDCD1'] = 'CD38- PD1+'
nCoV.integrated@meta.data[nCoV.integrated.CD38p.PD1n, 'CD38_PDCD1'] = 'CD38+ PD1-'
nCoV.integrated@meta.data[nCoV.integrated.CD38n.PD1n, 'CD38_PDCD1'] = 'CD38- PD1-'

Idents(nCoV.integrated) = 'CD38_PDCD1'

nCoV.integrated.Dual_Genes_Markers <- FindAllMarkers(nCoV.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")

write.csv(nCoV.integrated.Dual_Genes_Markers, file = paste(fig_names, '_', 'Dual_Genes', "_AllClusterMarkers.csv", sep = ''))

nCoV.integrated.Singlets.markers = read.csv2(paste(fig_names, '_', 'Dual_Genes', "_AllClusterMarkers.csv", sep = ''),sep = ',', row.names = 1)
top5 <- nCoV.integrated.Singlets.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
feat_hm <- top5$gene

nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, assay = "SCT", features = feat_hm)

dpi = 200
png(file=paste(fig_names, paste('CD38_PDCD1_cluster_heatmap_','.png', sep = ''), sep = ''), width = dpi*16, height = dpi*12, units = "px",res = dpi,type='cairo')
print(DoHeatmap(nCoV.integrated, features = feat_hm, assay = "SCT", slot = 'scale.data') + NoLegend())
dev.off()





###### GO ANALYSIS ##############################################################################################################
# BiocManager::install("topGO")
topGOterms = function( fg.genes = NULL,
                       bg.genes = NULL,
                       organism = "Human", 
                       ontology.use = "BP",
                       stats.use = "fisher",
                       algorithm.use = "weight01",
                       topnodes.print=20,
                       num.char=100){
  
  if (is.null(fg.genes) | is.null(bg.genes)){
    stop("Error : Both gene lists are empty")
  }
  
  require(topGO)
  if (organism == "Mouse"){
    mapping.use = "org.Mm.eg.db"
    library(org.Mm.eg.db)
  } else if (organism == "Human"){
    mapping.use = "org.Hs.eg.db"
    library(org.Hs.eg.db)
  } else {
    stop("Error : Organisms other than mouse not supported currently")
  }
  
  n = length(bg.genes)
  geneList = integer(n)
  names(geneList) = bg.genes
  geneList[intersect(names(geneList), fg.genes)]=1
  print(paste0("Total ", length(geneList), " genes. ", sum(geneList), " genes in the foreground"))
  geneList = factor(geneList)
  
  if (ontology.use %in% c("BP", "CC", "MF")){
    print(paste0("Using Ontology : ", ontology.use))
  } else {
    stop("Error: Ontology not available. Should be one of BP, CC or MF")
  }
  # Make GO object
  GOdata <- new("topGOdata",
                description = "GOanalysis",
                ontology = ontology.use,
                allGenes = geneList,
                annot = annFUN.org,
                mapping = mapping.use,
                ID = "SYMBOL",
                nodeSize = 10)
  print(paste0("Using the ", stats.use, " statistic with the ", algorithm.use, " algorithm"))
  res.result <- runTest(GOdata, statistic = stats.use, algorithm = algorithm.use)
  to.return = list()
  to.return$GOdata = GOdata
  to.return$res.table <- GenTable(GOdata, pval = res.result, topNodes = topnodes.print, numChar = num.char)
  return(to.return)
}

###############################################################################################################################

GO_results = list()

clusters = as.character(sort(as.numeric(unique(object_AXOSpatial_seurat@meta.data$SCT_snn_res.0.8))-1))


clusters = as.character(unique(nCoV.integrated@meta.data$CD38_PDCD1))

Idents(nCoV.integrated) = "CD38_PDCD1"


de_markers = FindAllMarkers(nCoV.integrated, only.pos = TRUE)



de_markers$cluster = as.character(de_markers$cluster)





for (cluster_ in clusters) {
  
  de_markers_cluster = dplyr::filter(de_markers, cluster == cluster_)$gene
  
  if (length(de_markers_cluster) != 0) {
    
    GOterms = topGOterms(fg.genes = de_markers_cluster, bg.genes = rownames(nCoV.integrated), organism = "Human")
    
    GO_results[[cluster_]] = GOterms$res.table  %>% top_n(n = 10, wt = pval)
  }
  
}



GO_results_df = bind_rows(GO_results, .id = "id")
GO_results_df$pval = as.numeric(GO_results_df$pval)
GO_results_df = filter(GO_results_df, pval < 0.05)

# GO_results_df$id = as.numeric(GO_results_df$id)
# 
# 
# GO_results_df = GO_results_df[order(GO_results_df$id), ]

############################# GO FIGURE #######################################################################################
GO_results_df$pval = as.numeric(GO_results_df$pval)



data1 <- GO_results_df

colnames(GO_results_df)[1] = 'Cluster ID'

S1 <- ggplot(data1, aes(x = id, y = Term , color = -log(pval), size = log(Significant))) + 
  geom_point() + xlab("Cluster") + #scale_x_continuous(breaks = 1:12) +
  ylab(NULL) + 
  labs(colour = "-Log10(P.adjust)", size = "Log10(Number of Genes)") +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.text=element_text(size=10,colour ="black"),
        axis.title=element_text(size=14,face="bold",colour ="black"),
        legend.text=element_text(size=12),
        legend.title = element_text(colour="black", size=14, face="bold"),
        axis.line = element_line(colour = "black", size =1),
        panel.grid.major = element_line(colour = "lightgray", size = 0.5), 
        # panel.grid.minor = element_line(colour = "lightgray", size = 0.5),
        panel.background = element_blank())

# S1 = S1+scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab", limit = c(min(data1$pval), max(data1$pval)), n.breaks = 5)
S1+scale_size(range = c(0, length(unique(data1$id))))
S1


####################################################################################################################




DoHeatmap(nCoV.integrated)



case_when(test_score_vector >= 90 ~ 'A'
          ,test_score_vector >= 80 ~ 'B'
          ,test_score_vector >= 70 ~ 'C'
          ,test_score_vector >= 60 ~ 'D'
          ,TRUE ~ 'F'
)











library(SeuratDisk)
library(SeuratData)

reference_directory = 'C:/Users/Max/Dropbox/Work Saar/Seurat_References/'

setwd(reference_directory)

Tabula_Sapiens_Blood <- LoadH5Seurat("TS_Blood.h5seurat",assays = "RNA")

Tabula_Sapiens_Blood <- SCTransform(Tabula_Sapiens_Blood, verbose = TRUE) 
Tabula_Sapiens_Blood <- RunPCA(Tabula_Sapiens_Blood)
Tabula_Sapiens_Blood <- RunUMAP(Tabula_Sapiens_Blood, dims = 1:20)

Tabula_Sapiens_Blood$free_annotation

Tabula_Sapiens_Bone_Marrow <- LoadH5Seurat("TS_Bone_Marrow.h5seurat",assays = "RNA")

Tabula_Sapiens_Bone_Marrow <- SCTransform(Tabula_Sapiens_Bone_Marrow, verbose = TRUE) 
Tabula_Sapiens_Bone_Marrow <- RunPCA(Tabula_Sapiens_Bone_Marrow)
Tabula_Sapiens_Bone_Marrow <- RunUMAP(Tabula_Sapiens_Bone_Marrow, dims = 1:20)


sample.tmp.seurat = readRDS("reference_SCT_BM.rds")

cell_type_ident_TabulaS_PBMC = TRUE
cell_type_ident_TabulaS_BM = TRUE

if (cell_type_ident_TabulaS_PBMC == TRUE) {
  
  if (is.na(reference_directory)) {
    
    reference_SCT_TabulaS_PBMC = readRDS('reference_SCT_TabulaS_PBMC.rds')
    
  } else {
    
    reference_SCT_TabulaS_PBMC = readRDS(paste(reference_directory,'reference_SCT_TabulaS_PBMC.rds', sep = ''))
    
  }
  
  
}

if (cell_type_ident_TabulaS_BM == TRUE) {
  
  if (is.na(reference_directory)) {
    
    reference_SCT_TabulaS_BM = readRDS('reference_SCT_TabulaS_BM.rds')
    
  } else {
    
    reference_SCT_TabulaS_BM = readRDS(paste(reference_directory,'reference_SCT_TabulaS_BM.rds', sep = ''))
    
  }
  
  
}


if (cell_type_ident_TabulaS_PBMC == TRUE) {
  
  anchors <- FindTransferAnchors(
    reference = reference_SCT_TabulaS_PBMC,
    query = sample.tmp.seurat,
    normalization.method = "SCT",
    reference.reduction = "pca",
    dims = 1:50
  )
  
  sample.tmp.seurat <- MapQuery(
    anchorset = anchors,
    query = sample.tmp.seurat,
    reference = reference_SCT_TabulaS_PBMC,
    refdata = list(
      TabulaS_Blood_cell_ontology_class = "cell_ontology_class",
      TabulaS_Blood_free_annotation = "free_annotation"
    ),
    reference.reduction = "pca"
  )
  
  
}
if (cell_type_ident_TabulaS_BM == TRUE) {
  
  anchors <- FindTransferAnchors(
    reference = reference_SCT_TabulaS_BM,
    query = sample.tmp.seurat,
    normalization.method = "SCT",
    reference.reduction = "pca",
    dims = 1:50
  )
  
  sample.tmp.seurat <- MapQuery(
    anchorset = anchors,
    query = sample.tmp.seurat,
    reference = reference_SCT_TabulaS_BM,
    refdata = list(
      TabulaS_BM_cell_ontology_class = "cell_ontology_class",
      TabulaS_BM_free_annotation = "free_annotation"
    ),
    reference.reduction = "pca"
  )
  
  
}








unique(sample.tmp.seurat$predicted.Tabula_Sapiens_Blood_cell_ontology_class)
unique(sample.tmp.seurat$predicted.Tabula_Sapiens_Blood_free_annotation)


saveRDS(Tabula_Sapiens_Blood, "reference_SCT_TabulaS_PBMC.rds")
saveRDS(Tabula_Sapiens_Bone_Marrow, "reference_SCT_TabulaS_BM.rds")








healthy_BM_aml = readRDS('Healthy.rds')

healthy_BM = readRDS('WTA_projected.rds')

AML_BMs = readRDS('AMLs_Scano_projected.rds')

healthy_BM <- SCTransform(healthy_BM, verbose = TRUE) 
AML_BMs <- SCTransform(AML_BMs, verbose = TRUE) 

reference_SCT_BM <- RunSPCA(reference_SCT_BM)
reference_SCT_BM <- RunUMAP(reference_SCT_BM, dims = 1:20)
saveRDS(reference_SCT_BM, 'reference_SCT_BM.rds')
reference_SCT_BM_AML <- RunSPCA(reference_SCT_BM_AML)
reference_SCT_BM_AML <- RunUMAP(reference_SCT_BM_AML, dims = 1:20)
saveRDS(reference_SCT_BM_AML, 'reference_SCT_BM_AML.rds')

unique(colnames(reference_SCT_BM_AML@meta.data))

DimPlot(reference_SCT_BM_AML,group.by = "ct")

                
library(DropletUtils)                
write10xCounts(x = reference_SCT_BM@assays$RNA@counts,path = 'Healthy_BM.h5',type = 'HDF5')


test_h5 = Read10X_h5('Healthy_BM.h5')



assay(reference_SCT_BM)


saveRDS(healthy_BM, 'reference_SCT_BM.rds')
saveRDS(AML_BMs, 'reference_SCT_BM_AML.rds')


colnames(healthy_BM@meta.data)


colnames(AML_BMs@meta.data)


DimPlot(healthy_BM,label = TRUE, repel = T,group.by = 'Prediction_Healthy')
DimPlot(healthy_BM,label = TRUE, repel = T,group.by = 'Prediction_HCA')
DimPlot(healthy_BM,label = TRUE, repel = T,group.by = 'Prediction_BM3')

DimPlot(AML_BMs,label = TRUE, repel = T,group.by = 'Prediction_Ind')



markers1 = cbind(rownames(markers), markers)
colnames(markers1)[1] = 'Gene'
rownames(markers1) = c(1:nrow(markers1))
markers1 = subset(markers1, p_val_adj < 0.05)
markers1 = cbind(c(1:nrow(markers1)), markers1)
markers1$`c(1:nrow(markers1))` = sub('\\d+','_', markers1$`c(1:nrow(markers1))`)
colnames(markers1)[1] = 'Original out p_adj < 0.05 >>'


markers2 = markers1
markers2$pct.1 = NULL
markers2$pct.2 = NULL
markers2$p_val = NULL
markers2 = markers2[order(-markers2$avg_log2FC), ]

colnames(markers2)[1] = 'Sorted logFC asc >>'


markers3 = markers2
markers3 = markers3[order(markers3$avg_log2FC), ]
colnames(markers3)[1] = 'Sorted logFC desc >>'

markers4 = cbind(markers1,markers2,markers3)


dat$V1 <- sub("^", "chr", dat$V1 )

df[order(-df$var1), ]
# genes = as.data.frame(rownames(nCoV.integrated))
# 
# VlnPlot(object = nCoV.integrated.Patient, features = 'CSF2RA', split.by = 'Day', group.by = 'predicted.celltype.l1',log = TRUE, adjust = 2)
# 
# 
# marker_genes = c('MKI67','IL2','TNF','IFNG','CSF2','IL6','IL10','TGFB1','TGFB2','TGFB3','IL2RA','TNFRSF1A','IFNGR1','CSF2RA','IL6R','IL10RA','IL10RB','TGFBR1','TGFBR2')
# 
# marker_genes = c('SLC2A5','SLC2A8','ALDOB','ALDOC','KHK','AKR1B1','SORD')
# 
# install.packages('SeuratWrappers')
# library(SeuratWrappers)
# 
# split_ = c("Day")
# 
# FeaturePlot(object = nCoV.integrated.Patient, features = marker_genes, label = TRUE ,pt.size = 1.5, raster=FALSE, repel = TRUE, order = TRUE)
# 
# nCoV.integrated.Patient.10 = readRDS('nCoV.integrated.Patient 10.rds')
# nCoV.integrated.Patient.5 = nCoV.integrated.Patients[['Patient 5']]
# 
# nCoV.integrated.Patient.10 = RunALRA(nCoV.integrated.Patient.10, assay = 'SCT')
# # 
# nCoV.integrated.Patient = nCoV.integrated.Patient.10
# 
# name = 'Patient 10'
# 
# rm(nCoV.integrated.Patients)
# rm(nCoV.integrated.Patient)
# gc()
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

