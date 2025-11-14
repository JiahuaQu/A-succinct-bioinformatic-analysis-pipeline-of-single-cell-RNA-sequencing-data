library(clustree)
library(ggplotify)
library(data.table)
library(tidyverse)
library(Seurat)
library(hdf5r)
library(harmony)
library(EnhancedVolcano)


### Dataset 1:
h5_files <- "sample1/outs/filtered_feature_bc_matrix.h5"
h5_read <- Read10X_h5(h5_files)
h5_seurat <- CreateSeuratObject(h5_read,project="sample1")
saveRDS(h5_seurat, "h5_seurat.rds", compress = T)

mydata <- readRDS("h5_seurat.rds")

# Percentage of mitochondrial genes
mydata <- PercentageFeatureSet(mydata, pattern = "^MT-", col.name = "percent_mito",assay = "RNA")
# Percentage of ribosomal genes
mydata <- PercentageFeatureSet(mydata, pattern = "^RP[SL]", col.name = "percent_ribo",assay = "RNA")
saveRDS(mydata, "mydata-before_filter.rds", compress = T)

feats1 <- c("nFeature_RNA", "nCount_RNA")
feats2 <- c("percent_mito", "percent_ribo")

Vlnplot1 <- VlnPlot(mydata, group.by = "orig.ident", features = feats1, pt.size = 0.000001, ncol = 2)
ggsave(filename="Vlnplot1.pdf",plot=Vlnplot1, width=5, height = 5)
ggsave(filename="Vlnplot1.png",plot=Vlnplot1, width=5, height = 5, units = "in", dpi = 300)

Vlnplot2 <- VlnPlot(mydata, group.by = "orig.ident", features = feats2, pt.size = 0.000001, ncol = 2) 
ggsave(filename="Vlnplot2.pdf",plot=Vlnplot2, width=5, height = 5)
ggsave(filename="Vlnplot2.png",plot=Vlnplot2, width=5, height = 5, units = "in", dpi = 300)


qc_cutoffs <- list(
  nCount_RNA   = quantile(mydata$nCount_RNA,   probs = c(0.10, 0.90)),
  nFeature_RNA = quantile(mydata$nFeature_RNA, probs = c(0.10, 0.90)),
  percent_mito = quantile(mydata$percent_mito, probs = 0.90),  # only upper cutoff
  percent_ribo = quantile(mydata$percent_ribo, probs = 0.90)   # only upper cutoff
)

qc_cutoffs

mydata <- subset(mydata, subset = nCount_RNA >= 4417.9  & 
                   nCount_RNA <= 20937.6  &
                   nFeature_RNA >= 1956.0  & 
                   nFeature_RNA <= 5014.4  &
                   percent_mito <= 11.18076  & 
                   percent_ribo <= 29.44844)

Vlnplot3 <- VlnPlot(mydata, group.by = "orig.ident", features = feats1, pt.size = 0.000001, ncol = 2) 
ggsave(filename="Vlnplot3.pdf",plot=Vlnplot3, width=5, height = 5)
ggsave(filename="Vlnplot3.png",plot=Vlnplot3, width=5, height = 5, units = "in", dpi = 300)

Vlnplot4 <- VlnPlot(mydata, group.by = "orig.ident", features = feats2, pt.size = 0.000001, ncol = 2) 
ggsave(filename="Vlnplot4.pdf",plot=Vlnplot4, width=5, height = 5)
ggsave(filename="Vlnplot4.png",plot=Vlnplot4, width=5, height = 5, units = "in", dpi = 300)

saveRDS(mydata, "mydata-after_filter.rds", compress = T)


### Dataset 2:
h5_files <- "sample2/outs/filtered_feature_bc_matrix.h5"
h5_read <- Read10X_h5(h5_files)
h5_seurat_2 <- CreateSeuratObject(h5_read,project="sample2")
saveRDS(h5_seurat_2, "h5_seurat_2.rds", compress = T)

# Filter
mydata <- h5_seurat_2
rm(h5_seurat_2)
# Percentage of mitochondrial genes
mydata <- PercentageFeatureSet(mydata, pattern = "^MT-", col.name = "percent_mito",assay = "RNA")
# Percentage of ribosomal genes
mydata <- PercentageFeatureSet(mydata, pattern = "^RP[SL]", col.name = "percent_ribo",assay = "RNA")
saveRDS(mydata, "mydata-before_filter.rds", compress = T)

feats1 <- c("nFeature_RNA", "nCount_RNA")
feats2 <- c("percent_mito", "percent_ribo")

Vlnplot5 <- VlnPlot(mydata, group.by = "orig.ident", features = feats1, pt.size = 0.000001, ncol = 2)
ggsave(filename="Vlnplot5.pdf",plot=Vlnplot5, width=5, height = 5)
ggsave(filename="Vlnplot5.png",plot=Vlnplot5, width=5, height = 5, units = "in", dpi = 300)

Vlnplot6 <- VlnPlot(mydata, group.by = "orig.ident", features = feats2, pt.size = 0.000001, ncol = 2) 
ggsave(filename="Vlnplot6.pdf",plot=Vlnplot6, width=5, height = 5)
ggsave(filename="Vlnplot6.png",plot=Vlnplot6, width=5, height = 5, units = "in", dpi = 300)

qc_cutoffs <- list(
  nCount_RNA   = quantile(mydata$nCount_RNA,   probs = c(0.10, 0.90)),
  nFeature_RNA = quantile(mydata$nFeature_RNA, probs = c(0.10, 0.90)),
  percent_mito = quantile(mydata$percent_mito, probs = 0.90),  # only upper cutoff
  percent_ribo = quantile(mydata$percent_ribo, probs = 0.90)   # only upper cutoff
)

qc_cutoffs
# $nCount_RNA
# 10%     90% 
#   7529.0 32994.6 
# 
# $nFeature_RNA
# 10%    90% 
#   2543.7 5935.0 
# 
# $percent_mito
# 90% 
# 4.363675 
# 
# $percent_ribo
# 90% 
# 38.06148 

mydata <- subset(mydata, subset = nCount_RNA >= 7529  & 
                   nCount_RNA <= 32994.6  &
                   nFeature_RNA >= 2543.7  & 
                   nFeature_RNA <= 5935.0  &
                   percent_mito <= 4.363675  & 
                   percent_ribo <= 38.06148)

Vlnplot7 <- VlnPlot(mydata, group.by = "orig.ident", features = feats1, pt.size = 0.000001, ncol = 2) 
ggsave(filename="Vlnplot7.pdf",plot=Vlnplot7, width=5, height = 5)
ggsave(filename="Vlnplot7.png",plot=Vlnplot7, width=5, height = 5, units = "in", dpi = 300)

Vlnplot8 <- VlnPlot(mydata, group.by = "orig.ident", features = feats2, pt.size = 0.000001, ncol = 2) 
ggsave(filename="Vlnplot8.pdf",plot=Vlnplot8, width=5, height = 5)
ggsave(filename="Vlnplot8.png",plot=Vlnplot8, width=5, height = 5, units = "in", dpi = 300)

saveRDS(mydata, "mydata-after_filter_2.rds", compress = T)
#mydata_2 <- mydata

### Merge
mydata_1 <- readRDS("mydata-after_filter.rds")   # The first dataset
mydata_2 <- readRDS("mydata-after_filter_2.rds")   # The second dataset
 
### Merge
cell_ids <- c("sample1", "sample2")
merged <- merge(
  x = mydata_1, 
  y = mydata_2, 
  add.cell.ids = cell_ids, 
  project = "Merge"
)

# Update orig.ident
unique(merged@meta.data$orig.ident)
saveRDS(merged, "mydata_merge.rds", compress = T)


### Join layers
mydata <- readRDS("mydata_merge.rds")

Layers(mydata[["RNA"]])
# "counts.sample1" "counts.sample2"

# JoinLayers
mydata[["RNA"]] <- JoinLayers(mydata[["RNA"]])

Layers(mydata[["RNA"]])
# "counts"

saveRDS(mydata, "mydata-merge_joinlayers.rds", compress = T)


### Downstream
mydata <- readRDS("mydata-merge_joinlayers.rds")

n.genes <- nrow(mydata)
feature.genes <- rownames(mydata)

mydata <- NormalizeData(object = mydata,normalization.method = "LogNormalize",scale.factor = 10000,margin = 1, verbose = FALSE)
mydata <- FindVariableFeatures(object = mydata, selection.method = "vst", nfeatures = n.genes)
mydata <- ScaleData(object = mydata, features = feature.genes)
mydata <- RunPCA(object = mydata, verbose = FALSE)
p1 <- ElbowPlot(object = mydata, ndims = 50,reduction="pca") 
ggsave(filename="ElbowPlot.pdf",plot=p1, width=5, height = 5)
ggsave(filename="ElbowPlot.png",plot=p1, width=5, height = 5, units = "in", dpi = 300)
p2 <- DimPlot(object = mydata, reduction = "pca", group.by="orig.ident")
ggsave(filename="PCA_by_orig.ident.pdf",plot=p2, width=5, height = 5)
ggsave(filename="PCA_by_orig.ident.png",plot=p2, width=5, height = 5, units = "in", dpi = 300)

saveRDS(mydata, "mydata-processed_after_pca.rds", compress = T)


###
mydata <- readRDS("mydata-processed_after_pca.rds")
### Harmony
library(harmony)
mydata <- RunHarmony(object = mydata, group.by.vars = "orig.ident", reduction = "pca",reduction.save = "harmony")

p3 <- ElbowPlot(object = mydata, ndims = 50, reduction="harmony") 
ggsave(filename="ElbowPlot-harmony.pdf",plot=p3, width=5, height = 5)
ggsave(filename="ElbowPlot-harmony.png",plot=p3, width=5, height = 5, units = "in", dpi = 300)
p4 <- DimPlot(object = mydata, reduction = "harmony", group.by="orig.ident")
ggsave(filename="PCA_by_orig.ident-harmony.pdf",plot=p4, width=5, height = 5)
ggsave(filename="PCA_by_orig.ident-harmony.png",plot=p4, width=5, height = 5, units = "in", dpi = 300)

saveRDS(mydata, "mydata-processed_after_pca-harmony.rds", compress = T)


###
mydata <- readRDS("mydata-processed_after_pca-harmony.rds")

harmony_dim <- 1:15

mydata <- RunUMAP(mydata, reduction = "harmony", dims = harmony_dim)

mydata <- RunTSNE(mydata, reduction = "harmony", dims = harmony_dim)

saveRDS(mydata, "mydata-processed_after_pca-harmony_UMAP_tSNE.rds", compress = T)


### clustering
#mydata <- readRDS("mydata-processed_after_pca-harmony_UMAP_tSNE.rds")
DefaultAssay(mydata) <- "RNA"
mydata <- FindNeighbors(object = mydata, reduction = "harmony", dims = harmony_dim)

res <- seq(0.1, 1, 0.1)
mydata <- FindClusters(object = mydata, resolution = res, verbose = FALSE)

p5 <- clustree(mydata)
p5 <- as.ggplot(p5)
ggsave(filename="clustree.pdf",plot=p5, width=5, height = 6)
ggsave(filename="clustree.png",plot=p5, width=5, height = 6, units = "in", dpi = 300)

# save after clustering
saveRDS(mydata, "mydata-processed_after_clustering.rds", compress = T)


##
mydata <- readRDS("mydata-processed_after_clustering.rds")
### Each resolution
subDir <- "resolution/"
if (!file.exists(subDir)){
  dir.create(subDir)
} 

colnames(mydata@meta.data)
res <- seq(0.1, 1, 0.1)

for (i in res) {
  ### Plot and color clusters in UMAP
  group_i <- paste0("RNA_snn_res.", i)
  Idents(object = mydata) <- group_i
  
  p6 <- DimPlot(object = mydata, reduction = "umap", group.by = group_i, label = T, repel = T) 
  ggsave(filename=paste0(subDir,group_i,'.pdf'), plot=p6, width = 7, height = 6) 
  ggsave(filename=paste0(subDir,group_i,'.png'), plot=p6, width = 7, height = 6, units = "in", dpi = 300) 
  
  ### biomarkers
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  markers <- FindAllMarkers(object = mydata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  fwrite(markers, paste0(subDir,group_i,'-markers.csv'))
  
  markers.top5 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  fwrite(markers.top5, paste0(subDir,group_i,'-markers.top5.csv'))
  
  features.top5 <- unique(markers.top5$gene)
  height_i <- (length(features.top5) + 2) * 0.15
  p7 <- DotPlot(object = mydata, features = features.top5, group.by=group_i) + coord_flip()
  ggsave(filename=paste0(subDir,group_i,'-markers.top5.pdf'), plot=p7, width = 10, height = height_i)
  ggsave(filename=paste0(subDir,group_i,'-markers.top5.png'), plot=p7, width = 10, height = height_i, units = "in", dpi = 300)
  
  ### heatmap
  p8 <- DoHeatmap(object = mydata, features = features.top5, group.by=group_i, label = F)
  ggsave(filename=paste0(subDir,group_i,'-heatmap.top5.pdf'), plot=p8, width = 10, height = height_i)
  ggsave(filename=paste0(subDir,group_i,'-heatmap.top5.png'), plot=p8, width = 10, height = height_i, units = "in", dpi = 300)
}

### Select resolution = 0.1
##
mydata <- readRDS("mydata-processed_after_clustering.rds")

group_i <- "RNA_snn_res.0.1"
unique(mydata@meta.data$RNA_snn_res.0.1)
mydata@meta.data$RNA_snn_res.0.1 <- factor(mydata@meta.data$RNA_snn_res.0.1,levels=0:8)

markers.top5 <- fread("resolution/RNA_snn_res.0.1-markers.top5.csv")
features.top5 <- markers.top5$gene
NROW(unique(features.top5))   # 45
example_9 <- features.top5[c(1,6,11,16,21,26,31,36,41)]
example_9
# "GZMK" "S100A8" "ANKRD55" "FGFBP2" "CD8B"           
# "VPREB3" "CTLA4" "TOX" "ENSG00000290592"


### Plots
Res <- "res.0.1/"
if (!file.exists(Res)){
  dir.create(Res)
}

p9 <- VlnPlot(object = mydata, group.by = group_i, features = example_9, pt.size = 0.000001, ncol = 3)
ggsave(filename=paste0(Res,"Vlnplot-markers.pdf"),plot=p9, width=15, height = 15)
ggsave(filename=paste0(Res,"Vlnplot-markers.png"),plot=p9, width=15, height = 15, units = "in", dpi = 300)

p10 <- FeaturePlot(object = mydata, reduction = "umap", features = example_9, pt.size = 0.000001, ncol = 3)
ggsave(filename=paste0(Res,"FeaturePlot-markers.pdf"),plot=p10, width=15, height = 15)
ggsave(filename=paste0(Res,"FeaturePlot-markers.png"),plot=p10, width=15, height = 15, units = "in", dpi = 300)


### AddModuleScores
modules <- list(
  module0=features.top5[1:5],
  module1=features.top5[6:10],
  module2=features.top5[11:15],
  module3=features.top5[16:20],
  module4=features.top5[21:25],
  module5=features.top5[26:30],
  module6=features.top5[31:35],
  module7=features.top5[36:40],
  module8=features.top5[41:45]
)

mydata <- Seurat::AddModuleScore(object = mydata,features = modules, name = "ModuleScore", ctrl = 100,seed = 123)

colnames(mydata@meta.data)
# New columns
# "ModuleScore1" "ModuleScore2" "ModuleScore3" "ModuleScore4" "ModuleScore5"
# "ModuleScore6"    "ModuleScore7"    "ModuleScore8"    "ModuleScore9"

# Rename
mydata@meta.data <- mydata@meta.data %>%
  rename(module0 = ModuleScore1, 
         module1 = ModuleScore2,
         module2 = ModuleScore3,
         module3 = ModuleScore4,
         module4 = ModuleScore5,
         module5 = ModuleScore6, 
         module6 = ModuleScore7,
         module7 = ModuleScore8,
         module8 = ModuleScore9)

colnames(mydata@meta.data)
# New name
# "module0" "module1" "module2" "module3" "module4"
# "module5"         "module6"         "module7"         "module8"

module_9 <- c("module0", "module1", "module2", "module3", "module4",
              "module5", "module6", "module7", "module8")

p11 <- VlnPlot(object = mydata, group.by = group_i, features = module_9, pt.size = 0.000001, ncol = 3)
ggsave(filename=paste0(Res,"Vlnplot-modules.pdf"),plot=p11, width=15, height = 15)
ggsave(filename=paste0(Res,"Vlnplot-modules.png"),plot=p11, width=15, height = 15, units = "in", dpi = 300)

p12 <- FeaturePlot(object = mydata, reduction = "umap", features = module_9, pt.size = 0.000001, ncol = 3)
ggsave(filename=paste0(Res,"FeaturePlot-modules.pdf"),plot=p12, width=15, height = 15)
ggsave(filename=paste0(Res,"FeaturePlot-modules.png"),plot=p12, width=15, height = 15, units = "in", dpi = 300)

# Assign cell annotation
mydata@meta.data <- mydata@meta.data %>%
  mutate(celltype = recode(RNA_snn_res.0.1,
                           `0` = "cluster0",
                           `1` = "cluster1",
                           `2` = "cluster2",
                           `3` = "cluster3",
                           `4` = "cluster4",
                           `5` = "cluster5",
                           `6` = "cluster6",
                           `7` = "cluster7",
                           `8` = "cluster8"
  ))
head(mydata@meta.data$celltype)

mydata@meta.data$celltype <- factor(mydata@meta.data$celltype,
                                    levels=c("cluster0","cluster1","cluster2",
                                             "cluster3","cluster4","cluster5",
                                             "cluster6","cluster7","cluster8"))  
head(mydata@meta.data$celltype)

saveRDS(mydata, "mydata-processed_after_score_celltype.rds", compress = T)


### DEG between two cell types or among all cell types
mydata <- readRDS("mydata-processed_after_score_celltype.rds.rds")

# Compare "cluster0" and "cluster1" as an example
DEG <- FindMarkers(
  object = mydata,
  ident.1 = "cluster0",  
  ident.2 = "cluster1", 
  group.by = "celltype", 
  test.use = "MAST" 
)

head(DEG,2)
DEG$gene <- rownames(DEG)
fwrite(DEG, paste0(Res,'DEG-cluster0_cluster1.csv'))
saveRDS(DEG, paste0(Res,'DEG-cluster0_cluster1.rds'), compress=T)

# Volcano plot
#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

p13 <- EnhancedVolcano(
  DEG,
  lab = NA,
  x = 'avg_log2FC',
  y = 'p_val_adj',
  ylab = bquote( ~ -Log[10] ~ "adjusted p-value"),
  pCutoff = 0.05,
  FCcutoff = 1.333,
  title = "cluster0 vs. cluster1",
  subtitle = NULL,
  legendLabels = c(
    'adjusted p-value >= 0.05 &\nabsolute Log2FC < 1.333',
    'adjusted p-value >= 0.05 &\nabsolute Log2FC >= 1.333',
    'adjusted p-value < 0.05 &\nabsolute Log2FC < 1.333',
    'adjusted p-value < 0.05 &\nabsolute Log2FC >= 1.333'
  ),
  legendPosition = 'right'
  )
p13

significant_down <- sum(DEG$avg_log2FC <= -1.333 & DEG$p_val_adj < 0.05, na.rm = TRUE)
significant_up <- sum(DEG$avg_log2FC >= 1.333 & DEG$p_val_adj < 0.05, na.rm = TRUE)

p14 <- p13 +
  annotate("text", x = -10, y = max(-log10(DEG$p_val_adj), na.rm = TRUE), 
           label = paste("Down:\n", significant_down), 
           hjust = 0, size = 5) +
  annotate("text", x = 7.5, y = max(-log10(DEG$p_val_adj), na.rm = TRUE), 
           label = paste("Up:\n", significant_up), 
           hjust = 0, size = 5) +
  guides(
    color = guide_legend(
      keyheight = unit(1.5, "cm"),  
      override.aes = list(size = 7) 
    )
    )
p14
ggsave(filename=paste0(Res,"volcano-cluster0_cluster1.pdf"),plot=p14, width=8, height = 5)
ggsave(filename=paste0(Res,"volcano-cluster0_cluster1.png"),plot=p14, width=8, height = 5, units = "in", dpi = 300)


# Export the package information 
sink("sessionInfo.txt")
sessionInfo()
sink()




