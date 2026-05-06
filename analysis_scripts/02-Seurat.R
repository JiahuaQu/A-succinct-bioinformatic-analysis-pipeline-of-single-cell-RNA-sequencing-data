
### Install R packages
# https://cran.r-project.org/web/packages/clustree/readme/README.html
install.packages("clustree")

# https://github.com/GuangchuangYu/ggplotify
install.packages("ggplotify")

# https://github.com/Rdatatable/data.table
install.packages("data.table")

# https://cran.r-project.org/web/packages/tidyverse/readme/README.html
install.packages("tidyverse")

# https://satijalab.org/seurat/articles/install_v5
install.packages("Seurat")

# https://github.com/hhoeflin/hdf5r
install.packages("hdf5r")

# https://cran.r-project.org/web/packages/harmony/vignettes/quickstart.html
install.packages('harmony')

install.packages("patchwork")

# https://github.com/kevinblighe/EnhancedVolcano
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')

# https://github.com/RGLab/MAST
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install("MAST")

### Load in packages
library(clustree)
#library(clustree,lib.loc="/home/jqu/R/workbench/4.3.2")
library(ggplotify)
library(data.table)
library(tidyverse)
library(Seurat)
library(hdf5r)
library(harmony)
library(EnhancedVolcano)
library(MAST)
library(patchwork)


### Read in and preprocess datasets
### Dataset 1:
# Path to folder/directory
# /mnt/iminfo/jqu/project/scRNAseq_pipeline/count/outs/filtered_feature_bc_matrix.h5

h5_files <- "../count/outs/filtered_feature_bc_matrix.h5"

h5_read <- Read10X_h5(h5_files)

# Create Seurat object
h5_seurat <- CreateSeuratObject(h5_read,project="sample1")
saveRDS(h5_seurat, "h5_seurat.rds")

h5_seurat

meta <- h5_seurat@meta.data

head(h5_seurat@meta.data)

# Save intermediate object
mydata <- readRDS("h5_seurat.rds")
#mydata <- readRDS("../R/h5_seurat.rds")

mydata <- h5_seurat

# For QC
# Percentage of mitochondrial genes
mydata <- PercentageFeatureSet(mydata, pattern = "^MT-", col.name = "percent_mito",assay = "RNA")
# Percentage of ribosomal genes
mydata <- PercentageFeatureSet(mydata, pattern = "^RP[SL]", col.name = "percent_ribo",assay = "RNA")

head(mydata@meta.data)

saveRDS(mydata, "mydata-before_filter.rds")

# mydata <- readRDS("mydata-before_filter.rds")

feats1 <- c("nFeature_RNA", "nCount_RNA")
feats2 <- c("percent_mito", "percent_ribo")

# Observe before filter
Vlnplot1 <- VlnPlot(mydata, group.by = "orig.ident", features = feats1, pt.size = 0.000001, ncol = 2)
Vlnplot1
ggsave(filename="Vlnplot1.pdf",plot=Vlnplot1, width=5, height = 5)
ggsave(filename="Vlnplot1.png",plot=Vlnplot1, width=5, height = 5, units = "in", dpi = 300)

Vlnplot2 <- VlnPlot(mydata, group.by = "orig.ident", features = feats2, pt.size = 0.000001, ncol = 2) 
Vlnplot2
ggsave(filename="Vlnplot2.pdf",plot=Vlnplot2, width=5, height = 5)
ggsave(filename="Vlnplot2.png",plot=Vlnplot2, width=5, height = 5, units = "in", dpi = 300)

qc_cutoffs <- list(
  nCount_RNA   = quantile(mydata$nCount_RNA,   probs = c(0.01, 0.99)),
  nFeature_RNA = quantile(mydata$nFeature_RNA, probs = c(0.01, 0.99)),
  percent_mito = quantile(mydata$percent_mito, probs = 0.99),  # only upper cutoff
  percent_ribo = quantile(mydata$percent_ribo, probs = 0.99)   # only upper cutoff
)

qc_cutoffs
# $nCount_RNA
# 1%      99% 
#   670.58 33033.83 
# 
# $nFeature_RNA
# 1%     99% 
#   31.29 6151.65 
# 
# $percent_mito
# 99% 
# 97.10126 
# 
# $percent_ribo
# 99% 
# 35.21121 

### Plots with dotted lines
make_vln_with_cutoff <- function(object, feature, group.by = "orig.ident") {
  
  vals <- object[[feature]][, 1]
  qs <- quantile(vals, probs = c(0.01, 0.99), na.rm = TRUE)
  
  VlnPlot(
    object,
    group.by = group.by,
    features = feature,
    pt.size = 0.000001
  ) +
    geom_hline(yintercept = qs[1], linetype = "dotted", linewidth = 0.5) +
    geom_hline(yintercept = qs[2], linetype = "dotted", linewidth = 0.5) +
    ggtitle(feature) +
    NoLegend()
}

plots1 <- lapply(feats1, function(x) {
  make_vln_with_cutoff(mydata, feature = x, group.by = "orig.ident")
})

Vlnplot1_2 <- wrap_plots(plots1, ncol = 2)

Vlnplot1_2

ggsave(filename="Vlnplot1_2.pdf",plot=Vlnplot1_2, width=5, height = 5)
ggsave(filename="Vlnplot1_2.png",plot=Vlnplot1_2, width=5, height = 5, units = "in", dpi = 300)

plots2 <- lapply(feats2, function(x) {
  make_vln_with_cutoff(mydata, feature = x, group.by = "orig.ident")
})

Vlnplot2_2 <- wrap_plots(plots2, ncol = 2)

Vlnplot2_2

ggsave(filename="Vlnplot2_2.pdf",plot=Vlnplot2_2, width=5, height = 5)
ggsave(filename="Vlnplot2_2.png",plot=Vlnplot2_2, width=5, height = 5, units = "in", dpi = 300)




q1 <- unname(qc_cutoffs$nCount_RNA[1]);   q2 <- unname(qc_cutoffs$nCount_RNA[2])
q3 <- unname(qc_cutoffs$nFeature_RNA[1]); q4 <- unname(qc_cutoffs$nFeature_RNA[2])

q5 <- min(as.numeric(qc_cutoffs$percent_mito), 10)
q6 <- min(as.numeric(qc_cutoffs$percent_ribo), 30)

# Filter
mydata <- subset(
  mydata,
  subset =
    nCount_RNA   >= q1 &
    nCount_RNA   <= q2 &
    nFeature_RNA >= q3 &
    nFeature_RNA <= q4 &
    percent_mito <= q5 &
    percent_ribo <= q6
)

head(mydata@meta.data)

# Observe after filter
Vlnplot3 <- VlnPlot(mydata, group.by = "orig.ident", features = feats1, pt.size = 0.000001, ncol = 2) 
Vlnplot3
ggsave(filename="Vlnplot3.pdf",plot=Vlnplot3, width=5, height = 5)
ggsave(filename="Vlnplot3.png",plot=Vlnplot3, width=5, height = 5, units = "in", dpi = 300)

Vlnplot4 <- VlnPlot(mydata, group.by = "orig.ident", features = feats2, pt.size = 0.000001, ncol = 2) 
Vlnplot4
ggsave(filename="Vlnplot4.pdf",plot=Vlnplot4, width=5, height = 5)
ggsave(filename="Vlnplot4.png",plot=Vlnplot4, width=5, height = 5, units = "in", dpi = 300)

saveRDS(mydata, "mydata-after_filter.rds")


### Dataset 2:
h5_files <- "../count-2-2/outs/filtered_feature_bc_matrix.h5"
h5_read <- Read10X_h5(h5_files)
h5_seurat_2 <- CreateSeuratObject(h5_read,project="sample2")
saveRDS(h5_seurat_2, "h5_seurat_2.rds")

# Filter
mydata <- h5_seurat_2
rm(h5_seurat_2)

#mydata <- readRDS("../R-2/h5_seurat_2.rds")

# Percentage of mitochondrial genes
mydata <- PercentageFeatureSet(mydata, pattern = "^MT-", col.name = "percent_mito",assay = "RNA")
# Percentage of ribosomal genes
mydata <- PercentageFeatureSet(mydata, pattern = "^RP[SL]", col.name = "percent_ribo",assay = "RNA")
saveRDS(mydata, "mydata-before_filter_2.rds")

feats1 <- c("nFeature_RNA", "nCount_RNA")
feats2 <- c("percent_mito", "percent_ribo")

Vlnplot5 <- VlnPlot(mydata, group.by = "orig.ident", features = feats1, pt.size = 0.000001, ncol = 2)
ggsave(filename="Vlnplot5.pdf",plot=Vlnplot5, width=5, height = 5)
ggsave(filename="Vlnplot5.png",plot=Vlnplot5, width=5, height = 5, units = "in", dpi = 300)

Vlnplot6 <- VlnPlot(mydata, group.by = "orig.ident", features = feats2, pt.size = 0.000001, ncol = 2) 
ggsave(filename="Vlnplot6.pdf",plot=Vlnplot6, width=5, height = 5)
ggsave(filename="Vlnplot6.png",plot=Vlnplot6, width=5, height = 5, units = "in", dpi = 300)

qc_cutoffs <- list(
  nCount_RNA   = quantile(mydata$nCount_RNA,   probs = c(0.01, 0.99)),
  nFeature_RNA = quantile(mydata$nFeature_RNA, probs = c(0.01, 0.99)),
  percent_mito = quantile(mydata$percent_mito, probs = 0.99),  # only upper cutoff
  percent_ribo = quantile(mydata$percent_ribo, probs = 0.99)   # only upper cutoff
)

qc_cutoffs
# $nCount_RNA
# 1%      99% 
#   907.14 66181.52 
# 
# $nFeature_RNA
# 1%     99% 
#   400.28 7919.59 
# 
# $percent_mito
# 99% 
# 35.44623 
# 
# $percent_ribo
# 99% 
# 42.73233 

q1 <- unname(qc_cutoffs$nCount_RNA[1]);   q2 <- unname(qc_cutoffs$nCount_RNA[2])
q3 <- unname(qc_cutoffs$nFeature_RNA[1]); q4 <- unname(qc_cutoffs$nFeature_RNA[2])

q5 <- min(as.numeric(qc_cutoffs$percent_mito), 10)
q6 <- min(as.numeric(qc_cutoffs$percent_ribo), 30)

# Filter
mydata <- subset(
  mydata,
  subset =
    nCount_RNA   >= q1 &
    nCount_RNA   <= q2 &
    nFeature_RNA >= q3 &
    nFeature_RNA <= q4 &
    percent_mito <= q5 &
    percent_ribo <= q6
)

Vlnplot7 <- VlnPlot(mydata, group.by = "orig.ident", features = feats1, pt.size = 0.000001, ncol = 2) 
ggsave(filename="Vlnplot7.pdf",plot=Vlnplot7, width=5, height = 5)
ggsave(filename="Vlnplot7.png",plot=Vlnplot7, width=5, height = 5, units = "in", dpi = 300)

Vlnplot8 <- VlnPlot(mydata, group.by = "orig.ident", features = feats2, pt.size = 0.000001, ncol = 2) 
ggsave(filename="Vlnplot8.pdf",plot=Vlnplot8, width=5, height = 5)
ggsave(filename="Vlnplot8.png",plot=Vlnplot8, width=5, height = 5, units = "in", dpi = 300)

saveRDS(mydata, "mydata-after_filter_2.rds")
#mydata_2 <- mydata

### Merge multiple datasets
# Read in each dataset
mydata_1 <- readRDS("mydata-after_filter.rds")   # The first dataset
mydata_2 <- readRDS("mydata-after_filter_2.rds")   # The second dataset
 
# Merge
cell_ids <- c("sample1", "sample2")
merged <- merge(
  x = mydata_1, 
  y = mydata_2, 
  add.cell.ids = cell_ids, 
  project = "Merge"
)

head(merged@meta.data)

# Update orig.ident
unique(merged@meta.data$orig.ident)
saveRDS(merged, "mydata_merge.rds")


### Join layers from different objects
# mydata <- merged
# rm(merged)
mydata <- readRDS("mydata_merge.rds")

mydata <- merged

Layers(mydata[["RNA"]])
# "counts.sample1" "counts.sample2"

# JoinLayers
mydata[["RNA"]] <- JoinLayers(mydata[["RNA"]])

Layers(mydata[["RNA"]])
# "counts"

saveRDS(mydata, "mydata-merge_joinlayers.rds")


### Downstream
mydata <- readRDS("mydata-merge_joinlayers.rds")

mydata <- NormalizeData(object = mydata,normalization.method = "LogNormalize",scale.factor = 10000,margin = 1, verbose = FALSE)
mydata <- FindVariableFeatures(object = mydata, selection.method = "vst")
mydata <- ScaleData(object = mydata)
mydata <- RunPCA(object = mydata, verbose = FALSE)
p1 <- ElbowPlot(object = mydata, ndims = 50,reduction="pca") 
p1
ggsave(filename="ElbowPlot.pdf",plot=p1, width=5, height = 5)
ggsave(filename="ElbowPlot.png",plot=p1, width=5, height = 5, units = "in", dpi = 300)

p2 <- DimPlot(object = mydata, reduction = "pca", group.by="orig.ident")
p2
ggsave(filename="PCA_by_orig.ident.pdf",plot=p2, width=5, height = 5)
ggsave(filename="PCA_by_orig.ident.png",plot=p2, width=5, height = 5, units = "in", dpi = 300)

saveRDS(mydata, "mydata-processed_after_pca.rds")

### Harmony to remove batch effect
mydata <- readRDS("mydata-processed_after_pca.rds")
library(harmony)
mydata <- RunHarmony(object = mydata, group.by.vars = "orig.ident", reduction = "pca",reduction.save = "harmony")

p3 <- ElbowPlot(object = mydata, ndims = 50, reduction="harmony") 
p3
ggsave(filename="ElbowPlot-harmony.pdf",plot=p3, width=5, height = 5)
ggsave(filename="ElbowPlot-harmony.png",plot=p3, width=5, height = 5, units = "in", dpi = 300)
p4 <- DimPlot(object = mydata, reduction = "harmony", group.by="orig.ident")
p4
ggsave(filename="PCA_by_orig.ident-harmony.pdf",plot=p4, width=5, height = 5)
ggsave(filename="PCA_by_orig.ident-harmony.png",plot=p4, width=5, height = 5, units = "in", dpi = 300)

p4 <- DimPlot(object = mydata, reduction = "harmony", group.by="orig.ident", shuffle = T)
p4


saveRDS(mydata, "mydata-processed_after_pca-harmony.rds")


### UMAP and t-SNE for visualization
mydata <- readRDS("mydata-processed_after_pca-harmony.rds")

harmony_dim <- 1:20

mydata <- RunUMAP(mydata, reduction = "harmony", dims = harmony_dim)
p4_2 <- DimPlot(object = mydata, reduction = "umap", group.by="orig.ident", shuffle=T)
p4_2

mydata <- RunTSNE(mydata, reduction = "harmony", dims = harmony_dim)

saveRDS(mydata, "mydata-processed_after_pca-harmony_UMAP_tSNE.rds")


### clustering
#mydata <- readRDS("mydata-processed_after_pca-harmony_UMAP_tSNE.rds")
DefaultAssay(mydata) <- "RNA"
mydata <- FindNeighbors(object = mydata, reduction = "harmony", dims = harmony_dim)

res <- seq(0.1, 1, 0.1)
mydata <- FindClusters(object = mydata, resolution = res, verbose = FALSE)

head(mydata@meta.data)

# Loop through each resolution and set factor levels
for (r in res) {
  colname <- paste0("RNA_snn_res.", r)
  if (colname %in% colnames(mydata@meta.data)) {
    vals <- mydata[[colname]][, 1]
    levs <- 0:(NROW(unique(vals))-1)
    mydata[[colname]] <- factor(vals, levels = levs)
  }
}

p5 <- clustree(mydata)
p5 <- as.ggplot(p5)
ggsave(filename="clustree.pdf",plot=p5, width=6, height = 8)
ggsave(filename="clustree.png",plot=p5, width=6, height = 8, units = "in", dpi = 300)

# save after clustering
saveRDS(mydata, "mydata-processed_after_clustering.rds")


### Each resolution
mydata <- readRDS("mydata-processed_after_clustering.rds")

subDir <- "resolution/"
if (!file.exists(subDir)){
  dir.create(subDir)
} 

# res <- seq(0.1, 1, 0.1)
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

### Select resolution = 0.2 as an example for downstream analysis
#mydata <- readRDS("mydata-processed_after_clustering.rds")

group_i <- "RNA_snn_res.0.2"

markers.top5 <- fread("resolution/RNA_snn_res.0.2-markers.top5.csv")
features.top5 <- markers.top5$gene
NROW(unique(features.top5))   # 
example <- features.top5[seq(1, length(features.top5), by = 5)]
example


### Plots at res=0.2 as an example
Res <- "res.0.2/"
if (!file.exists(Res)){
  dir.create(Res)
}

p9 <- VlnPlot(object = mydata, group.by = group_i, features = example, pt.size = 0.000001, ncol = 4)
ggsave(filename=paste0(Res,"Vlnplot-markers.pdf"),plot=p9, width=20, height = 15)
ggsave(filename=paste0(Res,"Vlnplot-markers.png"),plot=p9, width=20, height = 15, units = "in", dpi = 300)

p10 <- FeaturePlot(object = mydata, reduction = "umap", features = example, pt.size = 0.000001, ncol = 4)
ggsave(filename=paste0(Res,"FeaturePlot-markers.pdf"),plot=p10, width=20, height = 15)
ggsave(filename=paste0(Res,"FeaturePlot-markers.png"),plot=p10, width=20, height = 15, units = "in", dpi = 300)

# p10 <- FeaturePlot(object = mydata, reduction = "umap", features = "GATA3", pt.size = 0.000001)
# p10

### AddModuleScores
modules <- split(features.top5, ceiling(seq_along(features.top5) / 5))
names(modules) <- paste0("module", seq_along(modules) - 1)

mydata <- Seurat::AddModuleScore(object = mydata, features = modules, name = "ModuleScore", ctrl = 100, seed = 123)

# Old name
colnames(mydata@meta.data)


# Rename
mydata@meta.data <- mydata@meta.data %>%
  rename_with(
    ~ paste0("module", seq_along(.) - 1),
    starts_with("ModuleScore")
  )

# mydata@meta.data <- mydata@meta.data %>%
#   rename(module0 = ModuleScore1, 
#          module1 = ModuleScore2,
#          module2 = ModuleScore3,
#          module3 = ModuleScore4,
#          module4 = ModuleScore5,
#          module5 = ModuleScore6, 
#          module6 = ModuleScore7,
#          module7 = ModuleScore8,
#          module8 = ModuleScore9,
#          module9 = ModuleScore10, 
#          module10 = ModuleScore11,
#          module11 = ModuleScore12)

# New name
colnames(mydata@meta.data)

module_names <- paste0("module", seq_along(modules) - 1)

p11 <- VlnPlot(object = mydata, group.by = group_i, features = module_names, pt.size = 0.000001, ncol = 4)
ggsave(filename=paste0(Res,"Vlnplot-modules.pdf"),plot=p11, width=20, height = 15)
ggsave(filename=paste0(Res,"Vlnplot-modules.png"),plot=p11, width=20, height = 15, units = "in", dpi = 300)

p12 <- FeaturePlot(object = mydata, reduction = "umap", features = module_names, pt.size = 0.000001, ncol = 4)
ggsave(filename=paste0(Res,"FeaturePlot-modules.pdf"),plot=p12, width=20, height = 15)
ggsave(filename=paste0(Res,"FeaturePlot-modules.png"),plot=p12, width=20, height = 15, units = "in", dpi = 300)

# Assign cell annotation
mydata@meta.data <- mydata@meta.data %>%
  mutate(celltype = recode(RNA_snn_res.0.2,
                           `0` = "celltype0",
                           `1` = "celltype1",
                           `2` = "celltype2",
                           `3` = "celltype3",
                           `4` = "celltype4",
                           `5` = "celltype5",
                           `6` = "celltype6",
                           `7` = "celltype7",
                           `8` = "celltype8",
                           `9` = "celltype9",
                           `10` = "celltype10",
                           `11` = "celltype11",
  ))
head(mydata@meta.data$celltype)

mydata@meta.data$celltype <- factor(mydata@meta.data$celltype,
                                    levels=c("celltype0","celltype1","celltype2",
                                             "celltype3","celltype4","celltype5",
                                             "celltype6","celltype7","celltype8",
                                             "celltype9","celltype10","celltype11"
                                             )
                                    )  
head(mydata@meta.data$celltype)

saveRDS(mydata, "mydata-processed_after_score_celltype.rds")


### DEG between two cell types
#mydata <- readRDS("mydata-processed_after_score_celltype.rds.rds")

# Compare "celltype0" and "celltype1" as an example
DEG <- FindMarkers(
  object = mydata,
  ident.1 = "celltype0",  
  ident.2 = "celltype1", 
  group.by = "celltype", 
  test.use = "MAST" 
)

head(DEG,2)
DEG$gene <- rownames(DEG)
fwrite(DEG, paste0(Res,'DEG-celltype0_celltype1.csv'))
saveRDS(DEG, paste0(Res,'DEG-celltype0_celltype1.rds'))

### Volcano plot
p13 <- EnhancedVolcano(
  DEG,
  lab = NA,
  x = 'avg_log2FC',
  y = 'p_val_adj',
  ylab = bquote( ~ -Log[10] ~ "adjusted p-value"),
  pCutoff = 0.05,
  FCcutoff = 1.333,
  title = "celltype0 vs. celltype1",
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
ggsave(filename=paste0(Res,"volcano-celltype0_celltype1.pdf"),plot=p14, width=8, height = 5)
ggsave(filename=paste0(Res,"volcano-celltype0_celltype1.png"),plot=p14, width=8, height = 5, units = "in", dpi = 300)


### Export the package information 
sink("sessionInfo.txt")
sessionInfo()
sink()



