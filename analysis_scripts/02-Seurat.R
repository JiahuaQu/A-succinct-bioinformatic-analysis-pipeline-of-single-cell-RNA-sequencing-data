library(clustree)
library(ggplotify)
library(data.table)
library(tidyverse)
library(Seurat)
library(hdf5r)

# h5_files <- "/home/jqu/project/scRNAseq_pipeline/count/outs/filtered_feature_bc_matrix.h5"
# h5_read <- Read10X_h5(h5_files)
# h5_seurat <- CreateSeuratObject(h5_read,project="sample1")
# saveRDS(h5_seurat, "h5_seurat.rds", compress = T)

# mydata <- readRDS("h5_seurat.rds")
# mydata <- h5_seurat
# rm(h5_seurat)

### Subset
# # 获取所有细胞名称
all_cells <- colnames(mydata)

# 随机打乱细胞顺序
set.seed(123)  # 确保随机性可重复
shuffled_cells <- sample(all_cells)

# 确定子集数和分割索引
num_subsets <- 4  # 子集数量
subset_sizes <- floor(length(shuffled_cells) / num_subsets)
remainder <- length(shuffled_cells) %% num_subsets  # 计算多余细胞数

# 为每个子集分配细胞索引
indices <- rep(1:num_subsets, each = subset_sizes)
if (remainder > 0) {
  indices <- c(indices, sample(1:num_subsets, remainder))  # 随机分配多余细胞
}

# 将细胞分配到子集
split_cells <- split(shuffled_cells, indices)

# 创建 Seurat 对象的子集列表
subsets <- lapply(split_cells, function(cells) subset(mydata, cells = cells))

names(subsets) <- c("subset1", "subset2", "subset3", "subset4")

# Merge
# 定义子集名称
cell_ids <- c("subset1", "subset2", "subset3", "subset4")

# 将子集合并为一个 Seurat 对象
subset_merge <- merge(
  x = subsets[[1]], 
  y = subsets[2:length(subsets)], 
  add.cell.ids = cell_ids, 
  project = "Merge"
)

# 确保 orig.ident 被正确更新
subset_merge@meta.data$orig.ident <- gsub("_.*$", "", rownames(subset_merge@meta.data))

# 检查 orig.ident 是否正确
unique(subset_merge@meta.data$orig.ident)
saveRDS(subset_merge, "mydata_merge.rds", compress = T)


###
mydata <- readRDS("mydata_merge.rds")

# Need to join layers before performing any differential expression analysis
Layers(mydata[["RNA"]])
# [1] "counts.1"   "counts.2"   "counts.3"   "counts.4"

# JoinLayers
mydata[["RNA"]] <- JoinLayers(mydata[["RNA"]])
Layers(mydata[["RNA"]])
# [1] "counts"

saveRDS(mydata, "mydata-merge_joinlayers.rds", compress = T)


###
mydata <- readRDS("mydata-merge_joinlayers.rds")

mydata$type <- "all"
### Filter
mydata <- PercentageFeatureSet(object = mydata, pattern = "^MT-", col.name = "percent_mito",assay = "RNA")
mydata <- PercentageFeatureSet(object = mydata, pattern = "^RP[SL]", col.name = "percent_ribo",assay = "RNA")
mydata <- PercentageFeatureSet(object = mydata, pattern = "^HBA|^HBB", col.name = "percent_hb",assay = "RNA")

feats1 <- c("nFeature_RNA", "nCount_RNA")
Vlnplot1 <- VlnPlot(object = mydata, group.by = "orig.ident", features = feats1, pt.size = 0.000001, ncol = 2)
ggsave(filename="Vlnplot1.pdf",plot=Vlnplot1, width=5, height = 5)
ggsave(filename="Vlnplot1.png",plot=Vlnplot1, width=5, height = 5, units = "in", dpi = 300)

Vlnplot2 <- VlnPlot(object = mydata, group.by = "type", features = feats1, pt.size = 0.000001, ncol = 2) 
ggsave(filename="Vlnplot2.pdf",plot=Vlnplot2, width=4, height = 5)
ggsave(filename="Vlnplot2.png",plot=Vlnplot2, width=4, height = 5, units = "in", dpi = 300)

feats2 <- c("percent_mito", "percent_ribo", "percent_ribo", "percent_hb")
Vlnplot3 <- VlnPlot(object = mydata, group.by = "orig.ident", features = feats2, pt.size = 0.000001, ncol = 3) 
ggsave(filename="Vlnplot3.pdf",plot=Vlnplot3, width=7.5, height = 5)
ggsave(filename="Vlnplot3.png",plot=Vlnplot3, width=7.5, height = 5, units = "in", dpi = 300)

Vlnplot4 <- VlnPlot(object = mydata, group.by = "type", features = feats2, pt.size = 0.000001, ncol = 3) 
ggsave(filename="Vlnplot4.pdf",plot=Vlnplot4, width=5, height = 5)
ggsave(filename="Vlnplot4.png",plot=Vlnplot4, width=5, height = 5, units = "in", dpi = 300)

saveRDS(mydata, "mydata-processed_before_filter.rds", compress = T)


###
mydata <- readRDS("mydata-processed_before_filter.rds")

fivenum(mydata@meta.data$nCount_RNA)
# 514  7746  9975 13498 70018

fivenum(mydata@meta.data$nFeature_RNA)
# 15.0 2779.0 3285.5 3904.0 8132.0

fivenum(mydata@meta.data$percent_mito)
# 0.000000  5.677775  6.760239  8.303887 98.507463

fivenum(mydata@meta.data$percent_ribo)
# 0.00000 10.47143 16.72760 25.05176 42.11861

fivenum(mydata@meta.data$percent_hb)
# 0.00000  0.00000  0.00000  0.00000 97.27626

mydata <- subset(object = mydata, subset = nFeature_RNA >= 200 & 
                   nFeature_RNA <= 6000 &
                   percent_mito <= 10 & 
                   percent_ribo <= 40 & 
                   percent_hb < 20)

feats1 <- c("nFeature_RNA", "nCount_RNA")
Vlnplot5 <- VlnPlot(object = mydata, group.by = "orig.ident", features = feats1, pt.size = 0.000001, ncol = 2) 
ggsave(filename="Vlnplot5.pdf",plot=Vlnplot5, width=5, height = 5)
ggsave(filename="Vlnplot5.png",plot=Vlnplot5, width=5, height = 5, units = "in", dpi = 300)

Vlnplot6 <- VlnPlot(object = mydata, group.by = "type", features = feats1, pt.size = 0.000001, ncol = 2) 
ggsave(filename="Vlnplot6.pdf",plot=Vlnplot6, width=4, height = 5)
ggsave(filename="Vlnplot6.png",plot=Vlnplot6, width=4, height = 5, units = "in", dpi = 300)

feats2 <- c("percent_mito", "percent_ribo", "percent_ribo", "percent_hb")
Vlnplot7 <- VlnPlot(object = mydata, group.by = "orig.ident", features = feats2, pt.size = 0.000001, ncol = 3) 
ggsave(filename="Vlnplot7.pdf",plot=Vlnplot7, width=7.5, height = 5)
ggsave(filename="Vlnplot7.png",plot=Vlnplot7, width=7.5, height = 5, units = "in", dpi = 300)

Vlnplot8 <- VlnPlot(object = mydata, group.by = "type", features = feats2, pt.size = 0.000001, ncol = 3) 
ggsave(filename="Vlnplot8.pdf",plot=Vlnplot8, width=5, height = 5)
ggsave(filename="Vlnplot8.png",plot=Vlnplot8, width=5, height = 5, units = "in", dpi = 300)

saveRDS(mydata, "mydata-processed_after_filter.rds", compress = T)


###
mydata <- readRDS("mydata-processed_after_filter.rds")

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

saveRDS(mydata, "mydata-processed_after_filter-harmony.rds", compress = T)


###
mydata <- readRDS("mydata-processed_after_filter-harmony.rds")

harmony_dim <- 1:15
mydata <- RunUMAP(object = mydata, reduction = "harmony", dims = harmony_dim)

### clustering
#DefaultAssay(mydata) <- "RNA"
mydata <- FindNeighbors(object = mydata, reduction = "harmony", dims = harmony_dim)

res <- seq(0.1, 1, 0.1)
mydata <- FindClusters(object = mydata, resolution = res, verbose = FALSE)

p5 <- clustree(mydata)
p5 <- as.ggplot(p5)
ggsave(filename="clustree.pdf",plot=p5, width=7, height = 10)
ggsave(filename="clustree.png",plot=p5, width=7, height = 10, units = "in", dpi = 300)

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
mydata@meta.data$RNA_snn_res.0.1 <- factor(mydata@meta.data$RNA_snn_res.0.1,levels=0:4)

markers.top5 <- fread("resolution/RNA_snn_res.0.1-markers.top5.csv")
features.top5 <- markers.top5$gene
NROW(unique(features.top5))   # 25
example_5 <- features.top5[c(1,6,11,16,21)]
example_5
# "LEF1"    "S100A8"  "COL19A1" "SLC4A10" "SH2D1B" 


### Plots
Res <- "res.0.1/"
if (!file.exists(Res)){
  dir.create(Res)
}

p9 <- VlnPlot(object = mydata, group.by = group_i, features = example_5, pt.size = 0.000001, ncol = 2)
ggsave(filename=paste0(Res,"Vlnplot-markers.pdf"),plot=p9, width=15, height = 15)
ggsave(filename=paste0(Res,"Vlnplot-markers.png"),plot=p9, width=15, height = 15, units = "in", dpi = 300)

p10 <- FeaturePlot(object = mydata, reduction = "umap", features = example_5, pt.size = 0.000001, ncol = 2)
ggsave(filename=paste0(Res,"FeaturePlot-markers.pdf"),plot=p10, width=10, height = 15)
ggsave(filename=paste0(Res,"FeaturePlot-markers.png"),plot=p10, width=10, height = 15, units = "in", dpi = 300)


### AddModuleScores
modules <- list(
  module0=features.top5[1:5],
  module1=features.top5[6:10],
  module2=features.top5[11:15],
  module3=features.top5[16:20],
  module4=features.top5[21:25]
)

mydata <- Seurat::AddModuleScore(object = mydata,features = modules, name = "ModuleScore", ctrl = 100,seed = 123)

colnames(mydata@meta.data)
# New columns
# "ModuleScore1" "ModuleScore2" "ModuleScore3" "ModuleScore4" "ModuleScore5"

# Rename
mydata@meta.data <- mydata@meta.data %>%
  rename(module0 = ModuleScore1, 
         module1 = ModuleScore2,
         module2 = ModuleScore3,
         module3 = ModuleScore4,
         module4 = ModuleScore5)

colnames(mydata@meta.data)
# New name
# "module0" "module1" "module2" "module3" "module4"

module_5 <- c("module0", "module1", "module2", "module3", "module4")

p11 <- VlnPlot(object = mydata, group.by = group_i, features = module_5, pt.size = 0.000001, ncol = 2)
ggsave(filename=paste0(Res,"Vlnplot-modules.pdf"),plot=p11, width=15, height = 15)
ggsave(filename=paste0(Res,"Vlnplot-modules.png"),plot=p11, width=15, height = 15, units = "in", dpi = 300)

p12 <- FeaturePlot(object = mydata, reduction = "umap", features = module_5, pt.size = 0.000001, ncol = 2)
ggsave(filename=paste0(Res,"FeaturePlot-modules.pdf"),plot=p12, width=10, height = 15)
ggsave(filename=paste0(Res,"FeaturePlot-modules.png"),plot=p12, width=10, height = 15, units = "in", dpi = 300)

# Assign cell annotation
mydata@meta.data <- mydata@meta.data %>%
  mutate(celltype = recode(RNA_snn_res.0.1,
                           `0` = "cluster0",
                           `1` = "cluster1",
                           `2` = "cluster2",
                           `3` = "cluster3",
                           `4` = "cluster4"
  ))
head(mydata@meta.data$celltype)

mydata@meta.data$celltype <- factor(mydata@meta.data$celltype,
                                    levels=c("cluster0","cluster1","cluster2",
                                             "cluster3","cluster4"))  
head(mydata@meta.data$celltype)

saveRDS(mydata, "mydata-processed_after_score_celltype.rds", compress = T)


### DEG between two cell types or among all cell types
mydata <- readRDS("mydata-processed_after_score_celltype.rds.rds")

# Compare "cluster0" and "cluster1" as an example
DEG <- FindMarkers(
  object = mydata,
  ident.1 = "cluster0",  # 第一个比较组,分子
  ident.2 = "cluster1",  # 第二个比较组,分母
  group.by = "celltype",  # 按照 "celltype" 列分组
  test.use = "MAST"  # 默认使用 Wilcoxon 检验,改为MAST
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


# 计算显著增减点的数量
significant_down <- sum(DEG$avg_log2FC <= -1.333 & DEG$p_val_adj < 0.05, na.rm = TRUE)
significant_up <- sum(DEG$avg_log2FC >= 1.333 & DEG$p_val_adj < 0.05, na.rm = TRUE)

p14 <- p13 +
  annotate("text", x = -12.5, y = max(-log10(DEG$p_val_adj), na.rm = TRUE), 
           label = paste("Down:\n", significant_down), 
           hjust = 0, size = 5) +
  annotate("text", x = 7.5, y = max(-log10(DEG$p_val_adj), na.rm = TRUE), 
           label = paste("Up:\n", significant_up), 
           hjust = 0, size = 5) +
  guides(
    color = guide_legend(
      keyheight = unit(1.5, "cm"),  # 控制图例项之间的垂直间距
      override.aes = list(size = 7) # 调整图例中点的大小
    )
    )
ggsave(filename=paste0(Res,"volcano-cluster0_cluster1.pdf"),plot=p14, width=8, height = 5)
ggsave(filename=paste0(Res,"volcano-cluster0_cluster1.png"),plot=p14, width=8, height = 5, units = "in", dpi = 300)






