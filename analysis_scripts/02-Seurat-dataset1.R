library(data.table)
library(tidyverse)
library(Seurat)


h5_files <- "./count/outs/filtered_feature_bc_matrix.h5"
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

mydata <- readRDS("mydata-before_filter.rds")

# Compute 10th and 90th percentiles for QC metrics
qc_cutoffs <- list(
  nCount_RNA   = quantile(mydata$nCount_RNA,   probs = c(0.10, 0.90)),
  nFeature_RNA = quantile(mydata$nFeature_RNA, probs = c(0.10, 0.90)),
  percent_mito = quantile(mydata$percent_mito, probs = 0.90),  # only upper cutoff
  percent_ribo = quantile(mydata$percent_ribo, probs = 0.90)   # only upper cutoff
)

qc_cutoffs
# $nCount_RNA
#    10%     90% 
# 4417.9 20937.6 

# $nFeature_RNA
#    10%    90% 
# 1956.0 5014.4 

# $percent_mito
#      90% 
# 11.18076 

# $percent_ribo
#     90% 
# 29.44844 

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

#save the Seurat object
saveRDS(mydata, "mydata-after_filter.rds", compress = T)


