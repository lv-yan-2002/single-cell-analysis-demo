# 加载所需包
library(Seurat)
library(tidyverse)
library(ggplot2)

# 1. 读取数据（以 10x PBMC 3k 为例，替换成你的数据路径）
# 如果是自己的数据，替换成 Read10X() 读取mtx格式，或 ReadH5AD() 读取h5ad格式
pbmc.data <- Read10X(data.dir = "你的数据路径/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "PBMC3k", min.cells = 3, min.features = 200)

# 2. 计算质控指标
# 线粒体基因比例（关键质控指标）
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# 3. 质控可视化（查看分布，确定过滤阈值）
pdf("../results/figures/01_qc_violin.pdf", width = 10, height = 6)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# 4. 过滤低质量细胞（根据自己的数据调整阈值）
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# 5. 保存质控后的数据
saveRDS(pbmc, file = "../results/tables/01_qc_filtered_pbmc.rds")

# 输出质控结果
cat("质控完成！保留细胞数：", ncol(pbmc), "\n")
cat("保留基因数：", nrow(pbmc), "\n")