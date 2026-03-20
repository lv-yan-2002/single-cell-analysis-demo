# ====================== 1. 加载包 ======================
library(Seurat)
library(ggplot2)
library(patchwork)

# ====================== 2. 加载本地 pbmc3k 数据 ======================
pbmc3k <- Read10X(data.dir = getwd())
pbmc3k <- CreateSeuratObject(
  counts = pbmc3k,
  project = "pbmc3k",
  min.cells = 3,
  min.features = 200
)

# ====================== 3. 预处理 ======================
pbmc3k <- NormalizeData(pbmc3k)
pbmc3k <- FindVariableFeatures(pbmc3k)
pbmc3k <- ScaleData(pbmc3k, features = rownames(pbmc3k)) # 缩放所有基因，解决PTPRCAP警告

# ====================== 4. 降维+聚类 ======================
pbmc3k <- RunPCA(pbmc3k)
pbmc3k <- FindNeighbors(pbmc3k, dims = 1:10)
pbmc3k <- FindClusters(pbmc3k, resolution = 0.5)
pbmc3k <- RunUMAP(pbmc3k, dims = 1:10)

# ====================== 5. 细胞类型注释 ======================
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", 
                     "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc3k)
pbmc3k <- RenameIdents(pbmc3k, new.cluster.ids)

# ====================== 6. 可视化（所有示例） ======================
# 定义特征基因
features <- c("LYZ", "CCL5", "IL32", "PTPRCAP", "FCGR3A", "PF4")

# RidgePlot
RidgePlot(pbmc3k, features = features, ncol = 2)

# FeaturePlot（单基因）
FeaturePlot(pbmc3k, features = "MS4A1")

# FeaturePlot（双基因共表达）
FeaturePlot(pbmc3k, features = c("MS4A1", "CD79A"), blend = TRUE)

# 分组FeaturePlot（先创建随机分组）
pbmc3k$groups <- sample(c("group1", "group2"), size = ncol(pbmc3k), replace = TRUE)
FeaturePlot(pbmc3k, features = c("MS4A1", "CD79A"), split.by = "groups")

# 计算线粒体比例+分组小提琴图
pbmc3k <- PercentageFeatureSet(pbmc3k, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(pbmc3k, features = "percent.mt", split.by = "groups")