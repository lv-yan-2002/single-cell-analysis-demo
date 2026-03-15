# ====================== 02_标准化+聚类+差异分析.R ======================
# 作者：lvyan
# 日期：2026-03-15
# 用途：基于质控后的 pbmc 对象进行标准化、聚类和差异分析

# ---------------------- 1. 加载质控后的 pbmc 对象 ----------------------
pbmc <- readRDS("results/tables/01_qc_filtered_pbmc.rds")

# ---------------------- 2. 标准化 ----------------------
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# ---------------------- 3. 识别高变基因 ----------------------
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# ---------------------- 4. 数据缩放 ----------------------
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# ---------------------- 5. PCA 降维 ----------------------
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# ---------------------- 6. UMAP 聚类 ----------------------
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# ---------------------- 7. 绘制 UMAP 图 ----------------------
pdf("results/figures/02_umap_clusters.pdf", width = 8, height = 6)
print(DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE))
dev.off()

# ---------------------- 8. 差异表达分析 ----------------------
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(markers)

# ---------------------- 9. 标记基因可视化 ----------------------
FeaturePlot(pbmc, features = c("CD3D", "CD19", "CD14", "GNLY"), ncol = 2)