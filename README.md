# single-cell-analysis-demo
Demo for single-cell RNA-seq analysis (Seurat/Scanpy）
# 单细胞 RNA-seq 分析示例仓库
本仓库用于学习和记录单细胞测序数据分析流程，基于 R/Seurat 或 Python/Scanpy 实现。

## 分析流程
1. 数据质控（QC）：过滤低质量细胞、线粒体基因比例过高的细胞
2. 数据标准化与归一化
3. 细胞聚类与降维（PCA + UMAP）
4. 差异基因分析
5. 细胞类型注释

## 环境依赖
### R 环境（Seurat）
- R version ≥ 4.2.0
- Seurat ≥ 4.0.0
- tidyverse, ggplot2

### Python 环境（Scanpy）
- Python ≥ 3.9
- scanpy ≥ 1.9.0
- pandas, numpy, matplotlib

## 文件说明
- `scripts/`：分析脚本（分步骤）
- `data/`：数据说明（无原始数据）
- `results/`：分析结果（图片、表格）

## 数据来源
示例数据：PBMC 3k 数据集（10x Genomics 公开数据）