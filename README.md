# GeneLowPaper

<!-- badges: start -->
<!-- badges: end -->

**GeneLowPaper** 是一个专为 "单基因泛癌低分文章" 分析流程设计的 R 包。它将原始的 `Single_Gene_Low_paper_folw.R` 脚本中的全部功能模块化封装，提供了从**泛癌差异表达、预后分析、临床相关性、功能富集、免疫浸润到单细胞可视化**的一站式分析函数。

---

## 安装

```r
# 1. 安装依赖（部分包来自 Bioconductor / GitHub）
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "limma", "edgeR", "org.Hs.eg.db",
                       "clusterProfiler", "GSVA", "TCGAplot"))

# GitHub 包（如有需要）
# remotes::install_github("IOBR/IOBR")
# remotes::install_github("your/CyDataPro")
# remotes::install_github("your/SCNT")
# remotes::install_github("your/SeuratExtend")

# 2. 从本地安装 GeneLowPaper
devtools::install("path/to/GeneLowPaper")

# 3. 生成帮助文档（可选，如果 .Rd 文件尚未生成）
devtools::document("path/to/GeneLowPaper")
```

---

## 依赖包

| 类别 | 包名 |
|------|------|
| 差异分析 | `DESeq2`, `limma`, `edgeR` |
| 生存分析 | `survival`, `survminer`, `rms`, `ezcox` |
| 可视化 | `ggplot2`, `ggpubr`, `ggvenn`, `ggforce`, `ggbeeswarm`, `ggsignif`, `ggalluvial`, `patchwork`, `CyDataPro`, `RColorBrewer` |
| 富集分析 | `clusterProfiler`, `org.Hs.eg.db` |
| 免疫分析 | `IOBR`, `GSVA` |
| 单细胞 | `Seurat`, `SCNT`, `SeuratExtend`, `CellChat` |
| 数据处理 | `data.table`, `tidyverse`, `broom`, `forcats`, `igraph` |
| TCGA 数据 | `TCGAplot` |

---

## 快速开始

以下工作流完全对应原始脚本 `work_pipline.R` 的分析流程，以 **IL36RN** 在 **肺腺癌 (LUAD)** 中的分析为例。

```r
library(GeneLowPaper)

# 设置全局配色
color_dis <- color_dis_default

# 目标基因
gene <- "IL36RN"
TCGA_pro <- "LUAD"

# ========== 1. 泛癌差异表达 (Figure 1A) ==========
Fig1A_Pancer(
  gene = gene,
  plotWidth = 12, plotHeight = 6,
  outdir = "./",
  color_dis = color_dis
)

# ========== 2. 泛癌 KM 预后 (Figure 1B) ==========
Fig1B_PanKM(
  gene = gene,
  plotWidth = 6, plotHeight = 8,
  outdir = "./",
  color_dis = color_dis
)

# ========== 3. 差异癌种与预后癌种交集 Venn (Figure 1D) ==========
Fig1D_Venn(
  outdir = "./",
  plotWidth = 6, plotHeight = 6,
  show_element = TRUE,
  color_dis = color_dis
)

# ========== 4. GEO 验证集查询 ==========
BM_Gene_Find_datasets(
  TCGA_pro = TCGA_pro,
  gene = gene,
  data_dir = "F:/BioMed/预分析/data/GEO"   # 请修改为实际路径
)

# ========== 5. GEO 差异表达验证 (Figure 2A) ==========
Fig2A_GEO_box(
  GSE = "GSE19188",
  gene = gene,
  TCGA_pro = TCGA_pro,
  outdir = "./",
  plotWidth = 4, plotHeight = 4,
  color_dis = color_dis
)

# ========== 6. GEO 预后验证 (Figure 2B) ==========
Fig2A_GEO_KM(
  GSE = "GSE31210",
  gene = gene,
  TCGA_pro = TCGA_pro,
  plotWidth = 4, plotHeight = 4,
  outdir = "./",
  color_dis = color_dis
)

# ========== 7. TCGA 基因表达与临床相关性 (Figure 2D) ==========
Fig2D_TCGA_Clin(
  TCGA_pro = TCGA_pro,
  gene = gene,
  gene_type = "mRNA",
  outdir = "./",
  plotWidth = 8, plotHeight = 6,
  coxWidth = 8, coxHeight = 6,
  color_dis = color_dis
)

# ========== 8. GO / KEGG 富集分析 (Figure 3) ==========
Fig3_Enrich(
  TCGA_pro = TCGA_pro,
  gene = gene,
  gene_type = "mRNA",
  outdir = "./",
  plotHeight = 8,
  logFC_cut = 0.58,
  pvalue_cut = 0.05,
  color_dis = color_dis
)

# ========== 9. 免疫浸润相关性 (Figure 5A) ==========
Fig5_ImmFil(
  gene = gene,
  TCGA_pro = TCGA_pro,
  gene_type = "mrna",
  outdir = "./",
  plotWidth = 8, plotHeight = 6,
  corr_cut = 0.3,
  Imm_method = "mcpcounter",
  color_dis = color_dis
)

# ========== 10. ESTIMATE & IPS 免疫评分差异 (Figure 5B) ==========
Fig5_Est_IPS(
  gene = gene,
  TCGA_pro = TCGA_pro,
  gene_type = "mrna",
  outdir = "./",
  plotWidth = 6, plotHeight = 6.5,
  color_dis = color_dis
)

# ========== 11. 免疫检查点基因相关性 (Figure 5C) ==========
Fig5_ICC_cor(
  gene = gene,
  TCGA_pro = TCGA_pro,
  gene_type = "mrna",
  cor_method = "pearson",
  p_threshold = 0.05,
  r_threshold = 0,
  outdir = "./",
  plotWidth = 10, plotHeight = 10,
  color_dis = color_dis
)

# ========== 12. 单细胞数据可视化 (Figure 6) ==========
sc_basicPlot(
  scPath = "F:/Pritect_ana/115_IL36RN_肺腺癌/data/NSCLC_GSE148071_Single.Rds",
  gene = gene,
  reduction = "umap",
  outdir = "./",
  plotWidth = 8, plotHeight = 6,
  color_dis = color_dis
)

# ========== 13. 单基因细胞通讯网络图 (Figure 6B) ==========
gene_Cellchat1(
  scPath = "F:/Pritect_ana/115_IL36RN_肺腺癌/data/NSCLC_GSE148071_Single.Rds",
  gene = gene,
  Tar_celltype = "Malignant",
  outdir = "./",
  plotWidth = 8, plotHeight = 6,
  color_dis = color_dis
)

# ========== 14. 单基因通讯差异受配体 (Figure 6C) ==========
gene_Cellchat2(
  scPath = "F:/Pritect_ana/115_IL36RN_肺腺癌/data/NSCLC_GSE148071_Single.Rds",
  gene = gene,
  Tar_celltype = "Malignant",
  outdir = "./",
  plotWidth = 8, plotHeight = 6,
  color_dis = color_dis
)

# ========== 15. 单细胞基因分组差异与富集 (Figure 6D) ==========
sc_GeneEnrichr(
  scPath = "F:/Pritect_ana/115_IL36RN_肺腺癌/data/NSCLC_GSE148071_Single.Rds",
  gene = gene,
  Tar_celltype = "Malignant",
  pvlCut = 0.05, avgFC_cut = 0.3,
  organism = "hsa",
  kegg_data_path = "F:/BioMed/预分析/data/KEGG/KEGG_hsa_data.Rds",
  showterms = 6,
  outdir = "./",
  plotWidth = 6, plotHeight = 6,
  color_dis = color_dis
)

# ========== 16. Slingshot 拟时序轨迹分析 (Figure 6E) ==========
sc_slingshot(
  scPath = "F:/Pritect_ana/115_IL36RN_肺腺癌/data/NSCLC_GSE148071_Single.Rds",
  gene = gene,
  cell_start = NULL,
  cell_end = NULL,
  Reduction = "UMAP",
  outdir = "./",
  plotWidth = 6, plotHeight = 6,
  color_dis = color_dis
)

# ========== 17. Dorothea 转录因子活性分析 (Figure 6F) ==========
sc_dorotheaTF(
  scPath = "F:/Pritect_ana/115_IL36RN_肺腺癌/data/NSCLC_GSE148071_Single.Rds",
  gene = gene,
  Tar_celltype = "Malignant",
  topTF = 25,
  outdir = "./",
  plotWidth = 8, plotHeight = 4,
  color_dis = color_dis
)
```

---

## 数据准备

### 1. TCGA 数据

函数 `Fig2D_TCGA_Clin`、`Fig3_Enrich`、`Fig5_ImmFil`、`Fig5_Est_IPS`、`Fig5_ICC_cor` 需要读取 TCGA 表达矩阵。默认路径格式为：

```
{data_dir}/TCGA-{TCGA_pro}/TCGA-{TCGA_pro}_{gene_type}_expr_tpm.csv.gz
{data_dir}/TCGA-{TCGA_pro}/TCGA-{TCGA_pro}_clinical_OS.csv.gz
```

### 2. GEO 数据

函数 `BM_Gene_Find_datasets`、`Fig2A_GEO_box`、`Fig2A_GEO_KM` 需要预处理的 GEO 数据，目录结构如下：

```
{data_dir}/{TCGA_pro}/
├── 差异/
│   └── {GSE}_DEG.Rds
└── 预后/
    └── {GSE}_OS.Rds
```

每个 `.Rds` 文件应为包含 `exprmatix`（表达矩阵）和 `sample`（样本信息）的 `list` 对象。

### 3. 免疫检查点基因列表

`Fig5_ICC_cor` 默认读取：

```
F:/BioMed/预分析/data/IGG79Gene.csv
```

可通过 `icc_path` 参数修改。

### 4. 单细胞数据

`sc_basicPlot`、`gene_Cellchat1`、`gene_Cellchat2` 需要预处理好的 Seurat 对象 `.Rds` 文件，`meta.data` 中需包含 `celltype` 列。

---

## 核心函数概览

| 函数名 | 功能描述 | 对应图号 |
|--------|---------|---------|
| `BM_DEG_analysis()` | 差异表达分析集成（DESeq2 / limma / edgeR / Wilcoxon） | 通用 |
| `BM_Gene_Find_datasets()` | 扫描 GEO 数据集，筛选有差异和预后的验证集 | 表 1 |
| `Fig1A_Pancer()` | 单基因泛癌差异表达箱线图 | Fig 1A |
| `Fig1B_PanKM()` | 单基因泛癌 KM 预后（含森林图） | Fig 1B |
| `Fig1D_Venn()` | 差异癌种与预后癌种交集 Venn | Fig 1D |
| `Fig2A_GEO_box()` | GEO 差异表达验证（小提琴+箱线图） | Fig 2A |
| `Fig2A_GEO_KM()` | GEO 预后 KM 验证 | Fig 2B |
| `Fig2D_TCGA_Clin()` | TCGA 临床特征关联 + 单/多变量 Cox + 列线图 | Fig 2D |
| `Fig3_Enrich()` | GO / KEGG 富集分析 | Fig 3 |
| `Fig5_ImmFil()` | 免疫浸润相关性（棒棒糖图） | Fig 5A |
| `Fig5_Est_IPS()` | ESTIMATE & IPS 免疫评分差异 | Fig 5B |
| `Fig5_ICC_cor()` | 免疫检查点基因相关性 | Fig 5C |
| `sc_basicPlot()` | 单细胞注释 UMAP + 基因 DotPlot + FeaturePlot | Fig 6A |
| `gene_Cellchat1()` | 单基因细胞通讯网络图 | Fig 6B |
| `gene_Cellchat2()` | 基因高/低组通讯差异对比 | Fig 6C |
| `sc_GeneEnrichr()` | 单细胞基因分组差异 + 火山图 + GO/KEGG 富集 | Fig 6D |
| `sc_slingshot()` | Slingshot 拟时序轨迹分析 | Fig 6E |
| `sc_dorotheaTF()` | Dorothea 转录因子活性热图 | Fig 6F |

---

## 帮助文档

所有函数均使用 **roxygen2** 格式编写帮助文档。安装后可通过以下方式查看：

```r
?BM_DEG_analysis
?Fig1A_Pancer
?Fig2D_TCGA_Clin
?Fig5_ImmFil
# ... 等等
```

---

## 注意事项

1. **路径问题**：本包中的部分函数依赖外部预处理数据（TCGA、GEO、单细胞 RDS）。请将 `data_dir`、`scPath`、`icc_path` 等参数修改为您实际的数据路径。
2. **字体问题**：`mytheme` 使用了 `"serif"` 字体。在 Linux 服务器上若出现字体警告，可修改 `mytheme` 或安装相应字体。
3. **内存问题**：CellChat 分析（`gene_Cellchat1/2`）和免疫浸润分析（`Fig5_ImmFil`）在大型数据集上可能占用较多内存，建议在内存充足的环境中运行。
4. **包依赖**：部分依赖包（如 `IOBR`、`CyDataPro`、`SCNT`、`SeuratExtend`）可能需从 GitHub 安装，请提前确认。

---

## 作者

本包由原始分析脚本封装而来，保留了全部原始分析逻辑并增加了模块化参数控制。

---

## License

MIT License
