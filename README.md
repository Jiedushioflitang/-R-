# 小鼠卵巢基质细胞ECM重塑相关转录组特征分析

## 项目概述

本项目基于Niu & Spradling 2020 (PNAS)小鼠卵巢发育单细胞数据 (GEO: GSE136441)，专注于分析卵巢基质细胞在卵泡形成阶段的ECM重塑、TGF-β信号增强、基质重塑和细胞黏附相关转录组特征。
原始数据存放于美国国家生物技术信息中心（NCBI）的
基因表达综合数据库（Gene Expression Omnibus, GEO），
系列登录号为 GSE136441\footnote{数据集下载页面：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136441；
关联生物项目编号：PRJNA562536；
原始测序数据存放于 SRA 数据库，登录号：SRP219755}，
测序平台为 Illumina NextSeq 500（GEO 平台编号 GPL19057）。

### 研究背景

- 原文重点研究BPG/EPG前颗粒细胞分化路径
- 本项目从新角度出发，分析卵巢基质细胞在卵泡形成阶段的转录组变化
- 关注ECM重塑和TGF-β相关信号是否在该阶段增强

### 分析重点

1. **ECM_core**: 细胞外基质核心结构成分（胶原、纤连蛋白、基质蛋白聚糖等）
2. **TGFb_signaling**: TGF-β/BMP相关发育和组织重塑信号
3. **Matrix_remodeling**: ECM降解、修饰、交联和重排过程
4. **Adhesion_integrin**: 细胞与ECM的黏附、定位和相互作用

## 环境配置

### R版本要求
- R 4.0或以上版本

### 必需R包
```r
# 基础分析包
cran_packages <- c(
  "dplyr",
  "data.table",
  "tidyverse",
  "ggplot2",
  "Matrix"
)

# 单细胞分析包
seurat_packages <- c(
  "Seurat",
  "SeuratObject",
  "patchwork"
)

# 生物信息学包
bioconductor_packages <- c(
  "clusterProfiler",
  "org.Mm.eg.db",
  "enrichplot"
)

# 统计分析包
stats_packages <- c(
  "car",
  "FSA",
  "multcomp",
  "clustree"
)

# 降维和去批次包
reduction_packages <- c(
  "harmony",
  "FNN",
  "decontX"
)

# 双细胞检测包
doublet_packages <- c(
  "DoubletFinder"
)

# 安装命令
install.packages(cran_packages, dependencies = TRUE)
BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "enrichplot"))
BiocManager::install("harmony")
BiocManager::install("decontX")
install.packages("DoubletFinder")
```

### 系统要求
- 最低内存：16GB RAM（推荐32GB+）
- 硬盘空间：至少10GB可用空间
- 推荐操作系统：Windows 10+, macOS 10.14+, 或 Linux发行版

## 数据源

- **数据来源**: Niu & Spradling 2020 (PNAS) 小鼠卵巢发育单细胞数据
- **GEO编号**: GSE136441
- **发育时期**: E11.5 - P1
- **细胞类型**: 主要关注Mesenchymal基质细胞

## 分析流程

### 第一部分：单细胞数据预处理

1. **数据读取**
   - 读取10X单细胞表达矩阵
   - 创建Seurat对象
   - 设置质量控制参数

2. **质量控制**
   - 计算线粒体基因比例 (`percent.mt`)
   - 计算血红蛋白基因比例 (`percent.hb`)
   - 根据 `nFeature_RNA`、`percent.mt`、`percent.hb` 进行质控过滤
   - 生成质控前后对比图表

3. **数据标准化**
   - LogNormalize标准化
   - 筛选高变基因 (2000个)
   - 数据缩放与校正

4. **降维分析**
   - PCA主成分分析
   - 使用Harmony去除批次效应
   - t-SNE和UMAP降维

5. **双细胞检测与过滤**
   - 使用decontX工具进行去污染
   - 手动实现双细胞检测算法
   - 保留高质量单细胞

6. **细胞聚类与注释**
   - 基于Harmony结果进行聚类
   - 生成多种分辨率聚类结果
   - 细胞类型注释

### 第二部分：基质细胞ECM重塑分析

1. **提取Mesenchymal基质细胞**
   - 从总细胞中提取基质细胞亚群
   - 按发育时期分组 (E11.5-P1)

2. **构建基因集**
   - ECM_core: 细胞外基质核心成分
   - TGFb_signaling: TGF-β/BMP信号通路
   - Matrix_remodeling: 基质重塑相关基因
   - Adhesion_integrin: 细胞黏附相关基因

3. **模块评分分析**
   - 使用AddModuleScore计算单细胞模块分数
   - 比较不同发育时期的模块分数变化
   - 绘制模块评分小提琴图和趋势图

4. **关键基因表达分析**
   - 选择代表性ECM/TGF-β/黏附相关基因
   - 绘制热图和DotPlot

5. **差异表达分析**
   - 早期(E11.5-E14.5) vs 晚期(E16.5-P1) 基质细胞比较
   - Wilcoxon秩和检验
   - 绘制火山图

6. **GO富集分析**
   - 对晚期上调和早期上调基因分别进行GO BP富集
   - 生成富集图和结果表

### 第三部分：结果可视化与报告

1. **数据质量检查**
   - 生成QC图表
   - 双细胞分布图
   - 细胞类型标记验证

2. **结果可视化**
   - 细胞类型组成图
   - 模块评分图
   - 基因表达热图
   - 差异表达火山图
   - GO富集图

## 主要结果文件

- `fig1_cell_type_fraction_by_period.svg`: 各时期细胞类型组成图
- `fig2_mesenchymal_tsne_by_period.svg`: 基质细胞tSNE时期分布图
- `fig3_module_scores_violin.svg`: 四类基质相关基因程序的模块评分小提琴图
- `fig4_module_score_trends.svg`: 模块评分趋势图
- `fig5_key_gene_heatmap.svg`: 代表性ECM/TGF-β/黏附相关基因表达热图
- `fig6_key_gene_dotplot.svg`: 气泡图
- `fig7_late_vs_early_volcano.svg`: 早期和晚期基质细胞差异表达火山图
- `fig8_late_up_GO_BP.svg`: 晚期上调基因GO Biological Process富集图
- `fig8_early_up_GO_BP.svg`: 早期上调基因GO Biological Process富集图

## 运行说明

### 1. 环境准备
```bash
# 安装所需R包
Rscript install_dependencies.R
```

### 2. 数据准备
确保以下文件存在于指定路径：
- `D:/database/all/49/data/`: 10X单细胞数据
- `D:/database/all/49/seurat_object.RData`: 预处理的Seurat对象

### 3. 运行分析
```bash
# 主分析流程
Rscript f_code.R
```

### 4. 生成PPT报告
```bash
# 生成PPT演示文稿
python ppt.py
```

## 输出文件说明

### 表格文件
- `tables/cell_type_counts_by_period.csv`: 细胞类型计数统计表
- `tables/mesenchymal_counts_by_period.csv`: 基质细胞计数统计表
- `tables/curated_gene_sets_detected.csv`: 基因集检测情况表
- `tables/mesenchymal_late_vs_early_DE.csv`: 差异表达分析结果
- `tables/module_score_pairwise_wilcox.csv`: 模块评分统计检验表

### 图像文件
- `figures/`: 所有生成的SVG矢量图

## 注意事项

1. **数据处理耗时**: 单细胞数据处理需要较长时间，特别是双细胞检测和降维分析
2. **内存使用**: 大型单细胞数据集需要大量内存，建议使用高性能计算机
3. **随机性**: 某些步骤（如AddModuleScore）包含随机抽样，使用固定seed确保结果可重复
4. **结果解释**: 单细胞转录组分析主要提供表达相关性证据，结论应表述为转录组关联性观察而非因果关系

## 项目结构

```
project/
├── f_code.R          # 主分析脚本
├── ppt.py           # PPT生成脚本
├── README.md        # 项目说明文档
├── figures/         # 输出图像文件
├── tables/          # 输出统计表格
└── data/            # 输入数据
```

## 致谢

本项目基于Niu & Spradling 2020 (PNAS)的研究数据，感谢原作者提供高质量的单细胞数据资源。

## 许可证

仅供学术研究使用。
