

#1.数据准备
#remove list
#garbage collection
rm(list=ls())
gc()

load("D:/database/all/49/seurat_object.RData")
load(file = "D:/database/all/49/newall_objects.RData")
load(file = "D:/database/all/49/all_objects.RData")
library(dplyr)
library(data.table)
library(Seurat)


library(car)
library(FSA)
library(multcomp)

#BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "enrichplot"))
library(clusterProfiler)  
library(org.Mm.eg.db)    
library(enrichplot) 
library(devtools)
library(ggplot2)
library(tidyverse)
library(clustree)
library(decontX)
library(patchwork)
library(harmony)
library(FNN) 

setwd("D:/database/all/49")
set.seed(123)
path_10x <- "D:/database/all/49/data"
#list.files()
dir <- list.files(path_10x)
dir

seurat_list <- list()

for(i in 1:length(dir)){
  print(i)
  data.path <- file.path(path_10x,dir[i])
  seurat_data <- Read10X(data.dir = data.path)
  seurat=CreateSeuratObject(counts = seurat_data,
                            project = dir[i],
                            min.features = 200,    #只保留那些检测到至少200个不同基因的细胞。这能过滤掉空液滴或死细胞
                            min.cells = 3)   #只保留那些在至少3个细胞中表达的基因。这能过滤掉测序噪音或极低表达的基因
  seurat_list=append(seurat_list,seurat)
}
seurat_list

#合并seurat_list中的7个成员，变为一个完整的seurat_object
seurat_object <- merge(seurat_list[[1]],
                       y=seurat_list[-1],
                       add.cell.ids=dir)            #这一步把dir(即样本分组也就是标识)赋给了Idents

seurat_object <- JoinLayers(seurat_object)      #合并层，把counts1 counts2 counts3等全部counts合并成一个counts
View(seurat_object@meta.data)
#Idents(seurat_object)
table(Idents(seurat_object))

#seurat_object后面跟$是访问meta.data
seurat_object$sample=Idents(seurat_object)

#save(seurat_object,file="D:/database/all/49")




#identical(rownames(seurat_object@assays$RNA@layers$counts), Features(seurat_object))

#为了方便储存，matrix即原始矩阵counts，行名是1 2 3等，列名是v1 v2 v3等，现在给他本来对应的基因名和细胞名
rownames(seurat_object@assays$RNA@layers$counts) <- Features(seurat_object)
colnames(seurat_object@assays$RNA@layers$counts) <- Cells(seurat_object)


#2. 质控(QC)
view(seurat_object[[]])
#用线粒体作为质控指标
#实际操作后发现，样本可能混入少量血液，主要是被红细胞污染，所以质控指标需要加入红细胞
#卵巢好像就是容易混入血液(不确定)
#但原文其实是由血液相关细胞的，“Blood-related 细胞” 指卵巢微环境中的内皮细胞、巨噬细胞等（是卵巢组织的 “正常背景细胞”），和红细胞还是不一样
#免疫细胞集中在 P1/P5、少量在 E18.5、其他胚胎期无
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object,pattern="^mt-")
seurat_object[["percent.hb"]] <- PercentageFeatureSet(seurat_object,pattern = "^Hbb|^Hba")
#neutrophil_markers <- c("S100a8", "S100a9", "Ngp", "Lcn2", "Camp", "Cd177")
#seurat_object[["percent.neutrophil"]] <- PercentageFeatureSet(seurat_object, features = neutrophil_markers)

#画质控前图片，便于后续质控
p1=VlnPlot(seurat_object,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)
p1
#ggsave("01_QCbefore.pdf",width = 15,height=6,plot=p1)
#ggsave("01_QCbefore.tiff",width=15,height=6,dpi=300,plot=p1)


plot1 <- FeatureScatter(seurat_object,feature1 = "nCount_RNA",feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
p2=plot1+plot2
p2
#ggsave("02_FeatureScatter.pdf",width = 12,height = 6,plot = p2)
#ggsave("02_FeatureScatter.tiff",width = 12,height = 6,dpi=300,plot=p2)


#过滤 “低质量细胞” 并验证过滤效果
#后两项过滤掉的不是特别多
seurat_object <- subset(seurat_object,subset=nFeature_RNA>300&nFeature_RNA<4500&percent.mt<10)
seurat_object <- subset(seurat_object, subset = percent.hb < 1)
#seurat_object <- subset(seurat_object, subset = percent.neutrophil < 1)
table(seurat_object$orig.ident)
#绘制质控后的小提琴图
p3 <- VlnPlot(seurat_object,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3,pt.size = 0)
p3
#ggsave("03_QCafer.pdf",width = 15,height = 6,plot = p3)
#ggsave("03_QCafter.tiff",width = 15,height = 6,dpi = 300,plot = p3)


#3. 标准化与特征选择
#数据标准化，进行缩放
#消除不同细胞间 “测序深度差异” 对基因表达量的影响
seurat_object <- NormalizeData(seurat_object,normalization.method = "LogNormalize",scale.factor=10000)
seurat_object@assays$RNA@layers$counts[1:15,1:15]
seurat_object@assays$RNA@layers$data[1:15,1:15]


#筛选高变基因 vst是方差稳定转换法 结果seurat_object会按程度高低排序
seurat_object <- FindVariableFeatures(seurat_object,selection.method = "vst",nfeatures=2000)

#可视化高变基因的分布特征，并标记出变异程度最高的 10 个基因
top10 <- head(VariableFeatures(seurat_object),10)
plot1 <- VariableFeaturePlot(seurat_object)
plot2 <- LabelPoints(plot=plot1,points=top10,repel = TRUE,xnudge=0,ynudge=0)
combined_plot=plot1+plot2+plot_layout(ncol=2,guides="collect")&theme(legend.position="bottom")
print(combined_plot)


#4. 缩放与降维

#decontX工具：专门用于单细胞 RNA-seq 数据中 “环境 RNA 污染” 分析的工具
#decontX工具识别并过滤单细胞数据中可能存在 “环境 RNA 污染” 的细胞
#实际结果去污染应该比较必要，有的污染值已经到了0.9多了
decontX_results <- decontX::decontX(seurat_object@assays$RNA@layers$counts)
seurat_object$Contamination <- decontX_results$contamination
rt <- seurat_object@meta.data
rt <- rt[which(rt$Contamination<0.3),]
#替换为低污染细胞
seurat_object <- seurat_object[,row.names(rt)]

#按样本拆分 RNA 数据并进行数据标准化（缩放），不同样本操作不一样可能会有影响
#同时校正线粒体基因比例对表达数据的影响，为后续降维和聚类分析做准备
seurat_object[["RNA"]] <- split(seurat_object[["RNA"]],f=seurat_object$sample)
#RNA下面的counts和data变成了counts1 data1 counts2 data2
#数据缩放与校正
#ScaleData功能：对标准化后的表达数据（data层）进一步处理，包括中心化（均值为 0）和缩放（标准差为 1），使不同基因的表达量具有可比性，同时校正指定的干扰因素（这里是线粒体基因比例
seurat_object <- ScaleData(seurat_object,vars.to.regress = c("percent.mt"))
#处理后的数据存储在RNA assay 的scale.data层（每个样本单独生成scale.data），后续的降维（如 PCA）和聚类分析会基于这层数据，结果更可靠。

#PCA 是一种降维方法，能将高维度的基因表达数据（每个细胞有数千个基因的表达量）压缩到低维度空间（通常几十个主成分），同时保留数据中的主要变异信息。
#即主成分分析
#选择高变基因的原因：它们更能反映细胞间的生物学差异，用这些基因做 PCA 能更清晰地分离不同细胞群体。
seurat_object <- RunPCA(seurat_object,features = VariableFeatures(object = seurat_object))
DimPlot(object=seurat_object,reduction = "pca",pt.size = 0.1,group.by = "sample")

#5. 去除批次效应
library(harmony)
#用 Harmony 消除批次效应
seurat_object <- RunHarmony(seurat_object,group.by.vars="sample",plot_convergence=TRUE)
p6 <- DimPlot(object = seurat_object,reduction="harmony",pt.size = 0.1,group.by = "sample")
p6
#ggsave("06_Harmony.pdf",width=8,height = 6,plot = p6)
#ggsave("06_Harmony.tiff",width=8,height=6,dpi=300,plot=p6)
table(seurat_object$orig.ident)

#5.5  去除双细胞
# ==============================


# ==============================
# 前置准备：加载包 + 内存优化
# ==============================

 # 用于K近邻计算（需安装：install.packages("FNN")）
#install.packages("FNN")
# 内存优化
# if (.Platform$OS.type == "windows") {
#   memory.limit(size = 200000)  # 200GB
# } else {
#   options(mc.cores = parallel::detectCores())
# }

# 备份原始对象
seurat_object_original <- seurat_object
cat("过滤双细胞前总细胞数：", ncol(seurat_object), "\n")

# ==============================
# 核心：纯手动实现双细胞检测（等效DoubletFinder）
# ==============================
# 提取所有样本
samples <- unique(seurat_object$sample)
seurat_list_doublet_removed <- list()

# 双细胞检测核心参数
pN <- 0.25    # 人工混样比例
pK <- 0.09    # K近邻比例
PCs <- 1:30   # 使用的PCA维度

for (s in samples) {
  cat("\n==============================")
  cat("\n开始处理样本：", s, "\n")
  
  # 1. 子集化样本
  obj_sub <- subset(seurat_object, sample == s)
  n_cells <- ncol(obj_sub)
  cat("样本", s, "原始细胞数：", n_cells, "\n")
  
  # 2. 确保PCA已运行
  if (!"pca" %in% Reductions(obj_sub)) {
    obj_sub <- RunPCA(obj_sub, 
                      features = VariableFeatures(obj_sub), 
                      npcs = 30, 
                      verbose = FALSE)
  }
  
  # 3. 提取PCA矩阵（纯数值矩阵）
  pca_data <- as.matrix(obj_sub@reductions$pca@cell.embeddings[, PCs])
  rownames(pca_data) <- colnames(obj_sub)
  
  # 4. 计算K近邻数（基于pK比例）
  k <- round(n_cells * pK)
  k <- max(k, 5)  # 至少5个近邻
  cat("样本", s, "K近邻数：", k, "\n")
  
  # 5. 估算预期双细胞数（行业标准）
  if (n_cells <= 500) {
    n_doublets <- round(n_cells * 0.016)
  } else if (n_cells <= 1000) {
    n_doublets <- round(n_cells * 0.024)
  } else if (n_cells <= 2000) {
    n_doublets <- round(n_cells * 0.039)
  } else if (n_cells <= 5000) {
    n_doublets <- round(n_cells * 0.060)
  } else {
    n_doublets <- round(n_cells * 0.075)
  }
  cat("样本", s, "预期双细胞数：", n_doublets, "\n")
  
  # 6. 核心算法：计算双细胞概率（pANN）
  # 步骤1：生成人工混样
  n_artificial <- round(n_cells * pN)
  #set.seed(123)  # 固定随机种子，结果可重复
  cell_idx1 <- sample(1:n_cells, n_artificial, replace = TRUE)
  cell_idx2 <- sample(1:n_cells, n_artificial, replace = TRUE)
  artificial_cells <- (pca_data[cell_idx1, ] + pca_data[cell_idx2, ]) / 2
  
  # 步骤2：合并真实细胞+人工混样
  combined_data <- rbind(pca_data, artificial_cells)
  
  # 步骤3：计算K近邻
  knn_result <- get.knn(combined_data, k = k)
  
  # 步骤4：计算每个真实细胞的双细胞概率（pANN）
  real_knn <- knn_result$nn.index[1:n_cells, ]  # 仅真实细胞的近邻
  artificial_indicator <- ifelse(real_knn > n_cells, 1, 0)  # 近邻是否为人工混样
  pANN <- rowMeans(artificial_indicator)  # 双细胞概率
  
  # 7. 划分Singlet/Doublet（基于pANN阈值）
  obj_sub$pANN <- pANN
  # 取pANN最高的n_doublets个为双细胞
  doublet_threshold <- quantile(pANN, 1 - n_doublets / n_cells)
  obj_sub$doublet_status <- ifelse(obj_sub$pANN >= doublet_threshold, "Doublet", "Singlet")
  
  # 8. 过滤双细胞
  obj_sub_singlet <- subset(obj_sub, doublet_status == "Singlet")
  n_retained <- ncol(obj_sub_singlet)
  cat("样本", s, "过滤后保留细胞数：", n_retained, "（保留率：", round(n_retained/n_cells*100, 2), "%）\n")
  
  # 9. 恢复原始meta.data
  if (n_retained > 0) {
    obj_sub_singlet@meta.data <- cbind(
      obj_sub_singlet@meta.data,
      seurat_object_original@meta.data[colnames(obj_sub_singlet), 
                                       setdiff(colnames(seurat_object_original@meta.data), 
                                               colnames(obj_sub_singlet@meta.data))]
    )
  }
  
  # 10. 存储样本
  seurat_list_doublet_removed[[s]] <- obj_sub_singlet
}

# ==============================
# 合并样本 + 重新下游分析
# ==============================
# 过滤空样本
seurat_list_doublet_removed <- seurat_list_doublet_removed[sapply(seurat_list_doublet_removed, function(x) ncol(x) > 0)]

if (length(seurat_list_doublet_removed) > 0) {
  # 合并样本
  seurat_object_clean <- merge(seurat_list_doublet_removed[[1]],
                               y = seurat_list_doublet_removed[-1],
                               add.cell.ids = names(seurat_list_doublet_removed))
  seurat_object_clean <- JoinLayers(seurat_object_clean)
  
  # 重新执行关键下游步骤
  cat("\n==============================")
  cat("\n过滤双细胞后重新标准化/降维...\n")
  
  seurat_object_clean <- ScaleData(seurat_object_clean,
                                   vars.to.regress = c("percent.mt"),
                                   verbose = FALSE)
  seurat_object_clean <- RunPCA(seurat_object_clean,
                                features = VariableFeatures(seurat_object_clean),
                                npcs = 30,
                                verbose = FALSE)
  seurat_object_clean <- RunHarmony(seurat_object_clean,
                                    group.by.vars = "sample",
                                    plot_convergence = TRUE,
                                    verbose = FALSE)
  seurat_object_clean <- FindNeighbors(seurat_object_clean,
                                       reduction = "harmony",
                                       dims = 1:20,
                                       verbose = FALSE)
  seurat_object_clean <- FindClusters(seurat_object_clean,
                                      resolution = 0.8,
                                      verbose = FALSE)
  seurat_object_clean <- RunTSNE(seurat_object_clean,
                                 reduction = "harmony",
                                 dims = 1:20,
                                 verbose = FALSE)
  seurat_object_clean <- RunUMAP(seurat_object_clean,
                                 reduction = "harmony",
                                 dims = 1:20,
                                 verbose = FALSE)
  
  # ==============================
  # 结果验证 + 保存
  # ==============================
  cat("\n==============================")
  cat("\n✅ 双细胞过滤结果汇总：\n")
  total_before <- ncol(seurat_object_original)
  total_after <- ncol(seurat_object_clean)
  total_removed <- total_before - total_after
  
  cat("过滤前总细胞数：", total_before, "\n")
  cat("过滤后总细胞数：", total_after, "\n")
  cat("去除双细胞数：", total_removed, "\n")
  cat("整体保留率：", round(total_after/total_before*100, 2), "%\n")
  
  # 样本级详细统计
  cell_count_before <- as.data.frame(table(seurat_object_original$sample))
  colnames(cell_count_before) <- c("sample", "count_before")
  cell_count_after <- as.data.frame(table(seurat_object_clean$sample))
  colnames(cell_count_after) <- c("sample", "count_after")
  cell_count <- merge(cell_count_before, cell_count_after, by = "sample", all.x = TRUE)
  cell_count$count_after[is.na(cell_count$count_after)] <- 0
  cell_count$retained_rate <- round(cell_count$count_after / cell_count$count_before * 100, 2)
  
  cat("\n📊 样本级保留率：\n")
  print(cell_count)
  
  # 保存最终对象
  seurat_object <- seurat_object_clean
  save_path <- "D:/database/all/49/seurat_object_doublet_removed.RData"
  save(seurat_object, file = save_path)
  cat("\n💾 过滤后的对象已保存至：", save_path, "\n")
  
} else {
  cat("\n❌ 错误：所有样本过滤后无细胞保留，请检查参数！\n")
}







#6. 聚类与可视化
#合并之前按sample拆分的RNA（count data)
seurat_object <- JoinLayers(seurat_object)
#降维可视化，查看最合适的主成分个数
p7=ElbowPlot(seurat_object,ndims = 50,reduction = "pca")+theme_bw()
p7


#细胞聚类的核心步骤，用于将基因表达模式相似的细胞归为一类（即 “细胞聚类”），为后续细胞类型注释和差异分析奠定基础
#指定基于 “Harmony 校正后的降维结果” 计算近邻。因为 Harmony 已消除样本批次效应，基于它构建的近邻关系更能反映细胞的真实生物学相似性
seurat_object <- FindNeighbors(seurat_object,reduction = "harmony",dims = 1:20)
#FindClusters(...)：Seurat 中实现细胞聚类的核心函数，默认使用 Louvain 算法（一种社区发现算法），基于FindNeighbors构建的近邻网络，将连接紧密的细胞聚为一类
#在meta.data中生成seurat_clusters列
seurat_object <- FindClusters(seurat_object,resolution=seq(from=0.1,to=2.0,by=0.1))
p8=clustree(seurat_object)
p8
#ggsave("08_clustree.pdf",width = 12,height = 10,plot = p8)
#ggsave("08_clustree.tiff",width=12,height=10,dpi=300,plot=p8)

#这行代码是单独计算分辨率为 1时的聚类结果（之前已经用seq生成过 1 的结果，这一步也可以省略，直接使用已有的RNA_snn_res.1列）
#seurat_object <- FindClusters(seurat_object,resolution = 1)
#p8=clustree(seurat_object)
#p8

#确定并设置最终细胞聚类方案，但是这是临时的
#seurat_object$seurat_clusters列依然以聚类最后一个分辨率为准，所以后续会产生错误
Idents(seurat_object) <- "RNA_snn_res.1"


#之前少了这关键的一步
seurat_object$seurat_clusters <- seurat_object$RNA_snn_res.1

#聚类可视化
#通过 UMAP 和 t-SNE 两种算法将高维的细胞数据映射到二维平面，直观展示细胞聚类结果
seurat_object <- RunUMAP(seurat_object,reduction = "harmony",dims = 1:20)
seurat_object <- RunTSNE(seurat_object,reduction = "harmony",dims = 1:20)
p9 <- DimPlot(seurat_object,reduction = "umap",label = T)
p10 <- DimPlot(seurat_object,reduction = "tsne",label = T)
p9+p10
p10
#ggsave("09_UMAP.pdf",width = 8,height = 6,plot = p9)
#ggsave("09_UMAP.tiff",width = 8,height = 6,dpi=300,plot = p9)
#ggsave("10_TSNE.pdf",width = 8,height = 6,plot = p10)
#ggsave("10_TSNE.tiff",width = 8,height = 6,dpi=300,plot = p10)
ggsave("new10_TSNE.tiff",width = 8,height = 6,dpi=300,plot = p10)

pall <- DimPlot(
  seurat_object, 
  reduction = "tsne", 
  group.by = "orig.ident", 
  pt.size = 0.8,
  label = F, 
  alpha = 0.5
) 
pall

# 提取DimPlot的颜色映射
color_mapping <- ggplot_build(p10.5)$data[[1]] %>% 
  dplyr::distinct(group, colour) %>%  # group对应orig.ident，colour是颜色值
  dplyr::rename(orig.ident = group, color_hex = colour)

# 打印查看
print(color_mapping)
ggsave("p10.5.tiff",width = 10,height = 8,dpi=600,plot = p10.5)
ggsave("pall.tiff",width = 10,height = 8,dpi=600,plot = p10.5)
ggsave("newp10.5.tiff",width = 10,height = 8,dpi=600,plot = p10.5)
#按时期分
#p_tsne_stage <- DimPlot(
#  seurat_object, 
#  reduction = "tsne", 
#  group.by = "orig.ident", 
#  label = T, 
#  pt.size = 0.5
#) + 
#  ggtitle("tSNE by Developmental Stage") + 
#  theme_bw()

#p_tsne_stage
#ggsave("p_tsne_stage.tiff",width = 10,height = 8,dpi=600,plot = p_tsne_stage)

# germ_cells <- WhichCells(seurat_object, expression = Ddx4 > 1 | Dazl > 1)
# 
# epi_cells <- WhichCells(seurat_object, 
#                         expression = Upk3b > 1 | Krt19 > 1)
# 
# pregran_cells <- WhichCells(seurat_object, 
#                             expression = Wnt6 > 1 | Wnt4 > 1 | Kitl > 1 | Foxl2 > 1)
# 
# mesench_cells <- WhichCells(seurat_object, expression = Nr2f2 > 1 | Col1a1 > 1)
# 
# blood_cells <- WhichCells(seurat_object, expression = Cldn5 > 1 | Car2 > 1|Lcn2>1|Cx3cr1>1)



# p_germ <- DimPlot(seurat_object, cells.highlight = germ_cells, cols.highlight = "red", 
#                   pt.size = 0.5, reduction = "tsne",label = FALSE) +
#   ggtitle("Germ populations") +
#   annotate("text", x = 10, y = 50, label = "Ddx4\nDazl", color = "red") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 
# p_epi <- DimPlot(seurat_object, cells.highlight = epi_cells, cols.highlight = "blue", 
#                  pt.size = 0.5,reduction = "tsne") + ggtitle("Epithelial populations (Upk3b/Krt19)")
# 
# p_pregran <- DimPlot(seurat_object, cells.highlight = pregran_cells, cols.highlight = "red", 
#                      pt.size = 0.5,reduction = "tsne") + ggtitle("Pregranulosa populations (Wnt6/Wnt4/Kitl/Foxl2)")
# 
# 
# p_mesench <- DimPlot(seurat_object, cells.highlight = mesench_cells, cols.highlight = "red", 
#                      pt.size = 0.5, label = FALSE) +
#   ggtitle("Mesenchymal populations") +
#   annotate("text", x = 10, y = 50, label = "Nr2f2\nCol1a1", color = "red") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# p_blood <- DimPlot(seurat_object, cells.highlight = blood_cells, cols.highlight = "red", 
#                    pt.size = 0.5, label = FALSE) +
#   ggtitle("Blood-related populations") +
#   annotate("text", x = 10, y = 50, label = "Cldn5\nCar2", color = "red") +
#   theme(plot.title = element_text(hjust = 0.5))



#五个亚群
germ_Ddx4 <- WhichCells(seurat_object, expression = Ddx4 > 1)  
germ_Dazl <- WhichCells(seurat_object, expression = Dazl > 1) 

germ_cells_list <- list(Ddx4 = germ_Ddx4, Dazl = germ_Dazl)

epi_Upk3b <- WhichCells(seurat_object,expression = Upk3b>1)
epi_Krt19 <- WhichCells(seurat_object,expression = Krt19>1)

epi_cells_list <- list(Upk3b = epi_Upk3b, Krt19 = epi_Krt19)

pregran_Wnt6 <- WhichCells(seurat_object, expression = Wnt6 > 1)
pregran_Wnt4 <- WhichCells(seurat_object, expression = Wnt4 > 1)
pregran_Kitl <- WhichCells(seurat_object, expression = Kitl > 1)
pregran_Foxl2 <- WhichCells(seurat_object, expression = Foxl2 > 1)

pregran_cells_list <- list(Wnt6 = pregran_Wnt6, 
                           Wnt4 = pregran_Wnt4, 
                           Kitl = pregran_Kitl, 
                           Foxl2 = pregran_Foxl2)

mesench_Nr2f2 <- WhichCells(seurat_object, expression = Nr2f2 > 1)
mesench_Col1a1 <- WhichCells(seurat_object, expression = Col1a1 > 1)

mesench_cells_list <- list(Nr2f2 = mesench_Nr2f2, Col1a1 = mesench_Col1a1)

blood_Cldn5 <- WhichCells(seurat_object,expression = Cldn5>1)
blood_Car2 <- WhichCells(seurat_object,expression = Car2>1)
blood_Lcn2 <- WhichCells(seurat_object,expression = Lcn2>1)
blood_Cx3cr1 <- WhichCells(seurat_object,expression = Cx3cr1>1)

blood_cells_list <- list(Cldn5=blood_Cldn5,Car2=blood_Car2,Lcn2=blood_Lcn2,Cx3cr1=blood_Cx3cr1)

p_germ <- DimPlot(seurat_object, 
                  cells.highlight = germ_cells_list,  
                  cols.highlight = c("red", "blue"),  # Ddx4绿色，Dazl紫色
                  pt.size = 0.1, 
                  reduction = "tsne",  
                  label = FALSE) +
  ggtitle("Germ populations") +
  annotate("text", x = -50, y = 50, label = "Ddx4", color = "red", size = 5) +
  annotate("text", x = -50, y = 45, label = "Dazl", color = "blue", size = 5) +
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none")

p_epi <- DimPlot(seurat_object, 
                 cells.highlight = epi_cells_list,  
                 cols.highlight = c("red", "blue"),  # Upk3b蓝色，Krt19青色
                 pt.size = 0.1, 
                 reduction = "tsne") +  
  ggtitle("Epithelial populations (Upk3b/Krt19)") +
  annotate("text", x = -49, y = 50, label = "Upk3b", color = "red", size = 5) +
  annotate("text", x = -50, y = 45, label = "Krt19", color = "blue", size = 5) +
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none")

p_pregran <- DimPlot(seurat_object, 
                     cells.highlight = pregran_cells_list, 
                     cols.highlight = c("red", "orange", "blue", "cyan"),  
                     pt.size = 0.1, 
                     reduction = "tsne") +
  ggtitle("Pregranulosa populations") +
  annotate("text", x = -50, y = 50, label = "Wnt6", color = "red", size = 5) +
  annotate("text", x = -50, y = 45, label = "Wnt4", color = "orange", size = 5) +
  annotate("text", x = -50, y = 40, label = "Kitl", color = "blue", size = 5) +
  annotate("text", x = -50, y = 35, label = "Foxl2", color = "cyan", size = 5) +
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none")

p_mesench <- DimPlot(seurat_object, 
                     cells.highlight = mesench_cells_list, 
                     cols.highlight = c("red", "blue"), 
                     pt.size = 0.1, 
                     reduction = "tsne", 
                     label = FALSE) +
  ggtitle("Mesenchymal populations") +
  annotate("text", x = -50, y = 50, label = "Nr2f2", color = "red", size = 5) +
  annotate("text", x = -50, y = 45, label = "Col1a1", color = "blue", size = 5) +
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none")

p_blood <- DimPlot(seurat_object, 
                   cells.highlight = blood_cells_list, 
                   cols.highlight = c("red", "cyan", "green", "blue"), 
                   pt.size = 0.1, 
                   reduction = "tsne", 
                   label = FALSE) +
  ggtitle("Blood-related populations") +
  annotate("text", x = -50, y = 50, label = "Lcn2", color = "red", size = 5) +
  annotate("text", x = -50, y = 45, label = "Cx3cr1", color = "cyan", size = 5) +
  annotate("text", x = -50, y = 40, label = "Cldn5", color = "green", size = 5) +
  annotate("text", x = -50, y = 35, label = "Car2", color = "blue", size = 5) +
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none")

p_germ
p_epi
p_pregran
p_mesench
p_blood

ggsave("5群newp_germ.jpg",width=8,height = 6,dpi=600,plot=p_germ)
ggsave("5群newp_epi.jpg",width=8,height = 6,dpi=600,plot=p_epi)
ggsave("5群newp_pregran.jpg",width=8,height = 6,dpi=600,plot=p_pregran)
ggsave("5群newp_mesench.jpg",width=8,height = 6,dpi=600,p_mesench)
ggsave("5群newp_blood.jpg",width=8,height = 6,dpi=600,plot=p_blood)

#7.筛选标记基因，细胞注释
#FindAllMarkers() 是 Seurat 包中专门用于批量识别每个细胞簇（Cluster）的特异性标记基因的函数
#标记基因”，是指在某一细胞簇中高表达，而在其他细胞簇中低表达的基因，这些基因是后续 “细胞类型注释” 的核心依据
#FindAllMarkers函数默认基于 “当前Idents指定的聚类结果” 计算标记基因，而非直接使用seurat_object$seurat_clusters列
allmarkers=FindAllMarkers(seurat_object,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
write.csv(allmarkers,"newallmarkers.csv")
allmarkers=read.csv("newallmarkers.csv",row.names=1,check.names = F)
#按细胞簇分组，提取每组中表达差异最大的前 10 个基因
Top10.coarse=allmarkers %>%
  group_by(cluster) %>%
  slice_max(n=10,order_by = avg_log2FC)
write.csv(Top10.coarse,"newTop10.coarse.csv")

#细胞注释
#recode()是dplyr包中用于 “按规则替换向量值” 的函数，语法为recode(原始向量, "旧值1"="新值1", "旧值2"="新值2", ...)
#在panglaodb中查看类型
#Fibroblasts成纤维细胞可能是“衰老 / 分化更成熟的间充质干细胞MSCs
# seurat_object$cell.type=recode(seurat_object$seurat_clusters,
#                                "0"="Fibroblasts ",
#                                "1"="Neurons",
#                                "2"="Fibroblasts",   
#                                "3"="Oligodendrocyte.",
#                                "4"="Unknown/Neutrophils",
#                                "5"="Germ cells",
#                                "6"="nknown/Fibroblasts/Astrocyte",
#                                "7"="Unknown/Fibroblasts",
#                                "8"="Neurons",
#                                "9"="EC",
#                                "10"="Germ cells",
#                                "11"="Unknown，EC",
#                                "12"="Fibroblasts",
#                                "13"="Unknown/Astrocytes",
#                                "14"="unknown/Enterocytes",
#                                "15"="Neurons",
#                                "16"="unknown/Neutrophils",
#                                "17"="Germ cells",
#                                "18"="Neutrophils",
#                                "19"="Pericytes",
#                                "20"="EC",
#                                "21"="unknown，Enterocytes",
#                                "22"="EC",
#                                "23"="Smooth muscle c.",
#                                "24"="unknown,Neurons",
#                                "25"="Erythroid-like ."
# )

seurat_object$cell.type=recode(seurat_object$RNA_snn_res.1,
                               "3"="Pregranulosa",
                               "5"="Pregranulosa",
                               "6"="Pregranulosa",
                          
)
table(seurat_object$cell.type)
pregran_cells <- WhichCells(seurat_object, idents = c("3", "5", "6"))  # 按聚类编号筛选目标细胞

# 绘制TSNE图，高亮Pregranulosa（3、5、6簇）
p12 <- DimPlot(seurat_object,
                         reduction = "tsne",
                         cells.highlight = pregran_cells,  # 指定需高亮的细胞
                         cols.highlight = "red",  # 高亮颜色（如红色，可自定义）
                         #cols.background = "gray90",  # 其他细胞颜色（淡化，如浅灰）
                         group.by = "cell.type",  # 仍按细胞类型分组（显示Pregranulosa标签）
                         label = T,  # 显不显示细胞类型标签（仅Pregranulosa会显示，其他为原始簇号）
                         pt.size = 0.1) +
  theme(legend.position = "none")+
  ggtitle("pregranulosa populations")

p12
#p12=DimPlot(seurat_object,group.by = "cell.type",reduction = "tsne",label = T,pt.size=0.1)
#ggsave("p12.pdf",width=8,height=6,plot=p12)
#ggsave("p12.tiff",width=10,height=6,plot=p12,dpi=300)


save.image(file = "D:/database/all/49/newall_objects.RData")
save.image(file = "D:/database/all/49/all_objects.RData")
save(seurat_object,file = "D:/database/all/49/seurat_object.RData")

load(file = "D:/database/all/49/all_objects.RData")


# ============================================================
# ============================================================
# ===========================================#

# 清空当前 R 环境，避免之前运行留下的变量影响本次分析。
rm(list = ls())

# 主动触发垃圾回收，释放内存。单细胞对象通常较大，这一步可以减少内存压力。
gc()

# 加载分析需要的 R 包。
# suppressPackageStartupMessages 可以隐藏包加载时的大量提示，使运行日志更干净。
suppressPackageStartupMessages({
  # Seurat：单细胞对象处理、细胞提取、模块打分、差异表达和可视化。
  library(Seurat)
  # dplyr/tidyr：整理元数据和统计表。
  library(dplyr)
  library(tidyr)
  # ggplot2：自定义绘图。
  library(ggplot2)
  # clusterProfiler/org.Mm.eg.db：小鼠基因 ID 转换和 GO 富集分析。
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

# 固定随机种子，保证 AddModuleScore 等含随机抽样步骤的结果可重复。
set.seed(123)

# ============================================================
# 图形输出字体与格式设置
# ============================================================
# 本项目后续所有图都统一使用 Times New Roman 字体，并保存为 SVG 矢量图。
# SVG 的优点：
#   1. 放进 PPT 或 Word 后放大不容易糊；
#   2. 文字和线条保持矢量清晰；
#   3. 后续也可以用 Adobe Illustrator、Inkscape 等软件继续编辑。
#
# 注意：
#   这里使用 svglite::svglite() 输出 SVG。
#   svglite 会在 SVG 中保留字体族信息，适合后续插入 PPT/Word 或继续编辑。
plot_font <- "Times New Roman"

# 统一封装 SVG 保存函数。
# 参数：
#   filename：输出的 SVG 文件名；
#   plot：需要保存的 ggplot/Seurat 图对象；
#   width/height：图的宽和高，单位为英寸。
save_svg <- function(filename, plot, width, height) {
  ggsave(
    filename = file.path(fig_dir, filename),
    plot = plot,
    width = width,
    height = height,
    device = function(filename, width, height, ...) {
      svglite::svglite(
        file = filename,
        width = width,
        height = height,
        bg = "white",
        system_fonts = list(
          sans = plot_font,
          serif = plot_font,
          mono = plot_font
        )
      )
    }
  )
}

save_png <- function(filename, plot, width, height) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(
      filename = file.path(fig_dir, filename),
      plot = plot,
      width = width,
      height = height,
      dpi = 600,
      device = ragg::agg_png,
      bg = "white"
    )
  } else {
    ggsave(
      filename = file.path(fig_dir, filename),
      plot = plot,
      width = width,
      height = height,
      dpi = 600,
      bg = "white"
    )
  }
}

# 设置项目输出目录。所有图表和表格都集中输出到这个文件夹，便于作为作业附件提交。
project_dir <- "E:/R/final_stromal_ecm_project"
fig_dir <- file.path(project_dir, "figures")
table_dir <- file.path(project_dir, "tables")

# 如果输出文件夹不存在，则自动创建。
# recursive = TRUE 表示如果上级目录不存在，也会一起创建。
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

# 读取已经预处理好的 Seurat 对象。
# 该对象来自前期脚本，已经包含 cell_subtype、period_renamed 等元数据。
source_object <- "D:/database/all/49/uptodate4.5_seurat_object.RData"
load(source_object)
obj <- seurat_object

# 指定默认 assay 为 RNA。后续 FetchData、DotPlot、FindMarkers 等函数会默认使用这个 assay。
DefaultAssay(obj) <- "RNA"

# 统一发育时期的顺序。
# 如果不显式设置 factor levels，R 可能按字母顺序排列时期，导致图中顺序不符合发育时间。
period_levels <- c("E11.5", "E12.5", "E14.5", "E16.5", "E18.5", "P1")
obj$period_renamed <- factor(as.character(obj$period_renamed), levels = period_levels)

# ============================================================
# 1. 统计各时期细胞类型组成
# ============================================================
# 目的：
#   先从全局角度展示当前对象中各类细胞在不同发育时期的比例。
#   这一步可以说明为什么本项目选择 Mesenchymal 基质细胞作为主要分析对象：
#   它在多个时期均有细胞，数量也足够支持后续统计比较。

cell_count_table <- obj@meta.data %>%
  # 再次确保 period_renamed 按发育顺序排序。
  mutate(period_renamed = factor(as.character(period_renamed), levels = period_levels)) %>%
  # 按时期和细胞类型计数。
  count(period_renamed, cell_subtype, name = "n_cells") %>%
  # 在每个时期内部计算总细胞数和各细胞类型所占比例。
  group_by(period_renamed) %>%
  mutate(period_total = sum(n_cells), fraction = n_cells / period_total) %>%
  ungroup()

# 保存细胞组成统计表，报告中可以引用。
write.csv(cell_count_table, file.path(table_dir, "cell_type_counts_by_period.csv"), row.names = FALSE)

# 绘制堆叠柱状图：横轴为发育时期，纵轴为不同细胞类型比例。
cell_prop_plot <- ggplot(cell_count_table, aes(x = period_renamed, y = fraction, fill = cell_subtype)) +
  geom_col(width = 0.78, color = "white", linewidth = 0.15) +
  # 将小数比例显示为百分比。
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
  # 手动指定颜色，保证报告和 PPT 中颜色稳定一致。
  scale_fill_manual(values = c(
    "Epithelial" = "#4C78A8",
    "Mesenchymal" = "#59A14F",
    "Other" = "#BAB0AC",
    "Pregranulosa" = "#F28E2B"
  )) +
  labs(x = "Developmental stage", y = "Fraction of cells", fill = "Cell subtype") +
  theme_classic(base_size = 12, base_family = plot_font) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

# 保存细胞组成图，SVG 矢量格式。
save_svg("fig1_cell_type_fraction_by_period.svg", cell_prop_plot, width = 8, height = 5)

# ============================================================
# 2. 提取 Mesenchymal 基质细胞并设置早晚期分组
# ============================================================
# 目的：
#   本项目的科学问题集中在卵巢基质细胞微环境。
#   因此后续所有模块打分、差异表达和富集分析都只在 Mesenchymal 细胞中进行。

# 从总对象中提取 cell_subtype 为 Mesenchymal 的细胞。
# 同时限定时期为 period_levels 中已有的 E11.5-P1。
mes <- subset(obj, subset = cell_subtype == "Mesenchymal" & period_renamed %in% period_levels)
mes$period_renamed <- factor(as.character(mes$period_renamed), levels = period_levels)

# 构造早晚期分组，用于差异表达分析。
# 早期：E11.5、E12.5、E14.5，代表卵巢发育较早、基质细胞尚在建立阶段。
# 晚期：E16.5、E18.5、P1，代表卵泡形成和围产期卵巢结构成熟阶段。
mes$phase_group <- ifelse(mes$period_renamed %in% c("E11.5", "E12.5", "E14.5"),
                          "Early_E11.5_E14.5", "Late_E16.5_P1")
mes$phase_group <- factor(mes$phase_group, levels = c("Early_E11.5_E14.5", "Late_E16.5_P1"))

# 统计每个时期的基质细胞数量，用于报告中说明样本覆盖情况。
mes_count_table <- mes@meta.data %>%
  count(period_renamed, phase_group, name = "n_cells") %>%
  arrange(period_renamed)
write.csv(mes_count_table, file.path(table_dir, "mesenchymal_counts_by_period.csv"), row.names = FALSE)

# 如果对象中已经存在 tSNE 降维结果，则直接绘制基质细胞在 tSNE 空间中的时期分布。
# 这里不重新 RunTSNE，是为了复用前期分析结果，并避免引入新的随机差异。
if ("tsne" %in% Reductions(mes)) {
  tsne_plot <- DimPlot(mes, reduction = "tsne", group.by = "period_renamed", pt.size = 0.35) +
    labs(title = NULL, color = "Stage") +
    theme_classic(base_size = 12, base_family = plot_font)
  save_svg("fig2_mesenchymal_tsne_by_period.svg", tsne_plot, width = 6.6, height = 5.4)
}

# ============================================================
# 3. 构建 ECM/TGF-beta/基质重塑/黏附相关基因集
# ============================================================
# 目的：
#   与其只观察单个基因，不如把同一生物学过程中的多个基因合并成 gene set，
#   然后计算每个细胞的模块分数。这样可以更稳定地反映一个功能程序的整体活性。
#
# 四类基因集的含义：
#   ECM_core：
#     细胞外基质核心成分，如胶原、纤连蛋白、纤维蛋白、基质蛋白聚糖等。
#   TGFb_signaling：
#     TGF-beta/BMP 相关配体、受体、Smad 转录因子和潜在调控因子。
#   Matrix_remodeling：
#     基质降解、交联、修饰和平衡调控相关基因，如 Mmp、Timp、Lox、Plod 等。
#   Adhesion_integrin：
#     细胞和 ECM 之间的黏附、整合素受体和胶原受体相关基因。

gene_sets_raw <- list(
  ECM_core = c("Col1a1", "Col1a2", "Col3a1", "Col5a1", "Col5a2", "Col6a1", "Col6a2", "Col6a3",
               "Fn1", "Fbn1", "Fbn2", "Lum", "Dcn", "Bgn", "Sparc", "Postn", "Tnc", "Vcan",
               "Mgp", "Lama2", "Lamb1", "Nid1", "Nid2"),
  TGFb_signaling = c("Tgfb1", "Tgfb2", "Tgfb3", "Tgfbr1", "Tgfbr2", "Tgfbr3", "Smad2", "Smad3",
                     "Smad4", "Smad5", "Smad6", "Smad7", "Ltbp1", "Ltbp2", "Ltbp3", "Thbs1",
                     "Thbs2", "Bmp2", "Bmp4", "Bmpr1a", "Bmpr2", "Id1", "Id2", "Id3"),
  Matrix_remodeling = c("Mmp2", "Mmp14", "Mmp19", "Adamts1", "Adamts2", "Adamts5", "Adamts9",
                        "Timp1", "Timp2", "Timp3", "Plod1", "Plod2", "Lox", "Loxl1", "Loxl2"),
  Adhesion_integrin = c("Itga1", "Itga5", "Itga6", "Itga8", "Itgav", "Itgb1", "Itgb5", "Cd44",
                        "Sdc1", "Sdc4", "Cdh11", "Alcam", "Ddr1", "Ddr2")
)

# 只保留在当前 Seurat 对象中真实存在的基因。
# 这样可以避免 AddModuleScore 因部分基因不存在而报错，也能记录每个基因集实际参与分析的基因。
gene_sets <- lapply(gene_sets_raw, function(x) intersect(x, rownames(mes)))

# 生成基因集检测情况表：
#   n_input 表示人工输入的基因数；
#   n_detected 表示在当前数据中能找到的基因数；
#   detected_genes 列出实际用于打分的基因。
gene_set_table <- data.frame(
  gene_set = names(gene_sets),
  n_input = vapply(gene_sets_raw, length, integer(1)),
  n_detected = vapply(gene_sets, length, integer(1)),
  detected_genes = vapply(gene_sets, paste, character(1), collapse = ";")
)
write.csv(gene_set_table, file.path(table_dir, "curated_gene_sets_detected.csv"), row.names = FALSE)

# ============================================================
# 4. 使用 AddModuleScore 计算单细胞模块分数
# ============================================================
# AddModuleScore 的基本思想：
#   对一个目标基因集，计算这些基因在每个细胞中的平均表达；
#   再减去一组表达水平相近的随机对照基因的平均表达；
#   得到的分数越高，说明该细胞越倾向于表达这个功能基因程序。
#
# 参数说明：
#   features：传入基因集；
#   name：生成的元数据列名前缀；
#   ctrl = 50：每个目标基因匹配 50 个对照基因；
#   seed = 123：保证对照基因抽样可重复。
# 注意：
#   Seurat 会自动在 name 后面加数字 1，所以最终列名为 ECM_core_score1 等。

mes <- AddModuleScore(mes, features = gene_sets["ECM_core"], name = "ECM_core_score", ctrl = 50, seed = 123)
mes <- AddModuleScore(mes, features = gene_sets["TGFb_signaling"], name = "TGFb_signaling_score", ctrl = 50, seed = 123)
mes <- AddModuleScore(mes, features = gene_sets["Matrix_remodeling"], name = "Matrix_remodeling_score", ctrl = 50, seed = 123)
mes <- AddModuleScore(mes, features = gene_sets["Adhesion_integrin"], name = "Adhesion_integrin_score", ctrl = 50, seed = 123)

# 保存四个模块分数列名，方便后面批量统计和作图。
score_cols <- c("ECM_core_score1", "TGFb_signaling_score1", "Matrix_remodeling_score1", "Adhesion_integrin_score1")

# 给模块分数列名设置更适合图表展示的标签。
score_labels <- c(
  ECM_core_score1 = "ECM core",
  TGFb_signaling_score1 = "TGF-beta signaling",
  Matrix_remodeling_score1 = "Matrix remodeling",
  Adhesion_integrin_score1 = "Adhesion / integrin"
)

# 将宽格式元数据转换成长格式。
# 原始格式：每个细胞一行，四个模块分数分别占四列。
# 长格式：每个细胞的每个模块分数单独一行，适合 ggplot facet_wrap 作图和批量统计检验。
score_long <- mes@meta.data %>%
  dplyr::select(period_renamed, phase_group, all_of(score_cols)) %>%
  pivot_longer(cols = all_of(score_cols), names_to = "score_name", values_to = "score") %>%
  mutate(score_label = factor(score_labels[score_name], levels = unname(score_labels)))

# 按“模块分数类型”和“发育时期”汇总均值、中位数、标准差和标准误。
# mean_score 用于趋势图，median_score 用于理解小提琴图中的整体位置。
score_summary <- score_long %>%
  group_by(score_name, score_label, period_renamed) %>%
  summarise(
    n_cells = n(),
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    se_score = sd_score / sqrt(n_cells),
    .groups = "drop"
  )
write.csv(score_summary, file.path(table_dir, "module_score_summary_by_period.csv"), row.names = FALSE)

# ============================================================
# 5. 模块分数的统计检验
# ============================================================
# 由于单细胞模块分数通常不是严格正态分布，且不同阶段细胞数不均衡，
# 这里使用非参数检验更稳妥。
#
# Kruskal-Wallis 检验：
#   用于判断某个模块分数在所有发育时期之间是否存在总体差异。
# BH 校正：
#   对四个模块分数的总体检验 P 值进行多重检验校正，控制假阳性。

kw_table <- score_long %>%
  group_by(score_name, score_label) %>%
  summarise(
    kruskal_p = kruskal.test(score ~ period_renamed)$p.value,
    .groups = "drop"
  ) %>%
  mutate(kruskal_p_adj = p.adjust(kruskal_p, method = "BH"))
write.csv(kw_table, file.path(table_dir, "module_score_kruskal_by_period.csv"), row.names = FALSE)

# 两两 Wilcoxon 检验：
#   在 Kruskal-Wallis 显示总体有差异后，进一步比较任意两个时期之间是否显著不同。
# Holm 校正：
#   对两两比较的 P 值进行更严格的多重检验校正。
pairwise_table <- bind_rows(lapply(score_cols, function(sc) {
  # 取出当前模块分数的数据。
  dat <- score_long %>% filter(score_name == sc)
  # exact = FALSE 用于避免大量重复值或大样本时精确检验计算过慢/报错。
  pw <- pairwise.wilcox.test(dat$score, dat$period_renamed, p.adjust.method = "holm", exact = FALSE)
  # pairwise.wilcox.test 输出是矩阵，这里转换成三列表格：group1、group2、校正后 P 值。
  as.data.frame(as.table(pw$p.value)) %>%
    filter(!is.na(Freq)) %>%
    transmute(score_name = sc, score_label = score_labels[sc], group1 = Var1, group2 = Var2, p_adj = Freq)
}))
write.csv(pairwise_table, file.path(table_dir, "module_score_pairwise_wilcox.csv"), row.names = FALSE)

# ============================================================
# 6. 模块分数可视化
# ============================================================
# 小提琴图：
#   展示每个时期每个细胞的模块分数分布。
#   小提琴越向上，说明该时期更多细胞具有较高的功能程序分数。
#   中间叠加箱线图，用于显示中位数和四分位数。

score_violin <- ggplot(score_long, aes(x = period_renamed, y = score, fill = period_renamed)) +
  geom_violin(scale = "width", trim = TRUE, color = "grey35", linewidth = 0.2) +
  geom_boxplot(width = 0.14, outlier.size = 0.15, linewidth = 0.22, fill = "white", alpha = 0.88) +
  facet_wrap(~ score_label, ncol = 2, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Developmental stage", y = "Module score") +
  theme_classic(base_size = 12, base_family = plot_font) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "#F2F2F2", color = NA),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )
save_svg("fig3_module_scores_violin.svg", score_violin, width = 9, height = 6.8)

# 趋势图：
#   使用每个时期的平均模块分数和标准误，展示 ECM/TGF-beta 等功能程序随发育时间的变化趋势。
#   这张图适合 PPT 中讲“E16.5 是基质重塑增强节点”。
score_trend <- ggplot(score_summary, aes(x = period_renamed, y = mean_score, group = score_label, color = score_label)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.1) +
  geom_errorbar(aes(ymin = mean_score - se_score, ymax = mean_score + se_score), width = 0.15, linewidth = 0.4) +
  scale_color_manual(values = c("#4C78A8", "#F28E2B", "#59A14F", "#B07AA1")) +
  labs(x = "Developmental stage", y = "Mean module score", color = "Gene program") +
  theme_classic(base_size = 12, base_family = plot_font) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_svg("fig4_module_score_trends.svg", score_trend, width = 7.8, height = 5.4)

# ============================================================
# 7. 关键基因表达展示
# ============================================================
# 目的：
#   模块分数可以说明整体功能程序变化，但报告中还需要展示代表性基因。
#   因此选取 ECM、TGF-beta、基质重塑和细胞黏附中比较典型的基因，
#   计算其在每个时期基质细胞中的平均表达，并绘制热图和气泡图。
#
# 关键基因举例：
#   Col1a1/Col1a2/Col3a1：胶原，代表 ECM 主体成分；
#   Dcn/Lum/Sparc/Fbn1/Fn1：基质组织、纤维结构和细胞外支架相关；
#   Mmp2/Timp1/Timp2：基质降解和抑制平衡；
#   Tgfb1/Tgfb2/Tgfbr1/Tgfbr3/Smad3/Bmpr2：TGF-beta/BMP 信号轴；
#   Itga5/Itgav/Itgb1/Ddr1：细胞和 ECM 之间的黏附或胶原识别。

key_genes <- c("Col1a1", "Col1a2", "Col3a1", "Fn1", "Fbn1", "Lum", "Dcn", "Sparc",
               "Mmp2", "Timp1", "Timp2", "Tgfb1", "Tgfb2", "Tgfbr1", "Tgfbr3",
               "Smad3", "Bmpr2", "Itga5", "Itgav", "Itgb1", "Ddr1")

# 只保留当前数据中存在的基因，避免 FetchData 报错。
key_genes <- intersect(key_genes, rownames(mes))

# FetchData 从 Seurat 对象中提取元数据和基因表达。
# 这里提取 period_renamed 和关键基因的归一化表达值。
expr_data <- FetchData(mes, vars = c("period_renamed", key_genes))

# 计算每个时期的平均表达，并对每个基因在各时期之间做 Z-score 标准化。
# Z-score 的意义：
#   同一个基因在不同时间点的相对高低。
#   红色表示该基因在该时期相对更高，蓝色表示相对更低。
# 注意：
#   Z-score 不是绝对表达量，不能用于比较不同基因谁更高，只能比较同一基因在不同时期的变化。
avg_expr <- expr_data %>%
  group_by(period_renamed) %>%
  summarise(across(all_of(key_genes), mean), .groups = "drop") %>%
  pivot_longer(cols = all_of(key_genes), names_to = "gene", values_to = "avg_expression") %>%
  group_by(gene) %>%
  mutate(z_score = as.numeric(scale(avg_expression))) %>%
  ungroup()
write.csv(avg_expr, file.path(table_dir, "key_gene_average_expression_by_period.csv"), row.names = FALSE)

# 绘制关键基因热图。
# 横轴为发育时期，纵轴为代表性基因，颜色为每个基因跨时期标准化后的平均表达。
heatmap_plot <- ggplot(avg_expr, aes(x = period_renamed, y = factor(gene, levels = rev(key_genes)), fill = z_score)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_gradient2(low = "#3B6EA8", mid = "white", high = "#C43C39", midpoint = 0, name = "Scaled\nmean") +
  labs(x = "Developmental stage", y = "Gene") +
  theme_minimal(base_size = 12, base_family = plot_font) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )
save_svg("fig5_key_gene_heatmap.svg", heatmap_plot, width = 7.2, height = 6.2)

# 绘制气泡图。
# DotPlot 同时展示两个信息：
#   颜色：该基因在对应时期中的平均表达；
#   点大小：表达该基因的细胞比例。
# 这比单纯热图多了“有多少细胞表达该基因”的信息。
dot_plot <- DotPlot(mes, features = key_genes, group.by = "period_renamed") +
  scale_color_gradient(low = "#E8EEF4", high = "#9E2F2F") +
  labs(x = "Gene", y = "Developmental stage", color = "Average\nexpression", size = "Percent\nexpressed") +
  theme_classic(base_size = 11, base_family = plot_font) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_svg("fig6_key_gene_dotplot.svg", dot_plot, width = 9.5, height = 4.8)

# ============================================================
# 8. 早期 vs 晚期基质细胞差异表达分析
# ============================================================
# 目的：
#   模块分数和关键基因图说明了 ECM/TGF-beta 程序变化。
#   差异表达分析进一步从全转录组层面回答：
#     晚期基质细胞相对于早期基质细胞，哪些基因显著升高？
#     这些基因是否集中在 ECM、TGF-beta、迁移、发育等生物学过程？
#
# 分组：
#   ident.1 = Late_E16.5_P1 表示晚期；
#   ident.2 = Early_E11.5_E14.5 表示早期；
#   avg_log2FC > 0 表示晚期更高；
#   avg_log2FC < 0 表示早期更高。

# 将 Seurat 的身份标签设置为早晚期分组，FindMarkers 会根据 Idents 进行比较。
Idents(mes) <- mes$phase_group

# 使用 Wilcoxon 秩和检验做差异表达分析。
# min.pct = 0.10：至少在任一组 10% 细胞中表达的基因才参与检验，过滤过低表达基因。
# logfc.threshold = 0.10：先用较宽松阈值保留更多基因，后面再按 0.25 判断上调方向。
de <- FindMarkers(
  mes,
  ident.1 = "Late_E16.5_P1",
  ident.2 = "Early_E11.5_E14.5",
  test.use = "wilcox",
  min.pct = 0.10,
  logfc.threshold = 0.10
)

# 将行名中的基因名转成普通列，便于写出 CSV 和后续筛选。
de$gene <- rownames(de)

# 为差异表达结果增加几个辅助列：
#   direction：
#     Late_up 表示晚期显著上调；
#     Early_up 表示早期显著上调；
#     Not_significant 表示不满足阈值。
#   neg_log10_adj：
#     火山图纵轴，调整后 P 值越小，该值越大。
#   technical_gene：
#     标记核糖体、线粒体和血红蛋白相关基因。
#     这些基因常反映测序深度、线粒体比例或红细胞污染等技术因素，
#     所以做 GO 富集前单独过滤，避免生物学解释被技术信号主导。
de <- de %>%
  mutate(
    direction = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.25 ~ "Late_up",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Early_up",
      TRUE ~ "Not_significant"
    ),
    neg_log10_adj = -log10(pmax(p_val_adj, 1e-300)),
    technical_gene = grepl("^(Rpl|Rps|mt-|Hbb|Hba)", gene)
  ) %>%
  arrange(p_val_adj)

# 保存完整差异表达表，保留技术相关基因，方便追溯所有结果。
write.csv(de, file.path(table_dir, "mesenchymal_late_vs_early_DE.csv"), row.names = FALSE)

# 去除核糖体、线粒体、血红蛋白基因后的差异表达表。
# 后续 GO 富集使用这个版本，更适合解释基质细胞真实生物学变化。
de_bio <- de %>%
  filter(!technical_gene)
write.csv(de_bio, file.path(table_dir, "mesenchymal_late_vs_early_DE_without_ribo_mt_hb.csv"), row.names = FALSE)

# 选择一批和本项目主题相关的代表性基因，在火山图中标注。
# intersect 的作用是只保留差异表达结果中实际存在的基因。
volcano_labels <- intersect(c("Col1a1", "Col3a1", "Fn1", "Fbn1", "Sparc", "Mmp2", "Timp1",
                              "Timp2", "Tgfb2", "Tgfbr1", "Tgfbr3", "Smad3", "Itgb1",
                              "Itga5", "Ddr1", "Lum", "Dcn", "Mgp", "Thbs1"), de$gene)

# 绘制火山图。
# 横轴 avg_log2FC：
#   越靠右，表示晚期相对早期表达越高；
#   越靠左，表示早期相对晚期表达越高。
# 纵轴 -log10 adjusted P：
#   越高，表示统计显著性越强。
volcano_plot <- ggplot(de, aes(x = avg_log2FC, y = neg_log10_adj, color = direction)) +
  geom_point(alpha = 0.55, size = 0.85) +
  # 虚线标出 log2FC 阈值和 FDR = 0.05 阈值，帮助读者判断显著上调区域。
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", linewidth = 0.25, color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.25, color = "grey50") +
  geom_text(data = de %>% filter(gene %in% volcano_labels),
            aes(label = gene), size = 3, family = plot_font, check_overlap = TRUE, vjust = -0.6, color = "black") +
  scale_color_manual(values = c("Late_up" = "#C43C39", "Early_up" = "#3B6EA8", "Not_significant" = "#B8B8B8")) +
  labs(x = "avg_log2FC (Late E16.5-P1 vs Early E11.5-E14.5)", y = "-log10 adjusted P", color = "DE class") +
  theme_classic(base_size = 12, base_family = plot_font)
save_svg("fig7_late_vs_early_volcano.svg", volcano_plot, width = 7.3, height = 5.6)

# ============================================================
# 9. GO Biological Process 富集分析
# ============================================================
# 目的：
#   差异表达只能告诉我们“哪些基因变了”。
#   GO 富集进一步回答“这些变化的基因集中在哪些生物学过程”。
#
# 这里分别对：
#   late_up：晚期显著上调基因；
#   early_up：早期显著上调基因；
# 做 GO Biological Process 富集。

# 从去除技术基因后的差异表达表中提取晚期上调和早期上调基因。
late_up <- de_bio %>% filter(direction == "Late_up") %>% pull(gene)
early_up <- de_bio %>% filter(direction == "Early_up") %>% pull(gene)

# 设置 GO 富集的背景基因集 universe。
# 这里使用在基质细胞中至少 10 个细胞表达过的基因作为背景。
# 这样比使用全基因组更合理，因为富集分析的背景应该接近“本实验中有机会被检测到的基因”。
expressed_universe <- rownames(mes)[Matrix::rowSums(GetAssayData(mes, assay = "RNA", layer = "data") > 0) >= 10]

# 简化 GO Biological Process 条目名称，用于图8/图9的 y 轴标签。
# 完整 GO 名称仍然保存在 late_up_GO_BP.csv 和 early_up_GO_BP.csv 中。
# 这里的简化风格参考报告表格中的写法，尽量保留核心含义，同时减少左侧文字过长的问题。
simplify_go_label <- function(x) {
  label_map <- c(
    "response to transforming growth factor beta" = "TGF-beta response",
    "cellular response to transforming growth factor beta stimulus" = "cellular TGF-beta response",
    "extracellular matrix organization" = "ECM organization",
    "external encapsulating structure organization" = "encapsulating structure",
    "extracellular structure organization" = "extracellular structure",
    "skeletal system development" = "skeletal development",
    "response to wounding" = "wound response",
    "ameboidal-type cell migration" = "ameboidal migration",
    "positive regulation of cell migration" = "cell migration (+)",
    "gland development" = "gland development",
    "small GTPase-mediated signal transduction" = "small GTPase signaling",
    "cell surface receptor protein serine/threonine kinase signaling pathway" = "Ser/Thr kinase signaling",
    "RNA splicing" = "RNA splicing",
    "RNA splicing, via transesterification reactions" = "RNA splicing (transesterif.)",
    "RNA splicing, via transesterification reactions with bulged adenosine as nucleophile" = "RNA splicing (adenosine)",
    "mRNA splicing, via spliceosome" = "mRNA splicing (spliceosome)",
    "mRNA processing" = "mRNA processing",
    "chromosome segregation" = "chromosome segregation",
    "chromosome organization" = "chromosome organization",
    "sister chromatid segregation" = "sister chromatid seg.",
    "mitotic sister chromatid segregation" = "mitotic chromatid seg.",
    "mitotic nuclear division" = "mitotic nuclear division",
    "nuclear chromosome segregation" = "nuclear chromosome seg.",
    "ribonucleoprotein complex biogenesis" = "RNP complex biogenesis"
  )
  out <- unname(label_map[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

# 自定义一个 GO 富集函数，避免晚期和早期重复写两遍相同代码。
# 参数：
#   genes：需要做富集的基因符号向量；
#   out_prefix：输出文件名前缀，比如 late_up 或 early_up。
run_go <- function(genes, out_prefix) {
  # clusterProfiler 的 enrichGO 需要 Entrez ID。
  # bitr 用于把小鼠基因 SYMBOL 转换为 ENTREZID。
  ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  universe_ids <- bitr(expressed_universe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

  # 如果可转换的基因太少，则不做富集，直接返回空表。
  # 这可以避免 enrichGO 因输入基因数不足而报错。
  if (nrow(ids) < 5) {
    return(data.frame())
  }

  # enrichGO 参数说明：
  #   gene：待富集基因 Entrez ID；
  #   universe：背景基因；
  #   ont = "BP"：分析 Biological Process；
  #   pAdjustMethod = "BH"：Benjamini-Hochberg 多重检验校正；
  #   pvalueCutoff/qvalueCutoff：筛选显著富集条目；
  #   readable = TRUE：把结果中的 Entrez ID 转回更易读的基因符号。
  ego <- enrichGO(
    gene = unique(ids$ENTREZID),
    universe = unique(universe_ids$ENTREZID),
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )

  # 转换为普通 data.frame 并保存完整 GO 结果。
  ego_df <- as.data.frame(ego)
  write.csv(ego_df, file.path(table_dir, paste0(out_prefix, "_GO_BP.csv")), row.names = FALSE)

  # 如果存在显著 GO 条目，则绘制前 12 个条目的柱状图。
  # 横轴 Count 表示该 GO 条目命中的差异基因数量；
  # 颜色 -log10(FDR) 表示显著性，颜色越深越显著。
  if (nrow(ego_df) > 0) {
    top_df <- ego_df %>%
      arrange(p.adjust) %>%
      slice_head(n = 12) %>%
      mutate(Description_short = simplify_go_label(Description))
    write.csv(
      top_df[, c("Description", "Description_short", "p.adjust", "Count", "geneID")],
      file.path(table_dir, paste0(out_prefix, "_GO_BP_top12_plot_labels.csv")),
      row.names = FALSE
    )
    top_df <- top_df %>%
      mutate(Description_short = factor(Description_short, levels = rev(Description_short)))
    p <- ggplot(top_df, aes(x = Count, y = Description_short, fill = -log10(p.adjust))) +
      geom_col(width = 0.72) +
      scale_fill_gradient(low = "#9CC7D8", high = "#9E2F2F", name = "-log10\nFDR") +
      labs(x = "Gene count", y = NULL) +
      theme_classic(base_size = 11, base_family = plot_font) +
      theme(
        axis.text.y = element_text(size = 10.5),
        axis.text.x = element_text(size = 10.5),
        axis.title.x = element_text(size = 11),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9.5),
        plot.margin = margin(5.5, 8, 5.5, 5.5)
      )
    save_svg(paste0("fig8_", out_prefix, "_GO_BP.svg"), p, width = 8, height = 5.2)
    save_png(paste0("fig8_", out_prefix, "_GO_BP.png"), p, width = 8, height = 5.2)
  }

  # 返回 GO 结果，便于在 R 环境中继续查看。
  ego_df
}

# 分别运行晚期上调和早期上调基因的 GO 富集分析。
late_go <- run_go(late_up, "late_up")
early_go <- run_go(early_up, "early_up")

# ============================================================
# 10. 保存分析摘要、Seurat 对象和运行环境信息
# ============================================================
# 目的：
#   让报告具有可重复性。
#   analysis_summary.txt 记录核心数字；
#   mesenchymal_ecm_analysis_object.RData 保存分析后的基质细胞对象；
#   sessionInfo.txt 记录 R 版本和包版本，方便以后复现。

summary_lines <- c(
  paste0("source_object=", source_object),
  paste0("total_cells_after_preprocessing=", ncol(obj)),
  paste0("mesenchymal_cells=", ncol(mes)),
  paste0("mesenchymal_period_counts=", paste(paste(mes_count_table$period_renamed, mes_count_table$n_cells, sep = ":"), collapse = ";")),
  paste0("late_up_genes_nontechnical=", length(late_up)),
  paste0("early_up_genes_nontechnical=", length(early_up)),
  paste0("ECM_detected_genes=", paste(gene_sets$ECM_core, collapse = ";")),
  paste0("TGFb_detected_genes=", paste(gene_sets$TGFb_signaling, collapse = ";"))
)
writeLines(summary_lines, file.path(table_dir, "analysis_summary.txt"))

# 保存已经添加模块分数和早晚期分组的 Mesenchymal Seurat 对象。
# 后续如果想继续做基质细胞亚群分析或配体-受体通讯分析，可以从这个对象开始。
save(mes, file = file.path(project_dir, "mesenchymal_ecm_analysis_object.RData"))

# 保存 R 运行环境
session_info <- capture.output(sessionInfo())
writeLines(session_info, file.path(table_dir, "sessionInfo.txt"))

# 在控制台输出完成提示。
cat("Analysis completed.\n")
cat("Project directory:", project_dir, "\n")


# ============================================================

##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())
library(dplyr)
library(SeuratObject)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(DoubletFinder)
library(ggplot2)
set.seed(1)

setwd("F:\\Bioinformatics course\\data\\data")

# 定义函数
seurat_standard_normalize_and_scale <- function(colon, cluster = FALSE, cluster_resolution = NULL){
  colon <- NormalizeData(colon, normalization.method = "LogNormalize", scale.factor = 10000)
  colon <- FindVariableFeatures(colon, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(colon)
  colon <- ScaleData(colon, features = all.genes)
  colon <- RunPCA(colon, features = VariableFeatures(object = colon))
  if (cluster){
    colon <- FindNeighbors(colon, dims = 1:20)
    colon <- FindClusters(colon, resolution = cluster_resolution)
  }
  colon <- RunUMAP(colon, dims = 1:20)
  return(colon)
}

make_seurat_object_and_doublet_removal <- function(data_directory, project_name){
  if (!dir.exists(data_directory)) {
    stop(paste("Data directory does not exist:", data_directory))
  }
  
  colon.data <- Read10X(data.dir = data_directory)
  
  cat("Data dimensions:", dim(colon.data), "\n")
  cat("First few gene names:", head(rownames(colon.data)), "\n")
  
  currentSample <- CreateSeuratObject(counts = colon.data, project = project_name, min.cells = 3, min.features = 40)
  
  cat("\nSeurat object created:\n")
  print(currentSample)
  
  # 计算线粒体基因百分比（小鼠：^mt- 不区分大小写）
  mt.genes <- grep("^MT-", rownames(currentSample), value = TRUE, ignore.case = TRUE)
  
  if (length(mt.genes) > 0) {
    cat("Found mitochondrial genes:", length(mt.genes), "\n")
    currentSample[["percent.mt"]] <- PercentageFeatureSet(
      currentSample, 
      features = mt.genes
    )
  } else {
    mt.genes2 <- grep("^mt-", rownames(currentSample), value = TRUE, ignore.case = TRUE)
    if (length(mt.genes2) > 0) {
      cat("Found mitochondrial genes (alternative pattern):", length(mt.genes2), "\n")
      currentSample[["percent.mt"]] <- PercentageFeatureSet(
        currentSample, 
        features = mt.genes2
      )
    } else {
      warning("No mitochondrial genes found. Setting percent.mt to 0.")
      currentSample[["percent.mt"]] <- 0
    }
  }
  
  # qc plot - 过滤前
  pdf(paste0("./qc_plots_", project_name, "_prefiltered.pdf"))
  print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05))
  dev.off()
  
  pdf(paste0("./qc_plots_", project_name, "_prefiltered_no_points.pdf"))
  print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0))
  dev.off()
  
  # 过滤
  cat("\nFiltering cells...\n")
  cat("Before filtering:", ncol(currentSample), "cells\n")
  
  currentSample <- subset(currentSample, 
                          subset = nFeature_RNA > 300 & 
                            nCount_RNA > 1000 & 
                            percent.mt < 20)
  
  cat("After filtering:", ncol(currentSample), "cells\n")
  
  if (ncol(currentSample) < 100) {
    warning("Too few cells after filtering. Using relaxed criteria.")
    currentSample <- subset(currentSample, 
                            subset = nFeature_RNA > 200 & 
                              nCount_RNA > 500 & 
                              percent.mt < 30)
    cat("After relaxed filtering:", ncol(currentSample), "cells\n")
  }
  
  # ========== 新增：保存过滤后的QC图 ==========
  pdf(paste0("./qc_plots_", project_name, "_filtered.pdf"))
  print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05))
  dev.off()
  
  pdf(paste0("./qc_plots_", project_name, "_filtered_no_points.pdf"))
  print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0))
  dev.off()
  # ===========================================
  
  # Normalize and make UMAP
  currentSample <- seurat_standard_normalize_and_scale(currentSample, cluster = FALSE)
  
  # Run doublet finder
  if (ncol(currentSample) < 100) {
    warning("Too few cells for DoubletFinder. Skipping doublet detection.")
    currentSample$doublet.class <- "Singlet"
  } else {
    cat("\nRunning DoubletFinder...\n")
    
    sweep.res.list <- paramSweep(currentSample, PCs = 1:20, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    
    pK_value <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    cat("Optimal pK value:", pK_value, "\n")
    
    nExp_poi <- round(0.075 * ncol(currentSample))
    cat("Estimated doublets to remove:", nExp_poi, "\n")
    
    currentSample <- doubletFinder(
      currentSample, 
      PCs = 1:20, 
      pN = 0.25, 
      pK = pK_value, 
      nExp = nExp_poi, 
      reuse.pANN = NULL, 
      sct = FALSE
    )
    
    df_class_col <- grep("^DF.classifications", colnames(currentSample@meta.data), value = TRUE)[1]   # 取第一个
    pann_col <- grep("^pANN", colnames(currentSample@meta.data), value = TRUE)[1]
    
    if (!is.na(df_class_col)) {
      currentSample$doublet.class <- currentSample[[df_class_col]]
      currentSample[[df_class_col]] <- NULL
    } else {
      warning("No doublet classifications found.")
      currentSample$doublet.class <- "Singlet"
    }
    
    if (!is.na(pann_col)) {
      currentSample$pANN <- currentSample[[pann_col]]
      currentSample[[pann_col]] <- NULL
    }
    
    # ========== 新增：绘制双细胞分布图（过滤前） ==========
    pdf(paste0("./UMAP_doublet_distribution_", project_name, ".pdf"))
    print(DimPlot(currentSample, reduction = "umap", group.by = "doublet.class",
                  cols = c("Singlet" = "#D51F26", "Doublet" = "#377EB8")))
    dev.off()
    # ===================================================
    
    # 过滤双细胞
    if ("doublet.class" %in% colnames(currentSample@meta.data)) {
      singlet_count <- sum(currentSample$doublet.class == "Singlet")
      doublet_count <- sum(currentSample$doublet.class == "Doublet")
      cat("Singlets:", singlet_count, "Doublets:", doublet_count, "\n")
      
      currentSample <- subset(currentSample, subset = doublet.class == "Singlet")
    }
  }
  
  # plot results (过滤后UMAP)
  pdf(paste0("./UMAP_post_double_removal_", project_name, ".pdf"))
  print(DimPlot(currentSample, reduction = "umap", cols = c("#D51F26")))
  dev.off()
  
  # 清理对象
  currentSample <- DietSeurat(currentSample, 
                              counts = TRUE, 
                              data = TRUE, 
                              scale.data = FALSE,
                              assays = "RNA")
  
  if ("scale.data" %in% Layers(currentSample, assay = "RNA")) {
    currentSample[["RNA"]]$scale.data <- NULL
  }
  
  return(currentSample)
}

seurat_qc_plots <- function(colon, sample_name){
  if (ncol(colon) > 0) {
    pdf(paste0("./seurat_nFeature_plots_", sample_name, ".pdf"), width = 40, height = 15)
    print(VlnPlot(colon, features = c("nFeature_RNA"), ncol = 1, pt.size = 0.2))
    dev.off()
    
    pdf(paste0("./seurat_nCount_plots_", sample_name, ".pdf"), width = 40, height = 15)
    print(VlnPlot(colon, features = c("nCount_RNA"), ncol = 1, pt.size = 0.2))
    dev.off()
    
    if ("percent.mt" %in% colnames(colon@meta.data)) {
      pdf(paste0("./seurat_pMT_plots_", sample_name, ".pdf"), width = 40, height = 15)
      print(VlnPlot(colon, features = c("percent.mt"), ncol = 1, pt.size = 0.2))
      dev.off()
    }
  }
}

# 开始分析
data_directory <- c("F:\\Bioinformatics course\\data\\data/E11/", 
                    "F:\\Bioinformatics course\\data\\data/E12/",
                    "F:\\Bioinformatics course\\data\\data/E14/",
                    "F:\\Bioinformatics course\\data\\data/E16/",
                    "F:\\Bioinformatics course\\data\\data/E18/",
                    "F:\\Bioinformatics course\\data\\data/P1/")
project_name <- c("E11.5", "E12.5","E14.5","E16.5","E18.5","P1")
samples <- project_name

sample1 <- make_seurat_object_and_doublet_removal(data_directory[1], samples[1])

### 多个样本合并  
seu_list <- sample1
for (i in 2:length(samples)){
  sc.i = make_seurat_object_and_doublet_removal(data_directory[i], samples[i])
  seu_list=merge(seu_list,sc.i)
}
table(seu_list$orig.ident)
scRNA_harmony=seu_list
scRNA_harmony  <- NormalizeData(scRNA_harmony ) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
library(harmony)
scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")
###问题 harmony之后的数据在哪里？
###一定要指定harmony
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution =1)

scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:20)

DimPlot(scRNA_harmony , reduction = "umap",label = T) 
DimPlot(scRNA_harmony, reduction = "umap", split.by ='orig.ident')

DimPlot(scRNA_harmony, reduction = "umap", group.by='orig.ident')
table(scRNA_harmony$orig.ident)  

# 你已经有了成功的Harmony结果：scRNA_harmony
# 检查它的状态
print(scRNA_harmony)

# 如果还没有聚类，进行聚类
if (!"seurat_clusters" %in% colnames(scRNA_harmony@meta.data)) {
  scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:20)
  scRNA_harmony <- FindClusters(scRNA_harmony, resolution = 1.0)
}

# 如果还没有UMAP，运行UMAP
if (!"umap" %in% names(scRNA_harmony@reductions)) {
  scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:20)
}

# 可视化
DimPlot(scRNA_harmony, reduction = "umap", group.by = "orig.ident")
DimPlot(scRNA_harmony, reduction = "umap", label = TRUE)

 # 保存结果
saveRDS(scRNA_harmony, file = "final_harmony_integrated.rds")









# ==============================================================================
# 0. 载入核心单细胞与绘图包
# ==============================================================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)

# 设置工作目录（请根据您的实际路径修改）
# setwd("F:/Bioinformatics course/data/data")

# ==============================================================================
# 1. 读取整合后的 Harmony RDS 文件
# ==============================================================================
cat("\n[Step 1] Loading integrated Seurat object from RDS...\n")
scRNA_harmony <- readRDS("final_harmony_integrated.rds")

# ==============================================================================
# 2. 转换与修复数据格式 (解决 Seurat v5 专属的多样本分层问题)
# ==============================================================================
cat("\n[Step 2] Repairing Seurat v5 layers...\n")
DefaultAssay(scRNA_harmony) <- "RNA"
scRNA_harmony <- JoinLayers(scRNA_harmony) # 将千层饼矩阵合并，确保下游计算正常

# ==============================================================================
# 3. 配置文献（Niu & Spradling, 2020）定义的五大类标准 Marker 基因集
# ==============================================================================
# 严格采用小鼠标准大小写格式
paper_markers <- list(
  "Germ cells"          = c("Ddx4", "Dazl"),
  "Epithelial cells"    = c("Upk3b", "Krt19"),
  "Pregranulosa cells"  = c("Wnt4", "Wnt6", "Kitl", "Foxl2"), 
  "Mesenchymal cells"   = c("Nr2f2", "Col1a1"),
  "Blood-related cells" = c("Cldn5", "Car2", "Lcn2", "Cx3cr1")
)

# ==============================================================================
# 4. 极致顶刊主题配置 (【安全修复】：不指定 family 参数，彻底告别字体罢工)
# ==============================================================================
journal_theme <- theme_bw() + theme(
  panel.grid.major = element_blank(),                       # 移除主要网格线，突出数据点
  panel.grid.minor = element_blank(),                       # 移除次要网格线
  panel.border = element_rect(fill = NA, color = "black", linewidth = 1.2), # 加粗外边框
  plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 10)),
  axis.text = element_text(size = 11, color = "black"),    # 坐标轴标签纯黑，避免印刷发灰
  axis.title = element_text(size = 12, face = "bold", color = "black"),
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 11, face = "bold"),
  axis.ticks = element_line(color = "black", linewidth = 0.8)
)

# ==============================================================================
# 5. 绘制并保存文献 Marker 气泡图 (采用顶刊双色渐变色阶)
# ==============================================================================
cat("\n[Step 5] Generating reference marker DotPlot...\n")
flat_markers <- unique(unlist(paper_markers, use.names = FALSE))
flat_markers <- intersect(flat_markers, rownames(scRNA_harmony)) # 过滤矩阵未检测基因

pdf("Paper_Markers_Check_DotPlot.pdf", width = 9, height = 6)
p_dot <- DotPlot(scRNA_harmony, features = flat_markers, group.by = "seurat_clusters", assay = "RNA",
                 cols = c("#E6E6E6", "#B2182B"), dot.scale = 6) + # 优雅的低表达灰到高表达红
  coord_flip() + # 颠倒横纵轴，基因作为行，Cluster作为列，符合标准图表规范
  labs(title = "Lineage-Specific Reference Markers Expression", 
       x = "Marker Genes (Niu & Spradling, 2020)", 
       y = "Seurat Clusters") +
  journal_theme
print(p_dot)
dev.off()

# ==============================================================================
# 6. 注释策略：文献标准双轨映射 (全自动打分映射 + 预留精准手动微调接口)
# ==============================================================================
cat("\n[Step 6] Annotating cell types based on paper lineage signatures...\n")
cluster_averages <- AverageExpression(scRNA_harmony, features = flat_markers, group.by = "seurat_clusters", slot = "data")$RNA

score_matrix <- sapply(names(paper_markers), function(cell_type) {
  genes <- paper_markers[[cell_type]]
  genes_present <- intersect(genes, rownames(cluster_averages))
  if(length(genes_present) > 1) {
    colMeans(cluster_averages[genes_present, , drop = FALSE])
  } else if(length(genes_present) == 1) {
    cluster_averages[genes_present, ]
  } else {
    rep(0, ncol(cluster_averages))
  }
})

auto_matches <- colnames(score_matrix)[max.col(score_matrix)]
names(auto_matches) <- rownames(score_matrix)

# ---- 精准手动微调映射接口 ----
manual_cluster_map <- auto_matches  
# 💡 建议：打开生成的 Paper_Markers_Check_DotPlot.pdf 后，如有个别群错配，在此行下方人工修正：
# manual_cluster_map["5"] <- "Pregranulosa cells"

scRNA_harmony$paper_cell_type <- manual_cluster_map[as.character(scRNA_harmony$seurat_clusters)]

# ==============================================================================
# 7. 降维空间检查与联合计算 (确保同时具备 t-SNE 空间)
# ==============================================================================
if (!"tsne" %in% names(scRNA_harmony@reductions)) {
  cat("\n[Step 7] Running t-SNE based on harmony space...\n")
  scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:20)
}

# ==============================================================================
# 8. 降维可视化输出：输出高级学术配色、符合审稿规范的 PDF 图像
# ==============================================================================
cat("\n[Step 8] Exporting dimension reduction plots...\n")

# 顶刊离散型学术色板（高辨识度红、蓝、绿、橙、浅蓝组合）
pub_colors <- c("#3C5488FF", "#4DBBD5FF", "#00A087FF", "#F39B7FFF", "#8491B4FF")

# ---- 8.1 文献细胞分群 UMAP 图 ----
pdf("UMAP_Ovary_Paper_Annotation.pdf", width = 8, height = 6)
p1 <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "paper_cell_type", 
              label = TRUE, repel = TRUE, label.size = 4, cols = pub_colors) +
  labs(title = "Ovarian Cell Types (UMAP Space)", x = "UMAP_1", y = "UMAP_2") +
  journal_theme
print(p1)
dev.off()

# ---- 8.2 文献细胞分群 t-SNE 图 ----
pdf("TSNE_Ovary_Paper_Annotation.pdf", width = 8, height = 6)
p2 <- DimPlot(scRNA_harmony, reduction = "tsne", group.by = "paper_cell_type", 
              label = TRUE, repel = TRUE, label.size = 4, cols = pub_colors) +
  labs(title = "Ovarian Cell Types (t-SNE Space)", x = "t-SNE_1", y = "t-SNE_2") +
  journal_theme
print(p2)
dev.off()

# ==============================================================================
# 9. 保存带有纯正文献注释标签的终版 RDS 对象
# ==============================================================================
cat("\n[Step 9] Archiving final paper-annotated RDS object...\n")
saveRDS(scRNA_harmony, file = "final_harmony_integrated_paper_annotated.rds")
cat("\n======================================================================\n")
cat("🎉 Success! All pipelines completed with NO errors.\n")
cat("======================================================================\n")
