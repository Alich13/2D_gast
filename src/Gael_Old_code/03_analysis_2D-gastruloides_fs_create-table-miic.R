library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(mgsa)
library(topGO)

# Set up the directory name
directory <- "D:/Gael/sc_rna_seq"

# Load the datasets
dataset.fs <- readRDS(file = paste(directory,"/rds/2dgas_t36_fs.rds",sep=""))
data.raw.read <- read.table(paste(directory,'/raw_data/data_2dgas_fs/RawReadMatrix.txt',sep=""))
data.coverage <- read.table(paste(directory,'/raw_data/data_2dgas_fs/TranscriptEstimationMatrix-Coverage.txt',sep=""))
data.transcript.v2 <- read.table(paste(directory,'/raw_data/data_2dgas_fs/TranscriptEstimationMatrix-v2.txt',sep="")) 
data.transcript.v3 <- read.table(paste(directory,'/raw_data/data_2dgas_fs/CountEstimationMatrix-v3.txt',sep="")) 
data.transcript.v4 <- read.table(paste(directory,'/raw_data/data_2dgas_fs/CountEstimationMatrix-v4.txt',sep="")) 
data.transcript.v5 <- read.table(paste(directory,'/raw_data/data_2dgas_fs/TranscriptEstimationMatrix-v5.txt',sep="")) 
data.transcript.v6 <- read.table(paste(directory,'/raw_data/data_2dgas_fs/TranscriptEstimationMatrix-v6.txt',sep="")) 

cells.use <- substring(Cells(dataset.fs), first = regexpr("_",Cells(dataset.fs))+1)

# Plot the expression of one gene
gene <- "5430430B14Rik"
FeaturePlot(dataset.fs,gene,order = TRUE,label=TRUE)
VlnPlot(dataset.fs,gene)
FeaturePlot(dataset.fs,gene,order = TRUE,label=TRUE,split.by = "condition")
VlnPlot(dataset.fs,gene,group.by = "condition")

# Arbitrary list of TF 
list.gene.use.1 <- sort(c("NANOG","NKX1-2","OTX2","POU3F1","POU5F1","SOX2","UTF1",
                   "EOMES","MIXL1","T","SP5",
                   "FOXA1","FOXA2","GSC","HHEX","LHX1","SOX17","SOX21",
                   "FOXJ1","MAK","NOTO",
                   "FOXC1","FOXC2","MEIS2","SNAI1","TBX6","TWIST1",
                   "GATA4","GATA6","MESP1","MESP2","MSX1",
                   "FOXF1","HOXB1","ISL1","PITX2","PRDM6","TBX3",
                   "ETV2","GATA2","GATA3","HAND1","LMO2","MSX2","SMAD1","TBX20",
                   "CDX2","CITED1","NMI","PARP9","PRDM1","PRDM14","SMAD6","SMAD7",
                   "BMP4","ID1","ID2","ID3",
                   "AXIN2","LEF1","WNT3",
                   "LEFTY1","LEFTY2","NODAL","TDGF1"))

filename <- paste(directory,"/miic/2Dgas-t36-fs_65-TFs/2Dgas-t36-fs_65-TFs_matrix_all.csv",sep="")
matrix.use.1 <- as.data.frame(t(data.raw.read[list.gene.use.1,cells.use]))
matrix.use.1$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use.1,filename,row.names = FALSE)

filename <- paste(directory,"/miic/2Dgas-t36-fs_65-TFs/2Dgas-t36-fs_65-TFs_matrix_WNT.csv",sep="")
matrix.WNT <- matrix.use.1[dataset.fs$condition == "WNT",]
write.csv(matrix.WNT ,filename,row.names = FALSE)

filename <- paste(directory,"/miic/2Dgas-t36-fs_65-TFs/2Dgas-t36-fs_65-TFs_matrix_BMP.csv",sep="")
matrix.BMP <- matrix.use.1[dataset.fs$condition == "BMP",]
write.csv(matrix.BMP,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use.1),
                                 var_type = c(rep(1,ncol(matrix.use.1)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use.1)))
filename <- paste(directory,"/miic/2Dgas-t36-fs_65-TFs/2Dgas-t36-fs_65-TFs_category.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use.1),
                                 var_type = c(rep(1,ncol(matrix.use.1)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use.1)),
                                 is_contextual = c(rep(0,ncol(matrix.use.1)-1),1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_65-TFs/2Dgas-t36-fs_65-TFs_category_all_with-contextual.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

# Import GO data
library(org.Mm.eg.db)
x <- org.Mm.eg.db
goall <- data.frame(AnnotationDbi::select(x, keys(x), "GOALL"))
symbol <- AnnotationDbi::select(x, keys(x), "SYMBOL")

#  transcription regulator activity : GO:0140110
#  canonical Wnt signaling pathway : GO:0060070,
#  BMP signaling pathway : GO:0030509,
#  activin receptor signaling pathway : GO:0032924,
#  transforming growth factor beta receptor binding : GO:0005160
#  fibroblast growth factor receptor signaling pathway : GO:0008543 
#  vascular endothelial growth factor signaling pathway : GO:0038084
#  retinoic acid receptor signaling pathway : GO:0048384
#  smoothened signaling pathway : GO:0007224 
#  Notch signaling pathway : GO:0007219

#  extra-cellular region : GO:0005576

# Heatmap of top TF
GO.tf <- "GO:0140110"
entrezid <- unique(na.omit(goall[goall$GOALL==GO.tf,])$ENTREZID)
GO.genes.tf <- symbol[symbol$ENTREZID %in% entrezid,]$SYMBOL
gene.list.tf <- sort(intersect(toupper(GO.genes.tf),row.names(dataset.fs[['RNA']]@counts)))

markers.all.tf <- FindAllMarkers(dataset.fs,only.pos = TRUE,features = gene.list.tf, min.diff.pct = 0.2)

# Export data 
markers.all.tf <- markers.all.tf[markers.all.tf$p_val_adj<0.01,]
top.tf<- markers.all.tf %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
other.tf <- c("ID3","LEF1","POU5F1","PRDM1","PRDM14","T","FOXH1",
              "SMAD1","SMAD2","SMAD3","SMAD4","SMAD5","SMAD6","SMAD7",
              "AXIN2","BMP4","FGF5","FGF8","LEFTY1","LEFTY2","NODAL","TDGF1","TGFB1","WNT3")

# 222 TF
list.gene.use.2 <- sort(union(unique(top.tf$gene),other.tf))

filename <- paste(directory,"/miic/2Dgas-t36-fs_222-TFs/2Dgas-t36-fs_222-TFs_matrix_all.csv",sep="")
matrix.use.2 <- as.data.frame(t(data.raw.read[list.gene.use.2,cells.use]))
matrix.use.2$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use.2,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use.2),
                                 var_type = c(rep(1,ncol(matrix.use.2)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use.2)),
                                 is_contextual = c(rep(0,ncol(matrix.use.2)-1),1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_222-TFs/2Dgas-t36-fs_222-TFs_category_all_with-contextual.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

# 215 TF (with normalized data)
list.gene.use.2 <- sort(union(unique(top.tf$gene),other.tf))

filename <- paste(directory,"/miic/2Dgas-t36-fs_215-TFs/2Dgas-t36-fs_215-TFs_normalized_matrix_all.csv",sep="")
matrix.use.2 <- as.data.frame(t(dataset.fs@assays[["RNA"]]@data[list.gene.use.2,]))
matrix.use.2$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use.2,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use.2),
                                 var_type = c(rep(1,ncol(matrix.use.2)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use.2)),
                                 is_contextual = c(rep(0,ncol(matrix.use.2)-1),1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_215-TFs/2Dgas-t36-fs_215-TFs_normalized_category_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

# WNT only
filename <- paste(directory,"/miic/2Dgas-t36-fs_222-TFs/2Dgas-t36-fs_222-TFs_matrix_WNT.csv",sep="")
matrix.WNT <- matrix.use.2[dataset.fs$condition == "WNT",]
variables.to.change <- names(which(apply(matrix.WNT, 2, function(x) length(unique(x)))<3))
matrix.WNT[,variables.to.change] <- 0
write.csv(matrix.WNT ,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.WNT),
                                 var_type = c(rep(1,ncol(matrix.WNT)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.WNT)))
category.order.use[category.order.use$var_names %in% variables.to.change,]$var_type <- 0
filename <- paste(directory,"/miic/2Dgas-t36-fs_222-TFs/2Dgas-t36-fs_222-TFs_category_WNT.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

# BMP only
filename <- paste(directory,"/miic/2Dgas-t36-fs_222-TFs/2Dgas-t36-fs_222-TFs_matrix_BMP.csv",sep="")
matrix.BMP <- matrix.use.2[dataset.fs$condition == "BMP",]
variables.to.change <- names(which(apply(matrix.BMP, 2, function(x) length(unique(x)))<3))
matrix.BMP[,variables.to.change] <- 0
write.csv(matrix.BMP ,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.BMP),
                                 var_type = c(rep(1,ncol(matrix.BMP)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.BMP)))
category.order.use[category.order.use$var_names %in% variables.to.change,]$var_type <- 0
filename <- paste(directory,"/miic/2Dgas-t36-fs_222-TFs/2Dgas-t36-fs_222-TFs_category_BMP.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

# TF miic - 1 

list.gene.use.3 <- sort(union(list.gene.use.2,VariableFeatures(dataset.fs)))
list.tf.3 <- sort(union(list.gene.use.2, intersect(list.gene.use.3,gene.list.tf)))
is.consequence.3 <- c(rep(1,length(list.gene.use.3)),0)
is.consequence.3[list.gene.use.3 %in% list.tf.3] <- 0

filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-1/2Dgas-t36-fs_TF-miic-1_matrix_all.csv",sep="")
matrix.use.3 <- as.data.frame(t(data.raw.read[list.gene.use.3,cells.use]))
matrix.use.3$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use.3,filename,row.names = FALSE)

variables.to.change <- names(which(apply(matrix.use.3, 2, function(x) length(unique(x)))<3))

category.order.use <- data.frame(var_names=colnames(matrix.use.3),
                                 var_type = c(rep(1,ncol(matrix.use.3)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use.3)),
                                 is_contextual = c(rep(0,ncol(matrix.use.3)-1),1),
                                 is_consequence = is.consequence.3)
category.order.use[category.order.use$var_names %in% variables.to.change,]$var_type <- 0

filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-1/2Dgas-t36-fs_TF-miic-1_category_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)


# TF miic - 2

dataset.fs <- FindVariableFeatures(dataset.fs, nfeatures = 1000)
list.gene.use.4 <- sort(union(list.gene.use.1,VariableFeatures(dataset.fs)))
list.tf.4 <- sort(union(list.gene.use.1, intersect(list.gene.use.4,gene.list.tf)))
is.consequence.4 <- c(rep(1,length(list.gene.use.4)),0)
is.consequence.4[list.gene.use.4 %in% list.tf.4] <- 0

filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-2/2Dgas-t36-fs_TF-miic-2_matrix_all.csv",sep="")
matrix.use.4 <- as.data.frame(t(data.raw.read[list.gene.use.4,cells.use]))
matrix.use.4$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use.4,filename,row.names = FALSE)

variables.to.change <- names(which(apply(matrix.use.4, 2, function(x) length(unique(x)))<3))

category.order.use <- data.frame(var_names=colnames(matrix.use.4),
                                 var_type = c(rep(1,ncol(matrix.use.4)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use.4)),
                                 is_contextual = c(rep(0,ncol(matrix.use.4)-1),1),
                                 is_consequence = is.consequence.4)
category.order.use[category.order.use$var_names %in% variables.to.change,]$var_type <- 0

filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-2/2Dgas-t36-fs_TF-miic-2_category_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

### TF miic - 2 - TF only ###

#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-2_TF-only/2Dgas-t36-fs_TF-miic-2_TF-only_matrix_all.csv",sep="")
matrix.use.5 <- as.data.frame(t(data.raw.read[list.tf.4,cells.use]))
matrix.use.5$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use.5,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use.5),
                                 var_type = c(rep(1,ncol(matrix.use.5)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use.5)),
                                 is_contextual = c(rep(0,ncol(matrix.use.5)-1),1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-2_TF-only/2Dgas-t36-fs_TF-miic-2_TF-only_category_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

#All (with counts)

filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-2_TF-only/2Dgas-t36-fs_TF-miic-2_TF-only-with-counts_matrix_all.csv",sep="")
matrix.use.5 <- as.data.frame(t(data.raw.read[list.tf.4,cells.use]))
matrix.use.5$counts <- as.numeric(colSums(data.raw.read[,cells.use]))
matrix.use.5$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use.5,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use.5),
                                 var_type = c(rep(1,ncol(matrix.use.5)-2),1,0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use.5)),
                                 is_contextual = c(rep(0,ncol(matrix.use.5)-2),1,1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-2_TF-only/2Dgas-t36-fs_TF-miic-2_TF-only-with-counts_category_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

# WNT only
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-2_TF-only/2Dgas-t36-fs_TF-miic-2_TF-only_matrix_WNT.csv",sep="")
matrix.WNT <- matrix.use.5[dataset.fs$condition == "WNT",]
variables.to.change <- names(which(apply(matrix.WNT, 2, function(x) length(unique(x)))<3))
matrix.WNT[,variables.to.change] <- 0
write.csv(matrix.WNT ,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.WNT),
                                 var_type = c(rep(1,ncol(matrix.WNT)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.WNT)))
category.order.use[category.order.use$var_names %in% variables.to.change,]$var_type <- 0
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-2_TF-only/2Dgas-t36-fs_TF-miic-2_TF-only_category_WNT.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

# BMP only
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-2_TF-only/2Dgas-t36-fs_TF-miic-2_TF-only_matrix_BMP.csv",sep="")
matrix.BMP <- matrix.use.5[dataset.fs$condition == "BMP",]
variables.to.change <- names(which(apply(matrix.BMP, 2, function(x) length(unique(x)))<3))
matrix.BMP[,variables.to.change] <- 0
write.csv(matrix.BMP ,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.BMP),
                                 var_type = c(rep(1,ncol(matrix.BMP)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.BMP)))
category.order.use[category.order.use$var_names %in% variables.to.change,]$var_type <- 0
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-2_TF-only/2Dgas-t36-fs_TF-miic-2_TF-only_category_BMP.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

### TF miic - 3 (use MI to select genes)

## Using MI (from count) with condition and using raw data
MI_condition <- readRDS(paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/MI-condition.rds",sep=""))
names(MI_condition) %in% gene.list.tf
MI_condition_top_TF <- sort(names(MI_condition)[names(MI_condition) %in% gene.list.tf][1:200])

#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-condition_top200-TF_matrix_all.csv",sep="")
matrix.use.6 <- as.data.frame(t(data.raw.read[MI_condition_top_TF,cells.use]))
matrix.use.6$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use.6,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use.6),
                                 var_type = c(rep(1,ncol(matrix.use.6)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use.6)),
                                 is_contextual = c(rep(0,ncol(matrix.use.6)-1),1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-condition_top200-TF_category_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

## Using MI (from count) with condition, adding all counts and using raw data
MI_condition <- readRDS(paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/MI-condition.rds",sep=""))
names(MI_condition) %in% gene.list.tf
MI_condition_top_TF <- sort(names(MI_condition)[names(MI_condition) %in% gene.list.tf][1:200])

#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-condition_top200-TF-with-counts_matrix_all.csv",sep="")
matrix.use.6 <- as.data.frame(t(data.raw.read[MI_condition_top_TF,cells.use]))
matrix.use.6$counts <- as.numeric(colSums(data.raw.read[,cells.use]))
matrix.use.6$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use.6,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use.6),
                                 var_type = c(rep(1,ncol(matrix.use.6)-2),1,0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use.6)),
                                 is_contextual = c(rep(0,ncol(matrix.use.6)-2),1,1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-condition_top200-TF-with-counts_category_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)


## Using MI (from data) with condition and using raw data
MI_condition <- readRDS(paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/MI-data-condition.rds",sep=""))
names(MI_condition) %in% gene.list.tf
MI_condition_top_TF <- sort(names(MI_condition)[names(MI_condition) %in% gene.list.tf][1:200])

#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-data-condition_top200-TF_matrix_all.csv",sep="")
matrix.use.6 <- as.data.frame(t(data.raw.read[MI_condition_top_TF,cells.use]))
matrix.use.6$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use.6,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use.6),
                                 var_type = c(rep(1,ncol(matrix.use.6)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use.6)),
                                 is_contextual = c(rep(0,ncol(matrix.use.6)-1),1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-data-condition_top200-TF_category_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

## Using MI (from data) with condition and using normalized data
MI_condition <- readRDS(paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/MI-data-condition.rds",sep=""))
names(MI_condition) %in% gene.list.tf
MI_condition_top_TF <- sort(names(MI_condition)[names(MI_condition) %in% gene.list.tf][1:200])

#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-data-condition_top200-TF_matrix-normalized_all.csv",sep="")
matrix.use.6 <- as.data.frame(t(dataset.fs@assays[["RNA"]]@data[MI_condition_top_TF,]))
matrix.use.6$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use.6,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use.6),
                                 var_type = c(rep(1,ncol(matrix.use.6)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use.6)),
                                 is_contextual = c(rep(0,ncol(matrix.use.6)-1),1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-data-condition_top200-TF_category-normalized_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)


## Using MI (from data) with condition, adding nCount and using normalized data
MI_condition <- readRDS(paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/MI-data-condition.rds",sep=""))
names(MI_condition) %in% gene.list.tf
MI_condition_top_TF <- sort(names(MI_condition)[names(MI_condition) %in% gene.list.tf][1:200])

#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-data-condition_top200-TF-with-counts_matrix-normalized_all.csv",sep="")
matrix.use.6 <- as.data.frame(t(dataset.fs@assays[["RNA"]]@data[MI_condition_top_TF,]))
matrix.use.6$counts <- as.numeric(dataset.fs$nCount_RNA)
matrix.use.6$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use.6,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use.6),
                                 var_type = c(rep(1,ncol(matrix.use.6)-2),1,0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use.6)),
                                 is_contextual = c(rep(0,ncol(matrix.use.6)-2),1,1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-data-condition_top200-TF-with-counts_category-normalized_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)


## Using MI (from data) with T and using raw data
MI_condition <- readRDS(paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/MI-data-T.rds",sep=""))
MI_condition <- MI_condition[MI_condition>0]
names(MI_condition) %in% gene.list.tf
MI_condition_top_TF <- sort(c("T",names(MI_condition)[names(MI_condition) %in% gene.list.tf][1:199]))

#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-data-T_top200-TF_matrix_all.csv",sep="")
matrix.use.6 <- as.data.frame(t(data.raw.read[MI_condition_top_TF,cells.use]))
matrix.use.6$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use.6,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use.6),
                                 var_type = c(rep(1,ncol(matrix.use.6)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use.6)),
                                 is_contextual = c(rep(0,ncol(matrix.use.6)-1),1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-data-T_top200-TF_category_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)


## Using MI (from data) with T and using normalized data
MI_condition <- readRDS(paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/MI-data-T.rds",sep=""))
MI_condition <- MI_condition[MI_condition>0]
names(MI_condition) %in% gene.list.tf
MI_condition_top_TF <- sort(c("T",names(MI_condition)[names(MI_condition) %in% gene.list.tf][1:199]))

#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-data-T_top200-TF_matrix-normalized_all.csv",sep="")
matrix.use.6 <- as.data.frame(t(dataset.fs@assays[["RNA"]]@data[MI_condition_top_TF,]))
matrix.use.6$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use.6,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use.6),
                                 var_type = c(rep(1,ncol(matrix.use.6)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use.6)),
                                 is_contextual = c(rep(0,ncol(matrix.use.6)-1),1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-data-T_top200-TF_category-normalized_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

## Using MI (from data) with T and using normalized v2 data
MI_condition <- readRDS(paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/MI-data-T.rds",sep=""))
MI_condition <- MI_condition[MI_condition>0]
names(MI_condition) %in% gene.list.tf
MI_condition_top_TF <- sort(c("T",names(MI_condition)[names(MI_condition) %in% gene.list.tf][1:199]))

#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-data-T_top200-TF_matrix-normalized-v2_all.csv",sep="")
matrix.use.6 <- as.data.frame(t(data.transcript.v2[MI_condition_top_TF,cells.use]))
matrix.use.6$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use.6,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use.6),
                                 var_type = c(rep(1,ncol(matrix.use.6)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use.6)),
                                 is_contextual = c(rep(0,ncol(matrix.use.6)-1),1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-data-T_top200-TF_category-normalized-v2_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

## Using MI (from data) with T and using normalized v3 data
MI_condition <- readRDS(paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/MI-data-T.rds",sep=""))
MI_condition <- MI_condition[MI_condition>0]
names(MI_condition) %in% gene.list.tf
MI_condition_top_TF <- sort(c("T",names(MI_condition)[names(MI_condition) %in% gene.list.tf][1:199]))

#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-data-T_top200-TF_matrix-normalized-v3_all.csv",sep="")
matrix.use.6 <- as.data.frame(t(data.transcript.v3[MI_condition_top_TF,cells.use]))
matrix.use.6$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use.6,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use.6),
                                 var_type = c(rep(1,ncol(matrix.use.6)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use.6)),
                                 is_contextual = c(rep(0,ncol(matrix.use.6)-1),1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-data-T_top200-TF_category-normalized-v3_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

## Using MI (from data) with T and using normalized v4 data
MI_condition <- readRDS(paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/MI-data-T.rds",sep=""))
MI_condition <- MI_condition[MI_condition>0]
names(MI_condition) %in% gene.list.tf
MI_condition_top_TF <- sort(c("T",names(MI_condition)[names(MI_condition) %in% gene.list.tf][1:199]))

#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-data-T_top200-TF_matrix-normalized-v4_all.csv",sep="")
matrix.use.6 <- as.data.frame(t(data.transcript.v4[MI_condition_top_TF,cells.use]))
matrix.use.6$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use.6,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use.6),
                                 var_type = c(rep(1,ncol(matrix.use.6)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use.6)),
                                 is_contextual = c(rep(0,ncol(matrix.use.6)-1),1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-data-T_top200-TF_category-normalized-v4_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

## Using MI (from data) with T and using normalized v5 data
MI_condition <- readRDS(paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/MI-data-T.rds",sep=""))
MI_condition <- MI_condition[MI_condition>0]
names(MI_condition) %in% gene.list.tf
MI_condition_top_TF <- sort(c("T",names(MI_condition)[names(MI_condition) %in% gene.list.tf][1:199]))

#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-data-T_top200-TF_matrix-normalized-v5_all.csv",sep="")
matrix.use.6 <- as.data.frame(t(data.transcript.v5[MI_condition_top_TF,cells.use]))
matrix.use.6$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use.6,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use.6),
                                 var_type = c(rep(1,ncol(matrix.use.6)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use.6)),
                                 is_contextual = c(rep(0,ncol(matrix.use.6)-1),1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-data-T_top200-TF_category-normalized-v5_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)


## Using MI (from data) with T and using normalized v6 data
MI_condition <- readRDS(paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/MI-data-T.rds",sep=""))
MI_condition <- MI_condition[MI_condition>0]
names(MI_condition) %in% gene.list.tf
MI_condition_top_TF <- sort(c("T",names(MI_condition)[names(MI_condition) %in% gene.list.tf][1:199]))

#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-data-T_top200-TF_matrix-normalized-v6_all.csv",sep="")
matrix.use.6 <- as.data.frame(t(data.transcript.v6[MI_condition_top_TF,cells.use]))
matrix.use.6$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use.6,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use.6),
                                 var_type = c(rep(1,ncol(matrix.use.6)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use.6)),
                                 is_contextual = c(rep(0,ncol(matrix.use.6)-1),1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-data-T_top200-TF_category-normalized-v6_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)



## Using MI with T
MI_condition <- readRDS(paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/MI-T.rds",sep=""))
names(MI_condition) %in% gene.list.tf
MI_condition_top_TF <- sort(c("T",names(MI_condition)[names(MI_condition) %in% gene.list.tf][1:199]))

#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-T_top200-TF_matrix_all.csv",sep="")
matrix.use.6 <- as.data.frame(t(data.raw.read[MI_condition_top_TF,cells.use]))
matrix.use.6$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use.6,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use.6),
                                 var_type = c(rep(1,ncol(matrix.use.6)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use.6)),
                                 is_contextual = c(rep(0,ncol(matrix.use.6)-1),1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-T_top200-TF_category_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

## Using MI with POU5F1
MI_condition <- readRDS(paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/MI-POU5F1.rds",sep=""))
names(MI_condition) %in% gene.list.tf
MI_condition_top_TF <- sort(c("POU5F1",names(MI_condition)[names(MI_condition) %in% gene.list.tf][1:199]))

#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-POU5F1_top200-TF_matrix_all.csv",sep="")
matrix.use.6 <- as.data.frame(t(data.raw.read[MI_condition_top_TF,cells.use]))
matrix.use.6$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use.6,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use.6),
                                 var_type = c(rep(1,ncol(matrix.use.6)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use.6)),
                                 is_contextual = c(rep(0,ncol(matrix.use.6)-1),1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-POU5F1_top200-TF_category_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)



## With coverage and using MI of T (and counts)
MI_condition <- readRDS(paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/MI-T.rds",sep=""))
names(MI_condition) %in% gene.list.tf
MI_condition_top_TF <- sort(c("T",names(MI_condition)[names(MI_condition) %in% gene.list.tf][1:199]))
#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-T_top200-TF_coverage-with-counts_all_matrix-all.csv",sep="")
matrix.use <- as.data.frame(t(data.coverage[MI_condition_top_TF,cells.use]))
matrix.use$counts <- as.numeric(colSums(data.coverage[,cells.use]))
matrix.use$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use),
                                 var_type = c(rep(1,ncol(matrix.use)-2),1,0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use)),
                                 is_contextual = c(rep(0,ncol(matrix.use)-2),1,1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_TF-miic-3_MI/2Dgas-t36-fs_TF-miic-3_MI-T_top200-TF_coverage-with-counts_all_category_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)


## With coverage and using top TF (and counts)
list.gene.use <- sort(union(unique(top.tf$gene),other.tf))
#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_215-TFs/2Dgas-t36-fs_215-TFs_coverage-with-counts_matrix_all.csv",sep="")
matrix.use <- as.data.frame(t(data.coverage[list.gene.use,cells.use]))
matrix.use$counts <- as.numeric(colSums(data.coverage[,cells.use]))
matrix.use$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use),
                                 var_type = c(rep(1,ncol(matrix.use)-2),1,0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use)),
                                 is_contextual = c(rep(0,ncol(matrix.use)-2),1,1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_215-TFs/2Dgas-t36-fs_215-TFs_coverage-with-counts_category_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)


## With coverage (normalized) and using top TF
list.gene.use <- sort(union(unique(top.tf$gene),other.tf))
#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_215-TFs/2Dgas-t36-fs_215-TFs_coverage-normalized_matrix_all.csv",sep="")
matrix.use <- as.data.frame(t(data.transcript.v6[list.gene.use,cells.use]))
matrix.use$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use),
                                 var_type = c(rep(1,ncol(matrix.use)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use)),
                                 is_contextual = c(rep(0,ncol(matrix.use)-1),1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_215-TFs/2Dgas-t36-fs_215-TFs_coverage-normalized_category_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

## With raw counts and using top TF
list.gene.use <- sort(union(unique(top.tf$gene),other.tf))
#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_215-TFs/2Dgas-t36-fs_215-TFs_matrix_all.csv",sep="")
matrix.use <- as.data.frame(t(data.raw.read[list.gene.use,cells.use]))
matrix.use$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use),
                                 var_type = c(rep(1,ncol(matrix.use)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use)),
                                 is_contextual = c(rep(0,ncol(matrix.use)-1),1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_215-TFs/2Dgas-t36-fs_215-TFs_category_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)


## With raw counts and using top TF (and counts)
list.gene.use <- sort(union(unique(top.tf$gene),other.tf))
#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_215-TFs/2Dgas-t36-fs_215-TFs_with-counts_matrix_all.csv",sep="")
matrix.use <- as.data.frame(t(data.raw.read[list.gene.use,cells.use]))
matrix.use$counts <- as.numeric(colSums(data.coverage[,cells.use]))
matrix.use$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use),
                                 var_type = c(rep(1,ncol(matrix.use)-2),1,0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use)),
                                 is_contextual = c(rep(0,ncol(matrix.use)-2),1,1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_215-TFs/2Dgas-t36-fs_215-TFs_with-counts_category_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

## With normalized v4 and using top TF
list.gene.use <- sort(union(unique(top.tf$gene),other.tf))
#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_215-TFs/2Dgas-t36-fs_215-TFs_matrix-normalized-v4_all.csv",sep="")
matrix.use <- as.data.frame(t(data.transcript.v4[list.gene.use,cells.use]))
matrix.use$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use),
                                 var_type = c(rep(1,ncol(matrix.use)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use)),
                                 is_contextual = c(rep(0,ncol(matrix.use)-1),1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_215-TFs/2Dgas-t36-fs_215-TFs_category-normalized-v4_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)


## With normalized v5 and using top TF
list.gene.use <- sort(union(unique(top.tf$gene),other.tf))
#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_215-TFs/2Dgas-t36-fs_215-TFs_matrix-normalized-v5_all.csv",sep="")
matrix.use <- as.data.frame(t(data.transcript.v5[list.gene.use,cells.use]))
matrix.use$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use),
                                 var_type = c(rep(1,ncol(matrix.use)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use)),
                                 is_contextual = c(rep(0,ncol(matrix.use)-1),1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_215-TFs/2Dgas-t36-fs_215-TFs_category-normalized-v5_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)




## With raw counts and using 20 TF
list.gene.use <- c("UTF1","SOX2","POU3F1","DNMT3A","OTX2",
                   "NANOG","EOMES","MIXL1","T","SP5",
                   "MESP1","MESP2","SNAI1","GATA6","HAND1",
                   "HHEX","SOX17","FOXA2","GSC","LHX1")
#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_20-TFs/2Dgas-t36-fs_20-TFs_matrix_all.csv",sep="")
matrix.use <- as.data.frame(t(data.raw.read[list.gene.use,cells.use]))
matrix.use$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use),
                                 var_type = c(rep(1,ncol(matrix.use)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use)),
                                 is_contextual = c(rep(0,ncol(matrix.use)-1),1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_20-TFs/2Dgas-t36-fs_20-TFs_category_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

## With raw counts and using 20 TF (and counts)
#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_20-TFs/2Dgas-t36-fs_20-TFs_with-counts_matrix_all.csv",sep="")
matrix.use <- as.data.frame(t(data.raw.read[list.gene.use,cells.use]))
matrix.use$counts <- as.numeric(colSums(data.coverage[,cells.use]))
matrix.use$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use),
                                 var_type = c(rep(1,ncol(matrix.use)-2),1,0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use)),
                                 is_contextual = c(rep(0,ncol(matrix.use)-2),1,1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_20-TFs/2Dgas-t36-fs_20-TFs_with-counts_category_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

## With normalized v4 and using 20 TF
#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_20-TFs/2Dgas-t36-fs_20-TFs_matrix-normalized-v4_all.csv",sep="")
matrix.use <- as.data.frame(t(data.transcript.v4[list.gene.use,cells.use]))
matrix.use$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use),
                                 var_type = c(rep(1,ncol(matrix.use)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use)),
                                 is_contextual = c(rep(0,ncol(matrix.use)-1),1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_20-TFs/2Dgas-t36-fs_20-TFs_category-normalized-v4_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

## With normalized v5 and using 20 TF
#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_20-TFs/2Dgas-t36-fs_20-TFs_matrix-normalized-v5_all.csv",sep="")
matrix.use <- as.data.frame(t(data.transcript.v5[list.gene.use,cells.use]))
matrix.use$condition <- as.numeric(factor(dataset.fs@meta.data$condition))
write.csv(matrix.use,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use),
                                 var_type = c(rep(1,ncol(matrix.use)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use)),
                                 is_contextual = c(rep(0,ncol(matrix.use)-1),1))
filename <- paste(directory,"/miic/2Dgas-t36-fs_20-TFs/2Dgas-t36-fs_20-TFs_category-normalized-v5_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)




### Test with t48
dataset.gas.48h <- readRDS(file = paste(directory,"/rds/2dgas_t48.rds",sep=""))
list.gene.use <- sort(intersect(union(unique(top.tf$gene),other.tf),row.names(dataset.gas.48h@assays$RNA@data)))
#All
filename <- paste(directory,"/miic/2Dgas-t36-fs_215-TFs/2Dgas-t48_215-TFs_matrix-normalized-v1_all.csv",sep="")
matrix.use <- as.data.frame(t(dataset.gas.48h@assays$RNA@data[list.gene.use,]))
write.csv(matrix.use,filename,row.names = FALSE)

filename <- paste(directory,"/miic/2Dgas-t36-fs_215-TFs/2Dgas-t48_215-TFs_matrix-count_all.csv",sep="")
matrix.use <- as.data.frame(t(dataset.gas.48h@assays$RNA@counts[list.gene.use,]))
write.csv(matrix.use,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use),
                                 var_type = c(rep(1,ncol(matrix.use))),
                                 levels_increasing_order = rep(NA,ncol(matrix.use)),
                                 is_contextual = c(rep(0,ncol(matrix.use))))
filename <- paste(directory,"/miic/2Dgas-t36-fs_215-TFs/2Dgas-t48_215-TFs_category_all.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)


### Test with embryo (Pijuan-Sala2019)
dataset.emb <- readRDS(file = paste(directory,"/rds/pijuan-sala_all.rds",sep=""))
list.emb <- c("Epiblast","Primitive Streak","Nascent mesoderm",
              "Anterior Primitive Streak","Def. endoderm")
dataset.emb <- subset(dataset.emb, (celltype %in% list.emb))

list.gene.use <- c("UTF1","SOX2","POU3F1","DNMT3A","OTX2",
                   "NANOG","EOMES","MIXL1","T","SP5",
                   "MESP1","MESP2","SNAI1","GATA6","HAND1",
                   "HHEX","SOX17","FOXA2","GSC","LHX1")

#All
filename <- paste(directory,"/miic/emb-E7.0_Epi-PS-NasMes-aPS-DE_20-TFs_matrix-normalized-v1.csv",sep="")
matrix.use <- as.data.frame(t(dataset.emb@assays$RNA@data[list.gene.use,]))
write.csv(matrix.use,filename,row.names = FALSE)

filename <- paste(directory,"/miic/emb-E7.0_Epi-PS-NasMes-aPS-DE_20-TFs_matrix-count.csv",sep="")
matrix.use <- as.data.frame(t(dataset.emb@assays$RNA@counts[list.gene.use,]))
write.csv(matrix.use,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use),
                                 var_type = c(rep(1,ncol(matrix.use))),
                                 levels_increasing_order = rep(NA,ncol(matrix.use)),
                                 is_contextual = c(rep(0,ncol(matrix.use))))
filename <- paste(directory,"/miic/emb-E7.0_Epi-PS-NasMes-aPS-DE_20-TFs_category.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

#All with nCount
filename <- paste(directory,"/miic/emb-E7.0_Epi-PS-NasMes-aPS-DE_20-TFs_matrix-count_with-nCount.csv",sep="")
matrix.use <- as.data.frame(t(dataset.emb@assays$RNA@counts[list.gene.use,]))
matrix.use$counts <- as.numeric(dataset.emb$nCount_RNA)
write.csv(matrix.use,filename,row.names = FALSE)

category.order.use <- data.frame(var_names=colnames(matrix.use),
                                 var_type = c(rep(1,ncol(matrix.use)-1),0),
                                 levels_increasing_order = rep(NA,ncol(matrix.use)),
                                 is_contextual = c(rep(0,ncol(matrix.use)-1),1))
filename <- paste(directory,"/miic/emb-E7.0_Epi-PS-NasMes-aPS-DE_20-TFs_category_with-nCount.txt",sep="")
write.table(category.order.use,filename,row.names = FALSE)

