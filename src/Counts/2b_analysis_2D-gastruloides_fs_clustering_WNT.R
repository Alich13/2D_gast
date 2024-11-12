library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(mgsa)
library(topGO)
library(plotly)
library(grid)
library(gridExtra)
library(gtable)
library(Cairo)
library(ggpubr)
library(stringr) 
library(cowplot)
library(viridis)
library(ggpointdensity)
library(clustree)

# Set up the directory name
directory <- "D:/Gael/sc_rna_seq"

# # Open the pre-processed dataset
# dataset.fs <- readRDS(file = paste(directory,"/rds/2dgas_t36_fs_PreProcessing.rds",sep=""))
# 
# # Only keep the WNT condition
# dataset.fs.WNT <- subset(dataset.fs, condition == "WNT")
# 
# # Detection of variable genes across the single cells (the main parameter is y.cutoff, the higher the more
# # strict is the selection, the default value is 1)
# dataset.fs.WNT <- FindVariableFeatures(dataset.fs.WNT, selection.method = "vst", nfeatures = 2000)
# sort(dataset.fs.WNT@assays$RNA@var.features)
# 
# # Perform linear dimensional reduction
# dataset.fs.WNT <- RunPCA(dataset.fs.WNT, npcs = 30, verbose = FALSE)
# 
# # ProjectPCA scores each gene in the dataset.fs.WNT
# Cairo(file=paste(directory,"/figures/2dgas-fs/WNT-only/PCA/PCA-Plot-1-2.svg",sep=""),type="svg",width = 4,height = 4,units="cm",dpi=20)
# PCAPlot(dataset.fs.WNT)
# dev.off()
# 
# dataset.fs.WNT <- JackStraw(dataset.fs.WNT,dims = 30)
# dataset.fs.WNT <- ScoreJackStraw(dataset.fs.WNT, dims = 1:30)
# 
# Cairo(file=paste(directory,"/figures/2dgas-fs/WNT-only/PCA/PCA-JackStrawPlot-raw.svg",sep=""),type="svg",width = 6,height = 4,units="cm",dpi=20)
# JackStrawPlot(dataset.fs.WNT,dims = 1:30)
# dev.off()
# 
# jackstraw.pval <- -log10(dataset.fs.WNT@reductions$pca@jackstraw@overall.p.values[,"Score"])
# ndims <- 20
# Cairo(file=paste(directory,"/figures/2dgas-fs/WNT-only/PCA/PCA-JackStrawPlot-pval.svg",sep=""),type="svg",width = 4,height = 4,units="cm",dpi=20)
# print(ggplot(data = data.frame(dims = 1:ndims, stdev = jackstraw.pval[1:ndims])) +
#   geom_point(mapping = aes_string(x = "dims", y = "stdev")) +
#   labs(x = "PC", y = "-log10(p.value)") +
#   geom_hline(yintercept=-log10(0.05),col="red")+
#   theme_cowplot())
# dev.off()
# 
# Cairo(file=paste(directory,"/figures/2dgas-fs/WNT-only/PCA/PCA-ElbowPlot.svg",sep=""),type="svg",width = 4,height = 4,units="cm",dpi=20)
# ElbowPlot(dataset.fs.WNT)
# dev.off()
# 
# nb_pc <- 10
# k <- 10
# 
# ## Run UMAP
# dataset.fs.WNT <- RunUMAP(dataset.fs.WNT, reduction = "pca", dims = 1:nb_pc,n.neighbors = k,
#                        n.components=2)
# 
# ## Reverse umap1 and umap2 
# umap1 <- dataset.fs.WNT@reductions$umap@cell.embeddings[,1]
# umap2 <- dataset.fs.WNT@reductions$umap@cell.embeddings[,2]
# dataset.fs.WNT@reductions$umap@cell.embeddings[,1] = umap1
# dataset.fs.WNT@reductions$umap@cell.embeddings[,2] = -umap2
# 
# size=10
# Cairo(file=paste(directory,"/figures/2dgas-fs/WNT-only/GenePlots/UMAP-Condition.svg",sep=""),type="svg",width = 5,height = 4,units="cm",dpi=20)
# DimPlot(dataset.fs.WNT, reduction = "umap",dims=c(1,2),group.by = "condition")+
#   theme(axis.text = element_text(size=size*2),
#         axis.title = element_text(size=size*2),
#         plot.title = element_blank())
# dev.off()
# 
# # Cluster stability analysis
# dataset.fs.WNT <- FindNeighbors(dataset.fs.WNT, reduction= "pca", dims = 1:nb_pc, k.param=k)
# for (res in c(1:10)/10){
#   dataset.fs.WNT <- FindClusters(object = dataset.fs.WNT, reduction= "pca", resolution = res)
# }
# 
# size=10
# filename = paste(directory,"/figures/2dgas-fs/WNT-only/GenePlots/Cluster-stability-analysis.svg",sep="")
# Cairo(file=filename,type="svg",width = size,height = size*0.9,units="cm",dpi=20)
# clustree(dataset.fs.WNT, prefix = "RNA_snn_res.",node_text_size = size*0.75,edge_width = size*0.25)+
#   theme(legend.key.size = unit(size/12.5, 'cm'),
#         legend.text = element_text(size=size*2.5),
#         legend.title = element_text(size=size*3))
# dev.off()
# 
# # Cluster the cells
# resolution = 0.5
# k <- 10
# dataset.fs.WNT <- FindClusters(object = dataset.fs.WNT, reduction= "pca", resolution = resolution)
# 
# Cairo(file=paste(directory,"/figures/2dgas-fs/WNT-only/GenePlots/UMAP-Cluster-raw.svg",sep=""),type="svg",width = 4,height = 4,units="cm",dpi=20)
# DimPlot(dataset.fs.WNT, reduction = "umap",label=TRUE)+
#   theme(axis.text = element_text(size=size*2),
#         axis.title = element_text(size=size*2))+
#   NoLegend()
# dev.off()
# 
gene <- "LEFTY2"
FeaturePlot(dataset.fs.WNT,gene,order = TRUE,label=TRUE)
VlnPlot(dataset.fs.WNT,gene)
# 
# # Change the name of the clusters
# current.cluster.ids <- c(0:6)
# new.cluster.ids <- c("NasMes","aPS","Epi1","PS","DE","AxMes","Epi2")
# new.cluster.ids.order <- c("Epi1","Epi2","PS","aPS","AxMes","DE","NasMes")
# 
# dataset.fs.WNT@active.ident <- plyr::mapvalues(x = dataset.fs.WNT@active.ident, from = current.cluster.ids, to = new.cluster.ids)
# dataset.fs.WNT@active.ident  <- factor(dataset.fs.WNT@active.ident , levels = new.cluster.ids.order)
# dataset.fs.WNT$celltype <- dataset.fs.WNT@active.ident
# 
# Cairo(file=paste(directory,"/figures/2dgas-fs/WNT-only/GenePlots/UMAP-Cluster-annotated.svg",sep=""),type="svg",width = 4,height = 4,units="cm",dpi=20)
# DimPlot(dataset.fs.WNT, reduction = "umap",dims=c(1,2),label=TRUE,label.size = size)+
#   theme(axis.text = element_text(size=size*2),
#         axis.title = element_text(size=size*2))+
#   NoLegend()
# dev.off()
# 
# # Save the dataset
# saveRDS(dataset.fs.WNT,file = paste(directory,"/rds/2dgas_t36_fs_WNT-only.rds",sep=""))

# Load the dataset from the commented lines above
dataset.fs.WNT <- readRDS(file = paste(directory,"/rds/2dgas_t36_fs_WNT-only.rds",sep=""))

# Test diffusion map
library(SingleCellExperiment)
library(destiny)
library(ggplot2)

sce <- as.SingleCellExperiment(dataset.fs.WNT)
table(sce$ident)
dm <- DiffusionMap(sce, verbose = TRUE,n_pcs=10,k=10)

cellLabels <- sce$ident
tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  DC3 = eigenvectors(dm)[, 3],
                  DC4 = eigenvectors(dm)[, 4],
                  Samples = cellLabels)
ggplot(tmp, aes(x = DC1, y = DC2, colour = Samples)) +
  geom_point()  + 
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2")

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
gene.list.tf <- sort(intersect(toupper(GO.genes.tf),row.names(dataset.fs.WNT[['RNA']]@counts)))

markers.all.tf <- FindAllMarkers(dataset.fs.WNT,only.pos = TRUE,features = gene.list.tf, min.diff.pct = 0.2)
top10.tf<- markers.all.tf %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
DoHeatmap(dataset.fs.WNT,top10.tf$gene)

# Heatmap of signalling proteins 
GO.extra.reg <- "GO:0005576"
list.GO <- c("GO:0060070","GO:0030509","GO:0030509","GO:0005160","GO:0008543",
             "GO:0038084","GO:0048384","GO:0007224","GO:0007219")
entrezid <- unique(na.omit(goall[goall$GOALL==GO.extra.reg,])$ENTREZID)
GO.genes.extra.reg <- symbol[symbol$ENTREZID %in% entrezid,]$SYMBOL

GO.genes <- c()
for (GO in list.GO){
  entrezid <- unique(na.omit(goall[goall$GOALL==GO,])$ENTREZID)
  GO.genes <- union(GO.genes,symbol[symbol$ENTREZID %in% entrezid,]$SYMBOL)
}
GO.genes.signalling <- intersect(GO.genes,GO.genes.extra.reg)
gene.list.signalling <- sort(intersect(toupper(GO.genes.signalling ),row.names(dataset.fs.WNT[['RNA']]@counts)))

markers.all.signalling <- FindAllMarkers(dataset.fs.WNT,only.pos = TRUE, features = gene.list.signalling, min.diff.pct = 0.1)
top10.signalling <- markers.all.signalling %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

DoHeatmap(dataset.fs.WNT,top10.signalling$gene)

# Heatmap of selected genes
gene.heatmap <- c("SOX2","POU3F1","UTF1","CDH1","CLDN6",
                  "NANOG","T","EOMES","WNT3","NODAL",
                  "FOXA2","CER1","SOX17","HHEX","KRT8",
                  "NOTO","FOXJ1","SOX9","CHRD","NOG",
                  "MESP1","MESP2","SNAI1","GATA6","LEFTY2",
                  "FOXC1","FOXC2","TBX6","MEIS2","TWIST1")
DoHeatmap(dataset.fs.WNT,gene.heatmap)

g <- DotPlot(dataset.fs.WNT,features=gene.heatmap)
for (i in 1:6){
  g <- g + geom_vline(xintercept=5*i+0.5)
}
filename = paste(directory,"/figures/2dgas-fs/WNT-only/GenePlots/dotplot_list-genes.svg",sep="")
Cairo(file=filename,type="svg",width = size,height = size,units="cm",dpi=20)
print(g + geom_hline(yintercept=10.6) + 
        labs(x="",y="")+
        coord_flip()+
        theme(legend.key.size = unit(size/6, 'cm'),
              legend.text = element_text(size=size*2),
              legend.title = element_text(size=size*2),
              axis.text = element_text(size=size*2),
              axis.text.x = element_text(angle = 30, hjust = 1)))
dev.off()

# Plot of selected genes
genes.to.plot <- c("SOX2","POU3F1","UTF1","CDH1","OTX2","CLDN3","CLDN6","POU5F1",
                   "POU5F1","WNT3","NANOG","NODAL","FGF8","FST",
                   "T","EOMES","FST","FOXA2","CDH2",
                   "CER1","SOX17","KRT8","KRT18","HHEX","LEFTY1","GSC",
                   "NOTO","FOXJ1","SOX9","CHRD","NOG",
                   "MESP1","MESP2","GATA6","SNAI1","LEFTY2","TDGF1",
                   "FOXF1","HAND1","HAND2","GATA6","CDX2","KDR",
                   "FOXC1","FOXC2","TBX6","MEOX1","MEOX2",
                   "TFAP2A","KLF5","GATA3",
                   "PRDM1","PRDM14","DND1",
                   "ID1","ID2","ID3","AXIN2","LEF1","LEFTY1",
                   "SMAD1","SMAD2","SMAD3","SMAD4","SMAD5","SMAD7","SMAD9","CTNNB1",
                   "TAL1","LYL1","SOX7")
for (gene in genes.to.plot){
  Cairo(file=paste(directory,"/figures/2dgas-fs/WNT-only/GenePlots/UMAP-",gene,".svg",sep=""),type="svg",width = 4,height = 4,units="cm",dpi=20)
  print(FeaturePlot(dataset.fs.WNT,gene,pt.size = 0.5,label=FALSE,order = TRUE,label.size = 5)+
          theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.title.x=element_blank(),
                axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y=element_blank(),
                axis.line=element_blank()))
  dev.off()
  Cairo(file=paste(directory,"/figures/2dgas-fs/WNT-only/GenePlots/VLN-",gene,".svg",sep=""),type="svg",width = 4,height = 4,units="cm",dpi=20)
  print(VlnPlot(dataset.fs.WNT,gene,pt.size=0.1)+NoLegend())
  dev.off()
}

list.conditions <- c("t35-BMP","t36-BMP","t37-BMP","t35-WNT","t36-WNT","t37-WNT")
for (celltype.i in new.cluster.ids.order){
  dataset.cells <- subset(dataset.fs.WNT, subset = celltype ==celltype.i)
  table.celltype <- table(dataset.cells$cells)
  table.celltype[is.na(table.celltype)] <- 0

  table.celltype <- table.celltype[list.conditions]
  Cairo(file=paste(directory,"/figures/2dgas-fs/WNT-only/BarplotCondition/BarPlot-Condition-",celltype.i,".jpg",sep=""),type="svg",width = 5,height = 5,units="cm",dpi=20)
  barplot(100*table.celltype/sum(table.celltype),col = c("Pink","Pink","Pink","Cyan","Cyan","Cyan"),ylim=c(0,50),
          main = paste(celltype.i," (n=",sum(table.celltype),")",sep=""))
  dev.off()
}

# Plot of selected genes by cluster+stimulation
dataset.fs.WNT$celltype.condition <- paste(gsub('[0-9]+', '', dataset.fs.WNT$celltype),dataset.fs.WNT$condition,sep=".")

celltype <- unique(gsub('[0-9]+', '', new.cluster.ids.order))
condition <- c("BMP","WNT")
celltype.condition.order <- rep("",length(celltype)*length(condition))
for (i in 1:length(celltype)){
  celltype.i <- celltype[i]
  for (j in 1:length(condition)){
    condition.j <- condition[j]
    celltype.condition.order[(i-1)*length(condition)+j] <- paste(celltype.i,condition.j,sep=".")
  }
}

dataset.fs.WNT$celltype.condition  <- factor(dataset.fs.WNT$celltype.condition , levels = celltype.condition.order)

for (gene in genes.to.plot){
  Cairo(file=paste(directory,"/figures/2dgas-fs/WNT-only/GenePlots/VLN2-",gene,".svg",sep=""),type="svg",width = 5,height = 5,units="cm",dpi=20)
  print(VlnPlot(dataset.fs.WNT,gene,group.by = "celltype.condition",pt.size=0)+
          stat_summary(fun.y = mean, geom='point', size = 20, colour = "black", shape = 95)+
        NoLegend()+theme(axis.title.x = element_blank()))
  dev.off()
}

## Comparison with embryonic dataset
library(pdist)

time.gas <- "t36"
list.time.emb <- c("E6.5","E6.75","E7.0","E7.25","E7.5")

for (time.emb in list.time.emb){
  
  dataset.emb <- readRDS(file = paste(directory,"/rds/pijuan-sala_",time.emb,"_without-exe.rds",sep=""))
  dataset.emb$dataset <- "embryo"
  
  intersect.genes <- intersect(row.names(dataset.emb),row.names(dataset.fs.WNT))
  
  #####################################################
  #### Workflow with mutual neirest neighors (MNN) ####
  #####################################################
  
  nb.pc <- 20
  k.anchor <- 5
  
  intersect_genes <- intersect(row.names(dataset.emb[['RNA']]),row.names(dataset.fs.WNT[['RNA']]))
  # union_var <- union(VariableFeatures(dataset.fs.WNT),VariableFeatures(dataset.emb))
  union_var <- VariableFeatures(dataset.emb)
  feature <- intersect(intersect_genes,union_var)
  dataset.emb[['RNA']]@counts <- dataset.emb[['RNA']]@counts[intersect_genes,]
  dataset.fs.WNT[['RNA']]@counts <- dataset.fs.WNT[['RNA']]@counts[intersect_genes,]
  dataset.emb[['RNA']]@data<- dataset.emb[['RNA']]@data[intersect_genes,]
  dataset.fs.WNT[['RNA']]@data <- dataset.fs.WNT[['RNA']]@data[intersect_genes,]
  dataset.emb[['RNA']]@scale.data<- dataset.emb[['RNA']]@scale.data[intersect_genes,]
  dataset.fs.WNT[['RNA']]@scale.data <- dataset.fs.WNT[['RNA']]@scale.data[intersect_genes,]
  
  dataset.mnn <- RunFastMNN(object.list = c(dataset.emb,dataset.fs.WNT),features=feature,
                            d=nb.pc, k = k.anchor)
  
  # Run the standard workflow for visualization and clustering
  dataset.mnn <- RunUMAP(dataset.mnn, reduction = "mnn", dims = 1:nb.pc)
  
  dataset.mnn$celltype.emb.only <- dataset.mnn$celltype
  dataset.mnn$celltype.fs.only <- dataset.mnn$celltype
  dataset.mnn$celltype.emb.only[dataset.mnn$celltype %in% unique(dataset.fs.WNT$celltype)] <- time.gas
  dataset.mnn$celltype.fs.only[dataset.mnn$celltype %in% unique(dataset.emb$celltype)] <- time.emb
  
  DimPlot(dataset.mnn, reduction = "umap",group.by = "celltype.fs.only",label = TRUE) + NoLegend()
  DimPlot(dataset.mnn, reduction = "umap",group.by = "celltype.emb.only",label = TRUE) + NoLegend()
  
  ## Correlation with MNN reduction
  list.celltype.emb <- intersect(c("ExE endoderm","Epiblast","Primitive Streak","Anterior Primitive Streak","Def. endoderm","Notochord","Nascent mesoderm",
                         "ExE mesoderm","Haematoendothelial progenitors","Pharyngeal mesoderm","Paraxial mesoderm","PGC"),unique(dataset.mnn$celltype))
  
  list.celltype.gas <- rev(sort(unique(Idents(dataset.fs.WNT))))
  n.celltype.emb <- length(list.celltype.emb)
  n.celltype.gas <- length(list.celltype.gas)
  
  matrix.cor <- matrix(nrow = n.celltype.emb,ncol = n.celltype.gas)
  matrix.dist <- matrix(nrow = n.celltype.emb,ncol = n.celltype.gas)
  for (i in 1:n.celltype.emb){
    celltype.emb <- list.celltype.emb[i]
    dataset.emb.mnn <- subset(dataset.mnn,subset = celltype == celltype.emb)
    df.cor.emb <- cor(t(dataset.emb.mnn@reductions$mnn@cell.embeddings),t(dataset.emb.mnn@reductions$mnn@cell.embeddings))
    cor.emb <- mean(apply(df.cor.emb,2,mean))
    
    for (j in 1:n.celltype.gas){
      celltype.gas <- list.celltype.gas[j]
      dataset.fs.WNT.mnn <- subset(dataset.mnn,subset = celltype == celltype.gas)
      
      df.cor.gas <- cor(t(dataset.emb.mnn@reductions$mnn@cell.embeddings),t(dataset.fs.WNT.mnn@reductions$mnn@cell.embeddings))
      cor.gas <- mean(apply(df.cor.gas,2,mean))
      matrix.cor[i,j] <- cor.gas # /cor.emb
      
      df.dist.gas <- as.matrix(pdist(dataset.emb.mnn@reductions$mnn@cell.embeddings,dataset.fs.WNT.mnn@reductions$mnn@cell.embeddings))
      dist.gas <- mean(apply(df.dist.gas,2,mean))
      matrix.dist[i,j] <- dist.gas
    }
  }
  row.names(matrix.cor) <- list.celltype.emb
  colnames(matrix.cor) <- list.celltype.gas
  row.names(matrix.dist) <- list.celltype.emb
  colnames(matrix.dist) <- list.celltype.gas
  
  # Create the heatmap
  library(reshape2)
  library(ggplot2)
  
  matrix.melted <- melt(matrix.cor)
  min.value <- min(matrix.melted$value)
  max.value <- max(matrix.melted$value)
  
  
  Cairo(file=paste(directory,"/figures/2dgas-fs/WNT-only/Comparison/",time.gas,"vs",time.emb,"_mnn_heatmap_pearson-correlation.svg",sep="")
        ,type="svg",width = 3+n.celltype.emb,height = 2+n.celltype.gas,units="cm",dpi=20)
  print(ggplot(data = matrix.melted, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0.5, limit = c(0,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ 
    labs(y = "2Dgas clusters",x=paste("embryo clusters (",time.emb,")",sep="")) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1,size = 30, hjust = 1),
          axis.text.y = element_text(size = 30),
          legend.title = element_text(size = 30),
          legend.text = element_text(size = 25),
          axis.title = element_text(size = 45))+
    coord_fixed())
  dev.off()
  
  matrix.melted <- melt(matrix.dist)
  max.value <- max(matrix.melted$value,na.rm=T)
  min.value <- min(matrix.melted$value,na.rm=T)
  
  Cairo(file=paste(directory,"/figures/2dgas-fs/WNT-only/Comparison/",time.gas,"vs",time.emb,"_mnn_heatmap_euclidean-distance.svg",sep="")
        ,type="svg",width = 3+n.celltype.emb,height = 2+n.celltype.gas,units="cm",dpi=20)
  print(ggplot(data = matrix.melted, aes(x=Var1, y=Var2, fill=value)) + 
          geom_tile(color = "white")+
          scale_fill_gradient2(low = "red", high = "blue", mid = "white", 
                               midpoint = (0.7+0.1)/2, limit = c(0.1,0.7), space = "Lab", 
                               name="Euclidean\nDistance") +
          theme_minimal()+ 
          labs(y = "2Dgas clusters",x=paste("embryo clusters (",time.emb,")",sep="")) + 
          theme(axis.text.x = element_text(angle = 45, vjust = 1,size = 20, hjust = 1),
                axis.text.y = element_text(size = 20),
                legend.title = element_text(size = 20),
                legend.text = element_text(size = 15),
                axis.title = element_text(size = 20))+
          coord_fixed())
  dev.off()
}

## Comparison with ss3
library(pdist)

time.fs <- "FS"
time.ss3 <- "SS3"

dataset.fs.WNT <- readRDS(file = paste(directory,"/rds/2dgas_t36_ss3.rds",sep=""))

dataset.fs.WNT$dataset <- time.fs
dataset.fs.WNT$dataset <- time.ss3

intersect.genes <- intersect(row.names(dataset.fs.WNT),row.names(dataset.fs.WNT))

#####################################################
#### Workflow with mutual neirest neighors (MNN) ####
#####################################################

nb.pc <- 10
k.anchor <- 5

intersect_genes <- intersect(row.names(dataset.fs.WNT[['RNA']]),row.names(dataset.fs.WNT[['RNA']]))
# union_var <- union(VariableFeatures(dataset.fs.WNT),VariableFeatures(dataset.fs.WNT))
union_var <- union(VariableFeatures(dataset.fs.WNT),VariableFeatures(dataset.fs.WNT))
feature <- intersect(intersect_genes,union_var)
dataset.fs.WNT[['RNA']]@counts <- dataset.fs.WNT[['RNA']]@counts[intersect_genes,]
dataset.fs.WNT[['RNA']]@counts <- dataset.fs.WNT[['RNA']]@counts[intersect_genes,]
dataset.fs.WNT[['RNA']]@data<- dataset.fs.WNT[['RNA']]@data[intersect_genes,]
dataset.fs.WNT[['RNA']]@data <- dataset.fs.WNT[['RNA']]@data[intersect_genes,]
dataset.fs.WNT[['RNA']]@scale.data<- dataset.fs.WNT[['RNA']]@scale.data[intersect_genes,]
dataset.fs.WNT[['RNA']]@scale.data <- dataset.fs.WNT[['RNA']]@scale.data[intersect_genes,]

dataset.mnn <- RunFastMNN(object.list = c(dataset.fs.WNT,dataset.fs.WNT),features=feature,
                          d=nb.pc, k = k.anchor)

# Run the standard workflow for visualization and clustering
dataset.mnn <- RunUMAP(dataset.mnn, reduction = "mnn", dims = 1:nb.pc)

dataset.mnn$celltype.ss3.only <- dataset.mnn$celltype
dataset.mnn$celltype.fs.only <- dataset.mnn$celltype
dataset.mnn$celltype.ss3.only[(dataset.mnn$celltype %in% unique(dataset.fs.WNT$celltype))&
                              (dataset.mnn$dataset == time.fs)] <- time.fs
dataset.mnn$celltype.fs.only[(dataset.mnn$celltype %in% unique(dataset.fs.WNT$celltype))&
                               (dataset.mnn$dataset == time.ss3)] <- time.ss3

DimPlot(dataset.mnn, reduction = "umap",group.by = "celltype.fs.only",label = TRUE) + NoLegend()
DimPlot(dataset.mnn, reduction = "umap",group.by = "celltype.ss3.only",label = TRUE) + NoLegend()

## Correlation with MNN reduction
list.celltype.ss3 <- sort(unique(Idents(dataset.fs.WNT)))
list.celltype.fs <- rev(sort(unique(Idents(dataset.fs.WNT))))

n.celltype.ss3 <- length(list.celltype.ss3)
n.celltype.fs <- length(list.celltype.fs)

matrix.cor <- matrix(nrow = n.celltype.ss3,ncol = n.celltype.fs)
matrix.dist <- matrix(nrow = n.celltype.ss3,ncol = n.celltype.fs)
for (i in 1:n.celltype.ss3){
  celltype.ss3 <- list.celltype.ss3[i]
  dataset.fs.WNT.mnn <- subset(dataset.mnn,subset = celltype == celltype.ss3 & dataset == time.ss3)
  df.cor.ss3 <- cor(t(dataset.fs.WNT.mnn@reductions$mnn@cell.embeddings),t(dataset.fs.WNT.mnn@reductions$mnn@cell.embeddings))
  cor.ss3 <- mean(apply(df.cor.ss3,2,mean))
  
  for (j in 1:n.celltype.fs){
    celltype.fs <- list.celltype.fs[j]
    dataset.fs.WNT.mnn <- subset(dataset.mnn,subset = celltype == celltype.fs & dataset == time.fs)
    
    df.cor.fs <- cor(t(dataset.fs.WNT.mnn@reductions$mnn@cell.embeddings),t(dataset.fs.WNT.mnn@reductions$mnn@cell.embeddings))
    cor.fs <- mean(apply(df.cor.fs,2,mean))
    matrix.cor[i,j] <- cor.fs # /cor.ss3
    
    df.dist.fs <- as.matrix(pdist(dataset.fs.WNT.mnn@reductions$mnn@cell.embeddings,dataset.fs.WNT.mnn@reductions$mnn@cell.embeddings))
    dist.fs <- mean(apply(df.dist.fs,2,mean))
    matrix.dist[i,j] <- dist.fs
  }
}
row.names(matrix.cor) <- list.celltype.ss3
colnames(matrix.cor) <- list.celltype.fs
row.names(matrix.dist) <- list.celltype.ss3
colnames(matrix.dist) <- list.celltype.fs

# Create the heatmap
library(reshape2)
library(ggplot2)

matrix.melted <- melt(matrix.cor)
min.value <- min(matrix.melted$value)
max.value <- max(matrix.melted$value)

Cairo(file=paste(directory,"/figures/2dgas-fs/WNT-only/Comparison/",time.fs,"vs",time.ss3,"_mnn_heatmap_pearson-correlation.svg",sep="")
      ,type="svg",width = 3+n.celltype.ss3,height = 2+n.celltype.fs,units="cm",dpi=20)
print(ggplot(data = matrix.melted, aes(x=Var1, y=Var2, fill=value)) + 
        geom_tile(color = "white")+
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                             midpoint = 0.5, limit = c(-0.1,1), space = "Lab", 
                             name="Pearson\nCorrelation") +
        theme_minimal()+ 
        labs(y = "FS",x = "SS3") + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1,size = 30, hjust = 1),
              axis.text.y = element_text(size = 30),
              legend.title = element_text(size = 30),
              legend.text = element_text(size = 25),
              axis.title = element_text(size = 45))+
        coord_fixed())
dev.off()

matrix.melted <- melt(matrix.dist)
max.value <- max(matrix.melted$value,na.rm=T)
min.value <- min(matrix.melted$value,na.rm=T)

Cairo(file=paste(directory,"/figures/2dgas-fs/WNT-only/Comparison/",time.fs,"vs",time.ss3,"_mnn_heatmap_euclidean-distance.svg",sep="")
      ,type="svg",width = 3+n.celltype.ss3,height = 2+n.celltype.fs,units="cm",dpi=20)
print(ggplot(data = matrix.melted, aes(x=Var1, y=Var2, fill=value)) + 
        geom_tile(color = "white")+
        scale_fill_gradient2(low = "red", high = "blue", mid = "white", 
                             midpoint = (0.7+0.1)/2, limit = c(0.1,0.7), space = "Lab", 
                             name="Euclidean\nDistance") +
        theme_minimal()+ 
        labs(y = "FS",x = "SS3") + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1,size = 20, hjust = 1),
              axis.text.y = element_text(size = 20),
              legend.title = element_text(size = 20),
              legend.text = element_text(size = 15),
              axis.title = element_text(size = 20))+
        coord_fixed())
dev.off()


# Plot expression of PS markers in the PS cluster to check whether there
# are differences between the BMP4 and WNT3A condition
celltype.i <- "PS1"
dataset.PS <- subset(dataset.fs.WNT, subset = celltype==celltype.i)
genes.to.plot <- c("POU5F1","NANOG","T","EOMES","SP5","FST","NODAL","WNT3","FOXA2","TDGF1",
                   "FGF8","MESP1","MESP2","CDX2","PRDM1","CHRD","CER1",
                   "ID1","ID2","ID3","AXIN2","LEF1","LEFTY1","LEFTY2",
                   "SMAD1","SMAD2","SMAD3","SMAD4","SMAD5","SMAD7")

markers <- FindMarkers(dataset.PS,"BMP","WNT", group.by = "condition")
top10.tf<- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)


for (gene in genes.to.plot){
  Cairo(file=paste(directory,"/figures/2dgas-fs/WNT-only/ExpressionByCondition/VLNPlot-ByCondition-",celltype.i,"-",gene,".svg",sep=""),type="svg",width = 2,height = 2.5,units="cm",dpi=20)
  print(VlnPlot(dataset.PS,gene,group.by = "condition")+ 
          scale_y_continuous(limits = c(0,NA))+
          stat_summary(fun.y = mean, geom='point', size = 25, colour = "green", shape = 95)+
          theme(legend.position="none"))
  dev.off()
}

## Test using t35-WNT only
dataset.t35.WNT <- subset(dataset.fs.WNT, subset = condition == "WNT" & time == "t35")
dataset.t35.WNT <- RunPCA(dataset.t35.WNT, npcs = 30, verbose = FALSE)

dataset.t35.WNT <- JackStraw(dataset.t35.WNT,dims = 30)
dataset.t35.WNT <- ScoreJackStraw(dataset.t35.WNT, dims = 1:30)

JackStrawPlot(dataset.t35.WNT,dims = 1:30)

## Run UMAP
nb_pc <- 5
k <- 15
dataset.t35.WNT <- RunUMAP(dataset.t35.WNT, reduction = "pca", dims = 1:nb_pc,n.neighbors = k,
                       n.components=2)

DimPlot(dataset.t35.WNT, reduction = "umap",label=TRUE)

gene <- "NOTO"
FeaturePlot(dataset.t35.WNT,gene,order = TRUE,label=TRUE)

################################################
#### Test signalling pathway quantification ####
################################################

## Test MAYA (https://github.com/One-Biosciences/MAYA/)

library(MAYA)
count_mat <- Matrix::Matrix(as.matrix(dataset.fs.WNT@assays$RNA@counts))
meta <- dataset.fs.WNT@meta.data[,c("nCount_RNA","nFeature_RNA","celltype")]

activity_summary<-MAYA_pathway_analysis(expr_mat=count_mat,
                                        modules_list = "hallmark",
                                        is_logcpm=F)

plot_heatmap_activity_mat(activity_mat = scale_0_1(activity_summary$activity_matrix), 
                          meta = meta, 
                          annot_name = "celltype")

plot_umap_annot(umap=activity_summary$umap,
                labels = meta$celltype,
                title = "Author annotation - HALLMARK")

plot_umap_pathway_activity(umap=activity_summary$umap,
                           PCA_object = activity_summary$PCA_obj,
                           module = "KEGG_MAPK_SIGNALING_PATHWAY")
plot_umap_pathway_activity(umap=activity_summary$umap,
                           PCA_object = activity_summary$PCA_obj,
                           module = "KEGG_TGF_BETA_SIGNALING_PATHWAY")

logcpm<-logcpmNormalization(count_mat)
plot_heatmap_pathway_top_contrib_genes(expr_mat=logcpm,
                                       PCA_object = activity_summary$PCA_obj,
                                       module = "KEGG_MAPK_SIGNALING_PATHWAY",
                                       n = 10,
                                       meta = meta,
                                       annot_name = "celltype")
