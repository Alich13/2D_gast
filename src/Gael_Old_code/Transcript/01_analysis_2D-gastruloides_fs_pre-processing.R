library(Seurat)
library(SeuratWrappers)
library(ggpubr)
library(dplyr)
library(Cairo)
library(stringr) 
library(clustree)

## Set up the directory name
directory <- "D:/Gael/sc_rna_seq"

#############################
#### Import the raw data ####
#############################

data.fs.transcript <- as.matrix(read.table(paste(directory,"/raw_data/data_2dgas_fs/TranscriptEstimationMatrix.txt",sep="")))

data.fs.umi.1 <- as.matrix(Read10X(data.dir = paste(directory,"/raw_data/data_2dgas_fs/Plate1/10XlikeMatrix_umi",sep="")))
data.fs.umi.2.3 <- as.matrix(Read10X(data.dir = paste(directory,"/raw_data/data_2dgas_fs/Plate2-3/10XlikeMatrix_umi",sep="")))

list.genes <- intersect(row.names(data.fs.transcript),intersect(row.names(data.fs.umi.1),row.names(data.fs.umi.2.3)))
data.fs.umi <- cbind(data.fs.umi.1[list.genes,],data.fs.umi.2.3[list.genes,])
data.fs.transcript.test <- data.fs.transcript[list.genes,]

size=10

data.x1 <- log1p(colSums(data.fs.transcript.test))
data.y1 <- log1p(colSums(data.fs.umi))

dat1 <- data.frame(x = data.x1, y = data.y1)

p <- ggplot(data = dat1, mapping = aes(x = x, y = y)) +
  geom_point()+
  theme_bw()+
  labs(x="Transcript mean expression (log)",y="UMI mean expression (log)")+
  theme(axis.text = element_text(size=size*5),
        axis.title = element_text(size=size*5),
        legend.key.size = unit(0.25*size, 'cm'),
        legend.title= element_text(size=size*5),
        legend.text = element_text(size=size*5),
        legend.title.align=0.5,
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
p

Cairo(file=paste(directory,"/figures/2dgas-fs/PreProcessing/MeanTranscriptVsMeanUmi.svg",sep=""),type="svg",width = size,height = size,units="cm",dpi=20)
print(p)
dev.off()


list.genes <- c("Actb","Gapdh","Rpl13","T")
for (gene.i in list.genes){
  data.x2 <- log1p(data.fs.transcript.test[gene.i,])
  data.y2 <- log1p(data.fs.umi[gene.i,])
  
  dat2 <- data.frame(x = data.x2, y = data.y2)
  
  p <- ggplot(data = dat2, mapping = aes(x = x, y = y)) +
    geom_point()+
    theme_bw()+
    labs(x="Transcript",y="UMI",title=bquote(italic(.(str_to_title(gene.i)))))+
    theme(axis.text = element_text(size=size*5),
          axis.title = element_text(size=size*5),
          plot.title = element_text(size=size*5,hjust=0.5),
          legend.key.size = unit(0.25*size, 'cm'),
          legend.title= element_text(size=size*5),
          legend.text = element_text(size=size*5),
          legend.title.align=0.5,
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))
  p
  
  Cairo(file=paste(directory,"/figures/2dgas-fs/PreProcessing/TranscriptVsUMI-",gene.i,".svg",sep=""),type="svg",width = size,height = size,units="cm",dpi=20)
  print(p)
  dev.off()
}

sample.information <- read.table(paste(directory,"/raw_data/data_2dgas_fs/sample_information.txt",sep=""),
                                 header = TRUE)[colnames(data.fs.transcript),]

cells <- paste(sample.information$time,sample.information$condition,sep="-")

row.names(data.fs.transcript) <- toupper(row.names(data.fs.transcript))
colnames(data.fs.transcript) <- paste(cells,"_",colnames(data.fs.transcript ),sep="")

dataset.fs <- CreateSeuratObject(data.fs.transcript , project = "fs", min.cells = 1)
dataset.fs@meta.data$cells <- dataset.fs@active.ident
dataset.fs@meta.data$time <- as.factor(sample.information$time)
dataset.fs@meta.data$condition <- as.factor(sample.information$condition)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

## QC and selecting cells for further analysis
# Get the percentage of mitochondrial genes in each cell
mito.genes <- grep(pattern = "^MT-", x = rownames(x = GetAssayData(dataset.fs, slot = "counts")), value = TRUE)
percent.mito <- Matrix::colSums(GetAssayData(dataset.fs, slot = "counts")[mito.genes, ])/Matrix::colSums(GetAssayData(dataset.fs, slot = "counts"))
dataset.fs <- AddMetaData(object = dataset.fs, metadata = percent.mito, col.name = "percent.mito")

dataset.fs <- subset(dataset.fs, subset = nCount_RNA > 0 & percent.mito < 0.05)

# QC plots
sums <- data.frame(nb_transcript=dataset.fs$nCount_RNA)
Cairo(file=paste(directory,"/figures/2dgas-fs/PreProcessing/QC/Histogram-transcript.svg",sep=""),type="svg",width = 4,height = 4,units="cm",dpi=20)
ggplot(sums, aes(x=nb_transcript)) + geom_histogram(binwidth=10000)+
  geom_vline(xintercept = mean(sums$nb_transcript), color = "red", size=1.5)
dev.off()

sums <- data.frame(nb_transcript=dataset.fs$nFeature_RNA)
Cairo(file=paste(directory,"/figures/2dgas-fs/PreProcessing/QC/Histogram-nGenes.svg",sep=""),type="svg",width = 4,height = 4,units="cm",dpi=20)
ggplot(sums, aes(x=nb_transcript)) + geom_histogram(binwidth=1000)+
  geom_vline(xintercept = mean(sums$nb_transcript), color = "red", size=1.5)
dev.off()

for (to.plot in c("nFeature_RNA", "nCount_RNA", "percent.mito")){
  Cairo(file=paste(directory,"/figures/2dgas-fs/PreProcessing/QC/VLN-",to.plot,".svg",sep=""),type="svg",width = 6,height = 6,units="cm",dpi=20)
  print(VlnPlot(object = dataset.fs, features=to.plot,slot = "counts"))
  dev.off()
}

Cairo(file=paste(directory,"/figures/2dgas-fs/PreProcessing/QC/Scatter-CountVsFeature.svg",sep=""),type="svg",width = 6,height = 6,units="cm",dpi=20)
FeatureScatter(object = dataset.fs, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",slot = "counts")
dev.off()

Cairo(file=paste(directory,"/figures/2dgas-fs/PreProcessing/QC/Scatter-CountVsGapdh.svg",sep=""),type="svg",width = 6,height = 6,units="cm",dpi=20)
FeatureScatter(object = dataset.fs, feature1 = "nCount_RNA", feature2 = "GAPDH",
               slot = "counts")
dev.off()

## Scale the data to have 80000 transcript (the mean is 85000) in each cell and then, normalizing the data with log1p
mean(dataset.fs$nCount_RNA)
dataset.fs <- NormalizeData(object = dataset.fs, verbose = FALSE,scale.factor = 80000)

## Look at T, Eomes and Pou5f1 expression when scale with nCount or with Gapdh

for (gene in c("T","EOMES","POU5F1")){
  data_x <- dataset.fs@assays$RNA@data[gene,]
  data_y <- log1p(max(dataset.fs@assays$RNA@counts[gene,])*dataset.fs@assays$RNA@counts[gene,]/dataset.fs@assays$RNA@counts["GAPDH",])
  c <- dataset.fs$condition
  df <- data.frame(x=data_x,y=data_y,condition=c)
  
  Cairo(file=paste(directory,"/figures/2dgas-fs/PreProcessing/QC/Scatter-Scalling-",gene,".svg",sep=""),type="svg",width = 4,height = 4,units="cm",dpi=20)
  print(ggplot(data=df, aes(x=x, y=y)) +
    geom_smooth(method="lm") +
    geom_point() +
    stat_regline_equation(label.x=0.2*max(data_x), label.y=0.84*max(data_y)) +
    stat_cor(aes(label=..rr.label..), label.x=0.2*max(data_x), label.y=0.76*max(data_y))+
    xlab(substitute(paste("log(",italic(p), ") (Scale with nCounts)", sep=""), list(p=str_to_title(tolower(gene)))))+
    ylab(substitute(paste("log(",italic(p), ") (Scale with ",italic(q),")", sep=""), list(p=str_to_title(tolower(gene)),q="Gapdh"))))
  dev.off()
}

## Detection of variable genes across the single cells (the main parameter is y.cutoff, the higher the more
## strict is the selection, the default value is 1)
dataset.fs <- FindVariableFeatures(dataset.fs, selection.method = "vst", nfeatures = 2000)
sort(dataset.fs@assays$RNA@var.features)

# ## Perform cell cycle regression
# dataset.fs <- CellCycleScoring(dataset.fs, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# dataset.fs <- ScaleData(object = dataset.fs,vars.to.regress = c("ntranscript","percent.mito","S.Score", "G2M.Score"),
#                          features = rownames(dataset.fs))

dataset.fs <- ScaleData(object = dataset.fs,features = row.names(dataset.fs[['RNA']]@counts))

#################################
#### Pipeline for clustering ####
#################################

# Perform linear dimensional reduction
dataset.fs <- RunPCA(dataset.fs, npcs = 30, verbose = FALSE)

# List of top genes in PCx
PCx <- 2
top20.neg <- names(sort(dataset.fs@reductions$pca@feature.loadings[,PCx]))[1:20]

# ProjectPCA scores each gene in the dataset.fs
Cairo(file=paste(directory,"/figures/2dgas-fs/PreProcessing/PCA/PCA-Plot-1-2.svg",sep=""),type="svg",width = 4,height = 4,units="cm",dpi=20)
PCAPlot(dataset.fs)
dev.off()

dataset.fs <- JackStraw(dataset.fs,dims = 30)
dataset.fs <- ScoreJackStraw(dataset.fs, dims = 1:30)

Cairo(file=paste(directory,"/figures/2dgas-fs/PreProcessing/PCA/PCA-JackStrawPlot-raw.svg",sep=""),type="svg",width = 6,height = 4,units="cm",dpi=20)
JackStrawPlot(dataset.fs,dims = 1:30)
dev.off()

jackstraw.pval <- -log10(dataset.fs@reductions$pca@jackstraw@overall.p.values[,"Score"])
ndims <- 20
Cairo(file=paste(directory,"/figures/2dgas-fs/PreProcessing/PCA/PCA-JackStrawPlot-pval.svg",sep=""),type="svg",width = 4,height = 4,units="cm",dpi=20)
print(ggplot(data = data.frame(dims = 1:ndims, stdev = jackstraw.pval[1:ndims])) + 
  geom_point(mapping = aes_string(x = "dims", y = "stdev")) + 
  labs(x = "PC", y = "-log10(p.value)") + 
  geom_hline(yintercept=-log10(0.05),col="red")+
  theme_cowplot())
dev.off()


Cairo(file=paste(directory,"/figures/2dgas-fs/PreProcessing/PCA/PCA-ElbowPlot.svg",sep=""),type="svg",width = 4,height = 4,units="cm",dpi=20)
ElbowPlot(dataset.fs)
dev.off()

nb_pc <- 10
k <- 10

## Run UMAP
dataset.fs <- RunUMAP(dataset.fs, reduction = "pca", dims = 1:nb_pc,n.neighbors = k,
                       n.components=2)
Cairo(file=paste(directory,"/figures/2dgas-fs/PreProcessing/GenePlots/UMAP-Condition.svg",sep=""),type="svg",width = 5,height = 4,units="cm",dpi=20)
DimPlot(dataset.fs, reduction = "umap",dims=c(1,2),group.by = "condition")+
  theme(axis.text = element_text(size=size*2),
        axis.title = element_text(size=size*2),
        plot.title = element_blank())
dev.off()

## Cluster the cells
resolution = 0.1
k <- 10
dataset.fs <- FindNeighbors(dataset.fs, reduction= "pca", dims = 1:nb_pc, k.param=k)
dataset.fs <- FindClusters(object = dataset.fs, reduction= "pca", resolution = resolution)

Cairo(file=paste(directory,"/figures/2dgas-fs/PreProcessing/GenePlots/UMAP-Cluster-raw.svg",sep=""),type="svg",width = 4,height = 4,units="cm",dpi=20)
DimPlot(dataset.fs, reduction = "umap",label=TRUE)+
  theme(axis.text = element_text(size=size*2),
        axis.title = element_text(size=size*2))+
  NoLegend()
dev.off()

for (res in c(1:10)/10){
  dataset.fs <- FindClusters(object = dataset.fs, reduction= "pca", resolution = res)
}

size=10
filename = paste(directory,"/figures/2dgas-fs/PreProcessing/GenePlots/Cluster-stability-analysis.svg",sep="")
Cairo(file=filename,type="svg",width = size,height = size*0.9,units="cm",dpi=20)
clustree(dataset.fs, prefix = "RNA_snn_res.",node_text_size = size*0.75,edge_width = size*0.25)+
  theme(legend.key.size = unit(size/12.5, 'cm'),
        legend.text = element_text(size=size*2.5),
        legend.title = element_text(size=size*3))
dev.off()

Cairo(file=paste(directory,"/figures/2dgas-fs/PreProcessing/GenePlots/UMAP-nCount_RNA.svg",sep=""),type="svg",width = 6,height = 6,units="cm",dpi=20)
FeaturePlot(dataset.fs,"nCount_RNA",order = TRUE)+
  theme(axis.text = element_text(size=size*2),
        axis.title = element_text(size=size*2))
dev.off()

Cairo(file=paste(directory,"/figures/2dgas-fs/PreProcessing/GenePlots/VLN-nCount_RNA.svg",sep=""),type="svg",width = 6,height = 6,units="cm",dpi=20)
VlnPlot(dataset.fs,"nCount_RNA",group.by = "RNA_snn_res.0.1")+
  theme(axis.text = element_text(size=size*2),
        axis.title = element_text(size=size*2))+
  NoLegend()
dev.off()

# Filter out the cells
dataset.fs <- subset(dataset.fs, subset = RNA_snn_res.0.1 != 3)

# Save the dataset
saveRDS(dataset.fs,file = paste(directory,"/rds/2dgas_t36_fs_PreProcessing.rds",sep=""))
