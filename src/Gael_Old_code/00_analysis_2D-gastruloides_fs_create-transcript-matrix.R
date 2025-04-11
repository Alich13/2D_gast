library(Seurat)
library(dplyr)
library(grid)
library(gridExtra)
library(millefy)

## Set up the directory name
directory <- "D:/Gael/sc_rna_seq"

#############################
#### Import the raw data ####
#############################

data.fs.read.1 <- as.matrix(Read10X(data.dir = paste(directory,"/raw_data/data_2dgas_fs/Plate1/10XlikeMatrix_read",sep="")))
data.fs.read.2.3 <- as.matrix(Read10X(data.dir = paste(directory,"/raw_data/data_2dgas_fs/Plate2-3/10XlikeMatrix_read",sep="")))

list.genes <- intersect(row.names(data.fs.read.1),row.names(data.fs.read.2.3))
data.fs.read <- cbind(data.fs.read.1[list.genes,],data.fs.read.2.3[list.genes,])
colnames(data.fs.read) <- substring(colnames(data.fs.read), first = 1, last = regexpr("_",unique(colnames(data.fs.read)))[1]-1)

data.fs.read.save <- data.fs.read
row.names(data.fs.read.save) <- toupper(row.names(data.fs.read.save))
write.table(data.fs.read.save,paste(directory,'/raw_data/data_2dgas_fs/RawReadMatrix.txt',sep=""))

# Import a result from a rsem run to get gene length
rsem.table <- read.table(paste(directory,'/raw_data/data_2dgas_fs/D1439T380_Sorted.genes.txt',sep=""),
                         header = T)

require(org.Mm.eg.db)
genes_id <- rsem.table$gene_id
genes_name <- mapIds(org.Mm.eg.db,
                     keys = genes_id,
                     column = 'SYMBOL',
                     keytype = 'ENSEMBL')
names(genes_id) <- genes_name

genes_length <- rep(NA,length(genes_id))
names(genes_length) <- names(genes_id)

# Create the vector with gene length
for (gene.i in names(genes_length)){
  genes_length[gene.i] <- rsem.table[rsem.table$gene_id==genes_id[gene.i],"length"]
}

# Add genes not found previously
path_gtf = c(paste(directory,"/raw_data/data_2dgas_fs/mm10_genes.gtf",sep=""))
dt_gtf_exon <- gtfToDtExon(path_gtf)
dt_gtf_exon$gene_name <- substring(dt_gtf_exon$transcript_name,first=1, last = regexpr("-2",dt_gtf_exon$transcript_name)-1)

other.genes <- unique(dt_gtf_exon$gene_name)[!(unique(dt_gtf_exon$gene_name) %in% names(genes_length))]
for (gene.i in other.genes){
  gene.i.exon <- dt_gtf_exon[dt_gtf_exon$transcript_name==paste(gene.i,"-201",sep=""),]
  genes_length[gene.i] <- sum(gene.i.exon$end) - sum(gene.i.exon$start)
}

to_save <- genes_length[sort(names(genes_length))]
write.csv(to_save,paste(directory,"/MouseGeneLength.csv",sep=""),
          row.names = names(to_save))
names(to_save) <- toupper(names(to_save))
to_save2 <- to_save[MI_condition_top_TF]
write.csv(to_save2,paste(directory,"/MouseGeneLength-GenesOfInterest.csv",sep=""),
          row.names = names(to_save2))

# Create the matrix that estimates transcript counts (version 1)
data.fs.transcript <- matrix(NA, nrow = nrow(data.fs.read), ncol = ncol(data.fs.read))
colnames(data.fs.transcript) <- colnames(data.fs.read)
row.names(data.fs.transcript) <- row.names(data.fs.read)

for (gene.i in row.names(data.fs.transcript)){
  data.fs.transcript[gene.i,] <- 100*data.fs.read[gene.i,]/genes_length[gene.i]
}

data.fs.transcript <- na.omit(round(data.fs.transcript+0.33))

write.table(data.fs.transcript,paste(directory,'/raw_data/data_2dgas_fs/TranscriptEstimationMatrix.txt',sep=""))


# Create the matrix that estimates transcript counts (version 2)
data.fs.transcript <- matrix(NA, nrow = nrow(data.fs.read), ncol = ncol(data.fs.read))
colnames(data.fs.transcript) <- colnames(data.fs.read)
row.names(data.fs.transcript) <- row.names(data.fs.read)
for (gene.i in row.names(data.fs.transcript)){
  data.fs.transcript[gene.i,] <- ceiling(100*data.fs.read[gene.i,]/genes_length[gene.i])
}
data.fs.transcript <- na.omit(data.fs.transcript)

data.fs.transcript.norm <- matrix(NA, nrow = nrow(data.fs.transcript), ncol = ncol(data.fs.transcript))
colnames(data.fs.transcript.norm) <- colnames(data.fs.transcript)
row.names(data.fs.transcript.norm) <- row.names(data.fs.transcript)
for (cell.i in colnames(data.fs.transcript.norm)){
  data.fs.transcript.norm[,cell.i] <- data.fs.transcript[,cell.i]/sum(data.fs.transcript[,cell.i])
}
for (gene.i in row.names(data.fs.transcript.norm)){
  min.gene.i <- min(data.fs.transcript.norm[gene.i,data.fs.transcript.norm[gene.i,]>0])
  data.fs.transcript.norm[gene.i,] <- round(data.fs.transcript.norm[gene.i,]/min.gene.i)
}

row.names(data.fs.transcript.norm) <- toupper(row.names(data.fs.transcript.norm))

write.table(data.fs.transcript.norm,paste(directory,'/raw_data/data_2dgas_fs/TranscriptEstimationMatrix-v2.txt',sep=""))

# Create the matrix that estimates transcript counts (version 3)
data.fs.read.norm <- matrix(NA, nrow = nrow(data.fs.read), ncol = ncol(data.fs.read))
colnames(data.fs.read.norm) <- colnames(data.fs.read)
row.names(data.fs.read.norm) <- row.names(data.fs.read)
for (cell.i in colnames(data.fs.read.norm)){
  data.fs.read.norm[,cell.i] <- data.fs.read[,cell.i]/sum(data.fs.read[,cell.i])
}
for (gene.i in row.names(data.fs.read.norm)){
  min.gene.i <- min(data.fs.read.norm[gene.i,data.fs.read.norm[gene.i,]>0])
  data.fs.read.norm[gene.i,] <- round(data.fs.read.norm[gene.i,]/min.gene.i)
}

row.names(data.fs.read.norm) <- toupper(row.names(data.fs.read.norm))

write.table(data.fs.read.norm,paste(directory,'/raw_data/data_2dgas_fs/CountEstimationMatrix-v3.txt',sep=""))

# Create the matrix that estimates transcript counts (version 4)
data.fs.read.norm <- matrix(NA, nrow = nrow(data.fs.read), ncol = ncol(data.fs.read))
colnames(data.fs.read.norm) <- colnames(data.fs.read)
row.names(data.fs.read.norm) <- row.names(data.fs.read)
for (cell.i in colnames(data.fs.read.norm)){
  data.fs.read.norm[,cell.i] <- ceiling(1000000*data.fs.read[,cell.i]/sum(data.fs.read[,cell.i]))
}

row.names(data.fs.read.norm) <- toupper(row.names(data.fs.read.norm))

write.table(data.fs.read.norm,paste(directory,'/raw_data/data_2dgas_fs/CountEstimationMatrix-v4.txt',sep=""))

# Create the matrix that estimates transcript counts (version 5)
data.fs.read.norm <- matrix(NA, nrow = nrow(data.fs.read), ncol = ncol(data.fs.read))
colnames(data.fs.read.norm) <- colnames(data.fs.read)
row.names(data.fs.read.norm) <- row.names(data.fs.read)
for (cell.i in colnames(data.fs.read.norm)){
  data.fs.read.norm[,cell.i] <- ceiling(4000000*data.fs.read[,cell.i]/sum(data.fs.read[,cell.i]))
}

data.fs.transcript.norm <- matrix(NA, nrow = nrow(data.fs.read.norm), ncol = ncol(data.fs.read.norm))
colnames(data.fs.transcript.norm) <- colnames(data.fs.read.norm)
row.names(data.fs.transcript.norm) <- row.names(data.fs.read.norm)
for (gene.i in row.names(data.fs.transcript.norm)){
  data.fs.transcript.norm[gene.i,] <- ceiling(100*data.fs.read.norm[gene.i,]/genes_length[gene.i])
}
data.fs.transcript.norm <- na.omit(data.fs.transcript.norm)
mean(colSums(data.fs.transcript.norm))

row.names(data.fs.transcript.norm) <- toupper(row.names(data.fs.transcript.norm))

write.table(data.fs.transcript.norm,paste(directory,'/raw_data/data_2dgas_fs/TranscriptEstimationMatrix-v5.txt',sep=""))

# Create the matrix that estimates transcript counts (version xx)
data.fs.read.norm <- matrix(NA, nrow = nrow(data.fs.read), ncol = ncol(data.fs.read))
colnames(data.fs.read.norm) <- colnames(data.fs.read)
row.names(data.fs.read.norm) <- row.names(data.fs.read)
for (cell.i in colnames(data.fs.read.norm)){
  data.fs.read.norm[,cell.i] <- ceiling(10000000*data.fs.read[,cell.i]/sum(data.fs.read[,cell.i]))
}

gene.i <- "Nodal"
for (gene.i in list.genes){
  if (gene.i %in% row.names(data.fs.read.norm)){
    to_plot <- data.fs.read.norm[gene.i,cells.use]
    to_plot <- to_plot[to_plot>0]
    
    length.gene.i = genes_length[gene.i]/100
    n.i = ceiling(max(to_plot)/length.gene.i)
    to_plot <- to_plot[to_plot<quantile(to_plot,0.99)]
    
    test <- hist(to_plot,breaks = 200,main = gene.i)
    print(gene.i)
    print(genes_length[gene.i])
  }
}

for (j in 1:10){
  abline(v = j*length.gene.i,col="red")
}

100*sum(data.fs.read.norm[gene.i,])/genes_length[gene.i] 

data.fs.read[,"D1439T318"]

# Create the matrix that estimates transcript counts using coverage (version xx)
library(rtracklayer)
path_gtf = c(paste(directory,"/raw_data/data_2dgas_fs/mm10_genes.gtf",sep=""))
gtf_mm10 <- import(path_gtf)
gtf_mm10 <- gtf_mm10[gtf_mm10$type == "transcript",]
list.genes <- sort(unique(gtf_mm10$gene_name))

df.gene.ranges <-data.frame(chr=character(),min=integer(),max=integer())
# row.names(df.gene.ranges) <- list.genes
for (gene.i in list.genes){
  gtf_mm10.i <- gtf_mm10[gtf_mm10$gene_name == gene.i,]
  chr.i <- gtf_mm10.i@seqnames@values[1]
  min.i <- min(gtf_mm10.i@ranges@start)
  max.i <- max(gtf_mm10.i@ranges@start+gtf_mm10.i@ranges@width-1)
  df.gene.ranges[gene.i,"chr"] <- as.character(chr.i)
  df.gene.ranges[gene.i,c("min","max")] <- c(min.i,max.i)
}

list_files <- list.files(paste(directory,"/raw_data/data_2dgas_fs/coverage",sep=""))
data.fs.transcript.cov <- matrix(0, nrow = length(list.genes), ncol = length(list_files))
# data.fs.transcript.cov <- read.table(paste(directory,'/raw_data/data_2dgas_fs/TranscriptEstimationMatrix-v6.txt',sep=""))

for (j in 736:length(list_files)){
  file.j <- list_files[j]
  cell.j <- substring(file.j, first = 1, last = regexpr("_",file.j)-1)
  matrix <- read.table(paste(directory,"/raw_data/data_2dgas_fs/coverage/",file.j,sep=""),
                       col.names = c("chr","min","max","cov"))
  expression <- rep(0,length(list.genes))
  for (i in 1:length(list.genes)){
    gene.i <- list.genes[i]
    df.gene.ranges.gene.i <- df.gene.ranges[row.names(df.gene.ranges)==gene.i,]
    matrix.i <- matrix[(matrix$chr == df.gene.ranges.gene.i$chr)&
                       (matrix$min >= df.gene.ranges.gene.i$min)& 
                       (matrix$max <= df.gene.ranges.gene.i$max),]
    if (nrow(matrix.i) > 0){
      expression[i] <- max(matrix.i$cov)
    }
  }
  data.fs.transcript.cov[,j] <- expression
}

colnames(data.fs.transcript.cov) <- substring(list_files, first = 1, last = regexpr("_",list_files)-1)
row.names(data.fs.transcript.cov) <- toupper(list.genes)

write.table(data.fs.transcript.cov,paste(directory,'/raw_data/data_2dgas_fs/TranscriptEstimationMatrix-Coverage.txt',sep=""))

# Same matrix but normalized to have the same number of transcripts per cell
data.fs.transcript.cov <- read.table(paste(directory,'/raw_data/data_2dgas_fs/TranscriptEstimationMatrix-Coverage.txt',sep=""))
data.fs.transcript.cov.norm <- matrix(NA, nrow = nrow(data.fs.transcript.cov), ncol = ncol(data.fs.transcript.cov))
colnames(data.fs.transcript.cov.norm) <- colnames(data.fs.transcript.cov)
row.names(data.fs.transcript.cov.norm) <- row.names(data.fs.transcript.cov)

for (cell.i in colnames(data.fs.transcript.cov.norm)){
  data.fs.transcript.cov.norm[,cell.i] <- ceiling(500000*data.fs.transcript.cov[,cell.i]/sum(data.fs.transcript.cov[,cell.i]))
}

write.table(data.fs.transcript.cov.norm,paste(directory,'/raw_data/data_2dgas_fs/TranscriptEstimationMatrix-v6.txt',sep=""))
