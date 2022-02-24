###############################################
#Normalized data visualization MDS and PCA#####
#################################################
library(limma)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(corrplot)

set.seed(123456)

setwd("C:/Users/transbio/Desktop/MM_analysis/Mouse_filtered")

# Differential expression -------------------------------------------------

### LIMMA
### Uses alinear model. Needs to work with log transformed data, or if not you
### should use voom transformation. Very useful for multifactorial design.

# Norm counts table (Log2 transformed with VOOM)
table_norm <-  read.table("DATA/norm_batch_corrected.txt", sep="\t", dec=".", row.names=1, check.names = FALSE)


# Experimental design matrix
design_matrix <- read.csv('DATA/Col_Data.csv', sep=";", header=T, stringsAsFactors = TRUE)
colnames(design_matrix)[1] <- 'Sample'
summary(design_matrix)
rownames(design_matrix)<- design_matrix$Sample

design_matrix<- droplevels(design_matrix[colnames(table_norm),])

design_matrix$Combined <-  as.factor(design_matrix$Combined)
design_matrix$Sex <-  as.factor(design_matrix$Sex)
design_matrix$Stage <-  as.factor(design_matrix$Stage)
design_matrix$Model <-  as.factor(design_matrix$Model)
design_matrix$Batch <-  as.factor(design_matrix$Batch)
design_matrix$Color<-  factor(design_matrix$Color, levels=c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#8DD3C7","#FFFFB3","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#A6CEE3","#1F78B4"))
design_matrix$Color_Batch<-  as.factor(design_matrix$Color_Batch)
str(design_matrix)



### 1- Exploratory analysis PCA, with all the samples (done after batch correction and normalization)

##colors are modified in the design_matrix


tcounts<-t(table_norm)
d<- dist(tcounts)  #euclidean distance between samples
fit<-cmdscale(d, eig=TRUE, k=2)


#plot MDS with raw data, but taking out B_cells
x<- fit$points[,1]
y<- fit$points[,2]

mplot <- data.frame(Coordinate_1=x, Coordinate_2=y, Combination_type=design_matrix$Combined,Stage=design_matrix$Stage)
p1 <- ggplot(mplot, aes(x=Coordinate_1, y=Coordinate_2, color = Combination_type, label=rownames(mplot),shape=Stage)) +
  geom_point(size=4,fill=design_matrix$Color) +
  geom_point(size=4, color="black") +
  scale_shape_manual(values=c(22,23,21,24))+
  scale_color_manual(values=levels(design_matrix$Color)) +  
  labs(title="MDS - Normalized Count table",x = "Coordinate 1", y = "Coordinate 2")

ggsave(plot = p1, width = 9, height = 7, dpi = 600,filename = "RESULTS/normalized_count_table_mds_corrected.pdf")

##MDS no B_cells

samples_excluded<-c("YFP50_YFP_Cg1_GC_B_cell_S23","YFP51_YFP_Cg1_GC_B_cell_S23","YFP72_YFP_Cg1_GC_B_cell_S24")
count_table2 <- table_norm[,!colnames(table_norm)%in%samples_excluded]
design_matrix2<- droplevels(design_matrix[colnames(count_table2),])


tcounts<-t(count_table2)

d<- dist(tcounts)  #euclidean distance between samples
fit<-cmdscale(d, eig=TRUE, k=2)


#plot MDS with raw data
x<- fit$points[,1]
y<- fit$points[,2]

mplot <- data.frame(Coordinate_1=x, Coordinate_2=y, Combination_type=design_matrix2$Combined,Stage=design_matrix2$Stage)
p1 <- ggplot(mplot, aes(x=Coordinate_1, y=Coordinate_2, color = Combination_type, label=rownames(mplot),shape=Stage)) +
  geom_point(size=4,fill=design_matrix2$Color) +
  geom_point(size=4, color="black") +
  scale_shape_manual(values=c(23,21,24))+
  scale_color_manual(values=levels(design_matrix2$Color)) +  
  labs(title="MDS - Normalized Count table",x = "Coordinate 1", y = "Coordinate 2")

ggsave(plot = p1, width = 9, height = 7, dpi = 600,filename = "RESULTS/normalized_count_table_mds_corrected_no_bcell.pdf")


##PCA witb B_cells
pca <- prcomp(t(table_norm)) 

#
pr <- summary(pca)$importance[,1:5]
#
mplot <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Combination_type=design_matrix$Combined, Stage=design_matrix$Stage)

p1 <- ggplot(mplot, aes(x=PC1, y=PC2, color = Combination_type, label=rownames(mplot),shape=Stage)) +
  geom_point(size=4,fill=design_matrix$Color) +
  geom_point(size=4, color="black") +
  scale_shape_manual(values=c(23,21,24,22))+
  scale_color_manual(values=levels(design_matrix$Color)) + 
  labs(title="PCA - Normalized Count table",x = paste("PC1 (",round(pr[2,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pr[2,2]*100,0),"%)", sep=""))

ggsave(plot = p1,width = 9, height = 7, dpi = 600, filename = "RESULTS/PCA_normalized_corrected.pdf")

##PCA no B_cells
pca <- prcomp(t(count_table2)) 

#
pr <- summary(pca)$importance[,1:5]
#
mplot <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Combination_type=design_matrix2$Combined, Stage=design_matrix2$Stage)

p1 <- ggplot(mplot, aes(x=PC1, y=PC2, color = Combination_type, label=rownames(mplot),shape=Stage)) +
  geom_point(size=4,fill=design_matrix2$Color) +
  geom_point(size=4, color="black") +
  scale_shape_manual(values=c(23,21,24))+
  scale_color_manual(values=levels(design_matrix2$Color)) + 
  labs(title="PCA - Normalized Count table",x = paste("PC1 (",round(pr[2,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pr[2,2]*100,0),"%)", sep=""))

ggsave(plot = p1,width = 9, height = 7, dpi = 600, filename = "RESULTS/PCA_normalized_corrected_on_bcell.pdf")


###1.2 MDS and PCA only with b_cell control and MM samples, take out MGUS

design_MM<- droplevels(design_matrix[!(design_matrix$Stage=="MGUS"),])
table_MM<- table_norm[,rownames(design_MM)]

tcounts<-t(table_MM)
d<- dist(tcounts)  #euclidean distance between samples
fit<-cmdscale(d, eig=TRUE, k=2)


#plot MDS with raw data, but taking out B_cells
x<- fit$points[,1]
y<- fit$points[,2]

mplot <- data.frame(Coordinate_1=x, Coordinate_2=y, Combination_type=design_MM$Combined,Stage=design_MM$Stage)
p1 <- ggplot(mplot, aes(x=Coordinate_1, y=Coordinate_2, color = Combination_type, label=rownames(mplot),shape=Stage)) +
  geom_point(size=4,fill=design_MM$Color) +
  geom_point(size=4, color="black") +
  scale_shape_manual(values=c(22,23,24))+
  scale_color_manual(values=levels(design_MM$Color)) +  
  labs(title="MDS - Normalized Count table",x = "Coordinate 1", y = "Coordinate 2")

ggsave(plot = p1, width = 9, height = 7, dpi = 600,filename = "RESULTS/normalized_count_table_mds_MM.pdf")

##MDS no B_cells

samples_excluded<-c("YFP50_YFP_Cg1_GC_B_cell_S23","YFP51_YFP_Cg1_GC_B_cell_S23","YFP72_YFP_Cg1_GC_B_cell_S24")
count_table2 <- table_MM[,!colnames(table_MM)%in%samples_excluded]
design_matrix2<- droplevels(design_MM[colnames(count_table2),])


tcounts<-t(count_table2)

d<- dist(tcounts)  #euclidean distance between samples
fit<-cmdscale(d, eig=TRUE, k=2)


#plot MDS with raw data
x<- fit$points[,1]
y<- fit$points[,2]

mplot <- data.frame(Coordinate_1=x, Coordinate_2=y, Combination_type=design_matrix2$Combined,Stage=design_matrix2$Stage)
p1 <- ggplot(mplot, aes(x=Coordinate_1, y=Coordinate_2, color = Combination_type, label=rownames(mplot),shape=Stage)) +
  geom_point(size=4,fill=design_matrix2$Color) +
  geom_point(size=4, color="black") +
  scale_shape_manual(values=c(22,24))+
  scale_color_manual(values=levels(design_matrix2$Color)) +  
  labs(title="MDS - Normalized Count table",x = "Coordinate 1", y = "Coordinate 2")


ggsave(plot = p1, width = 9, height = 7, dpi = 600,filename = "RESULTS/normalized_count_table_mds_corrected_no_bcell_MM.pdf")


##PCA witb B_cells
pca <- prcomp(t(table_MM)) 

#
pr <- summary(pca)$importance[,1:5]
#
mplot <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Combination_type=design_MM$Combined, Stage=design_MM$Stage)

p1 <- ggplot(mplot, aes(x=PC1, y=PC2, color = Combination_type, label=rownames(mplot),shape=Stage)) +
  geom_point(size=4,fill=design_MM$Color) +
  geom_point(size=4, color="black") +
  scale_shape_manual(values=c(22,23,24))+
  scale_color_manual(values=levels(design_MM$Color)) + 
  labs(title="PCA - Normalized Count table",x = paste("PC1 (",round(pr[2,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pr[2,2]*100,0),"%)", sep=""))

ggsave(plot = p1,width = 9, height = 7, dpi = 600, filename = "RESULTS/PCA_normalized_corrected_MM.pdf")

##PCA no B_cells
pca <- prcomp(t(count_table2)) 

#
pr <- summary(pca)$importance[,1:5]
#
mplot <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Combination_type=design_matrix2$Combined, Stage=design_matrix2$Stage)

p1 <- ggplot(mplot, aes(x=PC1, y=PC2, color = Combination_type, label=rownames(mplot),shape=Stage)) +
  geom_point(size=4,fill=design_matrix2$Color) +
  geom_point(size=4, color="black") +
  scale_shape_manual(values=c(22,24))+
  scale_color_manual(values=levels(design_matrix2$Color)) + 
  labs(title="PCA - Normalized Count table",x = paste("PC1 (",round(pr[2,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pr[2,2]*100,0),"%)", sep=""))

ggsave(plot = p1,width = 9, height = 7, dpi = 600, filename = "RESULTS/PCA_normalized_corrected_on_bcell_MM.pdf")

###1.3 get the matrix with gene symbol

library(biomaRt)
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
anno2 = getBM(
  values = rownames(table_norm),
  filters = c("ensembl_gene_id"),
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  mart = mouse
)

counts_symbol <- merge(anno2,table_norm,by.x="ensembl_gene_id",by.y=0,all.y=TRUE)
counts_symbol_2<-counts_symbol[!((counts_symbol$ensembl_gene_id=="ENSMUSG00000115016") & (counts_symbol$mgi_symbol=="Gm33906")),]

#save normalized table with gene_symbols
write.table(counts_symbol_2,"DATA/normalized_gene_symbol.txt", sep="\t", dec=".", quote=FALSE, row.names = TRUE, col.names = TRUE)


#>>> MDS AND PCA MYC MODEL

design_myc<-droplevels(design_matrix[(design_matrix$Model=="MYC" | design_matrix$Model=="Control"),])
table_myc<-droplevels(table_norm[,rownames(design_myc)])

#MDS

tcounts<-t(table_myc)
d<- dist(tcounts)  #euclidean distance between samples
fit<-cmdscale(d, eig=TRUE, k=2)


#plot MDS with raw data, but taking out B_cells
x<- fit$points[,1]
y<- fit$points[,2]

mplot <- data.frame(Coordinate_1=x, Coordinate_2=y, Combination_type=design_myc$Combined,Stage=design_myc$Stage)
p1 <- ggplot(mplot, aes(x=Coordinate_1, y=Coordinate_2, color = Combination_type, label=rownames(mplot),shape=Stage)) +
  geom_point(size=4,fill=design_myc$Color) +
  geom_point(size=4, color="black") +
  scale_shape_manual(values=c(23,21,24))+
  scale_color_manual(values=levels(design_myc$Color)) +  
  labs(title="MDS - Normalized Count table",x = "Coordinate 1", y = "Coordinate 2")

ggsave(plot = p1, width = 9, height = 7, dpi = 600,filename = "RESULTS/MYC_MDS.pdf")

#PCA

pca <- prcomp(t(table_myc)) 

#
pr <- summary(pca)$importance[,1:5]
#
mplot <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Combination_type=design_myc$Combined, Stage=design_myc$Stage)

p1 <- ggplot(mplot, aes(x=PC1, y=PC2, color = Combination_type, label=rownames(mplot),shape=Stage)) +
  geom_point(size=4,fill=design_myc$Color) +
  geom_point(size=4, color="black") +
  scale_shape_manual(values=c(23,21,24))+
  scale_color_manual(values=levels(design_myc$Color)) + 
  labs(title="PCA - Normalized Count table",x = paste("PC1 (",round(pr[2,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pr[2,2]*100,0),"%)", sep=""))

ggsave(plot = p1,width = 9, height = 7, dpi = 600, filename = "RESULTS/PCA_MYC.pdf")

