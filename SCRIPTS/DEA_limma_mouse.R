## LIBRARIES TO USE
# Bioconductor
# source("http://bioconductor.org/biocLite.R")
# biocLite("maSigPro")

library(limma)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(corrplot)

set.seed(123456)

setwd("C:/Users/transbio/Desktop/MM_analysis/mouse_ensembl")

# Differential expression -------------------------------------------------

### LIMMA
### Uses alinear model. Needs to work with log transformed data, or if not you
### should use voom transformation. Very useful for multifactorial design.

# Norm counts table (Log2 transformed with VOOM)
table_norm <-  read.table("Data/norm_batch_corrected.txt", sep="\t", dec=".", row.names=1, check.names = FALSE)


# Experimental design matrix
design_matrix <- read.csv('Data/Col_Data.csv', sep=";", header=T, stringsAsFactors = TRUE)
colnames(design_matrix)[1] <- 'Sample'
summary(design_matrix)
rownames(design_matrix)<- design_matrix$Sample

design_matrix<- droplevels(design_matrix[colnames(table_norm),])

design_matrix$Sample <-  as.factor(design_matrix$Sample)
design_matrix$Stage<-as.factor(design_matrix$Stage)
design_matrix$Combined<-as.factor(design_matrix$Combined)
design_matrix$Model<-as.factor(design_matrix$Model)
str(design_matrix)



### 1- Exploratory analysis PCA, with all the samples

##modify the colors
library(RColorBrewer)
display.brewer.all() #view all colors 

# Labels colors
colors <- c(brewer.pal(n=8, name = 'Dark2'), brewer.pal(n=8, name='Set3'),brewer.pal(n=12, name = 'Paired'))

color_types <- colors[1:length(levels(design_matrix$Combined))]
names(color_types) <- levels(design_matrix$Combined)


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
  scale_shape_manual(values=c(25,23,21,24,22))+
  scale_color_manual(values=color_types) +  
  labs(title="MDS - Normalized Count table",x = "Coordinate 1", y = "Coordinate 2")

ggsave(plot = p1, width = 9, height = 7, dpi = 600,filename = "Results/MM/normalized_table_all_populations.pdf")

##MDS no B_cells

samples_excluded<-c("YFP50_YFP_Cg1_GC_B_cell_S23","YFP51_YFP_Cg1_GC_B_cell_S23","YFP72_YFP_Cg1_GC_B_cell_S24")
count_table2 <- table_norm[,!colnames(table_norm)%in%samples_excluded]
design_matrix2<- design_matrix[colnames(count_table2),]


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
  scale_shape_manual(values=c(23,21,24,22))+
  scale_color_manual(values=color_types) +  
  labs(title="MDS - Normalized Count table",x = "Coordinate 1", y = "Coordinate 2")

ggsave(plot = p1, width = 9, height = 7, dpi = 600,filename = "Results/MM/normalized_table_no_Bcells.pdf")



##PCA no B_cells
pca <- prcomp(t(count_table2)) 

#
pr <- summary(pca)$importance[,1:5]
#
mplot <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Combination_type=design_matrix2$Combined, Stage=design_matrix2$Stage)

p1 <- ggplot(mplot, aes(x=PC1, y=PC2, color = Combination_type, label=rownames(mplot),shape=Stage)) +
  geom_point(size=4,fill=design_matrix2$Color) +
  geom_point(size=4, color="black") +
  scale_shape_manual(values=c(23,21,24,22))+
  scale_color_manual(values=color_types) + 
  labs(title="PCA - Normalized Count table",x = paste("PC1 (",round(pr[2,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pr[2,2]*100,0),"%)", sep=""))

ggsave(plot = p1,width = 9, height = 7, dpi = 600, filename = "Results/MM/PCA_no_Bcells.pdf")


###1.2 get the matrix with gene symbol

library(biomaRt)
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
anno2 = getBM(
  values = rownames(table_norm),
  filters = c("ensembl_gene_id"),
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  mart = mouse
)

counts_symbol <- merge(anno2,table_norm,by.x="ensembl_gene_id",by.y=0,all.y=TRUE)




#### 2- Differenital Expression Analysis, by genotype

DE.design <- model.matrix(~0 + Combined, data= design_matrix)
colnames(DE.design) <- gsub("Combined","",colnames(DE.design))

# you have to use voom is your data is not log transfromed. In this case, we have
# log transformed our data during the normalization with the cqn package
# fit = lmFit(voom(rnaseqCQN), DE.design)
fit = lmFit(table_norm, DE.design)
toptable(fit)

# Depending on your questions you have to define the contrasts
# In this case, we are considering the diff expressed genes between control a treatment
# at a specific time point, but removing those genes diff expressed between the same
# time points in the control group.



#**************************#
# MM vs Control all models #
#**************************#
myContrMatrix = makeContrasts(MYC_C = (MYC_MM - Control),
                              Maf_MYC_C = (Maf_MYC_MM - Control),
                              B2IC_C = (B2IC_MM - Control),
                              B2IKC_C = (B2IKC_MM - Control),
                              pB2IC_C = (pB2IC_MM - Control),
                              CD1B2IC_C = (CD1B2IC_MM - Control),
                              MMsetB2IC_C =(MMsetB2IC_MM - Control),
                              cMafB2IC_C =(cMafB2IC_MM - Control),
                              levels= DE.design)
                              

head(myContrMatrix)
res <- contrasts.fit(fit, myContrMatrix)
res <- eBayes(res)

results <- matrix()
list_sig_genes <- vector()

for(i in 1:ncol(myContrMatrix)){
  res_cont <- topTable(res,coef=i,number=nrow(table_norm),sort.by="none")
  sig_genes_up <- which(res_cont$adj.P.Val<0.01 & res_cont$logFC>2)
  sig_genes_down <- which(res_cont$adj.P.Val<0.01 & res_cont$logFC< -2)
  cat(colnames(myContrMatrix)[i],"\n")
  cat(length(c(sig_genes_up,sig_genes_down)),"\n")
  cat("------------\n")
  res_cont$sig_genes <- rep(0,nrow(res_cont))
  res_cont$sig_genes[sig_genes_up] <- 1
  res_cont$sig_genes[sig_genes_down] <- -1
  res_cont$genes_tendency<- rep(0,nrow(res_cont))
  res_cont$genes_tendency[res_cont$logFC>0.58] <- 0.5
  res_cont$genes_tendency[res_cont$logFC< -0.58] <- -0.5
  res_cont$genes_tendency[sig_genes_up] <- 1
  res_cont$genes_tendency[sig_genes_down] <- -1
  colnames(res_cont) <- paste(colnames(myContrMatrix)[i],colnames(res_cont),sep=":")
  results <- cbind(results,res_cont)
  list_sig_genes <- c(list_sig_genes,rownames(res_cont)[c(sig_genes_up,sig_genes_down)])
}

results_MM <- results[,-1]
list_sig_genes_MM <- unique(list_sig_genes) #2234 significant genes

results_MM_symbol <- merge(anno2,results_MM,by.x="ensembl_gene_id",by.y=0,all.y=TRUE)
list_sig_genes<-as.data.frame(list_sig_genes_MM)
colnames(list_sig_genes)<-"ensembl_gene_id"
list_sig_genes_names<-merge(anno2,list_sig_genes,by.x="ensembl_gene_id",by.y="ensembl_gene_id",all.y=TRUE)


#save results
write.table(results_MM,"Results/dea_MM_CONTROL.txt", sep="\t", dec=".", quote=FALSE, row.names = TRUE, col.names = TRUE)
write.table(list_sig_genes_names,"Results/MM/list_sig_genes_MM_CONTROL_names.txt",sep="\t",dec=".",quote=FALSE)


# not take into account in the expression heatmap B-cell and sMM

# take out b_cells
#--REMOVE OUTLIERS AND PCAs
samples_excluded<-c("YFP50_YFP_Cg1_GC_B_cell_S23","YFP51_YFP_Cg1_GC_B_cell_S23","YFP72_YFP_Cg1_GC_B_cell_S24")
count_table2 <- table_norm[,!colnames(table_norm)%in%samples_excluded]
design_matrix2<- design_matrix[colnames(count_table2),]

#>>> Expression heatmap
data_to_plot <- count_table2[list_sig_genes_MM,]

data_to_plot2<- t(scale(t(data_to_plot)))

#heatmap annotation

my_palette1 <- brewer.pal(8, "Dark2")
my_palette2 <- brewer.pal(11, "Paired")

ha1 <- HeatmapAnnotation(
  Stage = design_matrix2$Stage,
  Model = design_matrix2$Model,
  col = list(Stage = c("Control"=my_palette1[2],"MGUS"=my_palette1[3],
                       "MM"=my_palette1[4]),
             Model = c("Control"=my_palette1[2],"B2IC"=my_palette2[1],
                       "B2IKC"=my_palette2[2],"CD1B2IC"=my_palette2[3],"cMafB2IC"=my_palette2[4],
                       "Maf_MYC"=my_palette2[6],"MYC"=my_palette2[7],
                       "MMsetB2IC"=my_palette2[8], "pB2IC"=my_palette2[9])),
  show_annotation_name = TRUE
)

#plot
pdf("Results/MM/expression_heatmap_genotype.pdf") 
jpeg("Results/MM/expression_heatmap_genotype.jpg", width=5000, height=5000, res=600)#

Heatmap(data_to_plot2, name = "DEGs between MM and Control",
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE)



dev.off()




#>>> Differential Expression Analysis "binary" heatmap
#--Taking into account the tendency
data_to_plot_MM <- results_MM[list_sig_genes_MM,grep("tendency",colnames(results_MM))]
colnames(data_to_plot_MM) <- gsub(":genes_tendency","",colnames(data_to_plot_MM))

#--Binary heatmap, sig genes no tendency
data_to_plot <- results_MM[list_sig_genes_MM,grep("sig",colnames(results_MM))]
colnames(data_to_plot) <- gsub(":sig_genes","",colnames(data_to_plot))

#heatmap annotation
my_palette2 <- brewer.pal(11, "Paired")

ha1 <- HeatmapAnnotation(
  Model = colnames(data_to_plot_MM),
  col = list(Model = c("B2IC_C"=my_palette2[1],
                       "B2IKC_C"=my_palette2[2],"CD1B2IC_C"=my_palette2[3],"cMafB2IC_C"=my_palette2[4],
                       "Maf_MYC_C"=my_palette2[6],"MYC_C"=my_palette2[7],
                       "MMsetB2IC_C"=my_palette2[8], "pB2IC_C"=my_palette2[9])),
  show_annotation_name = TRUE
)

#plot
jpeg("heatmap_binary_sig_genes_MM_Control.jpg", width=5000, height=5000, res=600)
pdf("Results/MM_all_genotype/heatmap_binary_sig_genes_MM_Control.pdf")
Heatmap(data_to_plot_MM, name = "DEGs between MM and Control",
        col = colorRamp2(c(-1, 0, 1), c("darkblue","gray","darkorange")),
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE)
dev.off()

library(dendextend)
row_dend = as.dendrogram(hclust(dist(data_to_plot_MM)))

#plot
png("heatmap_binary_sig_genes_MM_Control_row_dend.png", width=1000, height=1000)
Heatmap(data_to_plot_MM, name = "DEGs between MM and Control",
        col = colorRamp2(c(-1, 0, 1), c("darkblue","gray","darkorange")),
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE)
dev.off()

# PCA without b-cells and sMM

##modify the colors
library(RColorBrewer)
display.brewer.all() #view all colors 

# Labels colors
colors <- c(brewer.pal(n=8, name = 'Dark2'), brewer.pal(n=8, name='Set3'),brewer.pal(n=12, name = 'Paired'))

color_types <- colors[1:length(levels(design_matrix$Combined))]
names(color_types) <- levels(design_matrix$Combined)


pca <- prcomp(t(count_table2)) 

#
pr <- summary(pca)$importance[,1:5]
#
mplot <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Combination_type=design_matrix2$Combined, Stage=design_matrix2$Stage)

p1 <- ggplot(mplot, aes(x=PC1, y=PC2, color = Combination_type, label=rownames(mplot),shape=Stage)) +
  geom_point(size=4,fill=design_matrix2$Color) +
  geom_point(size=4, color="black") +
  scale_shape_manual(values=c(23,21,24,22))+
  scale_color_manual(values=color_types) + 
  labs(title="PCA - Normalized MM data",x = paste("PC1 (",round(pr[2,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pr[2,2]*100,0),"%)", sep=""))

ggsave(plot = p1,width = 9, height = 7, dpi = 600, filename = "Results/MM/PCA_MM_B2IC_No_B_cell.pdf")
ggsave(plot = p1,width = 9, height = 7, dpi = 600, filename = "Results/MM/PCA_MM_B2IC_No_B_cell.jpg")


##DEA B2IC models only taking into account MM stage

myContrMatrix = makeContrasts(B2IC_C = (B2IC_MM - Control),
                              B2IKC_C = (B2IKC_MM - Control),
                              pB2IC_C = (pB2IC_MM - Control),
                              CD1B2IC_C = (CD1B2IC_MM - Control),
                              MMsetB2IC_C =(MMsetB2IC_MM - Control),
                              cMafB2IC_C =(cMafB2IC_MM - Control),
                              levels= DE.design)


head(myContrMatrix)
res <- contrasts.fit(fit, myContrMatrix)
res <- eBayes(res)

results <- matrix()
list_sig_genes <- vector()

for(i in 1:ncol(myContrMatrix)){
  res_cont <- topTable(res,coef=i,number=nrow(table_norm),sort.by="none")
  sig_genes_up <- which(res_cont$adj.P.Val<0.01 & res_cont$logFC>2)
  sig_genes_down <- which(res_cont$adj.P.Val<0.01 & res_cont$logFC< -2)
  cat(colnames(myContrMatrix)[i],"\n")
  cat(length(c(sig_genes_up,sig_genes_down)),"\n")
  cat("------------\n")
  res_cont$sig_genes <- rep(0,nrow(res_cont))
  res_cont$sig_genes[sig_genes_up] <- 1
  res_cont$sig_genes[sig_genes_down] <- -1
  res_cont$genes_tendency<- rep(0,nrow(res_cont))
  res_cont$genes_tendency[res_cont$logFC>0.58] <- 0.5
  res_cont$genes_tendency[res_cont$logFC< -0.58] <- -0.5
  res_cont$genes_tendency[sig_genes_up] <- 1
  res_cont$genes_tendency[sig_genes_down] <- -1
  colnames(res_cont) <- paste(colnames(myContrMatrix)[i],colnames(res_cont),sep=":")
  results <- cbind(results,res_cont)
  list_sig_genes <- c(list_sig_genes,rownames(res_cont)[c(sig_genes_up,sig_genes_down)])
}

results_B2IC <- results[,-1]
list_sig_genes_B2IC <- unique(list_sig_genes) #2234 significant genes

#save results
write.table(results_MM,"Results/MM/dea_MM_CONTROL_B2IC_models.txt",sep="\t",dec=".",quote=FALSE, row.names = TRUE, col.names = TRUE)
write.table(list_sig_genes_MM,"Results/MM/list_sig_genes_MM_CONTROL_B2IC_models.txt",sep="\t",dec=".",quote=FALSE)


design_B2IC<-design_matrix[!(design_matrix$Model=="MIC" | design_matrix$Model=="Maf_MIC") ,]

#así tenemos el design_matrix

norm_B2IC<-final_norm[,rownames(design_B2IC)]

#>>> Expression heatmap
data_to_plot <- norm_B2IC[list_sig_genes_B2IC,]

#scaled by rows (z-score)
data_to_plot2 <-  data_to_plot - rowMeans(data_to_plot)

#heatmap annotation
my_palette1 <- brewer.pal(8, "Dark2")
my_palette2 <- brewer.pal(11, "Paired")

ha1 <- HeatmapAnnotation(
  Stage = design_B2IC$Stage,
  Model = design_B2IC$Model,
  col = list(Stage = c("B_Cell"=my_palette1[1],"Control"=my_palette1[2],"MGUS"=my_palette1[3],
                       "MM"=my_palette1[4],"sMM"=my_palette1[5]),
             Model = c("B_Cell"=my_palette1[1],"Control"=my_palette1[2],"B2IC"=my_palette2[1],
                       "B2IKC"=my_palette2[2],"CD1B2IC"=my_palette2[3],"cMafB2IC"=my_palette2[4],
                       "Maf_B2IKC"=my_palette2[5],
                       "MMsetB2IC"=my_palette2[8], "pB2IC"=my_palette2[9])),
  show_annotation_name = TRUE
)

#plot
pdf("Results/MM/heatmap_sig_genes_MM_Control_B2IC.pdf")
jpeg("Results/MM/heatmap_sig_genes_MGUS_Control_B2IC.jpg",width=5000, height=5000, res=600)

Heatmap(data_to_plot2, name = "DEGs between MM and Control",
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE)
dev.off()



#>>> Expression heatmap no B-cells

#--REMOVE OUTLIERS AND PCAs
samples_excluded<-c("YFP50_YFP_Cg1_GC_B_cell_S23","YFP51_YFP_Cg1_GC_B_cell_S23","YFP72_YFP_Cg1_GC_B_cell_S24")
count_table2 <- norm_B2IC[,!colnames(norm_B2IC)%in%samples_excluded]
design_matrix2<- design_B2IC[colnames(count_table2),]

data_to_plot <- count_table2[list_sig_genes_B2IC,]

#scaled by rows (z-score)
data_to_plot2 <-  data_to_plot - rowMeans(data_to_plot)

#heatmap annotation
my_palette1 <- brewer.pal(8, "Dark2")
my_palette2 <- brewer.pal(11, "Paired")

ha1 <- HeatmapAnnotation(
  Stage = design_matrix2$Stage,
  Model = design_matrix2$Model,
  col = list(Stage = c("Control"=my_palette1[2],"MGUS"=my_palette1[3],
                       "MM"=my_palette1[4],"sMM"=my_palette1[5]),
             Model = c("Control"=my_palette1[2],"B2IC"=my_palette2[1],
                       "B2IKC"=my_palette2[2],"CD1B2IC"=my_palette2[3],"cMafB2IC"=my_palette2[4],
                       "Maf_B2IKC"=my_palette2[5], "Maf_MIC"=my_palette2[6],"MIC"=my_palette2[7],
                       "MMsetB2IC"=my_palette2[8], "pB2IC"=my_palette2[9])),
  show_annotation_name = TRUE
)

#plot
pdf("Results/MM/heatmap_sig_genes_MM_Control_B2IC_No_B_cell.pdf")
jpeg("Results/MM/heatmap_sig_genes_MGUS_Control_B2IC_no_B_cell.jpg",width=5000, height=5000, res=600)

Heatmap(data_to_plot2, name = "DEGs between MM and Control",
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE)
dev.off()

##PCA with all the genes and DEG

# PCA with b-cells
pca <- prcomp(t(count_table2)) 

#
pr <- summary(pca)$importance[,1:5]
#
mplot <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Combination_type=design_matrix2$Combined, Stage=design_matrix2$Stage)

p1 <- ggplot(mplot, aes(x=PC1, y=PC2, color = Combination_type, label=rownames(mplot),shape=Stage)) +
  geom_point(size=4,fill=design_matrix2$Color) +
  geom_point(size=4, color="black") +
  scale_shape_manual(values=c(23,21,24,22))+
  scale_color_manual(values=color_types) + 
  labs(title="PCA - DEG MM",x = paste("PC1 (",round(pr[2,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pr[2,2]*100,0),"%)", sep=""))

ggsave(plot = p1,width = 9, height = 7, dpi = 600, filename = "Results/MM/PCA_MM_B2IC_No_B_cell.pdf")
ggsave(plot = p1,width = 9, height = 7, dpi = 600, filename = "Results/MM/PCA_MM_B2IC_No_B_cell.jpg")


##MYC models


myContrMatrix = makeContrasts(MIC_C = (MIC_MM - Control),
                              Maf_MIC_C = (Maf_MIC_MM - Control),
                              levels= DE.design)


head(myContrMatrix)
res <- contrasts.fit(fit, myContrMatrix)
res <- eBayes(res)

results <- matrix()
list_sig_genes <- vector()

for(i in 1:ncol(myContrMatrix)){
  res_cont <- topTable(res,coef=i,number=nrow(table_norm),sort.by="none")
  sig_genes_up <- which(res_cont$adj.P.Val<0.01 & res_cont$logFC>2)
  sig_genes_down <- which(res_cont$adj.P.Val<0.01 & res_cont$logFC< -2)
  cat(colnames(myContrMatrix)[i],"\n")
  cat(length(c(sig_genes_up,sig_genes_down)),"\n")
  cat("------------\n")
  res_cont$sig_genes <- rep(0,nrow(res_cont))
  res_cont$sig_genes[sig_genes_up] <- 1
  res_cont$sig_genes[sig_genes_down] <- -1
  res_cont$genes_tendency<- rep(0,nrow(res_cont))
  res_cont$genes_tendency[res_cont$logFC>0.58] <- 0.5
  res_cont$genes_tendency[res_cont$logFC< -0.58] <- -0.5
  res_cont$genes_tendency[sig_genes_up] <- 1
  res_cont$genes_tendency[sig_genes_down] <- -1
  colnames(res_cont) <- paste(colnames(myContrMatrix)[i],colnames(res_cont),sep=":")
  results <- cbind(results,res_cont)
  list_sig_genes <- c(list_sig_genes,rownames(res_cont)[c(sig_genes_up,sig_genes_down)])
}

results_MYC <- results[,-1]
list_sig_genes_MYC <- unique(list_sig_genes) #2234 significant genes

#save results
write.table(results_MM,"Results/MM/dea_MM_CONTROL_MYC.txt",sep="\t",dec=".",quote=FALSE, row.names = TRUE, col.names = TRUE)
write.table(list_sig_genes_MM,"Results/MM/list_sig_genes_MM_CONTROL_MYC.txt",sep="\t",dec=".",quote=FALSE)

design_myc<-design_matrix[ design_matrix$Model=="MIC" | design_matrix$Model=="Maf_MIC" |design_matrix$Model=="Control" |design_matrix$Model=="B_Cell" ,]

#así tenemos el design_matrix

norm_myc<-final_norm[,rownames(design_myc)]

#>>> Expression heatmap
data_to_plot <- norm_myc[list_sig_genes_MYC,]

#scaled by rows (z-score)
data_to_plot2 <-  data_to_plot - rowMeans(data_to_plot)

#heatmap annotation
my_palette1 <- brewer.pal(8, "Dark2")
my_palette2 <- brewer.pal(11, "Paired")

ha1 <- HeatmapAnnotation(
  Stage = design_myc$Stage,
  Model = design_myc$Model,
  col = list(Stage = c("B_Cell"=my_palette1[1],"Control"=my_palette1[2],"MGUS"=my_palette1[3],
                       "MM"=my_palette1[4]),
             Model = c("B_Cell"=my_palette1[1],"Control"=my_palette1[2],"MIC"=my_palette2[1],
                       "Maf_MIC"=my_palette2[2])),
  show_annotation_name = TRUE
)

#plot
pdf("Results/MM/heatmap_sig_genes_MM_Control_MYC.pdf")
jpeg("Results/MM/heatmap_sig_genes_MM_Control_MYC.jpg",width=5000, height=5000, res=600)

Heatmap(data_to_plot2, name = "DEGs between MM and Control",
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE)
dev.off()



#>>> Expression heatmap no B-cells

#--REMOVE OUTLIERS AND PCAs
samples_excluded<-c("YFP50_YFP_Cg1_GC_B_cell_S23","YFP51_YFP_Cg1_GC_B_cell_S23","YFP72_YFP_Cg1_GC_B_cell_S24")
count_table2 <- norm_myc[,!colnames(norm_myc)%in%samples_excluded]
design_matrix2<- design_myc[colnames(count_table2),]

data_to_plot <- count_table2[list_sig_genes_MYC,]

#scaled by rows (z-score)
data_to_plot2 <-  data_to_plot - rowMeans(data_to_plot)

#heatmap annotation
my_palette1 <- brewer.pal(8, "Dark2")
my_palette2 <- brewer.pal(11, "Paired")

ha1 <- HeatmapAnnotation(
  Stage = design_matrix2$Stage,
  Model = design_matrix2$Model,
  col = list(Stage = c("Control"=my_palette1[2],"MGUS"=my_palette1[3],
                       "MM"=my_palette1[4],"sMM"=my_palette1[5]),
             Model = c("Control"=my_palette1[2],"MIC"=my_palette2[1],
                       "Maf_MIC"=my_palette2[2])),
  show_annotation_name = TRUE
)

#plot
pdf("Results/MM/heatmap_sig_genes_MM_Control_MYC_No_B_cell.pdf")
jpeg("Results/MM/heatmap_sig_genes_MM_Control_MYC_no_B_cell.jpg",width=5000, height=5000, res=600)

Heatmap(data_to_plot2, name = "DEGs between MM and Control",
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE)
dev.off()

##PCA with all the genes and DEG

# PCA with b-cells
pca <- prcomp(t(count_table2[list_sig_genes_MYC,])) 

#
pr <- summary(pca)$importance[,1:5]
#
mplot <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Combination_type=design_matrix2$Combined, Stage=design_matrix2$Stage)

p1 <- ggplot(mplot, aes(x=PC1, y=PC2, color = Combination_type, label=rownames(mplot),shape=Stage)) +
  geom_point(size=4,fill=design_matrix2$Color) +
  geom_point(size=4, color="black") +
  scale_shape_manual(values=c(23,21,24,22))+
  scale_color_manual(values=color_types) + 
  labs(title="PCA - DEG MM",x = paste("PC1 (",round(pr[2,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pr[2,2]*100,0),"%)", sep=""))

ggsave(plot = p1,width = 9, height = 7, dpi = 600, filename = "Results/MM/PCA_MM_MYC_DEG_no_B_cell.pdf")
ggsave(plot = p1,width = 9, height = 7, dpi = 600, filename = "Results/MM/PCA_MM_MYC_DEG_no_B_cell.jpg")


#### 2- Differenital Expression Analysis, by stage

DE.design <- model.matrix(~0 + Stage, data= design_matrix)
colnames(DE.design) <- gsub("Stage","",colnames(DE.design))

# you have to use voom is your data is not log transfromed. In this case, we have
# log transformed our data during the normalization with the cqn package
# fit = lmFit(voom(rnaseqCQN), DE.design)
fit = lmFit(table_norm, DE.design)
toptable(fit)

# Depending on your questions you have to define the contrasts
# In this case, we are considering the diff expressed genes between control a treatment
# at a specific time point, but removing those genes diff expressed between the same
# time points in the control group.



#**************************#
# MM vs Control all models #
#**************************#
myContrMatrix = makeContrasts(MM_C =(MM - Control),
                              levels= DE.design)


head(myContrMatrix)
res <- contrasts.fit(fit, myContrMatrix)
res <- eBayes(res)

results <- matrix()
list_sig_genes <- vector()

for(i in 1:ncol(myContrMatrix)){
  res_cont <- topTable(res,coef=i,number=nrow(table_norm),sort.by="none")
  sig_genes_up <- which(res_cont$adj.P.Val<0.01 & res_cont$logFC>2)
  sig_genes_down <- which(res_cont$adj.P.Val<0.01 & res_cont$logFC< -2)
  cat(colnames(myContrMatrix)[i],"\n")
  cat(length(c(sig_genes_up,sig_genes_down)),"\n")
  cat("------------\n")
  res_cont$sig_genes <- rep(0,nrow(res_cont))
  res_cont$sig_genes[sig_genes_up] <- 1
  res_cont$sig_genes[sig_genes_down] <- -1
  res_cont$genes_tendency<- rep(0,nrow(res_cont))
  res_cont$genes_tendency[res_cont$logFC>0.58] <- 0.5
  res_cont$genes_tendency[res_cont$logFC< -0.58] <- -0.5
  res_cont$genes_tendency[sig_genes_up] <- 1
  res_cont$genes_tendency[sig_genes_down] <- -1
  colnames(res_cont) <- paste(colnames(myContrMatrix)[i],colnames(res_cont),sep=":")
  results <- cbind(results,res_cont)
  list_sig_genes <- c(list_sig_genes,rownames(res_cont)[c(sig_genes_up,sig_genes_down)])
}

results_MM_Stage <- results[,-1]
list_sig_genes_MM_Stage <- unique(list_sig_genes) #2234 significant genes

results_MM_Stage_symbol <- merge(anno2,results_MM_Stage,by.x="ensembl_gene_id",by.y=0,all.y=TRUE)

list_sig_genes<-as.data.frame(list_sig_genes_MM_Stage)
colnames(list_sig_genes)<-"ensembl_gene_id"
list_sig_genes_names<-merge(anno2,list_sig_genes,by.x="ensembl_gene_id",by.y="ensembl_gene_id",all.y=TRUE)


#save results
write.table(results_MM_Stage_symbol,"Results/MM_all_Stage/dea_MM_CONTROL.txt",sep="\t",dec=".",quote=FALSE, row.names = TRUE, col.names = TRUE)
write.table(different_genes,"Results/MM_all_Stage/list_different_genes_dea_genotype.txt",sep="\t",dec=".",quote=FALSE)


# not take into account in the expression heatmap B-cell and sMM

# take out b_cells
#--REMOVE OUTLIERS AND PCAs
samples_excluded<-c("YFP50_YFP_Cg1_GC_B_cell_S23","YFP51_YFP_Cg1_GC_B_cell_S23","YFP72_YFP_Cg1_GC_B_cell_S24","8322_pB2IC_MGUS_S20","8326_pB2IC_MGUS_S10",
                    "8328_pB2IC_MGUS_S29","8330_pB2IC_MGUS_S19","8331_pB2IC_MGUS_S12")
count_table2 <- table_norm[,!colnames(table_norm)%in%samples_excluded]
design_matrix2<- design_matrix[colnames(count_table2),]

#>>> Expression heatmap
data_to_plot <- count_table2[list_sig_genes_MM_Stage,]

data_to_plot2<- t(scale(t(data_to_plot)))

#heatmap annotation

my_palette1 <- brewer.pal(8, "Dark2")
my_palette2 <- brewer.pal(11, "Paired")

ha1 <- HeatmapAnnotation(
  Stage = design_matrix2$Stage,
  Model = design_matrix2$Model,
  col = list(Stage = c("Control"=my_palette1[2],"MGUS"=my_palette1[3],
                       "MM"=my_palette1[4]),
             Model = c("Control"=my_palette1[2],"B2IC"=my_palette2[1],
                       "B2IKC"=my_palette2[2],"CD1B2IC"=my_palette2[3],"cMafB2IC"=my_palette2[4],
                       "Maf_MYC"=my_palette2[6],"MYC"=my_palette2[7],
                       "MMsetB2IC"=my_palette2[8], "pB2IC"=my_palette2[9])),
  show_annotation_name = TRUE
)

#plot
pdf("Results/MM_all_Stage/expression_heatmap_Stage.pdf") 
jpeg("Results/MM/expression_heatmap_genotype.jpg", width=5000, height=5000, res=600)#

Heatmap(data_to_plot2, name = "DEGs between MM and Control",
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE)



dev.off()



#****************************#
# MGUS vs Control all models #
#****************************#
myContrMatrix = makeContrasts(MYC_MGUS_C = (MYC_MGUS - Control),
                              B2IC_MGUS_C = (B2IC_MGUS - Control),
                              B2IKC_MGUS_C = (B2IKC_MGUS - Control),
                              pB2IC_CMGUS_ = (pB2IC_MGUS - Control),
                              CD1B2IC_MGUS_C = (CD1B2IC_MGUS - Control),
                              MMsetB2IC_MGUS_C =(MMsetB2IC_MGUS - Control),
                              cMafB2IC_MGUS_C =(cMafB2IC_MGUS - Control),
                              levels= DE.design)


head(myContrMatrix)
res <- contrasts.fit(fit, myContrMatrix)
res <- eBayes(res)

results <- matrix()
list_sig_genes <- vector()

for(i in 1:ncol(myContrMatrix)){
  res_cont <- topTable(res,coef=i,number=nrow(table_norm),sort.by="none")
  sig_genes_up <- which(res_cont$adj.P.Val<0.01 & res_cont$logFC>2)
  sig_genes_down <- which(res_cont$adj.P.Val<0.01 & res_cont$logFC< -2)
  cat(colnames(myContrMatrix)[i],"\n")
  cat(length(c(sig_genes_up,sig_genes_down)),"\n")
  cat("------------\n")
  res_cont$sig_genes <- rep(0,nrow(res_cont))
  res_cont$sig_genes[sig_genes_up] <- 1
  res_cont$sig_genes[sig_genes_down] <- -1
  res_cont$genes_tendency<- rep(0,nrow(res_cont))
  res_cont$genes_tendency[res_cont$logFC>0.58] <- 0.5
  res_cont$genes_tendency[res_cont$logFC< -0.58] <- -0.5
  res_cont$genes_tendency[sig_genes_up] <- 1
  res_cont$genes_tendency[sig_genes_down] <- -1
  colnames(res_cont) <- paste(colnames(myContrMatrix)[i],colnames(res_cont),sep=":")
  results <- cbind(results,res_cont)
  list_sig_genes <- c(list_sig_genes,rownames(res_cont)[c(sig_genes_up,sig_genes_down)])
}

results_MGUS <- results[,-1]
list_sig_genes_MGUS <- unique(list_sig_genes) #702 significant genes


#save results
write.table(results_MGUS,"Results/dea_MGUS_CONTROL.txt",sep="\t",dec=".",quote=FALSE, row.names = TRUE, col.names = TRUE)
write.table(list_sig_genes_MGUS,"Results/list_sig_genes_MGUS_CONTROL.txt",sep="\t",dec=".",quote=FALSE)



#>>> Expression heatmap
data_to_plot <- table_norm[list_sig_genes_MGUS,]

#scaled by rows (z-score)
data_to_plot2 <-  data_to_plot - rowMeans(data_to_plot)

#heatmap annotation
my_palette1 <- brewer.pal(8, "Dark2")
my_palette2 <- brewer.pal(11, "Paired")

ha1 <- HeatmapAnnotation(
  Stage = design_matrix$Stage,
  Model = design_matrix$Model,
  col = list(Stage = c("B_Cell"=my_palette1[1],"Control"=my_palette1[2],"MGUS"=my_palette1[3],
                       "MM"=my_palette1[4],"sMM"=my_palette1[5]),
             Model = c("B_Cell"=my_palette1[1],"Control"=my_palette1[2],"B2IC"=my_palette2[1],
                       "B2IKC"=my_palette2[2],"CD1B2IC"=my_palette2[3],"cMafB2IC"=my_palette2[4],
                       "Maf_B2IKC"=my_palette2[5], "Maf_MIC"=my_palette2[6],"MIC"=my_palette2[7],
                       "MMsetB2IC"=my_palette2[8], "pB2IC"=my_palette2[9])),
  show_annotation_name = TRUE
)

#>>> Expression heatmap no B-cells

#--REMOVE OUTLIERS AND PCAs
samples_excluded<-c("YFP50_YFP_Cg1_GC_B_cell_S23","YFP51_YFP_Cg1_GC_B_cell_S23","YFP72_YFP_Cg1_GC_B_cell_S24")
count_table2 <- table_norm[,!colnames(table_norm)%in%samples_excluded]
design_matrix2<- design_matrix[colnames(count_table2),]

data_to_plot <- count_table2[list_sig_genes_MGUS,]

#scaled by rows (z-score)
data_to_plot2 <-  data_to_plot - rowMeans(data_to_plot)

#heatmap annotation
my_palette1 <- brewer.pal(8, "Dark2")
my_palette2 <- brewer.pal(11, "Paired")

ha1 <- HeatmapAnnotation(
  Stage = design_matrix2$Stage,
  Model = design_matrix2$Model,
  col = list(Stage = c("Control"=my_palette1[2],"MGUS"=my_palette1[3],
                       "MM"=my_palette1[4],"sMM"=my_palette1[5]),
             Model = c("Control"=my_palette1[2],"B2IC"=my_palette2[1],
                       "B2IKC"=my_palette2[2],"CD1B2IC"=my_palette2[3],"cMafB2IC"=my_palette2[4],
                       "Maf_B2IKC"=my_palette2[5], "Maf_MIC"=my_palette2[6],"MIC"=my_palette2[7],
                       "MMsetB2IC"=my_palette2[8], "pB2IC"=my_palette2[9])),
  show_annotation_name = TRUE
)

#plot
pdf("Results/MGUS/heatmap_sig_genes_MGUS_Control_No_B_cells.pdf")
jpeg("Results/MGUS/heatmap_sig_genes_MGUS_Control_No_B_cells.jpg",width=5000, height=5000, res=600)

Heatmap(data_to_plot2, name = "DEGs between MGUS and Control",
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE)
dev.off()



#>>> Differential Expression Analysis "binary" heatmap

#--Taking into account the tendency
data_to_plot_MGUS<- results_MGUS[list_sig_genes_MGUS,grep("tendency",colnames(results_MGUS))]
colnames(data_to_plot_MGUS) <- gsub(":genes_tendency","",colnames(data_to_plot_MGUS))
#-- Only 1 and 0
data_to_plot <- results_MGUS[list_sig_genes_MGUS,grep("sig",colnames(results_MGUS))]
colnames(data_to_plot) <- gsub(":sig_genes","",colnames(data_to_plot))

#heatmap annotation
my_palette2 <- brewer.pal(11, "Paired")

ha1 <- HeatmapAnnotation(
  Model = colnames(data_to_plot_MGUS),
  col = list(Model = c("MIC_MGUS_C"=my_palette2[1],
                       "B2IC_MGUS_C"=my_palette2[2],"B2IKC_MGUS_C"=my_palette2[3],"pB2IC_CMGUS_"=my_palette2[4],
                       "CD1B2IC_MGUS_C"=my_palette2[5], "MMsetB2IC_MGUS_C"=my_palette2[6],
                       "cMafB2IC_MGUS_C"=my_palette2[7])),
  show_annotation_name = TRUE
)

#plot
pdf("heatmap_binary_sig_genes_MGUS_Control_tendency_ordered.pdf")
Heatmap(data_to_plot, name = "DEGs between MGUS and Control",
        col = colorRamp2(c(-1, 0, 1), c("darkblue","gray","darkorange")),
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = row_dend, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE)

dev.off()

#--Take the sig genes of MM in MGUS, and see tendency
data_to_plot_merge <- results_MGUS[list_sig_genes_MM,grep("tendency",colnames(results_MGUS))]
colnames(data_to_plot_merge) <- gsub(":genes_tendency","",colnames(data_to_plot_merge))

#plot
png("heatmap_binary_sig_genes_2234_in_MGUS.png", width=1000, height=1000)
Heatmap(data_to_plot_merge, name = "DEGs between MM and Control",
        col = colorRamp2(c(-1, 0, 1), c("darkblue","gray","darkorange")),
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = row_dend, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE)
dev.off()

#--Take the intersection


#### 3- Similarities and differencies between MGUS and MM

##Common genes
length(intersect(list_sig_genes_MM,list_sig_genes_MGUS)) #570 -- need to check if the change direction is the same

##Correlation
# with all genes
FC_MM <- results_MM[,grep("logFC",colnames(results_MM))]
colnames(FC_MM) <- gsub(":logFC","",colnames(FC_MM))

FC_MGUS <- results_MGUS[,grep("logFC",colnames(results_MGUS))]
colnames(FC_MGUS) <- gsub(":logFC","",colnames(FC_MGUS))
correlation_MGUS<-cor(FC_MGUS, method="spearman")
FC_all <- merge(FC_MM,FC_MGUS, by.x=0,by.y=0, sort=FALSE)
rownames(FC_all) <- FC_all[,1]
FC_all <- FC_all[,-1]

correlation_FC_spearman<-cor(FC_all, method = "spearman")

pdf("cor_spearman_corrplot_MGUS.pdf")
corrplot(correlation_MGUS, method = "number", type = "upper", order = "original", tl.col = "black", tl.srt = 45,col=brewer.pal(n=8, name="PuOr"))
         
dev.off()
corrplot.mixed(correlation_FC_spearman,lower.col = "black", number.cex = .7,col=brewer.pal(n=8, name="PuOr"))
corrplot.mixed(correlation_FC_spearman, lower = "number", upper = "circle", tl.pos = "lt")                                                            
               

setwd("C:/Users/transbio/Desktop/MM_analysis/mouse_results")
save.image("DE_MM_C.Rdata")

##############################
# MM_Bcell by genotype #######
##############################

myContrMatrix = makeContrasts(MYC_C = (MIC_MM - B_Cell),
                              Maf_MYC_C = (Maf_MIC_MM - B_Cell),
                              B2IC_C = (B2IC_MM - B_Cell),
                              B2IKC_C = (B2IKC_MM - B_Cell),
                              pB2IC_C = (pB2IC_MM - B_Cell),
                              CD1B2IC_C = (CD1B2IC_MM - B_Cell),
                              MMsetB2IC_C =(MMsetB2IC_MM - B_Cell),
                              cMafB2IC_C =(cMafB2IC_MM - B_Cell),
                              levels= DE.design)


head(myContrMatrix)
res <- contrasts.fit(fit, myContrMatrix)
res <- eBayes(res)

results <- matrix()
list_sig_genes <- vector()

for(i in 1:ncol(myContrMatrix)){
  res_cont <- topTable(res,coef=i,number=nrow(table_norm),sort.by="none")
  sig_genes_up <- which(res_cont$adj.P.Val<0.01 & res_cont$logFC>2)
  sig_genes_down <- which(res_cont$adj.P.Val<0.01 & res_cont$logFC< -2)
  cat(colnames(myContrMatrix)[i],"\n")
  cat(length(c(sig_genes_up,sig_genes_down)),"\n")
  cat("------------\n")
  res_cont$sig_genes <- rep(0,nrow(res_cont))
  res_cont$sig_genes[sig_genes_up] <- 1
  res_cont$sig_genes[sig_genes_down] <- -1
  res_cont$genes_tendency<- rep(0,nrow(res_cont))
  res_cont$genes_tendency[res_cont$logFC>0.58] <- 0.5
  res_cont$genes_tendency[res_cont$logFC< -0.58] <- -0.5
  res_cont$genes_tendency[sig_genes_up] <- 1
  res_cont$genes_tendency[sig_genes_down] <- -1
  colnames(res_cont) <- paste(colnames(myContrMatrix)[i],colnames(res_cont),sep=":")
  results <- cbind(results,res_cont)
  list_sig_genes <- c(list_sig_genes,rownames(res_cont)[c(sig_genes_up,sig_genes_down)])
}

results_MM_Bcell <- results[,-1]
list_sig_genes_MM_Bcell <- unique(list_sig_genes) #2234 significant genes

results_MM_Bcell_symbol <- merge(anno2,results_MM_Bcell,by.x="ensembl_gene_id",by.y=0,all.y=TRUE)
list_sig_genes<-as.data.frame(list_sig_genes_MM_Bcell)
colnames(list_sig_genes)<-"ensembl_gene_id"
list_sig_genes_names<-merge(anno2,list_sig_genes,by.x="ensembl_gene_id",by.y="ensembl_gene_id",all.y=TRUE)


#save results
write.table(results_MM_Bcell2,"Results/MM_Bcell_genotype/dea_MM_Bcell.txt", sep="\t", dec=".", quote=FALSE, row.names = TRUE, col.names = TRUE)
write.table(list_sig_genes_names,"Results/MM_Bcell_genotype/list_sig_genes_MM_B_CELL_names.txt",sep="\t",dec=".",quote=FALSE)


##PLOT EXPRESSION HEATMAP

samples_excluded<-c("8322_pB2IC_MGUS_S20","8326_pB2IC_MGUS_S10",
                    "8328_pB2IC_MGUS_S29","8330_pB2IC_MGUS_S19","8331_pB2IC_MGUS_S12")
count_table2 <- table_norm[,!colnames(table_norm)%in%samples_excluded]
design_matrix2<- design_matrix[colnames(count_table2),]

#>>> Expression heatmap
data_to_plot <- count_table2[list_sig_genes_MM_Bcell,]

data_to_plot2<- t(scale(t(data_to_plot)))

#heatmap annotation

my_palette1 <- brewer.pal(8, "Dark2")
my_palette2 <- brewer.pal(11, "Paired")

ha1 <- HeatmapAnnotation(
  Stage = design_matrix2$Stage,
  Model = design_matrix2$Model,
  col = list(Stage = c("B_Cell"=my_palette1[1],"Control"=my_palette1[2],"MGUS"=my_palette1[3],
                       "MM"=my_palette1[4]),
             Model = c("B_Cell"=my_palette1[1], "Control"=my_palette1[2],"B2IC"=my_palette2[1],
                       "B2IKC"=my_palette2[2],"CD1B2IC"=my_palette2[3],"cMafB2IC"=my_palette2[4], "Maf_B2IKC"=my_palette2[5],
                       "Maf_MIC"=my_palette2[6],"MIC"=my_palette2[7],
                       "MMsetB2IC"=my_palette2[8], "pB2IC"=my_palette2[9])),
  show_annotation_name = TRUE
)

#plot
pdf("Results/MM_Bcell_genotype/expression_heatmap_genotype.pdf") 
jpeg("Results/MM_Bcell_genotype/expression_heatmap_genotype.jpg", width=5000, height=5000, res=600)#

Heatmap(data_to_plot2, name = "DEGs between MM and B_cell",
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE)



dev.off()
