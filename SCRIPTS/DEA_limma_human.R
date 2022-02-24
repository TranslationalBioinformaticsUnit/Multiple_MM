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
library(biomaRt)


setwd("C:/Users/transbio/Desktop/MM_new_dataset")

# Differential expression -------------------------------------------------

### LIMMA
### Uses alinear model. Needs to work with log transformed data, or if not you
### should use voom transformation. Very useful for multifactorial design.

# Norm counts table (Log2 transformed with VOOM)
table_norm <- read.table("DATA/norm_data_rm_IG.txt",  sep="\t", dec=".", row.names=1, check.names = FALSE)


# Experimental design matrix
design_matrix <- read.csv('C:/Users/transbio/Desktop/MM_analysis/huma_data/Samples_metadata_short.csv', sep="", header=T, stringsAsFactors = TRUE)
colnames(design_matrix)[1] <- 'Sample'
summary(design_matrix)
rownames(design_matrix)<- design_matrix$Coded_Name

design_matrix<- design_matrix[colnames(table_norm),]


#### 1- Exploration of normdata (MDS)  ##Hazlo con tu código y con los colorines tan bonitos que generaste!! :)
#Ejemplo de prueba con PCA
pca_table_norm <- prcomp(t(table_norm))
plot(pca_table_norm$x[,1],pca_table_norm$x[,2], col=design_matrix$Combined)



#### 2- Differenital Expression Analysis
DE.design <- model.matrix(~0 + Diagnosis, data= design_matrix)
colnames(DE.design) <- gsub("Diagnosis","",colnames(DE.design))

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
myContrMatrix = makeContrasts(MM_C = (MM - Control),
                              levels=DE.design)
                            


head(myContrMatrix)
res <- contrasts.fit(fit, myContrMatrix)
res <- eBayes(res)

results <- matrix()
list_sig_genes <- vector()

for(i in 1:ncol(myContrMatrix)){
  res_cont <- topTable(res,coef=1,number=nrow(table_norm),sort.by="none")
  sig_genes_up <- which(res_cont$adj.P.Val<0.05 & res_cont$logFC>0.58)
  sig_genes_down <- which(res_cont$adj.P.Val<0.05 & res_cont$logFC< -0.58)
  cat(colnames(myContrMatrix)[1],"\n")
  cat(length(c(sig_genes_up,sig_genes_down)),"\n")
  cat("------------\n")
  res_cont$sig_genes <- rep(0,nrow(res_cont))
  res_cont$sig_genes[sig_genes_up] <- 1
  res_cont$sig_genes[sig_genes_down] <- -1
  res_cont$genes_tendency<- rep(0,nrow(res_cont))
  res_cont$genes_tendency[res_cont$logFC>0] <- 0.5
  res_cont$genes_tendency[res_cont$logFC< 0] <- -0.5
  res_cont$genes_tendency[sig_genes_up] <- 1
  res_cont$genes_tendency[sig_genes_down] <- -1
  colnames(res_cont) <- paste(colnames(myContrMatrix)[1],colnames(res_cont),sep=":")
  results <- cbind(results,res_cont)
  list_sig_genes <- c(list_sig_genes,rownames(res_cont)[c(sig_genes_up,sig_genes_down)])
}

results_MM <- results[,-1]
list_sig_genes_MM <- unique(list_sig_genes) #1465 significant genes ensembl

#save results
write.table(results_MM,"TABLES/dea_MM_CONTROL_0.05_0.58.txt",sep="\t",dec=".",quote=FALSE, row.names = TRUE, col.names = TRUE)
write.table(list_sig_genes_MM,"TABLES/list_sig_genes_MM_CONTROL.txt",sep="\t",dec=".",quote=FALSE)


##read new table

#>>> Expression heatmap
data_to_plot <- table_norm[list_sig_genes_MM,]

#scaled by rows (z-score)
#data_to_plot2 <-  data_to_plot - rowMeans(data_to_plot)

#scale 
data_to_plot2 <- t(scale(t(data_to_plot)))

#heatmap annotation
my_palette1 <- brewer.pal(8, "Dark2")
my_palette2 <- brewer.pal(11, "Paired")

ha1 <- HeatmapAnnotation(
  Stage = design_matrix$Diagnosis,
  #Model = design_matrix$Model,
  col = list(Stage = c("Control"=my_palette1[1],"MGUS"=my_palette1[3],
                       "MM"=my_palette1[4])),
             
  show_annotation_name = TRUE
)

#plot
png("RESULTS/heatmap_sig_genes_MM_Control_rm_IG.png", width=1000, height=1000)
Heatmap(data_to_plot2, name = "DEGs between MM and Control",
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE)
dev.off()




#--Load the mouse data
list_large<-read.table("C:/Users/transbio/Desktop/MM_analysis/mouse_ensembl/Data/list_sig_genes_MM_CONTROL.txt",sep="\t")
head(list_large)
short_positive<-read.table("C:/Users/transbio/Desktop/MM_analysis/mouse_ensembl/Data/list_sig_genes_MM_short_positive.txt")
head(short_positive)
short_negative<-read.table("C:/Users/transbio/Desktop/MM_analysis/mouse_ensembl/Data/list_sig_MM_short_negative.txt")

#--Load names files
names<-read.table("C:/Users/transbio/Desktop/MM_analysis/Huma_data/mouse_sig_in_human.txt",sep="\t")
names<-unique(names)
names$mouse_names<-as.factor(names$mouse_names)
names$human_name<-as.factor(names$human_name)

#--Short positive and negative in human
short_positive_human<-names[names$mouse_names%in%short_positive$V1,]
short_negative_human<-names[names$mouse_names%in%short_negative$V1,]

#--Intersection
intersection<-list_sig_genes_MM[list_sig_genes_MM%in%names$human_name]
#219 genes in the intersection human and mouse
#>>> Expression heatmap of the intersection
data_to_plot <- table_norm[intersection,]


#scale 
data_to_plot2 <- t(scale(t(data_to_plot)))


#heatmap annotation
my_palette1 <- brewer.pal(8, "Dark2")
my_palette2 <- brewer.pal(11, "Paired")

ha1 <- HeatmapAnnotation(
  Stage = design_matrix_RLT$Diagnosis,
  #Model = design_matrix$Model,
  col = list(Stage = c("Control"=my_palette1[1],"MGUS"=my_palette1[3],
                       "MM"=my_palette1[4])),
  
  show_annotation_name = TRUE
)

#plot
png("heatmap_sig_genes_MM_Control_intersection_0,05.png", width=1000, height=1000)
Heatmap(data_to_plot2, name = "DEGs between MM and Control",
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE)
dev.off()
#--Load annotation information
load("C:/Users/transbio/Desktop/MM_analysis/human_rm_IG/Data/IG_genes.Rdata")
#--get gene symbol of intersection
intersection_anno<-myannot[myannot$ensembl_gene_id%in%intersection,]

#--intersection mouse_names
intersection_mouse<-names[names$human_name%in%intersection,]

#--Load mouse data

mouse_MM<-read.table("C:/Users/transbio/Desktop/MM_analysis/mouse_ensembl/Data/dea_MM_CONTROL.txt",sep="\t",row.names=1)

table_norm_mouse <-  read.table("C:/Users/transbio/Desktop/MM_analysis/mouse_ensembl/Data/final_norm_def.txt", sep="\t", dec=".", row.names=1, check.names = FALSE)

design_matrix_mouse<-design_matrix <- read.csv('C:/Users/transbio/Desktop/MM_analysis/mouse_ensembl/Data/Col_Data.csv', sep=";", header=T, stringsAsFactors = TRUE)
design_matrix_mouse<-design_matrix_mouse[colnames(table_norm_mouse),]

data_to_plot<-table_norm_mouse[intersection_mouse$mouse_names,]
#scale 
data_to_plot2 <- t(scale(t(data_to_plot)))
#heatmap annotation
my_palette1 <- brewer.pal(8, "Dark2")
my_palette2 <- brewer.pal(11, "Paired")

ha1 <- HeatmapAnnotation(
  Stage = design_matrix_mouse$Stage,
  Model = design_matrix_mouse$Model,
  col = list(Stage = c("Control"=my_palette1[1],"MGUS"=my_palette1[3],
                       "MM"=my_palette1[4]),
             Model = c("Control"=my_palette1[1],"B2IC"=my_palette2[1],
                       "B2IKC"=my_palette2[2],"CD1B2IC"=my_palette2[3],"cMafB2IC"=my_palette2[4],
                       "Maf_B2IKC"=my_palette2[5], "Maf_MIC"=my_palette2[6],"MIC"=my_palette2[7],
                       "MMsetB2IC"=my_palette2[8], "pB2IC"=my_palette2[9])),
  show_annotation_name = TRUE
)

#plot
png("heatmap_sig_87_mouse.png", width=1000, height=1000)
Heatmap(data_to_plot2, name = "DEGs between MM and Control",
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE)
dev.off()

#--CREATE THE NEW MATRIX WITH THE NAMES


#--PLOT THE DATA WITH MOUSE SIGNATURE
library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

#human gene symbol
library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
anno2 = getBM(
  values = rownames(table_norm),
  filters = c("ensembl_gene_id"),
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  mart = human
)

anno_human <- merge(anno2,results_MM,by.x="ensembl_gene_id",by.y=0,all.y=TRUE)
mouse_sig <- read.table("C:/Users/transbio/Desktop/MM_analysis/mouse_ensembl/Data/results_MM_ano_complete.txt",sep="\t",dec=".",header=TRUE) 
mouse_sig <- mouse_sig[mouse_sig$signature!=0,] #2356 x 71
list_ref_human_sig <- as.character(unique(mouse_sig$Gene.stable.ID.1[!is.na(mouse_sig$Gene.stable.ID.1)]))
common_ensembl <- which(anno_human$ensembl_gene_id%in%list_ref_human_sig)
ref_genes <- unique(anno_human$ensembl_gene_id[common_ensembl])


##plot data
data_to_plot<-table_norm[ref_genes,]

pos_annotation<-vector()
for (i in 1:nrow(data_to_plot)){
  a<-which(anno_human$ensembl_gene_id==rownames(data_to_plot)[i])
  if
  (length(a)>1){
    print(anno_human[a,1:2])
        a <- readline(prompt="Ensembl position: ")
    
  }
  pos_annotation[i]<-a
}
anno_plot <- anno_human[pos_annotation,]
rownames(data_to_plot) <- paste(rownames(data_to_plot),as.character(anno_plot$hgnc_symbol),sep=" : ")

#scaled by rows (z-score)
data_to_plot2 <- t(scale(t(data_to_plot)))

# >>> Y ahora defino la info del ratón

# >>> cómo puede ser que haya ensembl repetidos, para poder generar la columna he hecho un bulce que i encuentra ensembl con más de una entrada que te pregunte que dirección quieres, si lo ejecutas ya verás

direction_info <- vector()
for(i in 1:nrow(anno_plot)){
  dinfo <- mouse_sig$signature[mouse_sig$Gene.stable.ID.1%in%anno_plot$ensembl_gene_id[i]]
  if(length(dinfo)>1){
    cat(dinfo,"\n")
    dinfo <- readline(prompt="Gene direction: ")
  }
  direction_info[i] <- dinfo
}
names(direction_info) <- anno_plot$ensembl_gene_id

# >>> Y la pongo en el atributo rowAnnotation para incluirla en el complex heatmap

row_ha <- rowAnnotation(Mouse = as.factor(direction_info), col=list(Mouse=c("-1"="blue","0"="gray","1"="red")))
# >>> Y ploteo:
pdf("RESULTS/mouse_sig_in_human.pdf")
draw(Heatmap(data_to_plot2, name = "mouse signature in human data",
             top_annotation = ha1, right_annotation = row_ha,
             show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 10),
             cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE))
dev.off()

#--Do the same with the intersection human and mouse
common_ensembl <- which(anno_human$ensembl_gene_id%in%list_ref_human_sig)
ref_genes <- unique(anno_human$ensembl_gene_id[common_ensembl])
ref_genes_intersection<-ref_genes[(ref_genes%in%rownames(sig_genes))]


##plot data
data_to_plot<-table_norm[ref_genes_intersection,]

pos_annotation<-vector()
for (i in 1:nrow(data_to_plot)){
  a<-which(anno_human$ensembl_gene_id==rownames(data_to_plot)[i])
  if
  (length(a)>1){
    print(anno_human[a,1:2])
    a <- readline(prompt="Ensembl position: ")
    
  }
  pos_annotation[i]<-a
}
anno_plot <- anno_human[pos_annotation,]
rownames(data_to_plot) <- paste(rownames(data_to_plot),as.character(anno_plot$hgnc_symbol),sep=" : ")

#scaled by rows (z-score)
data_to_plot2 <- t(scale(t(data_to_plot)))

# >>> Y ahora defino la info del ratón

# >>> cómo puede ser que haya ensembl repetidos, para poder generar la columna he hecho un bulce que i encuentra ensembl con más de una entrada que te pregunte que dirección quieres, si lo ejecutas ya verás

direction_info <- vector()
for(i in 1:nrow(anno_plot)){
  dinfo <- mouse_sig$signature[mouse_sig$Gene.stable.ID.1%in%anno_plot$ensembl_gene_id[i]]
  if(length(dinfo)>1){
    cat(dinfo,"\n")
    dinfo <- readline(prompt="Gene direction: ")
  }
  direction_info[i] <- dinfo
}
names(direction_info) <- anno_plot$ensembl_gene_id

# >>> Y la pongo en el atributo rowAnnotation para incluirla en el complex heatmap

row_ha <- rowAnnotation(Mouse = as.factor(direction_info), col=list(Mouse=c("-1"="blue","0"="gray","1"="red")))
# >>> Y ploteo:
pdf("RESULTS/intersection mouse_sig_in_human_all.pdf")
draw(Heatmap(data_to_plot2, name = "int mouse sig in human data",
             top_annotation = ha1, right_annotation = row_ha,
             show_column_names = FALSE, show_row_names = TRUE, row_names_gp = gpar(fontsize = 4),
             cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE))
dev.off()

#human direction of 
human_direction<-results_MM[ref_genes_intersection,grep("sig",colnames(results_MM))]
colnames(human_direction) <- gsub(":sig_genes","",colnames(human_direction))


human_direc<-vector()
#take only the MM_C case
for (i in 1:nrow(human_direction)){
  human_direc[i]<-human_direction[i,1]
}
names(human_direc)<-rownames(human_direction)

human_direction_matrix<-as.data.frame(human_direc)
mouse_direction<-as.data.frame(direction_info)
directions<-merge(mouse_direction,human_direction_matrix,by="row.names",all = TRUE)
rownames(directions)<-directions$Row.names
directions<-directions[,2:3]

#--take those ones in the same directions
same_direction<-vector()
for (i in 1:nrow(directions)){
  if (directions[i,1]==directions[i,2]){
    same_direction<-c(same_direction,rownames(directions)[i])
 }}
##132 in same direction
##plot data
data_to_plot<-table_norm[same_direction,]

pos_annotation<-vector()
for (i in 1:nrow(data_to_plot)){
  a<-which(anno_human$ensembl_gene_id==rownames(data_to_plot)[i])
  if
  (length(a)>1){
    print(anno_human[a,1:2])
    a <- readline(prompt="Ensembl position: ")
    
  }
  pos_annotation[i]<-a
}
anno_plot <- anno_human[pos_annotation,]
rownames(data_to_plot) <- paste(rownames(data_to_plot),as.character(anno_plot$hgnc_symbol),sep=" : ")

#scaled by rows (z-score)
data_to_plot2 <- t(scale(t(data_to_plot)))

# >>> Y ahora defino la info del ratón

# >>> cómo puede ser que haya ensembl repetidos, para poder generar la columna he hecho un bulce que i encuentra ensembl con más de una entrada que te pregunte que dirección quieres, si lo ejecutas ya verás

direction_info <- vector()
for(i in 1:nrow(anno_plot)){
  dinfo <- mouse_sig$signature[mouse_sig$Gene.stable.ID.1%in%anno_plot$ensembl_gene_id[i]]
  if(length(dinfo)>1){
    cat(dinfo,"\n")
    dinfo <- readline(prompt="Gene direction: ")
  }
  direction_info[i] <- dinfo
}
names(direction_info) <- anno_plot$ensembl_gene_id

# >>> Y la pongo en el atributo rowAnnotation para incluirla en el complex heatmap

row_ha <- rowAnnotation(Mouse = as.factor(direction_info), col=list(Mouse=c("-1"="blue","0"="gray","1"="red")))
# >>> Y ploteo:
pdf("mouse_sig_in_human_same_direction_new.pdf")
draw(Heatmap(data_to_plot2, name = "mouse signature in human data",
             top_annotation = ha1, right_annotation = row_ha,
             show_column_names = FALSE, show_row_names = TRUE, row_names_gp = gpar(fontsize = 4),
             cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE))
dev.off()

##sig_genes from human to mouse

library(biomaRt)
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
human= useMart()

#human gene symbol
anno_biomart = getBM(
  values = list_sig_genes_MM,
  filters = c("ensembl_gene_id"),
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  mart = mouse
)


##plot list_sig_genes_MM with mouse direction
data_to_plot<-table_norm[ref_genes,]

pos_annotation<-vector()
for (i in 1:nrow(data_to_plot)){
  a<-which(anno_human$ensembl_gene_id==rownames(data_to_plot)[i])
  if
  (length(a)>1){
    print(anno_human[a,1:2])
    a <- readline(prompt="Ensembl position: ")
    
  }
  pos_annotation[i]<-a
}
anno_plot <- anno_human[pos_annotation,]
rownames(data_to_plot) <- paste(rownames(data_to_plot),as.character(anno_plot$hgnc_symbol),sep=" : ")

#scaled by rows (z-score)
data_to_plot2 <- t(scale(t(data_to_plot)))

# >>> Y ahora defino la info del ratón

# >>> cómo puede ser que haya ensembl repetidos, para poder generar la columna he hecho un bulce que i encuentra ensembl con más de una entrada que te pregunte que dirección quieres, si lo ejecutas ya verás

direction_info <- vector()
for(i in 1:nrow(anno_plot)){
  dinfo <- mouse_sig$signature[mouse_sig$Gene.stable.ID.1%in%anno_plot$ensembl_gene_id[i]]
  if(length(dinfo)>1){
    cat(dinfo,"\n")
    dinfo <- readline(prompt="Gene direction: ")
  }
  direction_info[i] <- dinfo
}
names(direction_info) <- anno_plot$ensembl_gene_id

# >>> Y la pongo en el atributo rowAnnotation para incluirla en el complex heatmap

row_ha <- rowAnnotation(Mouse = as.factor(direction_info), col=list(Mouse=c("-1"="blue","0"="gray","1"="red")))
# >>> Y ploteo:
pdf("mouse_sig_in_human_1436_def.pdf")
draw(Heatmap(data_to_plot2, name = "mouse signature in human data",
             top_annotation = ha1, right_annotation = row_ha,
             show_column_names = FALSE, show_row_names = TRUE, row_names_gp = gpar(fontsize = 4),
             cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE))
dev.off()

##sig_genes from human to mouse

library(biomaRt)
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
human= useMart()

#human gene symbol
anno_biomart = getBM(
  values = list_sig_genes_MM,
  filters = c("ensembl_gene_id"),
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  mart = mouse
)






