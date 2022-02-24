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

setwd("C:/Users/transbio/Desktop/MM_analysis/Mouse_filtered")

# Differential expression -------------------------------------------------

### LIMMA
### Uses alinear model. Needs to work with log transformed data, or if not you
### should use voom transformation. Very useful for multifactorial design.

# Norm counts table (Log2 transformed with VOOM)
table_norm <-  read.table("DATA/norm_batch_corrected.txt", sep="\t", dec=".", row.names=1, check.names = FALSE)


# Experimental design matrix
design_matrix <- read.csv('Col_Data.csv', sep=";", header=T, stringsAsFactors = TRUE)
colnames(design_matrix)[1] <- 'Sample'
summary(design_matrix)
rownames(design_matrix)<- design_matrix$Sample

design_matrix<- droplevels(design_matrix[colnames(table_norm),])

design_matrix$Sample <-  as.factor(design_matrix$Sample)
design_matrix$Stage<-as.factor(design_matrix$Stage)
design_matrix$Combined<-as.factor(design_matrix$Combined)
design_matrix$Model<-as.factor(design_matrix$Model)
str(design_matrix)

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



##############################
######## MM vs Control #######
##############################
myContrMatrix = makeContrasts(MYC_C = (MYC_MM - Control),
                              Maf_MYC_C = (Maf_MYC_MM - Control),
                              B2IC_C = (B2IC_MM - Control),
                              B2IKC_C = (B2IKC_MM - Control),
                              pB2IC_C = (pB2IC_MM - Control),
                              CD1B2IC_C = (CD1B2IC_MM - Control),
                              MMsetB2IC_C =(MMsetB2IC_MM - Control),
                              cMafB2IC_C =(cMafB2IC_MM - Control),
                              levels= DE.design)
                              
myContrMatrix = makeContrasts(MM_C = (pB2IC_MM - Control),
                              sMM_C = (pB2IC_sMM - Control),
                              MGUS_C = (pB2IC_MGUS - Control),
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



#save results
write.table(results_MM,"DEA/MM_C/100920_dea_MM_CONTROL.txt", sep="\t", dec=".", quote=FALSE, row.names = TRUE, col.names = TRUE)
write.table(list_sig_genes_MM,"DEA/MM_C/100920_list_sig_genes_MM_CONTROL.txt",sep="\t",dec=".",quote=FALSE)

#include gene symbol in our data

library(biomaRt)
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
anno2 = getBM(
  values = rownames(table_norm),
  filters = c("ensembl_gene_id"),
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  mart = mouse
)


results_MM_symbol <- merge(anno2,results_MM,by.x="ensembl_gene_id",by.y=0,all.y=TRUE)
results_MM_symbol2<-results_MM_symbol[!((results_MM_symbol$ensembl_gene_id=="ENSMUSG00000115016") & (results_MM_symbol$mgi_symbol=="Gm33906")),]

write.table(results_MM_symbol2,"DEA/MM_C/100920_dea_MM_CONTROL_names.txt", sep="\t", dec=".", quote=FALSE, row.names = TRUE, col.names = TRUE)

list_sig_genes<-as.data.frame(list_sig_genes_MM)
colnames(list_sig_genes)<-"ensembl_gene_id"
list_sig_genes_names<-merge(anno2,list_sig_genes,by.x="ensembl_gene_id",by.y="ensembl_gene_id",all.y=TRUE)

write.table(list_sig_genes_names,"DEA/MM_C/100920_list_sig_genes_MM_CONTROL_names.txt",sep="\t",dec=".",quote=FALSE)


#INTERSECTION, ONLY THOSE GENES THAT ARE SIG IN ALL MODELS

rownames(results_MM_symbol2)<-results_MM_symbol2$ensembl_gene_id
results_MM_sig <- results_MM_symbol2[,grep("sig_genes",colnames(results_MM_symbol2))]
core_genes <-results_MM_sig[apply(results_MM_sig,1,FUN=function(x){return(sum(x!=0)==8)}),] #181 genes sig en todos los modelos

core_genes_names<-list_sig_genes_names[list_sig_genes_names$ensembl_gene_id%in%rownames(core_genes),]

write.table(core_genes_names,"DEA/MM_C/100920_list_core_genes.txt",sep="\t",dec=".",quote=FALSE)


##check differences with the previous analysis
old_list<-read.table("C:/Users/transbio/Desktop/MM_analysis/mouse_ensembl/Data/list_sig_genes_MM_CONTROL.txt")

sum(list_sig_genes_MM%in%old_list$x) # 2096 genes are common in both analysis

new_genes<-list_sig_genes_names[!(list_sig_genes_names$ensembl_gene_id %in% old_list$x),]

write.table(new_genes,"DEA/MM_C/new_genes.txt",sep="\t",dec=".",quote=FALSE)

lost_genes<-list_sig_genes_names[!(old_list$x%in% list_sig_genes_names$ensembl_gene_id),]

write.table(lost_genes,"DEA/MM_C/lost_genes.txt",sep="\t",dec=".",quote=FALSE)


#----EXPRESSION HEATMAP

# take out b_cells

samples_excluded<-c("YFP50_YFP_Cg1_GC_B_cell_S23","YFP51_YFP_Cg1_GC_B_cell_S23","YFP72_YFP_Cg1_GC_B_cell_S24")
count_table2 <- table_norm[,!colnames(table_norm)%in%samples_excluded]
design_matrix2<- design_matrix[colnames(count_table2),]

#>>> Expression heatmap with the large gene set

data_to_plot <- table_norm[list_sig_genes_MM,]

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
ha1 <- HeatmapAnnotation(
  Stage = design_matrix$Stage,
  Model = design_matrix$Model,
  col = list(Stage = c("Control"=my_palette1[2],"MGUS"=my_palette1[3],
                       "MM"=my_palette1[4],"sMM"=my_palette1[5]),
             Model = c("Control"=my_palette1[2],"pB2IC"=my_palette2[9])),
  show_annotation_name = TRUE
)

#plot
pdf("DEA/MM_C/expression_heatmap.pdf") 
#jpeg("Results/MM/expression_heatmap_genotype.jpg", width=5000, height=5000, res=600)#

Heatmap(data_to_plot2, name = "DEGs between MM and Control",
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE)



dev.off()

#>>> Expression heatmap core geneset

data_to_plot <- count_table2[core_genes_names$ensembl_gene_id,]

data_to_plot2<- t(scale(t(data_to_plot)))

#plot
pdf("expression_heatmap.pdf") 
#jpeg("Results/MM/expression_heatmap_genotype.jpg", width=5000, height=5000, res=600)#

Heatmap(data_to_plot2, name = "DEGs between different Stages and Control",
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

library(dendextend)
row_dend = as.dendrogram(hclust(dist(data_to_plot_MM)))

#plot
#jpeg("DEA/MM_C/heatmap_binary_sig_genes_MM_Control.jpg", width=5000, height=5000, res=600)
pdf("DEA/MM_C/heatmap_binary_sig_genes_MM_Control.pdf", width = 8, height = 8)
Heatmap(data_to_plot, name = "DEGs between MM and Control",
        col = colorRamp2(c(-1, 0, 1), c("darkblue","gray","darkorange")),
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = row_dend, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE)
dev.off()

##############################
# MM_Bcell by genotype #######
##############################

myContrMatrix = makeContrasts(MYC_C = (MYC_MM - B_Cell),
                              Maf_MYC_C = (Maf_MYC_MM - B_Cell),
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

#save results
write.table(results_MM_Bcell,"DEA/MM_Bcell/100920_dea_MM_Bcell.txt", sep="\t", dec=".", quote=FALSE, row.names = TRUE, col.names = TRUE)
write.table(list_sig_genes_MM,"DEA/MM_Bcell/100920_list_sig_genes_MM_Bcell.txt",sep="\t",dec=".",quote=FALSE)


# get gene symbol

results_MM_symbol <- merge(anno2,results_MM_Bcell,by.x="ensembl_gene_id",by.y=0,all.y=TRUE)
results_MM_symbol2<-results_MM_symbol[!((results_MM_symbol$ensembl_gene_id=="ENSMUSG00000115016") & (results_MM_symbol$mgi_symbol=="Gm33906")),]

write.table(results_MM_symbol2,"DEA/MM_Bcell/100920_dea_MM_Bcell_names.txt", sep="\t", dec=".", quote=FALSE, row.names = TRUE, col.names = TRUE)

list_sig_genes<-as.data.frame(list_sig_genes_MM_Bcell)
colnames(list_sig_genes)<-"ensembl_gene_id"
list_sig_genes_Bcell_names<-merge(anno2,list_sig_genes,by.x="ensembl_gene_id",by.y="ensembl_gene_id",all.y=TRUE)
list_sig_genes_Bcell_names<-list_sig_genes_Bcell_names[!((list_sig_genes_Bcell_names$ensembl_gene_id=="ENSMUSG00000115016") & (list_sig_genes_Bcell_names$mgi_symbol=="Gm33906")),]
write.table(list_sig_genes_Bcell_names,"DEA/MM_Bcell/100920_list_sig_genes_MM_Bcell_names.txt",sep="\t",dec=".",quote=FALSE)




##PLOT EXPRESSION HEATMAP


#>>> Expression heatmap
data_to_plot <- table_norm[list_sig_genes_MM_Bcell,]

data_to_plot2<- t(scale(t(data_to_plot)))

#heatmap annotation

my_palette1 <- brewer.pal(8, "Dark2")
my_palette2 <- brewer.pal(11, "Paired")

ha1 <- HeatmapAnnotation(
  Stage = design_matrix$Stage,
  Model = design_matrix$Model,
  col = list(Stage = c("B_Cell"=my_palette1[1],"Control"=my_palette1[2],"MGUS"=my_palette1[3],
                       "MM"=my_palette1[4]),
             Model = c("B_Cell"=my_palette1[1], "Control"=my_palette1[2],"B2IC"=my_palette2[1],
                       "B2IKC"=my_palette2[2],"CD1B2IC"=my_palette2[3],"cMafB2IC"=my_palette2[4], "Maf_B2IKC"=my_palette2[5],
                       "Maf_MYC"=my_palette2[6],"MYC"=my_palette2[7],
                       "MMsetB2IC"=my_palette2[8], "pB2IC"=my_palette2[9])),
  show_annotation_name = TRUE
)

#plot
pdf("DEA/MM_Bcell/expression_heatmap.pdf") 
#jpeg("DEA/MM_Bcell/expression_heatmap.jpg", width=5000, height=5000, res=600)#

Heatmap(data_to_plot2, name = "DEGs between MM and B_cell",
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE)



dev.off()

#>>> Correlation matrix different models of MM, based on the logFC

FC_MM <- results_MM[,grep("logFC",colnames(results_MM))]
colnames(FC_MM) <- gsub(":logFC","",colnames(FC_MM))

correlation_MM<-cor(FC_MM, method="spearman")

pdf("DEA/MM_C/cor_spearman_corrplot_circle.pdf")
corrplot(correlation_MM, method = "pie", type = "upper", order = "original", tl.col = "black", tl.srt = 45,col=brewer.pal(n=8, name="PuOr"))

dev.off()



#****************************#
# MGUS vs Control all models #
#****************************#
myContrMatrix = makeContrasts(MYC_MGUS_C = (MYC_MGUS - Control),
                              B2IC_MGUS_C = (B2IC_MGUS - Control),
                              B2IKC_MGUS_C = (B2IKC_MGUS - Control),
                              pB2IC_MGUS_C = (pB2IC_MGUS - Control),
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
list_sig_genes_MGUS <- unique(list_sig_genes) #723 significant genes


#save results
write.table(results_MGUS,"DEA/MGUS_C/14092020_dea_MGUS_CONTROL.txt",sep="\t",dec=".",quote=FALSE, row.names = TRUE, col.names = TRUE)
write.table(list_sig_genes_MGUS,"DEA/MGUS_C/14092020_list_sig_genes_MGUS_CONTROL.txt",sep="\t",dec=".",quote=FALSE)

#get the results with gene symbol

results_MGUS_symbol <- merge(anno2,results_MGUS_MM,by.x="ensembl_gene_id",by.y=0,all.y=TRUE)
results_MGUS_symbol2<-results_MGUS_symbol[!((results_MGUS_symbol$ensembl_gene_id=="ENSMUSG00000115016") & (results_MGUS_symbol$mgi_symbol=="Gm33906")),]

write.table(results_MGUS_symbol2,"DEA/MGUS_MM/210920_dea_MGUS_MM_names.txt", sep="\t", dec=".", quote=FALSE, row.names = TRUE, col.names = TRUE)

list_sig_genes<-as.data.frame(list_sig_genes_MGUS_MM)
colnames(list_sig_genes)<-"ensembl_gene_id"
list_sig_genes_MGUS_names<-merge(anno2,list_sig_genes,by.x="ensembl_gene_id",by.y="ensembl_gene_id",all.y=TRUE)

write.table(list_sig_genes_MGUS_names,"DEA/MGUS_MM/210920_list_sig_genes_MGUS_MM_names.txt",sep="\t",dec=".",quote=FALSE)


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
                       "Maf_B2IKC"=my_palette2[5], "Maf_MYC"=my_palette2[6],"MYC"=my_palette2[7],
                       "MMsetB2IC"=my_palette2[8], "pB2IC"=my_palette2[9])),
  show_annotation_name = TRUE
)

#plot
pdf("DEA/MGUS_MM/heatmap_sig_genes_MGUS_MM.pdf")
#jpeg("Results/MGUS/heatmap_sig_genes_MGUS_Control_No_B_cells.jpg",width=5000, height=5000, res=600)

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
  col = list(Model = c("MYC_MGUS_C"=my_palette2[1],
                       "B2IC_MGUS_C"=my_palette2[2],"B2IKC_MGUS_C"=my_palette2[3],"pB2IC_MGUS_C"=my_palette2[4],
                       "CD1B2IC_MGUS_C"=my_palette2[5], "MMsetB2IC_MGUS_C"=my_palette2[6],
                       "cMafB2IC_MGUS_C"=my_palette2[7])),
  show_annotation_name = TRUE
)

#plot
pdf("DEA/MGUS_C/heatmap_binary_sig_genes_MGUS_Control_tendency.pdf",width = 8, height = 8)
Heatmap(data_to_plot_MGUS, name = "DEGs between MGUS and Control",
        col = colorRamp2(c(-1, 0, 1), c("darkblue","gray","darkorange")),
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE)

dev.off()



#--Take the sig genes of MM in MGUS, and see tendency

data_to_plot_merge <- results_MGUS[list_sig_genes_MM,grep("tendency",colnames(results_MGUS))]
colnames(data_to_plot_merge) <- gsub(":genes_tendency","",colnames(data_to_plot_merge))

#plot
pdf("DEA/MGUS_C/heatmap_binary_sig_genes_MM_Control_in_MGUS_2237.pdf",width = 8, height = 8)
Heatmap(data_to_plot_merge, name = "DEGs between MGUS and Control",
        col = colorRamp2(c(-1, 0, 1), c("darkblue","gray","darkorange")),
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = row_dend, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE)
dev.off()

#--Take the intersection


#### 3- Similarities and differencies between MGUS and MM

##Common genes
length(intersect(list_sig_genes_MM,list_sig_genes_MGUS)) #592 -- need to check if the change direction is the same

##Correlation
# with all genes
FC_MM <- results_MM[,grep("logFC",colnames(results_MM))]
colnames(FC_MM) <- gsub(":logFC","",colnames(FC_MM))

FC_MGUS <- results_MGUS[,grep("logFC",colnames(results_MGUS))]
colnames(FC_MGUS) <- gsub(":logFC","",colnames(FC_MGUS))
correlation_MGUS<-cor(FC_MGUS, method="spearman")

pdf("DEA/MGUS_C/cor_spearman_corrplot_MGUS.pdf")
corrplot(correlation_MGUS, method = "number", type = "upper", order = "original", tl.col = "black", tl.srt = 45,col=brewer.pal(n=8, name="PuOr"))

dev.off()

pdf("DEA/MGUS_C/cor_spearman_corrplot_MGUS_circle.pdf")
corrplot(correlation_MGUS, method = "circle", type = "upper", order = "original", tl.col = "black", tl.srt = 45,col=brewer.pal(n=8, name="PuOr"))

dev.off()

FC_all <- merge(FC_MM,FC_MGUS, by.x=0,by.y=0, sort=FALSE)
rownames(FC_all) <- FC_all[,1]
FC_all <- FC_all[,-1]

correlation_FC_pearman<-cor(FC_all, method = "spearman")

pdf("DEA/MM_C/cor_spearman_corrplot_all_number.pdf")
corrplot(correlation_FC_pearman, method = "number", type = "upper", order = "original", tl.col = "black", tl.srt = 45,col=brewer.pal(n=8, name="PuOr"),number.cex=0.70)

dev.off()

pdf("DEA/MM_C/cor_spearman_corrplot_all_circle.pdf")
corrplot(correlation_FC_pearman, method = "circle", type = "upper", order = "original", tl.col = "black", tl.srt = 45,col=brewer.pal(n=8, name="PuOr"))

dev.off()
                                                
               

save.image("DE_MM_C.Rdata")


#human homologous of list_sig_MM_CONTROL
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

genesV2 = getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", 
                 values = list_sig_genes_MM , mart = mouse, 
                 attributesL = c("ensembl_gene_id"), martL = human, uniqueRows=T)

write.table(genesV2,"DEA/MM_C/list_sig_genes_human_homologous.txt",sep="\t")
core_homologous<-genesV2[genesV2$Gene.stable.ID %in% rownames(core_genes),]
write.table(core_homologous,"DEA/MM_C/core_genes_human_homologous.txt",sep="\t")



#****************************#
# MGUS vs MM all models #
#****************************#
myContrMatrix = makeContrasts(MYC_MGUS_MM = (MYC_MGUS - MYC_MM),
                              B2IC_MGUS_MM = (B2IC_MGUS - B2IC_MM),
                              B2IKC_MGUS_MM = (B2IKC_MGUS - B2IKC_MM),
                              pB2IC_MGUS_MM = (pB2IC_MGUS - pB2IC_MM),
                              CD1B2IC_MGUS_MM = (CD1B2IC_MGUS - CD1B2IC_MM),
                              MMsetB2IC_MGUS_MM =(MMsetB2IC_MGUS - MMsetB2IC_MM),
                              cMafB2IC_MGUS_MM =(cMafB2IC_MGUS - cMafB2IC_MM),
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

results_MGUS_MM <- results[,-1]
list_sig_genes_MGUS_MM <- unique(list_sig_genes) #723 significant genes


#save results
write.table(results_MGUS_MM,"DEA/MGUS_MM/21092020_dea_MGUS_MM.txt",sep="\t",dec=".",quote=FALSE, row.names = TRUE, col.names = TRUE)
write.table(list_sig_genes_MGUS,"DEA/MGUS_MM/21092020_list_sig_genes_MGUS_MM.txt",sep="\t",dec=".",quote=FALSE)

#get the results with gene symbol

results_MGUS_symbol <- merge(anno2,results_MGUS,by.x="ensembl_gene_id",by.y=0,all.y=TRUE)
results_MGUS_symbol2<-results_MGUS_symbol[!((results_MGUS_symbol$ensembl_gene_id=="ENSMUSG00000115016") & (results_MGUS_symbol$mgi_symbol=="Gm33906")),]

write.table(results_MGUS_symbol2,"DEA/MGUS_C//140920_dea_MGUS_CONTROL_names.txt", sep="\t", dec=".", quote=FALSE, row.names = TRUE, col.names = TRUE)

list_sig_genes<-as.data.frame(list_sig_genes_MGUS)
colnames(list_sig_genes)<-"ensembl_gene_id"
list_sig_genes_MGUS_names<-merge(anno2,list_sig_genes,by.x="ensembl_gene_id",by.y="ensembl_gene_id",all.y=TRUE)

write.table(list_sig_genes_MGUS_names,"DEA/MGUS_C/140920_list_sig_genes_MGUS_CONTROL_names.txt",sep="\t",dec=".",quote=FALSE)


#>>> Expression heatmap no B-cells

#--REMOVE OUTLIERS AND PCAs
samples_excluded<-c("YFP50_YFP_Cg1_GC_B_cell_S23","YFP51_YFP_Cg1_GC_B_cell_S23","YFP72_YFP_Cg1_GC_B_cell_S24")
count_table2 <- table_norm[,!colnames(table_norm)%in%samples_excluded]
design_matrix2<- design_matrix[colnames(count_table2),]

data_to_plot <- count_table2[list_sig_genes_MGUS_MM,]

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
                       "Maf_B2IKC"=my_palette2[5], "Maf_MYC"=my_palette2[6],"MYC"=my_palette2[7],
                       "MMsetB2IC"=my_palette2[8], "pB2IC"=my_palette2[9])),
  show_annotation_name = TRUE
)

#plot
pdf("DEA/MGUS_MM/heatmap_sig_genes_MGUS_MM.pdf")
#jpeg("Results/MGUS/heatmap_sig_genes_MGUS_Control_No_B_cells.jpg",width=5000, height=5000, res=600)

Heatmap(data_to_plot2, name = "DEGs between MGUS and MM",
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE)
dev.off()

#>>> Differential Expression Analysis "binary" heatmap

#--Taking into account the tendency
data_to_plot_MGUS_MM<- results_MGUS_MM[list_sig_genes_MGUS_MM,grep("tendency",colnames(results_MGUS_MM))]
colnames(data_to_plot_MGUS_MM) <- gsub(":genes_tendency","",colnames(data_to_plot_MGUS_MM))
#-- Only 1 and 0
data_to_plot <- results_MGUS_MM[list_sig_genes_MGUS_MM,grep("sig",colnames(results_MGUS_MM))]
colnames(data_to_plot) <- gsub(":sig_genes","",colnames(data_to_plot))

#heatmap annotation
my_palette2 <- brewer.pal(11, "Paired")

ha1 <- HeatmapAnnotation(
  Model = colnames(data_to_plot_MGUS_MM),
  col = list(Model = c("MYC_MGUS_MM"=my_palette2[1],
                       "B2IC_MGUS_MM"=my_palette2[2],"B2IKC_MGUS_MM"=my_palette2[3],"pB2IC_MGUS_MM"=my_palette2[4],
                       "CD1B2IC_MGUS_MM"=my_palette2[5], "MMsetB2IC_MGUS_MM"=my_palette2[6],
                       "cMafB2IC_MGUS_MM"=my_palette2[7])),
  show_annotation_name = TRUE
)

#plot
pdf("DEA/MGUS_MM/heatmap_binary_sig_genes_MGUS_MM.pdf",width = 8, height = 8)
Heatmap(data_to_plot, name = "DEGs between MGUS and MM",
        col = colorRamp2(c(-1, 0, 1), c("darkblue","gray","darkorange")),
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE)

dev.off()
