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


setwd("C:/Users/transbio/Desktop/MM_analysis/mouse_ensembl/Data")

# Differential expression -------------------------------------------------

### LIMMA
### Uses alinear model. Needs to work with log transformed data, or if not you
### should use voom transformation. Very useful for multifactorial design.

# Norm counts table (Log2 transformed with VOOM)
table_norm <-  read.table("final_norm_def.txt", sep="\t", dec=".", row.names=1, check.names = FALSE)


# Experimental design matrix
design_matrix <- read.csv('Col_Data.csv', sep=";", header=T, stringsAsFactors = TRUE)
colnames(design_matrix)[1] <- 'Sample'
summary(design_matrix)
rownames(design_matrix)<- design_matrix$Sample

design_matrix<- design_matrix[colnames(table_norm),]


#### 1- Exploration of normdata (MDS)  ##Hazlo con tu código y con los colorines tan bonitos que generaste!! :)
#Ejemplo de prueba con PCA
pca_table_norm <- prcomp(t(table_norm))
plot(pca_table_norm$x[,1],pca_table_norm$x[,2], col=design_matrix$Combined)



#### 2- Differenital Expression Analysis
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
myContrMatrix = makeContrasts(MIC_C = (MIC_MM - Control),
                              Maf_MIC_C = (Maf_MIC_MM - Control),
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

#save results
write.table(results_MM,"dea_MM_CONTROL.txt",sep="\t",dec=".",quote=FALSE, row.names = TRUE, col.names = TRUE)
write.table(list_sig_genes_MM,"list_sig_genes_MM_CONTROL.txt",sep="\t",dec=".",quote=FALSE)



#>>> Expression heatmap
data_to_plot <- count_table2[mouse_names,]

data_to_plot2<-scale(data_to_plot, center = TRUE, scale = apply(data_to_plot, 1, sd, na.rm = TRUE)/apply(data_to_plot,1,sd,na.rm=T))

scale = apply(data_to_plot, 2, sd, na.rm = TRUE)/apply(data_to_plot,1,sd,na.rm=T)
#scaled by rows (z-score)
data_to_plot2 <-  (data_to_plot - rowMeans(data_to_plot))


data_to_plot2<- t(scale(t(data_to_plot)))

# take out b_cells
#--REMOVE OUTLIERS AND PCAs
samples_excluded<-c("YFP50_YFP_Cg1_GC_B_cell_S23","YFP51_YFP_Cg1_GC_B_cell_S23","YFP72_YFP_Cg1_GC_B_cell_S24")
count_table2 <- table_norm[,!colnames(table_norm)%in%samples_excluded]
design_matrix2<- design_matrix[colnames(count_table2),]

#heatmap annotation
my_palette1 <- brewer.pal(8, "Dark2")
my_palette2 <- brewer.pal(11, "Paired")

ha1 <- HeatmapAnnotation(
  Stage = design_matrix2$Stage,
  Model = design_matrix2$Model,
  col = list(Stage = c("Control"=my_palette1[1],"MGUS"=my_palette1[3],
                    "MM"=my_palette1[4]),
             Model = c("Control"=my_palette1[1],"B2IC"=my_palette2[1],
                       "B2IKC"=my_palette2[2],"CD1B2IC"=my_palette2[3],"cMafB2IC"=my_palette2[4],
                       "Maf_B2IKC"=my_palette2[5], "Maf_MIC"=my_palette2[6],"MIC"=my_palette2[7],
                       "MMsetB2IC"=my_palette2[8], "pB2IC"=my_palette2[9])),
  show_annotation_name = TRUE
)

#plot
png("heatmap_sig_human_in_mouse_scale_no_Bcell_no_cluster.png", width=1000, height=1000)
Heatmap(data_to_plot2, name = "DEGs between MM and Control",
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, cluster_columns = FALSE, row_gap = unit(3, "mm"), use_raster=TRUE)
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
                       "Maf_MIC_C"=my_palette2[6],"MIC_C"=my_palette2[7],
                       "MMsetB2IC_C"=my_palette2[8], "pB2IC_C"=my_palette2[9])),
  show_annotation_name = TRUE
)

#plot
png("20200506_heatmap_binary_sig_genes_MM_Control.png", width=1000, height=1000)
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




#****************************#
# MGUS vs Control all models #
#****************************#
myContrMatrix = makeContrasts(MIC_MGUS_C = (MIC_MGUS - Control),
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
write.table(results_MGUS,"dea_MGUS_CONTROL.txt",sep="\t",dec=".",quote=FALSE, row.names = TRUE, col.names = TRUE)
write.table(list_sig_genes_MGUS,"list_sig_genes_MGUS_CONTROL.txt",sep="\t",dec=".",quote=FALSE)



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
                       "MM"=my_palette1[4]),
             Model = c("B_Cell"=my_palette1[1],"Control"=my_palette1[2],"B2IC"=my_palette2[1],
                       "B2IKC"=my_palette2[2],"CD1B2IC"=my_palette2[3],"cMafB2IC"=my_palette2[4],
                       "Maf_B2IKC"=my_palette2[5], "Maf_MIC"=my_palette2[6],"MIC"=my_palette2[7],
                       "MMsetB2IC"=my_palette2[8], "pB2IC"=my_palette2[9])),
  show_annotation_name = TRUE
)

#plot
png("heatmap_sig_genes_MGUS_Control.png", width=1000, height=1000)
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
png("20200506_heatmap_binary_sig_genes_MGUS_Control_tendency_ordered.png", width=1000, height=1000)
Heatmap(data_to_plot_MGUS, name = "DEGs between MGUS and Control",
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
corrplot(correlation_MGUS, type = "upper", order = "original", 
         tl.col = "black", tl.srt = 45,col=brewer.pal(n=8, name="PuOr"))
dev.off()
corrplot.mixed(correlation_FC_spearman,lower.col = "black", number.cex = .7,col=brewer.pal(n=8, name="PuOr"))
corrplot.mixed(correlation_FC_spearman, lower = "number", upper = "circle", tl.pos = "lt")                                                            
               

setwd("C:/Users/transbio/Desktop/MM_analysis/mouse_results")
save.image("DE_MM_C.Rdata")
