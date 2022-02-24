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

design_matrix$Sample <-  as.factor(design_matrix$Sample)
design_matrix$Stage<-as.factor(design_matrix$Stage)
design_matrix$Combined<-as.factor(design_matrix$Combined)
design_matrix$Model<-as.factor(design_matrix$Model)
str(design_matrix)

#### 2- Differenital Expression Analysis, by genotype

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



##############################
######## MM vs Control #######
##############################
myContrMatrix = makeContrasts(MM_C = (MM - Control),
                              MGUS_C = (MGUS - Control),
                              MM_MGUS = (MM - MGUS),
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
list_sig_genes_MM <- unique(list_sig_genes)

write.table(results_MM,"DEA/MM_C/100920_dea_Stages.txt", sep="\t", dec=".", quote=FALSE, row.names = TRUE, col.names = TRUE)
write.table(list_sig_genes_MM,"DEA/MM_C/100920_list_sig_genes_Stages.txt",sep="\t",dec=".",quote=FALSE)


samples_excluded<-c("YFP50_YFP_Cg1_GC_B_cell_S23","YFP51_YFP_Cg1_GC_B_cell_S23","YFP72_YFP_Cg1_GC_B_cell_S24")
count_table2 <- table_norm[,!colnames(table_norm)%in%samples_excluded]
design_matrix2<- design_matrix[colnames(count_table2),]

#>>> Expression heatmap with the large gene set

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
pdf("DEA/MM_C/expression_heatmap_STAGE.pdf") 
#jpeg("Results/MM/expression_heatmap_genotype.jpg", width=5000, height=5000, res=600)#

Heatmap(data_to_plot2, name = "DEGs between MM and Control",
        top_annotation = ha1, 
        show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE)



dev.off()


