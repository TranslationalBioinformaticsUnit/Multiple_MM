################
#ORDER NAMES####
################

library(limma)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(corrplot)
library(biomaRt)


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
list_sig_genes_MM <- unique(list_sig_genes)
  
  
pos_sig <- grep("sig",colnames(results_MM))
results_MM$signature_description <- as.factor(paste(results_MM[,pos_sig[1]],results_MM[,pos_sig[2]],results_MM[,pos_sig[3]],results_MM[,pos_sig[4]],results_MM[,pos_sig[5]],results_MM[,pos_sig[6]],results_MM[,pos_sig[7]],results_MM[,pos_sig[8]], sep=""))
results_MM$signature <- as.factor(apply(as.matrix(results_MM[,pos_sig]),1,FUN=function(x){
  if(sum(x!=0)>0){
    up <- sum(x==1)
    dw <- sum(x==-1)
    if(up>dw) { ss <- 1}
    if(dw>up) { ss <- -1}
  }
  if(sum(x!=0)==0){ ss <- 0 }
  return(ss)
}))

library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")


#mouse gene symbol
anno_biomart = getBM(
  values = rownames(results_MM),
  filters = c("ensembl_gene_id"),
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  mart = mouse
)

merge1 <- merge(anno_biomart,results_MM,by.x="ensembl_gene_id",by.y=0,all.y=TRUE)

#el merge me permite unir el results_MM y anno_biomart en base al ensembl_gene_id

human_correspondance = getLDS(attributes = c("ensembl_gene_id","mgi_symbol"), filters = "ensembl_gene_id",
                              values = rownames(results_MM) , mart = mouse,
                              attributesL = c("ensembl_gene_id","hgnc_symbol"), martL = human)
merge2 <- merge(merge1,human_correspondance, by.x="ensembl_gene_id", by.y="Gene.stable.ID", all.x=TRUE)




