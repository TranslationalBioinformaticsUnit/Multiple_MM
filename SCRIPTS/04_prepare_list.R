##################################################
###get moouse signature to compare to human GSEA##
#################################################


library(limma)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(corrplot)
library(biomaRt)

set.seed(123456)

setwd("C:/Users/transbio/Desktop/MM_analysis/Mouse_filtered")

#>>> LOAD DEA RESULTS

results_MM<-read.table('DEA/MM_C/100920_dea_MM_CONTROL.txt',sep="\t", dec=".", row.names=1, check.names = FALSE)

HS_names<-read.table('DEA/MM_C/list_sig_genes_human_homologous.txt',sep="\t", dec=".", row.names=1, check.names = FALSE)

#>>> GET THE GENES THAT ARE UP AND DOWN IN ALL MODELS

results_MM_sig <- results_MM[,grep("sig_genes",colnames(results_MM))]

sig_all <- results_MM_sig[apply(results_MM_sig,1,FUN=function(x){return(sum(x!=0)==8)}),]

sig_all_up<-sig_all[which(rowSums(sig_all)==8),]

sig_all_up_human<-HS_names[HS_names$Gene.stable.ID%in% rownames(sig_all_up),]
write.table(sig_all_up_human,"DEA/MM_C/GENE_lists/list_sig_all_up.txt",sep="\t",dec=".",quote=FALSE)


sig_all_down<-sig_all[which(rowSums(sig_all)== -8),]

sig_all_down_human<-HS_names[HS_names$Gene.stable.ID%in% rownames(sig_all_down),]
write.table(sig_all_down_human,"DEA/MM_C/GENE_lists/list_sig_all_down.txt",sep="\t",dec=".",quote=FALSE)

#>>> GET THE GENES THAT ARE UP AND DOWN IN AT LEAST ONE MODEL

list_sig_genes_up<-read.table("DEA/MM_C/100920_list_upregulated_genes_MM_CONTROL.txt",sep="\t", dec=".", row.names=1, check.names = FALSE)
list_sig_genes_up_human<-HS_names[HS_names$Gene.stable.ID%in% list_sig_genes_up$ensembl_gene_id,]
write.table(list_sig_genes_up_human,"DEA/MM_C/GENE_lists/list_sig_up.txt",sep="\t",dec=".",quote=FALSE)


list_sig_genes_down<-read.table("DEA/MM_C/100920_list_downregulated_genes_MM_CONTROL.txt",sep="\t", dec=".", row.names=1, check.names = FALSE)
list_sig_genes_down_human<-HS_names[HS_names$Gene.stable.ID%in% list_sig_genes_down$ensembl_gene_id,]
write.table(list_sig_genes_down_human,"DEA/MM_C/GENE_lists/list_sig_down.txt",sep="\t",dec=".",quote=FALSE)

#>>> GET THE LIST FOR EACH MODEL (sig genes, sig up genes and sig down genes)

Contrast<-c('MYC_C','Maf_MYC_C','B2IC_C','B2IKC_C','pB2IC_C','cMafB2IC_C','CD1B2IC_C','MMsetB2IC_C')

for(i in 1:length(Contrast)){

  model<-results_MM[,grep(Contrast[i],colnames(results_MM))]
  colnames(model) <- gsub(paste0("", Contrast[i]),'', colnames(model))
  
  # all genes
  sig_genes<-rownames(model[model$`:sig_genes`!=0,])

  sig_genes_Human<-HS_names[HS_names$Gene.stable.ID %in% sig_genes,]

  write.table(sig_genes_Human,paste0("DEA/MM_C/GENE_lists/sig_",Contrast[i],".txt", sep=""),sep="\t",dec=".",quote=FALSE)
  
  #upregulated genes
  
  sig_up<-rownames(model[model$`:sig_genes`==1,])
  
  sig_genes_up_Human<-HS_names[HS_names$Gene.stable.ID %in% sig_up,]
  
  write.table(sig_genes_up_Human,paste0("DEA/MM_C/GENE_lists/sig_upregulated_",Contrast[i],".txt", sep=""),sep="\t",dec=".",quote=FALSE)
  
  #downregulated genes
  
  sig_down<-rownames(model[model$`:sig_genes`== -1,])
  
  sig_genes_down_Human<-HS_names[HS_names$Gene.stable.ID %in% sig_down,]
  
  write.table(sig_genes_down_Human,paste0("DEA/MM_C/GENE_lists/sig_downregulated_",Contrast[i],".txt", sep=""),sep="\t",dec=".",quote=FALSE)

}
  
#