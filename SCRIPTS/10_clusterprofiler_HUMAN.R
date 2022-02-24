################################################################
###### GENE ONTOLOGY ANALYSIS: CLUSTER PROFILER MOUSE DATA #####
################################################################
#Source: https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html#go-analysis


#>>> LOAD LIBRARIES

library('clusterProfiler')
library('org.Hs.eg.db') # Human organism
library('org.Mm.eg.db') # Mouse organism
library("GSEABase")
library('DOSE')
library('enrichplot')
library(clusterProfiler)
library(biomaRt)


#>>> WORK DIRECTORY AND DATA
set.seed(1234567)
setwd("C:/Users/transbio/Desktop/MM_human")

dea_res<-read.table("20200527_dea_anno_gse47552.txt", header = T, check.names = F, stringsAsFactors = F, row.names = 1)

#>>> Converting IDs using bitr
genes_ID = bitr(rownames(dea_res), fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb="org.Mm.eg.db")
head(genes_ID)





#>>> GO Gene Set Enrichment Analysis 

Contrast<-c('MM_C','MGUS_C','MM_MGUS')

for(i in 1:length(Contrast)){
  ##2.1-Prepare geneList. ordered by logFC
  
  model<-dea_res[,grep(Contrast[i],colnames(dea_res))]
  colnames(model) <- gsub(paste0("", Contrast[i]),'', colnames(model))
  ranked_genes<-model[,':logFC']
  names(ranked_genes)<-rownames(dea_res)
  ranked_genes <- sort(ranked_genes, decreasing = T)
  head(ranked_genes)
  tail(ranked_genes)
  
  ##2.2 Gene Set Enrichment Analysis
  gse<-gseGO(geneList = ranked_genes, OrgDb = org.Hs.eg.db,
             keyType = 'ENSEMBL',
             ont = "BP", nPerm = 1000,
             minGSSize = 20, maxGSSize = 500,
             pvalueCutoff = 1,
             verbose=FALSE)
  
  gse<-setReadable(gse, OrgDb = org.Hs.eg.db, keyType = "auto")
  
  write.table(gse,paste0("GO_",Contrast[i],".txt", sep=""),
              sep="\t",quote=F,row.names=F)
  
}

#>>> GMT Gene Set enrichment analysis (First HALLMARK)


gmtfile_H <- read.gmt('C:/Users/transbio/Desktop/MM_analysis/GSEA_files/h.all.v7.1.entrez.gmt')

##2.1-Ensembl first to human and then entrezid
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

anno2 = getBM(
  values = rownames(dea_res),
  filters = c("ensembl_gene_id"),
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  mart = human
)



for(i in 1:length(Contrast)){
  
  ##2.2-Prepare the list, ordered by logFC
  
  model<-dea_res[,grep(Contrast[i],colnames(dea_res))]
  colnames(model) <- gsub(paste0("", Contrast[i]),'', colnames(model))
  HS_Genes <- merge(model, anno2, by.x='row.names', by.y='ensembl_gene_id')
  ranked_genes<-HS_Genes[,':logFC']
  names(ranked_genes)<-HS_Genes$entrezgene_id
  ranked_genes <- sort(ranked_genes, decreasing = T)
  head(ranked_genes)
  tail(ranked_genes)
  
  ##2.3-Gene Set Enrichment Analysis
  egmt <- GSEA(ranked_genes, TERM2GENE=gmtfile_H, verbose=T, minGSSize = 10, pvalueCutoff = 1)
  egmt<-setReadable(egmt, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
  
  #write.table(egmt,paste0("MOUSE_GSEA_",Contrast[i],".txt", sep=""),
  #sep="\t",quote=F,row.names=F)
  
  ##2.4- plot different ways to visualize the results
  
  ID <- as.data.frame(egmt)
  ID <- egmt[, c('Description')]
  #pdf(paste0(Contrast[i],'_GSEA_ALL.pdf'), width = 15, height = 7)
  #for (i in ID){
  #plot <- gseaplot2(egmt, geneSetID = i , title = i, pvalue_table = T) 
  #print(plot)}
  #dev.off()
  
  #pdf(paste0(Contrast[i],'_dotplot.pdf'))
  
  #dotpot(egmt, showCategory =20, col='pval')
  
  #dev.off()
  
}

setwd("C:/Users/transbio/Desktop/MM_human/GSEA")

#>>> GSEA with mouse DEA results


gmtfile_mouse <- read.gmt('C:/Users/transbio/Desktop/MM_analysis/GSEA_files/Mouse_filtered_complet_geneset.gmt.txt')


for(i in 1:length(Contrast)){
  
  ##2.2-Prepare the list, ordered by logFC
  
  model<-dea_res[,grep(Contrast[i],colnames(dea_res))]
  colnames(model) <- gsub(paste0("", Contrast[i]),'', colnames(model))
  ranked_genes<-model[,":logFC"]
  names(ranked_genes)<-dea_res$ensembl
  ranked_genes <- sort(ranked_genes, decreasing = T)
  head(ranked_genes)
  tail(ranked_genes)
  
  ##2.3-Gene Set Enrichment Analysis
  egmt <- GSEA(ranked_genes, TERM2GENE=gmtfile_mouse, verbose=T, minGSSize = 1, maxGSSize = 10000, pvalueCutoff = 1)
  #egmt<-setReadable(egmt, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
  
  write.table(egmt,paste0("Mouse_comp_",Contrast[i],".txt", sep=""),
              sep="\t",quote=F,row.names=F)
  
  ##2.4- plot different ways to visualize the results
  
  ID <- as.data.frame(egmt)
  ID <- egmt[, c('Description')]
  pdf(paste0(Contrast[i],'_GSEA_plot_all_mouse.pdf'), width = 15, height = 7)
  for (i in ID){
  plot <- gseaplot2(egmt, geneSetID = i , title = i, pvalue_table = T) 
  print(plot)}
  dev.off()
  
  pdf(paste0(Contrast[i],'_dotplot_mouse.pdf'))
  
  dotplot(egmt, showCategory =20, col='pval')
  
  dev.off()
  
  pdf(paste0(Contrast[i],'_ridgeplot_mouse.pdf'))
  
  ridgeplot(egmt)
  
  dev.off()
  
}

setwd("C:/Users/transbio/Desktop/MM_human")
#>>> GSEA with metanalysis results

metadata<-read.table("metadata_updated.txt",sep="\t")
res_integrated<-merge(dea_res,metadata, by.x=0, by.y=0)

  
  ##2.2-Prepare the list, ordered by logFC
  
  model<-res_integrated[,grep(Contrast[i],colnames(dea_res))]
  colnames(model) <- gsub(paste0("", Contrast[i]),'', colnames(model))
  ranked_genes<-res_integrated[,"meta_adj"]
  names(ranked_genes)<-res_integrated$Row.names
  ranked_genes <- sort(ranked_genes, decreasing = T)
  head(ranked_genes)
  tail(ranked_genes)
  
  ##2.3-Gene Set Enrichment Analysis
  egmt <- GSEA(ranked_genes, TERM2GENE=gmtfile_mouse, verbose=T, minGSSize = 1, maxGSSize = 10000, pvalueCutoff = 1)
  #egmt<-setReadable(egmt, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
  
  write.table(egmt,paste0("Integrated_vs_mouse.txt", sep=""),
              sep="\t",quote=F,row.names=F)
  
  ##2.4- plot different ways to visualize the results
  
  ID <- as.data.frame(egmt)
  ID <- egmt[, c('Description')]
  pdf('GSEA/integrated_vs_mouse_GSEA_plot_all.pdf', width = 15, height = 7)
  for (i in ID){
    plot <- gseaplot2(egmt, geneSetID = i , title = i, pvalue_table = T) 
    print(plot)}
  dev.off()
  
  pdf('GSEA/dotplot_integrated.pdf')
  
  dotplot(egmt, showCategory =20, col='pval')
  
  dev.off()
  
  pdf('GSEA/ridgeplot_integrated.pdf')
  
  ridgeplot(egmt)
  
  dev.off()

