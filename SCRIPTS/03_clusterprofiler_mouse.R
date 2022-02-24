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
setwd("C:/Users/transbio/Desktop/MM_analysis/Mouse_filtered")

  dea_res<-read.table("DEA/MM_C/100920_dea_MM_CONTROL.txt", header = T, check.names = F, stringsAsFactors = F, row.names = 1)

#>>> Converting IDs using bitr
genes_ID = bitr(rownames(dea_res), fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb="org.Mm.eg.db")
head(genes_ID)





#>>> GO Gene Set Enrichment Analysis 

Contrast<-c('MYC_C','Maf_MYC_C','B2IC_C','B2IKC_C','pB2IC_C','cMafB2IC_C','CD1B2IC_C','MMsetB2IC_C')



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
    gse<-gseGO(geneList = ranked_genes, OrgDb = org.Mm.eg.db,
               keyType = 'ENSEMBL',
               ont = "BP", nPerm = 1000,
               minGSSize = 20, maxGSSize = 500,
               pvalueCutoff = 1,
               verbose=FALSE)
    
    gse<-setReadable(gse, OrgDb = org.Mm.eg.db, keyType = "auto")
    
    write.table(gse,paste0("GO_",Contrast[i],".txt", sep=""),
                sep="\t",quote=F,row.names=F)
    
}


#>>> KEGG pathway

for(i in 1:length(Contrast)){

model<-dea_res[,grep(Contrast[i],colnames(dea_res))]
colnames(model) <- gsub(paste0("", Contrast[i]),'', colnames(model))
HS_Genes <- merge(model, genesV2, by.x='row.names', by.y='Gene.stable.ID')
ranked_genes<-HS_Genes[,':logFC']
names(ranked_genes)<-HS_Genes$NCBI.gene..formerly.Entrezgene..ID
ranked_genes <- sort(ranked_genes, decreasing = T)
head(ranked_genes)
tail(ranked_genes)

KEGG<-gseKEGG(ranked_genes, organism = "human", keyType = "kegg", nPerm = 1000, minGSSize = 20, maxGSSize = 500, pvalueCutoff = 1)

}
#>>> GMT Gene Set enrichment analysis (First HALLMARK)


gmtfile_H <- read.gmt('C:/Users/transbio/Desktop/MM_analysis/GSEA_files/h.all.v7.1.entrez.gmt') #HALLMARK GENESET

gmtfile_C6 <-read.gmt('C:/Users/transbio/Desktop/MM_analysis/GSEA_files/c6.all.v7.1.entrez.gmt')#Oncogenic gene signatures (C6)

gmtfile_C2 <-read.gmt('C:/Users/transbio/Desktop/MM_analysis/GSEA_files/c2.cgp.v7.2.entrez.gmt')#Oncogenic gene signatures (C6)

##2.1-Ensembl first to human and then entrezid
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

genesV2 = getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", 
                 values = rownames(dea_res) , mart = mouse, 
                 attributesL = c("ensembl_gene_id","entrezgene_id"), martL = human, uniqueRows=T)

for(i in 1:length(Contrast)){
  
  ##2.2-Prepare the list, ordered by logFC
  
  model<-dea_res[,grep(Contrast[i],colnames(dea_res))]
  colnames(model) <- gsub(paste0("", Contrast[i]),'', colnames(model))
  HS_Genes <- merge(model, genesV2, by.x='row.names', by.y='Gene.stable.ID')
  ranked_genes<-HS_Genes[,':logFC']
  names(ranked_genes)<-HS_Genes$NCBI.gene..formerly.Entrezgene..ID
  ranked_genes <- sort(ranked_genes, decreasing = T)
  head(ranked_genes)
  tail(ranked_genes)
  
  ##2.3-Gene Set Enrichment Analysis
  egmt <- GSEA(ranked_genes, TERM2GENE=gmtfile_C2, verbose=T, minGSSize = 10, pvalueCutoff = 1)
  #egmt<-setReadable(egmt, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")h
  
  #write.table(egmt,paste0("HALLMARK_",Contrast[i],".txt", sep=""),
              #sep="\t",quote=F,row.names=F)
  
  ##2.4- plot different ways to visualize the results
  
  #ID <- as.data.frame(egmt)
  #ID <- egmt[, c('Description')]
    #pdf(paste0(Contrast[i],'_GSEA_ALL.pdf'), width = 15, height = 7)
    #for (i in ID){
    #plot <- gseaplot2(egmt, geneSetID = i , title = i, pvalue_table = T) 
      #print(plot)}
    #dev.off()
  
  pdf(paste0(Contrast[i],'_dotplot.pdf'))
  
  dotplot(egmt, showCategory =20, col='pval')
  
  dev.off()
  
  pdf(paste0(Contrast[i],'_ridgeplot.pdf'),width = 15, height = 7 )
  
  ridgeplot(egmt)
  
  dev.off()
  
}

#>>> GSEA with human DEA results

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

genesV2 = getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", 
                 values = rownames(dea_res) , mart = mouse, 
                 attributesL = c("ensembl_gene_id","entrezgene_id"), martL = human, uniqueRows=T)

gmtfile_human <- read.gmt('C:/Users/transbio/Desktop/MM_analysis/GSEA_files/MM_C_human_def.txt')


for(i in 1:length(Contrast)){
  
  ##2.2-Prepare the list, ordered by logFC
  
  model<-dea_res[,grep(Contrast[i],colnames(dea_res))]
  colnames(model) <- gsub(paste0("", Contrast[i]),'', colnames(model))
  HS_Genes <- merge(model, genesV2, by.x='row.names', by.y='Gene.stable.ID')
  ranked_genes<-HS_Genes[,':logFC']
  names(ranked_genes)<-HS_Genes$Gene.stable.ID.1
  ranked_genes <- sort(ranked_genes, decreasing = T)
  head(ranked_genes)
  tail(ranked_genes)
  
  ##2.3-Gene Set Enrichment Analysis
  egmt <- GSEA(ranked_genes, TERM2GENE=gmtfile_human, verbose=T, minGSSize = 1, maxGSSize = 10000, pvalueCutoff = 1)
  #egmt<-setReadable(egmt, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
  
  write.table(egmt,paste0("GSEA/Human_comp_new",Contrast[i],".txt", sep=""),
  sep="\t",quote=F,row.names=F)
  
  ##2.4- plot different ways to visualize the results
  
  ID <- as.data.frame(egmt)
  ID <- egmt[, c('Description')]
  pdf(paste0(Contrast[i],'_GSEA_plot_all_human_new.pdf'), width = 15, height = 7)
  for (i in ID){
  plot <- gseaplot2(egmt, geneSetID = i , title = i, pvalue_table = T) 
  print(plot)}
  dev.off()
  
  #pdf(paste0(Contrast[i],'_dotplot_new.pdf'))
  
  #dotpot(egmt, showCategory =20, col='pval')
  
  #dev.off()
  
  }
