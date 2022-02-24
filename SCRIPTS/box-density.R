library(ggplot2)
library(hrbrthemes)
library(RColorBrewer)

setwd("C:/Users/transbio/Desktop/MM/MM_analysis/Mouse_filtered/Data")


table_norm<-read.table("norm_batch_corrected.txt", sep="\t", dec=".", row.names=1, check.names = FALSE)
colnames(table_norm) <- gsub('X', '', colnames(table_norm))
colnames(table_norm)

target_genes <- "SIX5"

stage<-design_matrix$Stage

for(i in 1:length(target_genes)){
  gene <- 'ENSMUSG00000006494'
  gene_data <- t(table_norm[gene,])
  
  data <- data.frame(gene=gene_data, stage=stage)
  
  #HISTOGRAM  
  p_hist <- ggplot(data, aes(x=gene_data, group=stage, fill=stage)) +
    geom_density(adjust=1.5, alpha=.6) +
    ggtitle(paste(target_genes[i]," - normalized gene expression", sep="")) +
    #xlim(min(data$)-2,max(data$CD9)+2)+
    scale_fill_manual(values=my_palette1)
    #theme_ipsum()
  
  ggsave(paste(target_genes[i],"_rna_hist.png",sep=""))
  
  
  #BOXPLOT
  p_box <- ggplot(data, aes(x=stage,y=SIX5, fill=stage)) + 
    geom_boxplot() +
    ggtitle(paste(target_genes," - normalized gene expression", sep="")) +
    scale_fill_manual(values=my_palette1)+
    #theme_ipsum()+
    theme(axis.text.x = element_text(size=10, angle = 45))
  
  ggsave(paste(target_genes,"_rna_boxplot.pdf",sep=""))
  
} 


#NURIA

stage <- factor(design_matrix$Stage, levels=c("Control","MGUS","sMM","MM"))

target_genes<-"Myc"


library(biomaRt)
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
anno2 = getBM(
  values = rownames(table_norm),
  filters = c("ensembl_gene_id"),
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  mart = mouse
)


for(i in 1:length(target_genes)){
  gene <- anno2$ensembl_gene_id[anno2$mgi_symbol==target_genes[i]]
  gene_data <- table_norm[gene,]
  
  data <- data.frame(gene=gene_data, stage=stage)

  #HISTOGRAM  
  p_hist <- ggplot(data, aes(x=gene, group=stage, fill=stage)) +
    geom_density(adjust=1.5, alpha=.6) +
    ggtitle(paste(target_genes[i]," - normalized gene expression", sep="")) +
    #xlim(min(data$gene)-2,max(data$gene)+2)+
    scale_fill_manual(values=my_palette1)+
    #theme_ipsum()
  
  ggsave(paste(target_genes[i],"_rna_hist.png",sep=""))
  
  
  #BOXPLOT
  p_box <- ggplot(data, aes(x=stage, y=ENSMUSG00000006494, fill=stage)) + 
    geom_boxplot() +
    ggtitle(paste(target_genes[i]," - normalized gene expression", sep="")) +
    scale_fill_manual(values=my_palette1)+
    #theme_ipsum()+
    theme(axis.text.x = element_text(size=10, angle = 45))
  
  ggsave(paste(target_genes[i],"_rna_boxplot.png",sep=""))
}

p_box <- ggplot(data, aes(x=cell_type, y=gene, fill=cell_type)) + 
  geom_boxplot() +
  ggtitle(paste(target_genes[i]," - normalized gene expression", sep="")) +
  scale_fill_manual(values=my_palette1)+
  theme_ipsum()+
  theme(axis.text.x = element_text(size=10, angle = 45))

ggsave(paste(target_genes[i],"_rna_boxplot.png",sep=""))


norm_data<-merge(table_norm, anno2, by.x=0, by.y="ensembl_gene_id")

num_columns=1
num_rows=1
num_pages=1
val=num_columns*num_rows*num_pages

vplots=list()

target_genes<-c("Six5","Tmem198", "Pdk1", "Bhlha15", "Esr2", "Fkbp2", "Itm2c")
for(i in 1:length(target_genes)){
  
  gene <- anno2$ensembl_gene_id[anno2$mgi_symbol==target_genes[i]]
  gene_data <- t(table_norm[gene,])
  data <- data.frame(gene=gene_data, stage=stage)
  
  colnames(data)<-c("gene", "stage")

  p<- ggplot(data, aes(x=stage, y=gene, fill=stage)) + 
    geom_boxplot() +
    scale_fill_manual(values=my_palette1)+
    labs(x = "Cell_type", y = target_genes[i])+
    ggtitle(paste(target_genes[i]," - normalized gene expression", sep=""))
  
  
  p
  vplots[[i]]=p 
  
}



pdf("Six5_related_genes_boxplot.pdf")
marrangeGrob(vplots, nrow = 1, ncol = 1)
dev.off()
