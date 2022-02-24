###########################
#COMPARE THREE DATASETS####
###########################

#LOAD DATA

iruña_MM<-read.table("C:/Users/transbio/Desktop/MM_analysis/human_rm_IG/Results/results_MM_annotated.txt",row.names=1,sep="\t")
results_6477<-read.table("C:/Users/transbio/Desktop/MM_analysis/human_rm_IG/Results/20200525_dea_anno_gse6477.txt",row.names=1,sep="\t")
results_47552<-read.table("C:/Users/transbio/Desktop/MM_analysis/human_rm_IG/Results/20200527_dea_anno_gse47552.txt",row.names=1,sep="\t")

#load list sig genes MM
#coger los significativos en cada uno de ellos
mouse_sig <- read.table("C:/Users/transbio/Desktop/MM_analysis/mouse_ensembl/Data/results_MM_ano_complete.txt",sep="\t",dec=".",header=TRUE) 
mouse_sig <- mouse_sig[mouse_sig$signature!=0,] #2356 x 71
list_ref_human_sig <- as.character(unique(mouse_sig$Gene.stable.ID.1[!is.na(mouse_sig$Gene.stable.ID.1)]))
common_ensembl <- which(iruña_MM$ensembl_gene_id%in%list_ref_human_sig)
ref_genes_iruña <- as.character(unique(iruña_MM$ensembl_gene_id[common_ensembl])) #1436 ref genes

common_ensembl_6477<-which(results_6477$ensembl%in%list_ref_human_sig)
ref_genes_6477 <- as.character(results_6477$ensembl[common_ensembl_6477]) #1596

common_ensembl_47552<-which(results_47552$ensembl%in%list_ref_human_sig)
ref_genes_47552 <- as.character(results_47552$ensembl[common_ensembl_47552]) #1375

#concatenate and take unique ones 
merge<-c(ref_genes_iruña,ref_genes_6477,ref_genes_47552)
merge<-unique(merge) #1586

#create the matrix with the directions
iruña_MM_sig<-iruña_MM[iruña_MM$ensembl_gene_id%in%merge,c("MM_C.sig_genes","ensembl_gene_id")]

results_6477_sig<-results_6477[results_6477$ensembl%in%merge,c("MM_NPC.sig_genes","ensembl")]

results_47552_sig<-results_47552[results_47552$ensembl%in%merge,c("MM_NPC.sig_genes","ensembl")]

sig_all<-merge(iruña_MM_sig,results_6477_sig,by.x="ensembl_gene_id",by.y="ensembl",all.x=T,all.y=T)
summary(sig_all)
significant_all<-merge(sig_all,results_47552_sig,by.x="ensembl_gene_id",by.y="ensembl",all.x=T,all.y=T)
#calcular 0s
no_sig_pos<-apply(significant_all[,-1], 1, FUN=function(x){return(sum(x==0)==3)})
sum(!is.na(no_sig_pos))

#filtrar significant all
significant_all_sets<-significant_all[no_sig_pos!=TRUE,]
significant_all_sets<-significant_all_sets[!is.na(significant_all_sets$ensembl_gene_id),]
significant_def<-unique(significant_all_sets)
#significant_def<-significant_def[unique(significant_def$ensembl_gene_id),]

# quietar el duplicado con significancia diferentes
significant_def<-significant_def[-10,]
significant_def<-significant_def[-16,]

#create an annotation matrix for these genes
library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#human gene symbol
anno_biomart = getBM(
  values = significant_def$ensembl_gene_id,
  filters = c("ensembl_gene_id"),
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  mart = human
)



#plotear la tabla normalizada con esos genes

data_to_plot<-significant_def


anno_plot <- anno_biomart[anno_biomart$ensembl_gene_id%in%data_to_plot$ensembl_gene_id,]
rownames(data_to_plot) <- paste(data_to_plot$ensembl_gene_id,as.character(anno_plot$hgnc_symbol),sep=" : ")

data_to_plot<-data_to_plot[,2:4]
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
pdf("merge_datasets.pdf")
draw(Heatmap(data_to_plot, name = "merge_datasets",
             #top_annotation = ha1,
             right_annotation = row_ha,
             show_column_names = TRUE, show_row_names = TRUE, row_names_gp = gpar(fontsize = 6),
             cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE))
dev.off()

##do the same, but with tendencys

iruña_tendency<-iruña_MM[iruña_MM$ensembl_gene_id%in%significant_def$ensembl_gene_id,c("MM_C.genes_tendency","ensembl_gene_id")]
tendency_6477<-unique(results_6477[results_6477$ensembl%in%significant_def$ensembl_gene_id,c("MM_NPC.genes_tendency","ensembl")])
tendency_47552<-unique(results_47552[results_47552$ensembl%in%significant_def$ensembl_gene_id,c("MM_NPC.genes_tendency","ensembl")])

#merge the datasets

merge<-merge(iruña_tendency,tendency_6477,by.x="ensembl_gene_id",by.y="ensembl",all.x=T,all.y=T)

tendency_all<-merge(merge,tendency_47552,by.x="ensembl_gene_id",by.y="ensembl",all.x=T,all.y=T)
tendency_def<-unique(tendency_all)



colnames(tendency_def)<-c('ensembl_gene_id','pamplona','set_6477','set_47552')

#plotear la tabla normalizada con esos genes

data_to_plot<-tendency_good


anno_plot <- anno_biomart[anno_biomart$ensembl_gene_id%in%data_to_plot$ensembl_gene_id,]
rownames(data_to_plot) <- paste(data_to_plot$ensembl_gene_id,as.character(anno_plot$hgnc_symbol),sep=" : ")

data_to_plot<-data_to_plot[,2:4]
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
pdf("merge_datasets_tendency.pdf")
draw(Heatmap(data_to_plot, name = "merge_datasets",
             #top_annotation = ha1,
             #right_annotation = row_ha,
             show_column_names = TRUE, show_row_names = TRUE, row_names_gp = gpar(fontsize = 6),
             cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE))
dev.off()
