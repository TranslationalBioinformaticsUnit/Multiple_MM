#####################
##metanalisis########
######################
library(dplyr)
setwd("C:/Users/transbio/Desktop/MM_human")
#--Load the three datasets

gse47552 <- read.table("C:/Users/transbio/Desktop/MM_analysis/human_rm_IG/Results/20200527_dea_anno_gse47552.txt", sep="\t",dec=".",header=TRUE)
gse6477 <- read.table("C:/Users/transbio/Desktop/MM_analysis/human_rm_IG/Results/20200525_dea_anno_gse6477.txt",row.names=1,sep="\t",dec=".",header=TRUE)
iruna <- read.table("C:/Users/transbio/Desktop/MM_human/dea_MM_control_anno.txt", sep="\t",dec=".", header=TRUE,row.names=1)
#dataset_4<-read.table("C:/Users/transbio/Desktop/MM_new_dataset/dataset_4.txt",sep="\t",dec=".", header=TRUE,row.names=1)

#take the most significant in case of duplicates
gse47552_filt<-gse47552 %>% group_by(ensembl) %>% filter((MM_NPC.P.Value)==min(MM_NPC.P.Value))
gse6477_filt<-gse6477 %>% group_by(ensembl) %>% filter((MM_NPC.P.Value)==min(MM_NPC.P.Value))
iruna_filt<-iruna[!(iruna$ensembl_gene_id=="ENSG00000276085" & iruna$hgnc_symbol=="CCL3L3"),]
#dataset4_filt<-dataset_4[!(dataset_4$ensembl_gene_id=="ENSG00000187510" & dataset_4$hgnc_symbol=="C12orf74"),] 


#--take the pvalues and the values of interest

iruna_pval<-iruna_filt[,c("MM_C.P.Value","ensembl_gene_id")]
gse47552_pval<-gse47552_filt[,c("ensembl","MM_NPC.P.Value")]
gse6477_pval<-gse6477_filt[,c("ensembl","MM_NPC.P.Value")]
#dataset4_pval<-dataset4_filt[,c("ensembl_gene_id","P.Value")]

#--merge all the pvalues

data_merged <- merge(gse47552_pval,gse6477_pval,by.x="ensembl",by.y="ensembl", sort=FALSE, all.x=TRUE, all.y=TRUE) #2,064 x 7
pvalues <- merge(data_merged,iruna_pval, by.x="ensembl",by.y="ensembl_gene_id", sort=FALSE, all.x=TRUE, all.y=TRUE) #2,456 x 9
#pvalues<- merge(merged2,dataset4_pval,by.x="ensembl",by.y="ensembl_gene_id", sort=FALSE, all.x=TRUE, all.y=TRUE)

#--quitar los duplicados, quedarnos con los más significativos

#ordenarlos por ensembl and pvalue

pvalues =  pvalues[order(pvalues$ensembl),]
#pvalues<-pvalues[,]
rownames(pvalues)<-pvalues$ensembl

pvalues_simple<- pvalues %>% group_by(ensembl) %>% filter((MM_NPC.P.Value.x)==min(MM_NPC.P.Value.x))  
pvalues_def<- pvalues_simple %>% group_by(ensembl) %>% filter((MM_NPC.P.Value.y)==min(MM_NPC.P.Value.y))                                                                                                       


#pvalues_def2<-pvalues[which(pvalues$MM_NPC.P.Value.x<0.05 | pvalues$MM_NPC.P.Value.y<0.05 |pvalues$MM_C.P.Value<0.05),]

#rownames(pvalues_def2)<-pvalues_def2$ensembl
##metadata
#install.packages("compute.es")
require(compute.es)

# Meta-Analysis with Mean Differences:
 #MAd package: http://CRAN.R-project.org/package=MAd
require(metafor)
require("MAd")
# http://www.creative-wisdom.com/teaching/WBI/es.shtml
##### COMBINING P-VALUES

pvalues_meta <- matrix(NA,nrow(pvalues),4) 
rownames(pvalues_meta) <- pvalues$ensembl

pvalues_meta[,1]<-pvalues$MM_NPC.P.Value.x
pvalues_meta[,2]<-pvalues$MM_NPC.P.Value.y
pvalues_meta[,3]<-pvalues$MM_C.P.Value
#pvalues_meta[,4]<-pvalues$P.Value

for(i in 1:nrow(pvalues_meta))
{#i<-1
  effect_iruna <- as.numeric(t_to_d(t=iruna_filt[iruna_filt$ensembl_gene_id%in%rownames(pvalues)[i],"MM_C.t"],n.1=43, n.2=7)[1])#, n.1=51, n.2=16
  effect_gse47552 <- as.numeric(t_to_d(t=gse47552_filt[gse47552_filt$ensembl%in%rownames(pvalues)[i],"MM_NPC.t"],n.1=41, n.2=5)[1]) # n.1=12, n.2=12)[1])
  effect_gse6477 <- as.numeric(t_to_d(t=gse6477_filt[gse6477_filt$ensembl%in%rownames(pvalues)[i],"MM_NPC.t"],n.1=75, n.2=15)[1])
  #effect_dataset4<- as.numeric(t_to_d(t=dataset4_filt[dataset4_filt$ensembl_gene_id%in%rownames(pvalues)[i],"t"],n.1=36, n.2=3)[1])
  var_iruna<- as.numeric(t_to_d(t=iruna_filt[iruna_filt$ensembl_gene_id%in%rownames(pvalues)[i],"MM_C.t"], n.1=43, n.2=7)[2])
  var_gse47552<- as.numeric(t_to_d(t=gse47552_filt[gse47552_filt$ensembl%in%rownames(pvalues)[i],"MM_NPC.t"],n.1=41, n.2=5)[2]) #  n.1=12, n.2=12)[2])
  var_gse6477<- as.numeric(t_to_d(t=gse6477_filt[gse6477_filt$ensembl%in%rownames(pvalues)[i],"MM_NPC.t"],n.1=75, n.2=15)[2])
  #var_dataset4<- as.numeric(t_to_d(t=dataset4_filt[dataset4_filt$ensembl_gene_id%in%rownames(pvalues)[i],"t"],n.1=36, n.2=3)[2])
  
  pvalues_meta[i,4] <- rma(yi=c(effect_iruna,effect_gse47552,effect_gse6477),
                            vi=c(var_iruna,var_gse47552,var_gse6477),weighted=F,
                            method="REML")$pval
}

metadata<-as.data.frame(pvalues_meta)
colnames(metadata)<-c('gse47552','gse6477','iruña','meta')
metadata$meta_adj<-p.adjust(metadata$meta,method="BH")



# Mouse signature ------------------------------------------------------------------

mouse_sig <- read.table("C:/Users/transbio/Desktop/MM_analysis/mouse_ensembl/Data/results_MM_ano_complete.txt",sep="\t",dec=".",header=TRUE) #1734
mouse_sig <- mouse_sig[mouse_sig$signature!=0,] #2356 x 71

list_ref_human_sig <- as.character(unique(mouse_sig$Gene.stable.ID.1[!is.na(mouse_sig$Gene.stable.ID.1)])) #1,698


metadata_common<-metadata[rownames(metadata)%in%list_ref_human_sig,]


meta_sig<-metadata_common[metadata_common$meta_adj<0.05,]

common_genes<-c(rownames(meta_sig),iruna_filtered$ensembl_gene_id)

## Differential expression results -------------------------------------------------

gse47552 <- read.table("C:/Users/transbio/Desktop/MM_analysis/human_rm_IG/Results/20200527_dea_anno_gse47552.txt", sep="\t",dec=".",header=TRUE)
gse6477 <- read.table("C:/Users/transbio/Desktop/MM_analysis/human_rm_IG/Results/20200525_dea_anno_gse6477.txt",row.names=1,sep="\t",dec=".",header=TRUE)
iruna <- read.table("C:/Users/transbio/Desktop/MM_human/dea_MM_control_anno.txt", sep="\t",dec=".", header=TRUE,row.names=1)


# Generate variables to plot

##>>>>> gse47552
gse47552_filt$to_plot <- gse47552_filt$MM_NPC.sig_genes
gse47552_filt$to_plot[gse47552_filt$MM_NPC.logFC>0 & gse47552_filt$MM_NPC.sig_genes==0] <- 0.5
gse47552_filt$to_plot[gse47552_filt$MM_NPC.logFC<0 & gse47552_filt$MM_NPC.sig_genes==0] <- -0.5

##>>>>> gse6477
gse6477_filt$to_plot <- gse6477_filt$MM_NPC.sig_genes
gse6477_filt$to_plot[gse6477_filt$MM_NPC.logFC>0 & gse6477_filt$MM_NPC.sig_genes==0] <- 0.5
gse6477_filt$to_plot[gse6477_filt$MM_NPC.logFC<0 & gse6477_filt$MM_NPC.sig_genes==0] <- -0.5

##>>>>>Iruña to plot
iruna_filt$to_plot <- rep(0,nrow(iruna_filt))
iruna_filt$to_plot[iruna_filt$MM_C.logFC>0.58 & iruna_filt$MM_C.adj.P.Val<0.05] <- 1
iruna_filt$to_plot[iruna_filt$MM_C.logFC< -0.58 & iruna_filt$MM_C.adj.P.Val<0.05] <- -1
iruna_filt$to_plot[iruna_filt$MM_C.logFC>0 & iruna_filt$to_plot==0] <- 0.5
iruna_filt$to_plot[iruna_filt$MM_C.logFC<0 & iruna_filt$to_plot==0] <- -0.5



# Filter data by mouse signature
gse47552_filtered <- gse47552_filt[gse47552_filt$ensembl%in%list_ref_human_sig,c("probe","symbol","ensembl","to_plot")] #1,375 x 4
gse6477_filtered <- gse6477_filt[gse6477_filt$ensembl%in%list_ref_human_sig,c("probe","symbol","ensembl","to_plot")] #1,596 x 4
iruna_filtered <- iruna_filt[iruna_filt$ensembl_gene_id%in%list_ref_human_sig,c("ensembl_gene_id","to_plot")] #1,436 x 24
dataset4_filtered<-dataset4_filt[dataset4_filt$ensembl_gene_id%in%list_ref_human_sig,c("ensembl_gene_id","to_plot")]

iruna_sig<-iruna[abs(iruna$to_plot)==1,]
gse47552_sig<-gse47552[abs(gse47552$to_plot)==1,]
gse6477_sig<-gse6477[abs(gse6477$to_plot)==1,]
dataset4_sig<-dataset_4[abs(dataset_4$to_plot)==1,]

m<-list(metanalysis=rownames(metadata),iruna=iruna_sig$ensembl_gene_id,gse47552=gse47552_sig$ensembl,gse6477=gse6477_sig$ensembl,dataset4=dataset4_sig$ensembl_gene_id)
venn.diagram(m,filename="venn_complete_with_metadata.png",fill=c("skyblue", "pink1", "mediumorchid", "orange","red"))

# plot heatmap with sig genes of metanalysis

iruna_meta<-iruna_filtered[iruna_filtered$ensembl_gene_id%in%rownames(meta_sig),]
gse47552_meta<-gse47552_filtered[gse47552_filtered$ensembl%in%rownames(meta_sig),]
gse6477_meta<-gse6477_filtered[gse6477_filtered$ensembl%in%rownames(meta_sig),]


data_merged <- merge(gse47552_meta,gse6477_meta,by.x="ensembl",by.y="ensembl", sort=FALSE, all.x=TRUE, all.y=TRUE) #2,064 x 7
data_merged2 <- merge(data_merged,iruna_meta, by.x="ensembl",by.y="ensembl_gene_id", sort=FALSE, all.x=TRUE, all.y=TRUE) #2,456 x 9


sel_significant <- vector()
for(i in 1:nrow(data_merged2)){
  sel_significant[i] <- (abs(data_merged2[i,"to_plot.x"])==1 | abs(data_merged2[i,"to_plot.y"])==1 | abs(data_merged2[i,"to_plot"])==1)
}

data_merged_significant <- data_merged2[which(sel_significant==TRUE),] #439 x 9
data_merged_significant$ensembl <- as.character(data_merged_significant$ensembl)
data_meta<-data_merged2[data_merged2$ensembl%in%rownames(meta_sig),]



##Plotting

library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

target_ensembl_ids <- unique(data_merged_significant$ensembl) #301 unique IDs

aa <- table(data_merged_significant$ensembl)
repited_probes <- names(aa)[aa>1]

to_remove <- list()
for(i in 1:length(repited_probes)){
  print(data_meta[which(as.character(data_meta$ensembl)%in%repited_probes[i]),])
  print(which(as.character(data_meta$ensembl)%in%repited_probes[i]))
  to_remove[[i]] <- as.integer(strsplit(readline(prompt="Positions to remove: "), " ")[[1]])
}

to_remove2 <- unlist(to_remove)

#36  37  38  83  14  15  16  17  18  40 129  86  87  88   8   9 132 133  3  74  56  91  22  25  94  76  77 122 123 124 125 113 114  32  20  46  67 108   5  48  61  62  52 137  69  70  99 118 119 121 110

#144    100    101    102    300    122    203    307    308    309
#45     47     48     49     50     77    305 219220    221    171
#256     79     80     81     36     37     16    184    232    141
#225     58     67     97    292    280    281    229    175      7
# 94    187    188     11    288    289    290    291    249      1
# 90     92     52    303    130    131    180    181    182    405
#   274    115    120    403    117    152    166     23     25     26
#  27     27     28     29     30     31      4     21    264     72
# 416    294    295    296    297     63    117    325    326    137
# 206     69     70     32     33     18     13     14    238    261
# 262    124    126    128    158    159    149    150    113    191
# 192    194    240    208    209    210    111    332    243    216
#217    218    234    313    315    271    333    334    335    358
# 266    267    268    223    409    277    278    279

kentzeko<-c("ENSG00000211445:GPX3","ENSG00000075218:GTSE1","ENSG00000153094:BCL2L11")

data_merged_significant2 <- data_merged_significant[-to_remove2,] #118 
data_merged_significant_def <- data_merged_significant2[,] #301
data_merged_significant2 <- data_meta[-to_remove2,]

#>>> Differential Expression Analysis "binary" heatmap
data_to_plot <- data_merged2[,c("to_plot.x","to_plot.y","to_plot")]
colnames(data_to_plot) <- c("gse47552","gse6477","Iruña")
rownames(data_to_plot) <- paste(data_merged2$ensembl,data_merged2$symbol.x,sep=":")
genes_our_data<-paste(data_merged_significant$ensembl,data_merged_significant$symbol.x,sep=":")
#heatmap annotation
my_palette2 <- brewer.pal(11, "Paired")

ha1 <- HeatmapAnnotation(
  Dataset = colnames(data_to_plot),
  col = list(Dataset = c("gse47552"=my_palette2[1],
                         "gse6477"=my_palette2[2],"Iruña"=my_palette2[3])),
  show_annotation_name = TRUE
)


direction_info <- vector()
for(i in 1:nrow(data_to_plot)){
  dinfo <- mouse_sig$signature[mouse_sig$Gene.stable.ID.1%in%data_merged2$ensembl[i]]
  if(length(dinfo)>1){
    cat(dinfo,"\n")
    dinfo <- readline(prompt="Gene direction: ")
  }
  direction_info[i] <- dinfo
}
names(direction_info) <- data_merged2$ensembl

row_ha <- rowAnnotation(Mouse = as.factor(direction_info), col=list(Mouse=c("-1"="blue","0"="gray","1"="red")))


data_to_plot[is.na(data_to_plot)] <- 0

#plot
pdf("meta_all.pdf")
draw(Heatmap(data_to_plot, name = "DEGs MM vs Control",
             col = colorRamp2(c(-1, 0, 1), c("darkblue","gray","red")),
             top_annotation = ha1, right_annotation = row_ha,
             show_column_names = FALSE, show_row_names = TRUE, row_names_gp = gpar(fontsize = 1), 
             cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE))
dev.off()
