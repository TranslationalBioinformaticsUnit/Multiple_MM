############################################
##  Mouse signature in human              ##
############################################




setwd("C:/Users/transbio/Desktop/MM_analysis/human_rm_IG/Results")

# Differential expression results -------------------------------------------------

gse47552 <- read.table("C:/Users/transbio/Desktop/MM_analysis/human_rm_IG/Results/20200527_dea_anno_gse47552.txt", sep="\t",dec=".",header=TRUE)
gse6477 <- read.table("C:/Users/transbio/Desktop/MM_analysis/human_rm_IG/Results/20200525_dea_anno_gse6477.txt",row.names=1,sep="\t",dec=".",header=TRUE)
iruna <- read.table("C:/Users/transbio/Desktop/MM_human/dea_MM_control_anno.txt", sep="\t",dec=".", header=TRUE,row.names=1)
dataset_4<-read.table("C:/Users/transbio/Desktop/MM_new_dataset/dataset_4.txt",sep="\t",dec=".", header=TRUE,row.names=1)

# Mouse signature ------------------------------------------------------------------

mouse_sig <- read.table("C:/Users/transbio/Desktop/MM_analysis/mouse_ensembl/Data/results_MM_ano_complete.txt",sep="\t",dec=".",header=TRUE) #1734
mouse_sig <- mouse_sig[mouse_sig$signature!=0,] #2356 x 71

list_ref_human_sig <- as.character(unique(mouse_sig$Gene.stable.ID.1[!is.na(mouse_sig$Gene.stable.ID.1)])) #1,698


# Generate variables to plot

##>>>>> gse47552
gse47552$to_plot <- gse47552$MM_NPC.sig_genes
gse47552$to_plot[gse47552$MM_NPC.logFC>0 & gse47552$MM_NPC.sig_genes==0] <- 0.5
gse47552$to_plot[gse47552$MM_NPC.logFC<0 & gse47552$MM_NPC.sig_genes==0] <- -0.5

##>>>>> gse6477
gse6477$to_plot <- gse6477$MM_NPC.sig_genes
gse6477$to_plot[gse6477$MM_NPC.logFC>0 & gse6477$MM_NPC.sig_genes==0] <- 0.5
gse6477$to_plot[gse6477$MM_NPC.logFC<0 & gse6477$MM_NPC.sig_genes==0] <- -0.5

##>>>>>Iruña to plot
iruna$to_plot <- rep(0,nrow(iruna))
iruna$to_plot[iruna$MM_C.logFC>0.58 & iruna$MM_C.adj.P.Val<0.05] <- 1
iruna$to_plot[iruna$MM_C.logFC< -0.58 & iruna$MM_C.adj.P.Val<0.05] <- -1
iruna$to_plot[iruna$MM_C.logFC>0 & iruna$to_plot==0] <- 0.5
iruna$to_plot[iruna$MM_C.logFC<0 & iruna$to_plot==0] <- -0.5

##>>>>>Dataset4 to plot
dataset_4$to_plot <- rep(0,nrow(dataset_4))
dataset_4$to_plot[dataset_4$MM_C.logFC>0.58 & dataset_4$MM_C.adj.P.Val<0.05] <- 1
dataset_4$to_plot[dataset_4$MM_C.logFC< -0.58 & dataset_4$MM_C.adj.P.Val<0.05] <- -1
dataset_4$to_plot[dataset_4$MM_C.logFC>0 & dataset_4$to_plot==0] <- 0.5
dataset_4$to_plot[dataset_4$MM_C.logFC<0 & dataset_4$to_plot==0] <- -0.5

# Filter data by mouse signature
gse47552_filtered <- gse47552[gse47552$ensembl%in%list_ref_human_sig,c("probe","symbol","ensembl","to_plot")] #1,375 x 4
gse6477_filtered <- gse6477[gse6477$ensembl%in%list_ref_human_sig,c("probe","symbol","ensembl","to_plot")] #1,596 x 4
iruna_filtered <- iruna[iruna$ensembl_gene_id%in%list_ref_human_sig,c("ensembl_gene_id","to_plot")] #1,436 x 24
dataset4_filtered <- dataset_4[dataset_4$ensembl_gene_id%in%list_ref_human_sig,c("ensembl_gene_id","to_plot")] #1,436 x 24

data_merged <- merge(gse47552_filtered,gse6477_filtered,by.x="ensembl",by.y="ensembl", sort=FALSE, all.x=TRUE, all.y=TRUE) #2,064 x 7
data_merged2 <- merge(data_merged,iruna_filtered, by.x="ensembl",by.y="ensembl_gene_id", sort=FALSE, all.x=TRUE, all.y=TRUE) #2,456 x 9


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

target_ensembl_ids <- unique(data_meta$ensembl) #301 unique IDs

aa <- table(data_meta$ensembl)
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
data_to_plot <- data_meta[,c("to_plot.x","to_plot.y","to_plot")]
colnames(data_to_plot) <- c("gse47552","gse6477","Iruña")
rownames(data_to_plot) <- paste(data_merged_significant2$ensembl,data_merged_significant2$symbol.x,sep=":")
genes_our_data<-paste(data_merged_significant2$ensembl,data_merged_significant2$symbol.x,sep=":")
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
  dinfo <- mouse_sig$signature[mouse_sig$Gene.stable.ID.1%in%data_merged_significant2$ensembl[i]]
  if(length(dinfo)>1){
    cat(dinfo,"\n")
    dinfo <- readline(prompt="Gene direction: ")
  }
  direction_info[i] <- dinfo
}
names(direction_info) <- data_merged_significant2$ensembl

row_ha <- rowAnnotation(Mouse = as.factor(direction_info), col=list(Mouse=c("-1"="blue","0"="gray","1"="red")))


data_to_plot[is.na(data_to_plot)] <- 0

#plot
pdf("20200601_heatmap_summary_mouse_signature_in_human_final.pdf")
draw(Heatmap(data_to_plot, name = "DEGs MM vs Control",
        col = colorRamp2(c(-1, 0, 1), c("darkblue","gray","red")),
        top_annotation = ha1, #right_annotation = row_ha,
        show_column_names = FALSE, show_row_names = TRUE, row_names_gp = gpar(fontsize = 4), 
        cluster_rows = TRUE, cluster_columns = TRUE, row_gap = unit(3, "mm"), use_raster=TRUE))
dev.off()


##solo teniendco en cuenta pamplona

iruna<-read.table("C:/Users/transbio/Desktop/MM_human/dea_MM_control_anno.txt", sep="\t",dec=".", header=TRUE,row.names=1)
common_ensembl <- which(iruna$ensembl_gene_id%in%list_ref_human_sig)
ref_genes <- unique(anno_human$ensembl_gene_id[common_ensembl])
ref_genes_intersection<-iruna[which(ref_genes%in%list_sig_genes_MM),]
