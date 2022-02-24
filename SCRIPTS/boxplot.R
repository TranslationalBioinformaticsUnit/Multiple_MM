#################################################
###Boxplot up and down genes of mouse in human###
#################################################

#--Set workdirectory
setwd("C:/Users/transbio/Desktop/MM_analysis/mouse_ensembl/Data")

#--Load the list of genes of the mouse data

list_large<-read.table("list_sig_genes_MM_CONTROL.txt",sep="\t")
head(list_large)
results_MM<-read.table("dea_MM_CONTROL.txt",sep="\t",row.names=1)
short_positive<-read.table("list_sig_genes_MM_short_positive.txt")
head(short_positive)
short_negative<-read.table("list_sig_MM_short_negative.txt")

#--Select those genes in human data

names<-read.table("C:/Users/transbio/Desktop/MM_analysis/Huma_data/mouse_sig_in_human.txt",sep="\t")
names<-unique(names)
names$mouse_names<-as.factor(names$mouse_names)
names$human_name<-as.factor(names$human_name)

short_positive_human<-names[names$mouse_names%in%short_positive$V1,]
short_negative_human<-names[names$mouse_names%in%short_negative$V1,]
short_all_human<-names[names$mouse_names%in%human_short_all,]
#--Load human data

MM_human<-read.table("C:/Users/transbio/Desktop/MM_analysis/Human_rm_IG/Results/dea_MM_CONTROL.txt")
head(MM_human)

human_short_neg<-MM_human[short_negative_human$human_name,grep("logFC",colnames(MM_human))]
colnames(human_short_neg) <- gsub(".logFC","",colnames(human_short_neg))

human_short_pos<-MM_human[short_positive_human$human_name,grep("logFC",colnames(MM_human))]
colnames(human_short_pos) <- gsub(".logFC","",colnames(human_short_pos))

human_short_all<-MM_human[short_all_human$human_name,grep("logFC",colnames(MM_human))]
colnames(human_short_all) <- gsub(".logFC","",colnames(human_short_all))

human_all<-MM_human[,grep("logFC",colnames(MM_human))]
colnames(human_all) <- gsub(".logFC","",colnames(human_all))
#--plot the boxplot for each case

p_box <- ggplot(human_short_neg, aes(x="", y=MM_C),fill="#E69F00")+
  
  geom_boxplot() +
  ggtitle(paste(" logFC normalized gene expression", sep="")) +
  #scale_fill_manual(values=my_palette1)+
  #theme_ipsum()+
  theme(axis.text.x = element_text(size=10, angle = 45))

ggsave("short_negative.png")

#intersection
intersection<-list_sig_genes_MM[list_sig_genes_MM%in%names$human_name]
intersection<-MM_human[intersection,]
