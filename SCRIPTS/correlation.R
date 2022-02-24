####################
#CORRELATION OF FCs#
####################

#-Load human and mouse data and select the logFC column

MM_human<-read.table("C:/Users/transbio/Desktop/MM_analysis/Human_results/20200506_dea_MM_CONTROL.txt",sep="\t",row.names=1)
MM_mouse<-read.table("C:/Users/transbio/Desktop/MM_analysis/mouse_ensembl/Data/dea_MM_CONTROL.txt",sep="\t",row.names=1)

human_FC<-MM_human[,grep("logFC",colnames(MM_human))]
colnames(human_FC) <- gsub(".logFC","",colnames(human_FC))

mouse_FC<-MM_mouse[,grep("logFC",colnames(MM_mouse))]
colnames(mouse_FC)<-gsub(".logFC","",colnames(mouse_FC))

#--Load the names and get only the ones that have in both samples

names<-read.table("C:/Users/transbio/Desktop/MM_analysis/Huma_data/mouse_sig_in_human.txt",sep="\t")
names_short_positive<-read.table()
human_FC2<-human_FC[names$human_name,grep("MM_C",colnames(MM_mouse))]
mouse_FC<-mouse_FC[names$mouse_names,]

#--Calculate the correlation

for(i in 1:ncol(mouse_FC)){
  correlation[i]<-cor(human_FC,mouse_FC[,i],method="spearman")
}
MIC<-cor(human_FC,mouse_FC[1,1],method="spearman")
