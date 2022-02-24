#http://www.bioconductor.org/packages/release/data/experiment/vignettes/gskb/inst/doc/gskb.pdf
#source("https://bioconductor.org/biocLite.R")
#biocLite("gskb")
require(gskb)
library(biomaRt)
library(gskb)

## SELECTING GO
data(mm_GO)
list_members <- lapply(mm_GO,length)
list_members_S <- list_members[list_members > 20 &
                                 list_members < 200 ]
list_members_S<-names(list_members_S[
  grep("_BP_",names(list_members_S))])
mm_GO_S<-mm_GO[names(mm_GO) %in% list_members_S]
## SELECTING GO
data(mm_metabolic)
list_members <- lapply(mm_metabolic,length)
list_members_S <- list_members[list_members > 20 &
                                 list_members < 200 ]
#list_members_S<-names(list_members_S[
grep("_BP_",names(list_members_S))
mm_metabolic_S<-mm_metabolic[names(mm_metabolic) %in% names(list_members_S)]


## SELECTING pathway
data(mm_pathway)
rm(list_members,list_members_S)
list_members <- lapply(mm_pathway,length)
list_members_S <- list_members[list_members > 20 &
                                 list_members < 200 ]
mm_pathway_S<-mm_pathway[names(mm_pathway) %in% names(list_members_S)]

## SELECTING TF
data(mm_TF)
list_members <- lapply(mm_TF,length)
list_members_S <- list_members[list_members > 20 &
                                 list_members < 200 ]
mm_TF_S<-mm_TF[names(mm_TF) %in% names(list_members_S)]


list_mm<-list(GO=mm_GO_S,
              pathways=mm_pathway_S,
              metabolic=mm_metabolic_S,
              TF=mm_TF_S)

#>>> LOAD OUR DATA

set.seed(1234567)
setwd("C:/Users/transbio/Desktop/MM_human")

dea_res<-read.table("dea_MM_control.txt", header = T, check.names = F, stringsAsFactors = F, row.names = 1)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

anno2 = getBM(
  values = rownames(dea_res),
  filters = c("ensembl_gene_id"),
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  mart = human
)


dea_anno<-merge(dea_res,anno2,by.x='row.names', by.y='ensembl_gene_id')
dea_sig<-dea_anno[dea_anno$`MM_C:sig_genes`!=0,]


anno_metadata = getBM(
  values = rownames(metadata),
  filters = c("ensembl_gene_id"),
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  mart = human
)

metadata<-read.table("metadata_updated.txt", header = T, check.names = F, stringsAsFactors = F, row.names = 1)
metadata_ano<-merge(metadata,anno_metadata,by.x='row.names', by.y='ensembl_gene_id')


metadata_sig<-subset(metadata_ano, metadata_ano[,6] < 0.05)

######## MM_C contrast
version <- "_v1_"
sig_names <- my_dataset[,"hgnc_symbol"]
list_setsgo <- list(sig_names)
gostart<-0

genessymb_ALL <- dea_anno$hgnc_symbol



for(i in 1:length(list_mm))
{#i<-1
  genesetgo<-list_mm[[i]]
  for(j in 1:length(list_setsgo))
  {#j<-1
    
    positive<-matrix(NA,length(genesetgo),10)
    list_setsgo1 <- unique(toupper(list_setsgo[[j]]))
    genessymb<-list_setsgo1[list_setsgo1 %in% genessymb_ALL]
    altgenessymb<-genessymb_ALL[!(genessymb_ALL %in% genessymb)]
    
    for(h in 1:length(genesetgo))
    {#h<-1
      genesetgo_use<-genesetgo[[h]][genesetgo[[h]] %in% genessymb_ALL]
      
      list_and_pathway <- sum(genessymb %in% genesetgo_use)
      list_and_NOTpathway <- sum(!(genessymb %in% genesetgo_use))
      NOTlist_and_pathway <- sum(altgenessymb %in% genesetgo_use)
      NOTlist_and_NOTpathway <- sum(!(altgenessymb %in% genesetgo_use))     
      
      mycontingencytable = matrix(c(list_and_pathway,list_and_NOTpathway,
                                    NOTlist_and_pathway,NOTlist_and_NOTpathway),2,2, byrow = TRUE)
      positive[h,1] <- names(genesetgo)[h]
      positive[h,2] <- fisher.test(x = mycontingencytable, alternative = "greater")$p.value
      positive[h,3] <- phyper(q = mycontingencytable[1,1]-1, m = rowSums(mycontingencytable)[1], 
                              n = rowSums(mycontingencytable)[2], k = colSums(mycontingencytable)[1], 
                              lower.tail = FALSE)
      positive[h,4] <-  list_and_pathway
      positive[h,5] <-  NOTlist_and_pathway
      
    }
    positive[,6]<-p.adjust(as.numeric(positive[,2]),method="BH")
    positive[,7]<-p.adjust(as.numeric(positive[,3]),method="BH")
    
    colnames(positive) <- c("NAME","pval_Fisher","pval_HyperG",
                            "elements","pathway",
                            "p.adjust_Fisher","p.adjust_HyperG",
                            "list","Gene-Set Type","Cluster_involved")
    positive<-positive[order(as.numeric(positive[,2])),]
    #  SOL[,8]<-list_setsgo[j]
    positive[,9]<-names(list_mm)[i]
    positive[,10] <-version
    # if(gostart==0)
    {
      # SOLALL<-IRF8_positive
      # gostart<-1
      # }else{
      #SOLALL<-rbind(SOLALL,IRF8_positive)
    }
    
    write.table(positive, file= "GSEA/GSA_my_dataset_only_genes.txt", sep=",")
    
  }
}

