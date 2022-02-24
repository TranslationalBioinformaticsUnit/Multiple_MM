#####################################################
###### GENE ONTOLOGY ANALYSIS: CLUSTER PROFILER #####
#####################################################
#Source: https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html#go-analysis

library('clusterProfiler')
library('org.Hs.eg.db') # Human
library('org.Mm.eg.db') # Mouse
library("GSEABase")
library('DOSE')
library('enrichplot')
library(clusterProfiler)
library('ReactomePA')
library(biomaRt)
library(pheatmap)

##### 0. Reading contrasts matrix ----------
setwd("C:/Users/transbio/Desktop/MM_analysis/mouse_results/MM")

limma.res <- read.table('20200506_dea_MM_CONTROL.txt', header = T, 
                        check.names = F, stringsAsFactors = F, row.names = 1)
colnames(limma.res)

# Filter by contrast
Contrast <- "MIC_C"
limma.res.filt <- limma.res[, grep(Contrast, colnames(limma.res))]
colnames(limma.res.filt)
colnames(limma.res.filt) <- gsub(paste0(".", Contrast), '', colnames(limma.res.filt))
colnames(limma.res.filt)

# If there is only one contrast
colnames(limma.res) <- c('A', 'Coef', 't', 'p.value', 'p.value.adj', 'F', 'F.p.value', 'Res')
limma.res.filt <- limma.res
colnames(limma.res.filt)

##### 1. Biological Id TranslatoR ----------
# ClusterProfiler provides bitr and bitr_kegg for converting ID types.

### Converting IDs using bitr
genes_ID = bitr(rownames(limma.res), fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")
head(genes_ID)

### Converting biological IDs using KEGG API (entrez ID is the same as kegg ID)
#genes_kegg_ID <- bitr_kegg(genes_ID$ENTREZID, fromType='ncbi-geneid', toType='kegg', organism='mmu')
#head(genes_kegg_ID)

##### 2. GO/KEGG Enrichment classification ----------------
## For Enrichment classification you need to define a set of DEGs based on a FC and p.value threshold

DEG <- subset(limma.res.filt, abs(Coef) >1 & p.value.adj < 0.05)
DEG <- subset(limma.res.filt, p.value.adj < 0.05)
DEG <- subset(limma.res.filt, p.value.adj_MM_C  < 0.05 & p.value.adj_MGUS_C < 0.05)


DEG_FC <- DEG[,'Coef']
names(DEG_FC) <- as.character(rownames(DEG))

# Run GO Enrichment
ggo <- groupGO(gene = rownames(DEG),
               OrgDb = org.Mm.eg.db,
               ont = "BP", # Options: BP, CC, MF 
               level = 2,
               keyType = 'SYMBOL') # if True, entrezIDs are converted to gene symbols

View(as.data.frame(ggo))
barplot(ggo, drop=T, showCategory=10)
cnetplot(ggo, foldChange = DEG_FC)
heatplot(ggo, foldChange = DEG_FC)

# Compare GO enrichment between Up/Downreg genes

mydf <- data.frame('Gene_IDs'=rownames(DEG), 'FC' = DEG$Coef)
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"
mydf <- merge(mydf, genes_ID[,c(1:2)], by.x='Gene_IDs', by.y='SYMBOL')

head(mydf)

# Multiple Groups

mydf <- data.frame('Gene_IDs'=rownames(DEG), 
                   'FC_MM_C' = DEG$Coef_MM_C, 'FC_MGUS_C'=DEG$Coef_MGUS_C)
mydf$group_MMvsC <- "up"
mydf$group_MMvsC[mydf$FC_MM_C < 0] <- "down"
#mydf$group_MMvsC[abs(mydf$FC_MM_C) < 1] <- "ns"
mydf$group_MGUSvsC <- "up"
mydf$group_MGUSvsC[mydf$FC_MGUS_C < 0] <- "down"
#mydf$group_MGUSvsC[abs(mydf$FC_MGUS_C) < 1] <- "ns"
mydf$group_MMvsMGUS <- "up"
mydf$group_MMvsMGUS[mydf$FC_MM_MGUS < 0] <- "down"
#mydf$group_MMvsMGUS[abs(mydf$FC_MM_MGUS) < 1] <- "ns"
head(mydf)

mydf <- merge(mydf, genes_ID[,c(1:2)], by.x='Gene_IDs', by.y='SYMBOL')

head(mydf)

# GO enrichment
formula_res_GO <- compareCluster(Gene_IDs~group, data=mydf, 
                                 fun="enrichGO", 
                                 OrgDb = org.Mm.eg.db, ont = "BP", keyType = 'SYMBOL')

formula_res_GO <- compareCluster(Gene_IDs~group_MMvsC+group_MGUSvsC, data=mydf, 
                                 fun="enrichGO", 
                                 OrgDb = org.Mm.eg.db, ont = "BP", keyType = 'SYMBOL')

# KEGG enrichment
formula_res_KEGG <- compareCluster(ENTREZID~group, data=mydf, 
                                   fun="enrichKEGG", 
                                   organism = 'mmu', pvalueCutoff=0.05)

formula_res_KEGG <- compareCluster(ENTREZID~group_MMvsC+group_MGUSvsC, data=mydf, 
                                   fun="enrichKEGG", 
                                   organism = 'mmu', pvalueCutoff=0.01)

# Disease Ontology Enrichment (few enriched categories!)
formula_res_DO <- compareCluster(ENTREZID~group, data=mydf, 
                                 fun="enrichDO")

# Pathwat enrichment
#require(ReactomePA)
formula_res_PATH <- compareCluster(ENTREZID~group, data=mydf, 
                                   fun="enrichPathway", organism='mouse', readable=T)

formula_res_PATH <- compareCluster(ENTREZID~group_MMvsC+group_MGUSvsC, data=mydf, 
                                   fun="enrichPathway", organism='mouse', readable=T)

View(as.data.frame(formula_res_PATH))

# Plot Results
pdf('RESULTS/B2IC_Model/GO_Up_Down_padj0.05.pdf', width = 9)
clusterProfiler::dotplot(formula_res_GO, x=~group, showCategory=6, color = 'p.adjust', by='geneRatio')
dev.off()

clusterProfiler::dotplot(formula_res_PATH, x=~Cluster, showCategory=5, color = 'p.adjust', by='geneRatio')

clusterProfiler::dotplot(formula_res_PATH, x=~group_MMvsC, showCategory=5, color = 'p.adjust', by='geneRatio') + ggplot2::facet_grid(~group_MGUSvsC)
clusterProfiler::dotplot(formula_res, x=~group_MMvsC, showCategory=5, color = 'p.adjust', by='geneRatio') 

#More visualizations http://bioconductor.org/packages/release/bioc/vignettes/enrichplot/inst/doc/enrichplot.html

##### 3. GO Gene Set Enrichment analysis --------------
## For Gene Set Enrichment analysis (GSEA) you need to provide a ranked list of genes (ALL) based on a FC or p.value

### Create a decreasing sorted vector according to FC values
ranked_genes <- limma.res.filt[ ,"Maf:logFC"]
names(ranked_genes) <- as.character(rownames(limma.res.filt))
ranked_genes <- sort(ranked_genes, decreasing = T)
head(ranked_genes)
tail(ranked_genes)

### Run GO analysis
ego <- gseGO(geneList = ranked_genes,
             OrgDb = org.Mm.eg.db,
             ont = 'BP',    # Options: BP, CC, MF 
             nPerm = 1000,
             minGSSize = 20,
             maxGSSize = 500,
             pvalueCutoff = 1,
             verbose = F,
             keyType = 'SYMBOL')
View(as.data.frame(ego))

barplot(ego, drop=T, showCategory=10)
cnetplot(ego, foldChange = DEG_FC)
heatplot(ego, foldChange = DEG_FC)

##### 4. GMT categories enrichment --------------------
## For Enrichment classification you need to define a set of DEGs based on a FC and p.value threshold

DEG <- subset(limma.res.filt, abs(Coef) >1 & p.value.adj < 0.05)
DEG <- subset(limma.res.filt, p.value.adj < 0.05)
DEG <- subset(limma.res.filt, p.value.adj < 0.1)

DEG_FC <- DEG[,'Coef']
names(DEG_FC) <- as.character(rownames(DEG))

##### Convert mouse IDs to homolog human IDs with BIOMART #####
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

genesV2 = getLDS(attributes = c("external_gene_name"), filters = "external_gene_name", 
                 values = rownames(limma.res.filt) , mart = mouse, 
                 attributesL = c("external_gene_name"), martL = human, uniqueRows=T)

HS_Genes <- merge(DEG, genesV2, by.x='row.names', by.y='Gene.name')

### Read .gtm files based on ENTREZ_IDs
# http://software.broadinstitute.org/gsea/downloads.jsp
### 1. All gene sets
gmtfile <- read.gmt('~/Documents/COMP_ANALYSIS/GSEA/msigdb.v6.2.symbols.gmt.txt')
### 2. Hallmark gene sets (H)
gmtfile <- read.gmt('~/Documents/COMP_ANALYSIS/GSEA/h.all.v6.2.symbols.gmt.txt')
### 3. Chemical and Genetic perturbations, curated gene sets (C2)
gmtfile <- read.gmt('~/Documents/COMP_ANALYSIS/GSEA/c2.cgp.v6.2.symbols.gmt.txt')
### 4. Canonical pathways, curated gene sets (C2)
gmtfile <- read.gmt('~/Documents/COMP_ANALYSIS/GSEA/c2.cp.v6.2.symbols.gmt.txt')
### 5. Transcription factor targets, motif gene sets (C3)
gmtfile <- read.gmt('~/Documents/COMP_ANALYSIS/GSEA/c3.tft.v6.2.symbols.gmt.txt')
### 6. KEGG Gene Sets (C2)
gmtfile <- read.gmt('~/Documents/COMP_ANALYSIS/GSEA/c2.cp.kegg.v6.2.symbols.gmt.txt')
### 7. Gene ontology (GO) biological processes (C5)
gmtfile <- read.gmt('~/Documents/COMP_ANALYSIS/GSEA/c5.bp.v6.2.symbols.gmt.txt')
### 8. BioCarta gene sets (C2)
gmtfile <- read.gmt('~/Documents/COMP_ANALYSIS/GSEA/c2.cp.biocarta.v6.2.symbols.gmt.txt')
### 9. Oncogenic gene signatures (C6)
gmtfile <- read.gmt('~/Documents/COMP_ANALYSIS/GSEA/c6.all.v6.2.symbols.gmt.txt')


gmt_enrich <- enricher(HS_Genes$Gene.name.1, TERM2GENE = gmtfile, pvalueCutoff = 0.9, minGSSize = 5)
View(as.data.frame(gmt_enrich))
dim(as.data.frame(gmt_enrich))

# Define FC vector for colors in plots
DEG_FC <- HS_Genes[,'Coef']
names(DEG_FC) <- as.character(HS_Genes$Hs.Symbol)

barplot(gmt_enrich, drop=T, showCategory=10)
cnetplot(gmt_enrich, foldChange = DEG_FC)
heatplot(gmt_enrich, foldChange = DEG_FC)

##### 4. GMT Gene Set enrichment analysis --------------
## For Gene Set Enrichment analysis (GSEA) you need to provide a ranked list of genes (ALL) based on a FC or p.value

# Filter by contrast
Contrast <- "MIC_MM_C"
limma.res.filt <- limma.res[, grep(Contrast, colnames(limma.res))]
colnames(limma.res.filt)
colnames(limma.res.filt) <- gsub(paste0(".", Contrast), '', colnames(limma.res.filt))
colnames(limma.res.filt)

##### Convert mouse IDs to homolog human IDs with BIOMART #####
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

genesV2 = getLDS(attributes = c("external_gene_name"), filters = "external_gene_name", 
                 values = rownames(limma.res.filt) , mart = mouse, 
                 attributesL = c("external_gene_name"), martL = human, uniqueRows=T)

HS_Genes <- merge(limma.res.filt, genesV2, by.x='row.names', by.y='Gene.name')

### Create a decreasing sorted vector according to FC values
ranked_genes <- limma.res.filt[ ,'MM_C.logFC']
names(ranked_genes) <- rownames(limma.res.filt)
ranked_genes <- sort(ranked_genes, decreasing = T)
head(ranked_genes)
tail(ranked_genes)

### Read .gtm files based on ENTREZ_IDs
# http://software.broadinstitute.org/gsea/downloads.jsp
### 1. All gene sets
gmtfile <- read.gmt('~/Documents/COMP_ANALYSIS/GSEA/msigdb.v6.2.symbols.gmt.txt')
### 2. Hallmark gene sets (H)
gmtfile <- read.gmt('~/Documents/COMP_ANALYSIS/GSEA/h.all.v6.2.symbols.gmt.txt')
### 3. Chemical and Genetic perturbations, curated gene sets (C2)
gmtfile <- read.gmt('~/Documents/COMP_ANALYSIS/GSEA/c2.cgp.v6.2.symbols.gmt.txt')
### 4. Canonical pathways, curated gene sets (C2)
gmtfile <- read.gmt('~/Documents/COMP_ANALYSIS/GSEA/c2.cp.v6.2.symbols.gmt.txt')
### 5. Transcription factor targets, motif gene sets (C3)
gmtfile <- read.gmt('~/Documents/COMP_ANALYSIS/GSEA/c3.tft.v6.2.symbols.gmt.txt')
### 6. KEGG Gene Sets (C2)
gmtfile <- read.gmt('~/Documents/COMP_ANALYSIS/GSEA/c2.cp.kegg.v6.2.symbols.gmt.txt')
### 7. Gene ontology (GO) biological processes (C5)
gmtfile <- read.gmt('~/Documents/COMP_ANALYSIS/GSEA/c5.bp.v6.2.symbols.gmt.txt')
### 8. BioCarta gene sets (C2)
gmtfile <- read.gmt('~/Documents/COMP_ANALYSIS/GSEA/c2.cp.biocarta.v6.2.symbols.gmt.txt')
### 9. Oncogenic gene signatures (C6)
gmtfile <- read.gmt('~/Documents/COMP_ANALYSIS/GSEA/c6.all.v6.2.symbols.gmt.txt')
### 10. My dataset
gmtfile <- read.gmt('C:/Users/transbio/Desktop/MM_analysis/GSEA_files/MM_mouse_geneset.gmt.txt')

### Calculate GSEA
egmt <- GSEA(ranked_genes, TERM2GENE=gmtfile, verbose=T, minGSSize = 10, pvalueCutoff = 1)
dim(as.data.frame(egmt))
View(as.data.frame(egmt))

pB2IC <- as.data.frame(egmt)

write.csv(as.data.frame(egmt), 'RESULTS/05.pB2IC/GSEA_Hallmarks_pB2IC_vs_B2IC.csv')
write.csv(as.data.frame(egmt), 'RESULTS/05.pB2IC/GSEA_TF_Targets_pB2IC_vs_B2IC.csv')
write.csv(as.data.frame(egmt), 'RESULTS/05.pB2IC/GSEA_ChemicalGeneticPerturbations.csv')
write.csv(as.data.frame(egmt), 'RESULTS/05.pB2IC/GSEA_Custom_MM_datasets.csv')

# Plot individual GSEA
ID <- c('DELPUECH_FOXO3_TARGETS_DN')
ID <- c('MENSSEN_MYC_TARGETS', 'YU_MYC_TARGETS_UP', 'COLLER_MYC_TARGETS_UP',
        'DANG_MYC_TARGETS_UP', 'KIM_MYC_AMPLIFICATION_TARGETS_UP')

ID <- as.data.frame(egmt)
ID <- egmt[, c('Description')]
#ID <- egmt[grep('RLT', ID$Description), c('Description')]

pdf(paste0(ID,'_GSEA.pdf'), width = 6, height = 5)
gseaplot2(egmt, geneSetID = ID , title = ID, pvalue_table = F)
dev.off()

pdf('GSEA_MM_signature_ALL.pdf', width = 6, height = 5)
for (i in ID){
  plot <- gseaplot2(egmt, geneSetID = i , title = i, pvalue_table = F) 
  print(plot)}
dev.off()

as.data.frame(egmt[ID, 'core_enrichment'])

# Define FC vector for colors in plots
HS_Genes_FC <- HS_Genes[,':logFC']
names(HS_Genes_FC) <- as.character(HS_Genes$NCBI.gene..formerly.Entrezgene..ID)

cnetplot(egmt, foldChange = HS_Genes_FC)
heatplot(egmt, foldChange = HS_Genes_FC)

# Ridgeplot of top10 enriched/depleted categories
# Expression distributions (Log2FC) of core enriched genes for GSEA enriched categories
a <- egmt@result
a <- a[order(a$enrichmentScore), ]
a <- a[c(1:10, c((nrow(a)-9):nrow(a))), ]
a <- a[grep('MYELOMA', a$ID), ]
a <- a[grep('MYC', a$ID), ]
a <- a[grep('MENSSEN_MYC_TARGETS|YU_MYC_TARGETS_UP|COLLER_MYC_TARGETS_UP|DANG_MYC_TARGETS_UP|KIM_MYC_AMPLIFICATION_TARGETS_UP',
            a$ID), ]
egmt_selected <- egmt
egmt_selected@result <- a
View(as.data.frame(egmt_selected))

cnetplot(egmt_selected, foldChange = HS_Genes_FC, circular=T, colorEdge=T)
heatplot(egmt_selected, foldChange = HS_Genes_FC)


ridgeplot(egmt_selected, showCategory = 20, fill = 'pvalue', core_enrichment = T )

dotplot(egmt, color = 'pvalue')

pdf('Ridgeplot_MYELOMA_GSEA.pdf', width = 10, height = 6)
ridgeplot(egmt, showCategory = 20, fill = 'p.adjust', core_enrichment = F )
dev.off()

##### Heatmap GSEA results #####
# Save GSEA results for Hallmarks
B2IC <- as.data.frame(egmt)

# Order tables
models <- c('MIC', 'Maf_MIC', 'B2IC', 'pB2IC', 'B2IKC')
B2IKC <- B2IKC[order(rownames(B2IKC)), ]


B2IKC_MOD <- B2IKC
for (i in rownames(B2IKC_MOD)){
  if (B2IKC_MOD[i,'p.adjust'] > 0.05){
    B2IKC_MOD[i, 'NES'] = 0
  }
}


merged_GSEA_MOD <- cbind(MIC_MOD[ ,c(1,5,7)], Maf_MIC_MOD[ ,c(5,7)], B2IC_MOD[ ,c(5,7)],
                         pB2IC_MOD[ ,c(5,7)], B2IKC_MOD[ ,c(5,7)])
colnames(merged_GSEA_MOD)

colnames(merged_GSEA_MOD) <- c('ID', 'MIC_NES', 'MIC_p.adj',
                               'Maf_MIC_NES', 'Maf_MIC_p.adj',
                               'B2IC_NES', 'B2IC_p.adj',
                               'pB2IC_NES', 'pB2IC_p.adj',
                               'B2IKC_NES', 'B2IKC_p.adj')
rownames(merged_GSEA_MOD) <- gsub('HALLMARK_', '', rownames(merged_GSEA_MOD))

hmcol <- colorRampPalette(c('navy', 'skyblue4', 'skyblue', 'white', 'white', 'indianred2', 'indianred3', 'indianred'))(100)
hmcol <- colorRampPalette(c('navy', 'skyblue', 'white', 'white', 'pink', 'indianred3', 'indianred3', 'darkred'))(100)


# Select those categories enriched in at least one model
select_IDs <- merged_GSEA_MOD[, grep('p.adj', colnames(merged_GSEA_MOD))]
select_IDs$min <- apply(select_IDs, 1, FUN=min)
select_IDs <- select_IDs[select_IDs$min < 0.05, ]

# Or in all models
select_IDs <- merged_GSEA[, grep('p.adj', colnames(merged_GSEA))]
select_IDs$max <- apply(select_IDs, 1, FUN=max)
select_IDs <- select_IDs[select_IDs$max < 0.1, ]

# Create a matrix of NES
A <- merged_GSEA_MOD[, grep('NES', colnames(merged_GSEA_MOD))]
colnames(A) <- gsub('_NES', '', colnames(A))

# Plot Heatmap
pdf('RESULTS/GSEA_Combined_padj0.05.pdf')
pheatmap(A[rownames(select_IDs), ],
         treeheight_row=0,
         treeheight_col=0,
         annotation_legend=T,
         annotation_names_col=F,
         color = hmcol)
dev.off()

B <- merged_GSEA[rownames(select_IDs), ]

##### Heatmap GSEA results --- filtered by qvalue #####
# Save GSEA results for Hallmarks
B2IC <- as.data.frame(egmt)

# Order tables
models <- c('MIC', 'Maf_MIC', 'B2IC', 'pB2IC', 'B2IKC')
MIC <- MIC[order(rownames(MIC)), ]
Maf_MIC <- Maf_MIC[order(rownames(Maf_MIC)), ]
B2IC <- B2IC[order(rownames(B2IC)), ]
pB2IC <- pB2IC[order(rownames(pB2IC)), ]
B2IKC <- B2IKC[order(rownames(B2IKC)), ]

# Change NES of non significant GSEAs to == 0
MIC_M <- MIC
for (i in rownames(MIC_M)){
  if (MIC_M[i,'qvalues'] > 0.05){
    MIC_M[i, 'NES'] = 0
  }
}

Maf_MIC_M <- Maf_MIC
for (i in rownames(Maf_MIC_M)){
  if (Maf_MIC_M[i,'qvalues'] > 0.05){
    Maf_MIC_M[i, 'NES'] = 0
  }
}

B2IC_M <- B2IC
for (i in rownames(B2IC_M)){
  if (B2IC_M[i,'qvalues'] > 0.05){
    B2IC_M[i, 'NES'] = 0
  }
}

pB2IC_M <- pB2IC
for (i in rownames(pB2IC_M)){
  if (pB2IC_M[i,'qvalues'] > 0.05){
    pB2IC_M[i, 'NES'] = 0
  }
}

B2IKC_M <- B2IKC
for (i in rownames(B2IKC_M)){
  if (B2IKC_M[i,'qvalues'] > 0.05){
    B2IKC_M[i, 'NES'] = 0
  }
}

# Merge tables
merged_GSEA_M <- cbind(MIC_M[ ,c(1,5,8)], Maf_MIC_M[ ,c(5,8)], B2IC_M[ ,c(5,8)],
                       pB2IC_M[ ,c(5,8)], B2IKC_M[ ,c(5,8)])
colnames(merged_GSEA_M)
colnames(merged_GSEA_M) <- c('ID', 'MIC_NES', 'MIC_q.val',
                             'Maf_MIC_NES', 'Maf_MIC_q.val',
                             'B2IC_NES', 'B2IC_q.val',
                             'pB2IC_NES', 'pB2IC_q.val',
                             'B2IKC_NES', 'B2IKC_q.val')
#rownames(merged_GSEA_M) <- gsub('HALLMARK_', '', rownames(merged_GSEA_M))

# Select those categories enriched in at least one model
select_IDs <- merged_GSEA_M[, grep('q.val', colnames(merged_GSEA_M))]
select_IDs$min <- apply(select_IDs, 1, FUN=min)
select_IDs <- select_IDs[select_IDs$min < 0.05, ]

# Or in all models
select_IDs <- merged_GSEA[, grep('p.adj', colnames(merged_GSEA))]
select_IDs$max <- apply(select_IDs, 1, FUN=max)
select_IDs <- select_IDs[select_IDs$max < 0.1, ]

# Create a matrix of NES
A <- merged_GSEA_M[, grep('NES', colnames(merged_GSEA_M))]
colnames(A) <- gsub('_NES', '', colnames(A))

# Plot Heatmap
hmcol <- colorRampPalette(c('navy', 'skyblue4', 'skyblue', 'white', 'white', 'indianred2', 'indianred3', 'indianred'))(100)
hmcol <- colorRampPalette(c('navy', 'skyblue', 'white', 'white', 'pink', 'indianred3', 'indianred3', 'darkred'))(100)
hmcol <- colorRampPalette(c('navy', 'skyblue', 'white', 'white', 'pink', 'darkred'))(100)


pdf('RESULTS/GSEA_MMsignatures_qval0.05.pdf')
pheatmap(A[rownames(select_IDs), ], 
         treeheight_row=0,
         treeheight_col=0,
         annotation_legend=T,
         annotation_names_col=F,
         color = hmcol)
dev.off()

B <- merged_GSEA[rownames(select_IDs), ]