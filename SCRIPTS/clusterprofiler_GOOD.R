############################
#CLUSTERPROFILER GSEA HUMAN#
############################

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
setwd("C:/Users/transbio/Desktop/MM_analysis/mouse_ensembl/")

limma.res <- read.table('Results/MM_C_GSEA/results_MM_symbol.txt',sep="\t")
colnames(limma.res)

rownames(limma.res)<-limma.res$ensembl_gene_id
# Filter by contrast
Contrast <- "MafB2IC_C:"
limma.res.filt <- limma.res[, grep(Contrast, colnames(limma.res))]
colnames(limma.res.filt)
colnames(limma.res.filt) <- gsub(paste0(".", Contrast), '', colnames(limma.res.filt))
colnames(limma.res.filt)

##### Convert mouse IDs to homolog human IDs with BIOMART #####
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")


genesV2 = getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", 
                 values = rownames(limma.res.filt) , mart = mouse, 
                 attributesL = c("ensembl_gene_id"), martL = human, uniqueRows = T)

HS_Genes <- merge(limma.res.filt, genesV2, by.x='row.names', by.y='Gene.stable.ID')

### Create a decreasing sorted vector according to FC values
ranked_genes <- HS_Genes[ ,'logFC']
names(ranked_genes) <- HS_Genes$Gene.stable.ID.1
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
gmtfile <- read.gmt('C:/Users/transbio/Desktop/MM_analysis/MM_mouse_geneset.gmt.txt')

### Calculate GSEA
egmt <- GSEA(ranked_genes, TERM2GENE=gmtfile, verbose=FALSE, pvalueCutoff = 1, maxGSSize = 5000)
dim(as.data.frame(egmt))
View(as.data.frame(egmt))

pB2IC <- as.data.frame(egmt)

write.csv(as.data.frame(egmt), 'GSEA_B2IKC_miren.csv')
write.csv(as.data.frame(egmt), 'RESULTS/05.pB2IC/GSEA_TF_Targets_pB2IC_vs_B2IC.csv')
write.csv(as.data.frame(egmt), 'RESULTS/05.pB2IC/GSEA_ChemicalGeneticPerturbations.csv')
write.csv(as.data.frame(egmt), 'RESULTS/05.pB2IC/GSEA_Custom_MM_datasets.csv')

# Plot individual GSEA
ID <- c('Sig_Down_0.05')
ID <- c('MENSSEN_MYC_TARGETS', 'YU_MYC_TARGETS_UP', 'COLLER_MYC_TARGETS_UP',
        'DANG_MYC_TARGETS_UP', 'KIM_MYC_AMPLIFICATION_TARGETS_UP')

ID <- as.data.frame(egmt)
ID <- egmt[, c('Description')]
#ID <- egmt[grep('RLT', ID$Description), c('Description')]
pdf("MIC_sig_up.pdf")
pdf(paste0('', ID,'_GSEA.pdf'), width = 6, height = 5)
gseaplot2(egmt, geneSetID = ID , title = ID, pvalue_table = F)
dev.off()

pdf('RESULTS/05.pB2IC/GSEA_MM_signature_ALL.pdf', width = 6, height = 5)
for (i in ID){
  plot <- gseaplot2(egmt, geneSetID = i , title = i, pvalue_table = F) 
  print(plot)}
dev.off()

as.data.frame(egmt[ID, 'core_enrichment'])

# Define FC vector for colors in plots
HS_Genes_FC <- HS_Genes[,'Coef']
names(HS_Genes_FC) <- as.character(HS_Genes$Gene.name.1)

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


ridgeplot(egmt, showCategory = 20, fill = 'pvalue', core_enrichment = T )

dotplot(egmt)

pdf('RESULTS/B2IC_Model/Ridgeplot_MYELOMA_GSEA.pdf', width = 10, height = 6)
ridgeplot(egmt_selected, showCategory = 20, fill = 'pvalue', core_enrichment = T )
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