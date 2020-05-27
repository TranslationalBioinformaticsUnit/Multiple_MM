library(biomaRt)


counts_table <- read.table('C:/Users/transbio/Desktop/MM_analysis/Huma_data/norm_data_rownames.txt', header=T, sep='\t', row.names = 1)
head(counts_table)

colnames(counts_table) <- gsub('X', '', colnames(counts_table))
colnames(counts_table)

biomartMM = useMart("ensembl", dataset = "hsapiens_gene_ensembl", ensemblRedirect = F)
atributos = listAttributes(biomartMM)    

atributos[grep("hgnc_id", atributos$name, ignore.case = TRUE),]

myannot = getBM(attributes = c("ensembl_gene_id", "external_gene_name","percentage_gene_gc_content", "gene_biotype", 'chromosome_name', 'hgnc_id'),
                filters = "ensembl_gene_id", values=rownames(counts_table), mart=biomartMM)


IG_genes <- myannot[grep('IG', myannot$gene_biotype), ]
counts_table_filt <- counts_table[!(rownames(counts_table) %in% IG_genes$ensembl_gene_id), ]
myannot_filt <- myannot[!(myannot$ensembl_gene_id %in% IG_genes$ensembl_gene_id), ]
