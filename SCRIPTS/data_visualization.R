###########################
#DATA VISUALIZATION HUMAN##
###########################

###DATA VISUALIZATION 

#-----LOAD LIBRARYS AND WORK DIRECTORY
library(ggplot2)
setwd("C:/Users/transbio/Desktop/MM_analysis/mouse_ensembl")


#-----LOAD COUNT TABLE AND DESIGN MATRIX


counts_table <- read.table('final_norm_def.txt', header=T, sep='', row.names = 1)
head(counts_table)

colnames(counts_table) <- gsub('X', '', colnames(counts_table))
colnames(counts_table)

##read the metadata of the data

design_matrix <- read.csv('Data/Col_Data.csv', sep=";", header=T, stringsAsFactors = F)
colnames(design_matrix)[1] <- 'Sample'
head(design_matrix)
str(design_matrix)
rownames(design_matrix)<- design_matrix$Coded_Name
rownames(design_matrix)<- design_matrix$Sample

design_matrix$Sample <-  as.factor(design_matrix$Sample)
design_matrix$Stage<-as.factor(design_matrix$Stage)
design_matrix$Combined<-as.factor(design_matrix$Combined)
design_matrix$Model<-as.factor(design_matrix$Model)
#design_matrix$Coded_Name <-  as.factor(design_matrix$Coded_Name)
#design_matrix$Diagnosis <-  as.factor(design_matrix$Diagnosis)
#design_matrix$Cell_Type <-  as.factor(design_matrix$Cell_Type)
#design_matrix$Run <-  as.factor(design_matrix$Run)
#design_matrix$Buffer<-  as.factor(design_matrix$Buffer)
#design_matrix$Hosp<-  as.factor(design_matrix$Hosp)
str(design_matrix)


#------PLOT THE MDS

##modify the colors
library(RColorBrewer)
#display.brewer.all() #view all colors 

# Labels colors
colors <- c(brewer.pal(n=8, name = 'Dark2'), brewer.pal(n=8, name='Set3'),brewer.pal(n=12, name = 'Paired'))

color_types <- colors[1:length(levels(as.factor(design_matrix$Model)))]
names(color_types) <- levels(design_matrix$Model)

##plot the mds
tcounts<-t(counts_table)

d<- dist(tcounts)  #euclidean distance between samples
fit<-cmdscale(d, eig=TRUE, k=2)


#plot MDS with raw data
x<- fit$points[,1]
y<- fit$points[,2]

mplot <- data.frame(Coordinate_1=x, Coordinate_2=y, Combination_type=design_matrix$Model,shape=ddesign_matrix$Stage)
p1 <- ggplot(mplot, aes(x=Coordinate_1, y=Coordinate_2, color = Combination_type, label=rownames(mplot))) +
  geom_point(size=3) + 
  scale_color_manual(values=color_types) + 
  labs(title="MDS - RAW Count table",x = "Coordinate 1", y = "Coordinate 2")

ggsave(plot = p1, width = 9, height = 7, dpi = 100,filename = "RESULTS/raw_count_table_complete.pdf")
ggsave(plot = p1+geom_text(size=2, aes(label=rownames(mplot))), width = 9, height = 7, dpi = 100, filename = "raw_count_table__completed_annotated.pdf")





pca_rna <- prcomp(t(counts_table), 
                  center = T, scale = F)

summary(pca_rna)

# Plot PCs
# Prepare data for the plot
df <- data.frame(PC1=pca_rna$x[,1], PC2=pca_rna$x[,2],
                 design_matrix_filt[,
                c('Diagnosis','Run', 'Buffer', 'Hosp')])


# Plot according to different factors
sum=summary(pca_rna)
imp=sum$importance
perc1=imp[2,1]*100
perc2=imp[2,2]*100

var <- c('Diagnosis', 'Run', 'Buffer', 'Hosp')
pdf('PCA_raw_filt.pdf', width = 7, height = 5)
for (v in var) {
  print(ggplot(df) +
          geom_point(aes(x=PC1,y=PC2,color=factor(df[ , v])),size=6,shape=20) +
          labs(x= paste('PC1 (', perc1, ')', sep=''), y= paste('PC2 (', perc2, ')', sep='')) +
          ggtitle(paste0('Principal Components for RNA data by ', v)))
}
dev.off()
#--Exclude mgus 3 and 4
samples_excluded<-c("MGUS_3","MGUS_4","MGUS_2","BMPC_5_2","BMPC_8_3","MGUS_10","MGUS_6")
count_table2 <- counts_table_short[,!colnames(counts_table_short)%in%samples_excluded]
design_matrix2<- design_matrix[colnames(count_table2),]
norm_data2<-norm_data[,!colnames(counts_table)%in%samples_excluded]

#--PCAs WITHOUT OUTLIERS

color_types <- colors[1:length(levels(as.factor(design_matrix$Diagnosis)))]
names(color_types) <- levels(design_matrix$Diagnosis)

##plot the mds
tcounts<-t(count_table2)

d<- dist(tcounts)  #euclidean distance between samples
fit<-cmdscale(d, eig=TRUE, k=2)


#plot solution
x<- fit$points[,1]
y<- fit$points[,2]

mplot <- data.frame(Coordinate_1=x, Coordinate_2=y, Combination_type=design_matrix2$Diagnosis)
p1 <- ggplot(mplot, aes(x=Coordinate_1, y=Coordinate_2, color = Combination_type, label=rownames(mplot))) +
  geom_point(size=3) + 
  scale_color_manual(values=color_types) + 
  labs(title="MDS - raw Count table",x = "Coordinate 1", y = "Coordinate 2")

ggsave(plot = p1, width = 9, height = 7, dpi = 100, filename = "raw_filt_short.pdf")
ggsave(plot = p1+geom_text(size=2, aes(label=rownames(mplot))), width = 9, height = 7, dpi = 100, filename = "raw_filt_anno_short.pdf")


require(edgeR)

#--NORMALIZATION
#create DGEList object
d0 <- DGEList(counts=counts_table)

#filter: remove rows that consistently have zero or very low counts
Diagnosis <- design_matrix$Diagnosis
buffer <- design_matrix$Buffer
design <- model.matrix(~0 + Diagnosis )

keep <- filterByExpr(d0, design, group=Diagnosis, min.count = 1) #21.894

d1 <- d0[keep,,keep.lib.sizes=FALSE]

#calculate normalization factors (TMM)
d1_norm <- calcNormFactors(d1)

#voom transformation
v <- voom(d1_norm, design, plot=TRUE)

norm_data <- v$E

design_matrix_ordered <- design_matrix[colnames(norm_data),]  
write.table(norm_data,"DATA/norm_data_not_TPC.txt", sep="\t", dec=".", quote=FALSE, row.names=TRUE, col.names=TRUE)

# Labels colors
colors <- c(brewer.pal(n=8, name = 'Dark2'), brewer.pal(n=8, name='Set3'),brewer.pal(n=12, name = 'Paired'))

color_types <- colors[1:length(levels(as.factor(design_matrix$Diagnosis)))]
names(color_types) <- levels(design_matrix$Diagnosis)

##plot the mds
tcounts<-t(norm_data)

d<- dist(tcounts)  #euclidean distance between samples
fit<-cmdscale(d, eig=TRUE, k=2)


#plot solution
x<- fit$points[,1]
y<- fit$points[,2]

mplot <- data.frame(Coordinate_1=x, Coordinate_2=y, Combination_type=design_matrix$Diagnosis)
p1 <- ggplot(mplot, aes(x=Coordinate_1, y=Coordinate_2, color = Combination_type, label=rownames(mplot))) +
  geom_point(size=3) + 
  scale_color_manual(values=color_types) + 
  labs(title="MDS - normalized Count table",x = "Coordinate 1", y = "Coordinate 2")

ggsave(plot = p1, filename = "RESULTS/normalized_data_no_TPC.pdf")
ggsave(plot = p1+geom_text(size=2, aes(label=rownames(mplot))), width = 9, height = 7, dpi = 100, filename = "RESULTS/normalized_annotated_no_TPC.pdf")
pca_rna <- prcomp(t(norm_data), 
                  center = T, scale = F)

summary(pca_rna)

# Plot PCs
# Prepare data for the plot
df <- data.frame(PC1=pca_rna$x[,1], PC2=pca_rna$x[,2],
                 design_matrix[colnames(norm_data),
                               c('Diagnosis', 'Buffer', 'Hosp')])


# Plot according to different factors
sum=summary(pca_rna)
imp=sum$importance
perc1=imp[2,1]*100
perc2=imp[2,2]*100

var <- c('Diagnosis', 'Buffer', 'Hosp')
pdf('norm_data.pdf', width = 7, height = 5)
for (v in var) {
  print(ggplot(df) +
          geom_point(aes(x=PC1,y=PC2,color=factor(df[ , v])),size=6,shape=20) +
          labs(x= paste('PC1 (', perc1, ')', sep=''), y= paste('PC2 (', perc2, ')', sep='')) +
          ggtitle(paste0('Principal Components for RNA data by ', v)))
}
dev.off()

#--DENSITY PLOT AND HISTOGRAM
library(scales)
sampleNames = vector()
intensities = vector()
diagnosis= vector()

counts<-counts_table+1 # convert all the 0s in 1 to the log10

for (i in 1:ncol(counts)){
  sampleNames = c(sampleNames,rep(colnames(counts)[i],nrow(counts)))
  intensities = c(intensities,counts[,i])
  diagnosis<-c(diagnosis,rep(design_matrix$Diagnosis[i],nrow(counts)))
}
arrayData <- data.frame(intensities,sampleNames,diagnosis)
g = ggplot(arrayData, aes(intensities, color=diagnosis ,fill=diagnosis))

g + geom_histogram(position="identity",binwidth=1) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  
  labs(title='Gene Expression Density Curve',
       x='Counts',
       y = 'Density')

#--REMOVE OUTLIERS AND PCAs
samples_excluded<-c("BMPC_5_2","BMPC_8_3","MGUS_10","MGUS_6")
count_table2 <- counts_table2[,!colnames(counts_table2)%in%samples_excluded]
design_matrix2<- design_matrix[colnames(norm_data2),]
norm_data2<-norm_data[,!colnames(count_table2)%in%samples_excluded]
#--PCAs WITHOUT OUTLIERS

color_types <- colors[1:length(levels(design_matrix2$Diagnosis))]
names(color_types) <- levels(design_matrix2$Diagnosis)

##plot the mds
tcounts<-t(norm_data2)

d<- dist(tcounts)  #euclidean distance between samples
fit<-cmdscale(d, eig=TRUE, k=2)


#plot solution
x<- fit$points[,1]
y<- fit$points[,2]

mplot <- data.frame(Coordinate_1=x, Coordinate_2=y, Combination_type=design_matrix2$Diagnosis)
p1 <- ggplot(mplot, aes(x=Coordinate_1, y=Coordinate_2, color = Combination_type, label=rownames(mplot))) +
  geom_point(size=3) + 
  scale_color_manual(values=color_types) + 
  labs(title="MDS - normalized Count table",x = "Coordinate 1", y = "Coordinate 2")

ggsave(plot = p1, width = 9, height = 7, dpi = 100, filename = "norm_diag_short_filt.pdf")
ggsave(plot = p1+geom_text(size=2, aes(label=rownames(mplot))), width = 9, height = 7, dpi = 100, filename = "norm_diag_short_filt_anno.pdf")


##mouse_data visualization
rownames(design_matrix)<- design_matrix$Sample
design_matrix<- design_matrix[colnames(counts_table),]
#------PLOT THE MDS

##modify the colors
library(RColorBrewer)
display.brewer.all() #view all colors 

# Labels colors
colors <- c(brewer.pal(n=8, name = 'Dark2'), brewer.pal(n=8, name='Set3'),brewer.pal(n=8, name = 'Paired2'))

color_types <- colors[1:length(levels(design_matrix$Combined))]
names(color_types) <- levels(design_matrix$Combined)

##plot the mds
tcounts<-t(counts_table)

d<- dist(tcounts)  #euclidean distance between samples
fit<-cmdscale(d, eig=TRUE, k=2)


#plot solution
x<- fit$points[,1]
y<- fit$points[,2]

mplot <- data.frame(Coordinate_1=x, Coordinate_2=y, Combination_type=design_matrix$Combined , Stage=design_matrix$Stage)
p1 <- ggplot(mplot, aes(x=Coordinate_1, y=Coordinate_2, color = Combination_type, label=rownames(mplot),shape=Stage)) +
  geom_point(size=3, fill=design_matrix$Color) +
  geom_point(size=3, color="black") +
  scale_shape_manual(values=c(25,23,21,24,22))+
  scale_color_manual(values=color_types) + 
  labs(title="MDS - Normalized Count table",x = "Coordinate 1", y = "Coordinate 2")

ggsave(plot = p1, filename = "Results/norm_count_table_mds.pdf")

#-------EXCLUDE THE TWO OUTLIERS AND THE B-CELLS


samples_excluded<-c("YFP50_YFP_Cg1_GC_B_cell_S23","YFP51_YFP_Cg1_GC_B_cell_S23","YFP72_YFP_Cg1_GC_B_cell_S24")

count_table2 <- counts_table[,!colnames(counts_table)%in%samples_excluded]


design_matrix2<- design_matrix[colnames(count_table2),]

colors <- c(brewer.pal(n=8, name = 'Dark2'), brewer.pal(n=8, name='Set3'),brewer.pal(n=12, name = 'Paired'))

color_types <- colors[1:length(levels(design_matrix2$Combined))]
names(color_types) <- levels(design_matrix2$Combined)

##plot the mds
tcount2<-t(count_table2)

d<- dist(tcount2)  #euclidean distance between samples
fit<-cmdscale(d, eig=TRUE, k=2)


#plot solution
x<- fit$points[,1]
y<- fit$points[,2]


mplot <- data.frame(Coordinate_1=x, Coordinate_2=y, Combination_type=design_matrix2$Combined, Stage=design_matrix2$Stage)
#TSNE1=raw_tsne$Y[,1], TSNE2=raw_tsne$Y[,2], 
p1 <- ggplot(mplot, aes(x=Coordinate_1, y=Coordinate_2, color = Combination_type, label=rownames(mplot),shape=Stage)) +
  geom_point(size=3, fill=design_matrix2$Color) +
  geom_point(size=3, color="black") +
  scale_shape_manual(values=c(23,21,24,22))+
  scale_color_manual(values=color_types) + 
  labs(title="MDS - Normalized Count table, NO B-cells",x = "Coordinate 1", y = "Coordinate 2")

ggsave(plot = p1, width = 9, height = 7, dpi = 100, filename = "Results/norm_count_table_mds_NO_Bcells.pdf")




