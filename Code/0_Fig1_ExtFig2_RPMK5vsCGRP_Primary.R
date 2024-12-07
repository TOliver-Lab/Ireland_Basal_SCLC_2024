# RPM K5 Cre vs CGRP Cre Primary Tumor Analysis
# Ireland et al., BioRxiv, 2024

setwd("/Users/abbieireland/Desktop/scRNAseq")


# Load necessary packages
#load necessary packages
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(SingleCellExperiment)
  library(ggplot2)
  library(cowplot)
  library(Seurat)
  library(SeuratObject)
  library(CellTagR)
  library(viridis)
  library(ggpubr)
  library(SeuratData)
  library(SeuratDisk)
  library(zellkonverter)
})


###############################################################################################################
############## Seurat conversion and analysis using leiden clustering from scanpy ##########
########################################################################################################

######### Now, bring into Seurat, from  anndata object processed in scanpy/scvi ##########################################################

# Read in adata object as SCE
getwd()
adata<-readH5AD("09.04.24_RPM_K5vCGRP_adata6.h5ad")
counts(adata)
adata
table(adata$Cre)
table(adata$leiden_scVI_1.6)
# 0    1    2    3    4    5 
# 8407 4379 3730 1244 1243  925 

table(adata$UnID)
# RPM_CGRP1    RPM_CGRP2    RPM_CGRP3    RPM_CGRP4    RPM_CGRP5 RPM_K5_total 
# 1236          217           82          122          240         6546 
# RPM_K5c      RPM_K5t 
# 1946         9539 

#Convert SCE to seurat
RPM_K5vCGRP <- CreateSeuratObject(counts = counts(adata), meta.data = as.data.frame(colData(adata)))
RPM_K5vCGRP

# An object of class Seurat 
# 55491 features across 19928 samples within 1 assay 
# Active assay: RNA (55491 features, 0 variable features)
# 1 layer present: counts

DefaultAssay(RPM_K5vCGRP)
adata

#Create an assay of normalized gene expression
norm_assay <- CreateAssayObject(counts = assay(adata,"norm"))
RPM_K5vCGRP[["norm"]] <- norm_assay

# Set assay norm as default assay for seurat object
DefaultAssay(RPM_K5vCGRP)<-'norm'
table(RPM_K5vCGRP@meta.data$UnID)

# Add embeddings of umap, and X_scVI_1.4 to assay norm
adata
reducedDim(adata, "X_umap")
test<-reducedDim(adata, "X_umap")
colnames(test)<-c("UMAP_1","UMAP_2")
rownames(test)<-colnames(RPM_K5vCGRP)
dim(test)
RPM_K5vCGRP[['umap']] <- CreateDimReducObject(test, key="UMAP_", assay = "RNA")
RPM_K5vCGRP[['umap']]@cell.embeddings

# Fig 1h
unid_cols<-c("darkorchid4","orange")
DimPlot(RPM_K5vCGRP,split.by='Cre',group.by='Cre',cols=unid_cols,reduction='umap',shuffle=TRUE)+NoAxes()
DimPlot(RPM_K5vCGRP,group.by='Cre',cols=unid_cols,reduction='umap',shuffle=TRUE)+NoAxes()

# Fig 1i
colors<-c('#F39B7F99','brown1','#ff7f0e', 'turquoise', '#2ca02c', 'deepskyblue4', '#e7298a',  '#0aa6d8',  '#A0522D', '#a1c299', 'gold', '#984ea3', '#3C548899', '#377eb8',  '#d62728','navyblue', '#a1c299', 'purple','gray','maroon','gray20')
DimPlot(RPM_K5vCGRP,group.by='leiden_scVI_1.6',cols=colors, reduction='umap',label=TRUE,label.size=6)&NoAxes()


# For visualization of scores, logarithmize and normalize raw counts in Seurat
DefaultAssay(RPM_K5vCGRP)<-'RNA'
RPM_K5vCGRP<-NormalizeData(RPM_K5vCGRP)
RPM_K5vCGRP

# Now can apply signatures
# Read in this table to convert between mouse and human homologs
mouse94<-read.csv("Signatures/mouse94.csv")


### Look at distribution of leiden clusters in tumors by Cre (Fig. 1i)
library(ggplot2)
library(ggpubr)

##### What % of cells occupy what cluster for all types ####
Idents(RPM_K5vCGRP)<-'Cre'
Idents(RPM_K5vCGRP)

x<-table(Idents(RPM_K5vCGRP),RPM_K5vCGRP@meta.data$leiden_scVI_1.6)
x
proportions <- as.data.frame(100*prop.table(x, margin = 1))
proportions

colnames(proportions)<-c("Cluster", "Sample", "Frequency")

ggbarplot(proportions, x="Sample", y="Frequency", fill = "Sample", group = "Sample", ylab = "Frequency (percent)", xlab="Phase", palette =c("darkorchid4","orange"))+ theme_bw()+ facet_wrap(facets = "Cluster", scales="free_y", ncol =4)+ theme(axis.text.y = element_text(size=12)) +rotate_x_text(angle = 45)

# Stacked
p<-ggplot(proportions, aes(fill=Sample, y=Frequency, x=Cluster)) +
  geom_bar(position="stack", stat="identity")

p + scale_fill_manual(values=colors)+ theme_bw()+ theme(axis.text.y = element_text(size=20), axis.text.x=element_text(size=20), axis.title.x =element_text(size=20), axis.title.y = element_text(size=20), legend.text = element_text(size=20), legend.title = element_text(size=20))


# Visualize A, N, P expression by split violin plot (Fig. 1j)
library(viridis)
VlnPlot(RPM_K5vCGRP,features = c("Ascl1","Neurod1","Pou2f3"), same.y.lims=FALSE,group.by=c("Genotype"),split.by=c("Cre"),split.plot=TRUE,cols=c("darkorchid4","orange"),pt.size=.01, ncol=3)


############### APPLY VARIOUS SIGNATURES ######################
## Archetype signature analysis (Extended Data Fig. 2c) ##

arch<-read.csv("Signatures/Archetype_Sigs_Maddox.csv")
a<-arch$SCLC.A
a2<-arch$SCLC.A2
n<-arch$SCLC.N
p<-arch$SCLC.P
y<-arch$SCLC.Y

a<-a[1:987]
a_sc<-subset(mouse94, mouse94$human_homolog %in% a)
a_sc_mouse<-(a_sc$gene_name)
a_sc_mouse

n_sc<-subset(mouse94, mouse94$human_homolog %in% n)
n_sc_mouse<-(n_sc$gene_name)
n_sc_mouse

a2
a2_sc<-subset(mouse94, mouse94$human_homolog %in% a2)
a2_sc_mouse<-(a2_sc$gene_name)
a2_sc_mouse

p
p_sc<-subset(mouse94, mouse94$human_homolog %in% p)
p_sc_mouse<-(p_sc$gene_name)
p_sc_mouse

y
y_sc<-subset(mouse94, mouse94$human_homolog %in% y)
y_sc_mouse<-(y_sc$gene_name)
y_sc_mouse


RPM_K5vCGRP<-AddModuleScore(
  object = RPM_K5vCGRP,
  features = list(a_sc_mouse),
  name = 'A_Archetype')

RPM_K5vCGRP<-AddModuleScore(
  object = RPM_K5vCGRP,
  features = list(a2_sc_mouse),
  name = 'A2_Archetype')

RPM_K5vCGRP<-AddModuleScore(
  object = RPM_K5vCGRP,
  features = list(n_sc_mouse),
  name = 'N_Archetype')

RPM_K5vCGRP<-AddModuleScore(
  object = RPM_K5vCGRP,
  features = list(p_sc_mouse),
  name = 'P_Archetype')

RPM_K5vCGRP<-AddModuleScore(
  object = RPM_K5vCGRP,
  features = list(y_sc_mouse),
  name = 'Y_Archetype')

library(viridis)

#VlnPlot(RPM_K5vCGRP,features = c("A_Archetype1","A2_Archetype1","N_Archetype1","P_Archetype1"), group.by=c("Genotype"),same.y.lims=FALSE,split.by=c("Cre"),split.plot=TRUE,cols=c("darkorchid4","orange"),pt.size=.01, ncol=5)
VlnPlot(RPM_K5vCGRP,features = c("A_Archetype1","N_Archetype1","P_Archetype1"), group.by=c("Genotype"),same.y.lims=FALSE,split.by=c("Cre"),split.plot=TRUE,cols=c("darkorchid4","orange"),pt.size=.01, ncol=3)


#############################################################################
# Apply human A N P scRNA seq sigs from Joe Chan. (Extended Data Fig. 2d)
sc_sclc_sigs<-read.csv("Signatures/hSCLC_Chan_sigs.csv")
sc_sclc_sigs

sc_sclc_sigs$hSCLC_A
a_sc<-sc_sclc_sigs$hSCLC_A[1:67]
a_sc<-subset(mouse94, mouse94$human_homolog %in% a_sc)
a_sc_mouse<-(a_sc$gene_name)
a_sc_mouse

sc_sclc_sigs$hSCLC_N

n_sc<-sc_sclc_sigs$hSCLC_N[1:73]
n_sc<-subset(mouse94, mouse94$human_homolog %in% n_sc)
n_sc_mouse<-(n_sc$gene_name)
n_sc_mouse

sc_sclc_sigs$hSCLC_P

p_sc<-sc_sclc_sigs$hSCLC_P
p_sc<-subset(mouse94, mouse94$human_homolog %in% p_sc)
p_sc_mouse<-(p_sc$gene_name)
p_sc_mouse


RPM_K5vCGRP<-AddModuleScore(
  object = RPM_K5vCGRP,
  features = list(a_sc_mouse),
  name = 'hSCLC_A')

RPM_K5vCGRP<-AddModuleScore(
  object = RPM_K5vCGRP,
  features = list(n_sc_mouse),
  name = 'hSCLC_N')

RPM_K5vCGRP<-AddModuleScore(
  object = RPM_K5vCGRP,
  features = list(p_sc_mouse),
  name = 'hSCLC_P')

library(viridis)
library(ggplot2)
library(gridExtra)

VlnPlot(RPM_K5vCGRP,features = c("hSCLC_A1","hSCLC_N1","hSCLC_P1"), same.y.lims=FALSE,group.by=c("Genotype"),split.by=c("Cre"),split.plot=TRUE,cols=c("darkorchid4","orange"),pt.size=.01, ncol=3)

###############################################
## Look at ND1 vs ASCl1 ChIP targets as scores (Extended Data Fig. 1e)

chip<-read.csv("Signatures/ASCL1_NEUROD1_POU2F3_ChIP_Targets.csv")
achip<-chip$Borromeo_ASCL1_Targets
nchip<-chip$Borromeo_Oliver_NEUROD1_Targets
pchip<-chip$POU2F3_Targets
pchip<-subset(mouse94, mouse94$human_homolog %in% pchip)
pchip_m<-(pchip$gene_name)
pchip_m

RPM_K5vCGRP<-AddModuleScore(
  object = RPM_K5vCGRP,
  features = list(achip),
  name = 'ASCL1_Targets')

RPM_K5vCGRP<-AddModuleScore(
  object = RPM_K5vCGRP,
  features = list(nchip),
  name = 'NEUROD1_Targets')

RPM_K5vCGRP<-AddModuleScore(
  object = RPM_K5vCGRP,
  features = list(pchip_m),
  name = 'POU2F3_Targets')


VlnPlot(RPM_K5vCGRP,features = c("ASCL1_Targets1","NEUROD1_Targets1","POU2F3_Targets1"), same.y.lims=FALSE,group.by=c("Genotype"),split.by=c("Cre"),
        split.plot=TRUE,cols=c("darkorchid4","orange"),pt.size=.01, ncol=3)



######## Add NE Score ##########

### For NE score assignemnt

#### Converting to SCE to Add NE Score and create Violin plots
# Convert seurat object to SCE
library(SingleCellExperiment)
library(SummarizedExperiment)

sce<-as.SingleCellExperiment(RPM_K5vCGRP)

# import the NE/Non-NE data from the literature
nedat_mouse<-read.csv("Signatures/nedat_mouse.csv")
nedat_mouse$gene_name
gene_idx <- rownames(sce) %in% nedat_mouse$gene_name
table(gene_idx)
X <- as.matrix(assay(sce, "logcounts")[gene_idx,])
X
# assure the gene order matches up
X <- X[nedat_mouse$gene_name,]
head(X)
dim(X)

# define the scoring function based on methods from Zhang et al TLCR 2018
ne_score <- function(cell, ne, notne, ...) { 
  # report missing value is no variation 
  if (sd(cell) == 0) { 
    NA
  } else { 
    (cor(x = cell, y = ne, ...) - cor(x = cell, y = notne, ...))/2
  }
}

# get the NE scores. note that 4 cells have missing values b/c 
# they had zero variation across the 50 genes. 

sce$NE_spearman <- apply(X = X, 2, ne_score,  
                         ne = nedat_mouse$NE, 
                         notne = nedat_mouse$NonNE, 
                         method = "spearman")


#Add NE score as metadata to Seurat object
ne_score<-sce$NE_spearman
RPM_K5vCGRP@meta.data$NE_spearman<-ne_score

VlnPlot(RPM_K5vCGRP,features = c("NE_spearman"), same.y.lims=FALSE,group.by=c("Genotype"),split.by=c("Cre"),split.plot=TRUE,cols=c("darkorchid4","orange"),pt.size=.01, ncol=1)



# Apply Student's unpaired t-tests to each split violin comparison using stat_compare_means
# For example, to compare Cgrp to K5 for SCLC-A archetype...
my_comparisons=list(c("CGRP","K5"))
t<-VlnPlot(RPM_K5vCGRP, features = c("A_Archetype1"), group.by=c("Cre"),pt.size=.01, y.max=0.4)

# Add pairwise comparisons p-value
t + stat_compare_means(comparisons = my_comparisons, method="t.test", paired=FALSE) 
# Add global p-value
t + stat_compare_means(comparisons = my_comparisons, method="t.test", paired=FALSE,label = "p.signif") # Add pairwise comparisons p-value

# Repeat above for all desired comparisons


# Save resulting seurat object with signatures
saveRDS(RPM_K5vCGRP,"110224_RPMK5vCGRP_Norm_Allsigs.rds")

