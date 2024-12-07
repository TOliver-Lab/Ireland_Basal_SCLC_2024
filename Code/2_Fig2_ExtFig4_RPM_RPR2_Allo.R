# RPM vs RPR2 TBO Allograft script
# 09/06/24 #

######## MOVING INTO SEURAT NOW ############

#### Pull in adata object from Scanpy 

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
})

############## Seurat conversion and analysis using xenograft leiden clustering from scanpy ##########

########################################################################################################
########################################################################################################
######### Now, bring into Seurat, from my anndata object processed according to onb paper methods ########################################################################################################
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(zellkonverter)

# Read in adata object as SCE
adata<-readH5AD("090824_adata_3_RPM_RPR2_TBOAllo.h5ad")
counts(adata)
table(adata$Cre)
table(adata$leiden_scVI_1.3)
# 0    1    2    3    4    5    6    7    8    9   10 
# 3010 2593 2112 1540 1519 1423 1199  698  506  317   44 
table(adata$UnID)
# RPM_Allo RPR2_Allo 
# 6440      8521 

#Convert SCE to seurat
RPMvRPR2 <- CreateSeuratObject(counts = counts(adata), meta.data = as.data.frame(colData(adata)))
RPMvRPR2
# 55491 features across 14961 samples within 1 assay 
# Active assay: RNA (55491 features, 0 variable features)
# 1 layer present: counts
DefaultAssay(RPMvRPR2)

#Create an assay of normalized gene expression
norm_assay <- CreateAssayObject(counts = assay(adata,"norm"))

RPMvRPR2[["norm"]] <- norm_assay


# Set assay norm as default assay for seurat object
DefaultAssay(RPMvRPR2)<-'norm'

adata
# Add embeddings of umap, and X_scVI_1.4 to assay norm
adata
reducedDim(adata, "X_umap")
test<-reducedDim(adata, "X_umap")
colnames(test)<-c("UMAP_1","UMAP_2")
rownames(test)<-colnames(RPMvRPR2)
dim(test)
RPMvRPR2[['umap']] <- CreateDimReducObject(test, key="UMAP_", assay = "RNA")
RPMvRPR2[['umap']]@cell.embeddings

colors<-c('#F39B7F99','brown1','#ff7f0e', 'turquoise', '#2ca02c', 'deepskyblue4', '#e7298a',  '#4ba4d3',  '#A0522D', '#a1c299', 'gold', '#984ea3', '#3C548899', '#377eb8',  '#d62728','navyblue', '#a1c299', 'purple','gray','maroon','gray20')
DimPlot(RPMvRPR2,group.by='leiden_scVI_1.3',cols=colors, reduction='umap',label=FALSE,label.size=6)&NoAxes()

unid_cols<-c("turquoise3","maroon")
DimPlot(RPMvRPR2,group.by='Genotype',cols=unid_cols,reduction='umap',shuffle=TRUE)+NoAxes()


# For visualization and signatures, apply log normalization to counts assay
DefaultAssay(RPMvRPR2)<-'RNA'
RPMvRPR2<-NormalizeData(RPMvRPR2)

## Cell Cycle scoring
mouse94<-read.csv("Signatures/mouse94.csv")
s.genes<-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
s.genes.list<-subset(mouse94, mouse94$human_homolog %in% s.genes)
s.genes.mouse<-(s.genes.list$gene_name)
s.genes.mouse

g2m.genes.list<-subset(mouse94, mouse94$human_homolog %in% g2m.genes)
g2m.genes.mouse<-(g2m.genes.list$gene_name)
g2m.genes.mouse

RPMvRPR2 <- CellCycleScoring(RPMvRPR2, s.features = s.genes.mouse, g2m.features = g2m.genes.mouse)
RPMvRPR2@meta.data$Phase<-factor(RPMvRPR2@meta.data$Phase,levels=c("G1","S","G2M"))

DimPlot(RPMvRPR2, reduction = "umap", group.by = "Phase",shuffle=TRUE,label=FALSE,pt.size=.3,cols=c("indianred3","green3","royalblue4"))+NoAxes()


### Look at distribution of labeled tumors by cluster
library(ggplot2)
library(ggpubr)

##### What % of cells occupy what cluster for all types (For Fig. 2m) ####

x<-table(Idents(RPMvRPR2),RPMvRPR2@meta.data$leiden_scVI_1.3)
x
proportions <- as.data.frame(100*prop.table(x, margin = 1))
proportions

colnames(proportions)<-c("Cluster", "Sample", "Frequency")

ggbarplot(proportions, x="Sample", y="Frequency", fill = "Sample", group = "Sample", ylab = "Frequency (percent)", xlab="Phase", palette =c("indianred3","green3","royalblue4"))+ theme_bw()+ facet_wrap(facets = "Cluster", scales="free_y", ncol =4)+ theme(axis.text.y = element_text(size=12)) +rotate_x_text(angle = 45)

# Stacked
p<-ggplot(proportions, aes(fill=Sample, y=Frequency, x=Cluster)) +
  geom_bar(position="stack", stat="identity")

p + scale_fill_manual(values=colors)+ theme_bw()+ theme(axis.text.y = element_text(size=20), axis.text.x=element_text(size=20), axis.title.x =element_text(size=20), axis.title.y = element_text(size=20), legend.text = element_text(size=20), legend.title = element_text(size=20))


############### APPLY VARIOUS SIGNATURES ######################

#######################################################
## Look at archetype sigs as scores (Ext Data Fig 4i)

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


RPMvRPR2<-AddModuleScore(
  object = RPMvRPR2,
  features = list(a_sc_mouse),
  name = 'A_Archetype')

RPMvRPR2<-AddModuleScore(
  object = RPMvRPR2,
  features = list(a2_sc_mouse),
  name = 'A2_Archetype')

RPMvRPR2<-AddModuleScore(
  object = RPMvRPR2,
  features = list(n_sc_mouse),
  name = 'N_Archetype')

RPMvRPR2<-AddModuleScore(
  object = RPMvRPR2,
  features = list(p_sc_mouse),
  name = 'P_Archetype')

RPMvRPR2<-AddModuleScore(
  object = RPMvRPR2,
  features = list(y_sc_mouse),
  name = 'Y_Archetype')


library(viridis)

VlnPlot(RPMvRPR2, features = c("A_Archetype1","N_Archetype1","P_Archetype1"),group.by=c("Genotype"), cols=c("turquoise3","maroon"),alpha=0.2,ncol=4)


#################################################################
## Look at ND1 vs ASCl1 ChIP targets as scores (Ext Data Fig 4j)

chip<-read.csv("Signatures/ASCL1_NEUROD1_POU2F3_ChIP_Targets.csv")
achip<-chip$Borromeo_ASCL1_Targets
nchip<-chip$Borromeo_Oliver_NEUROD1_Targets
pchip<-chip$POU2F3_Targets
pchip<-subset(mouse94, mouse94$human_homolog %in% pchip)
pchip_m<-(pchip$gene_name)
pchip_m

RPMvRPR2<-AddModuleScore(
  object = RPMvRPR2,
  features = list(achip),
  name = 'ASCL1_Targets')

RPMvRPR2<-AddModuleScore(
  object = RPMvRPR2,
  features = list(nchip),
  name = 'NEUROD1_Targets')

RPMvRPR2<-AddModuleScore(
  object = RPMvRPR2,
  features = list(pchip_m),
  name = 'POU2F3_Targets')

# Also add MYC ChIP data
MYC_targets<-read.csv("Signatures/50_Conserved_MYC_Targets.csv")
MYC_targets<-MYC_targets$Gene
MYC_targets

RPMvRPR2<-AddModuleScore(
  object = RPMvRPR2,
  features = list(MYC_targets),
  name = 'MYC_Targets')

## Ext Data Fig 4j ##
VlnPlot(RPMvRPR2, features = c("ASCL1_Targets1","NEUROD1_Targets1","POU2F3_Targets1","MYC_Targets1"),group.by=c("Genotype"), cols=c("turquoise3","maroon"),alpha=.2, ncol=4)


# Violin plots for Ext Data Fig. 4h and t-test method
my_comparisons <- list( c("RPR2", "RPM"))

t<-VlnPlot(RPMvRPR2, features = c("Mycl"), group.by=c("Genotype"), cols=c("turquoise3","maroon"),pt.size=.01, y.max=100)
VlnPlot(RPMvRPR2, features = c("Myc", "Mycl"), group.by=c("Genotype"), cols=c("turquoise3","maroon"),pt.size=.01)
VlnPlot(RPMvRPR2, features = c("Ascl1","Neurod1","Pou2f3"), group.by=c("Genotype"), cols=c("turquoise3","maroon"),pt.size=.01)


t + stat_compare_means(comparisons = my_comparisons, method="t.test", paired=FALSE) # Add pairwise comparisons p-value

# Save rds object
saveRDS(RPMvRPR2,"090824_RPMvRPR2_TBOAllo_Seurat.rds")
RPMvRPR2<-readRDS("090824_RPMvRPR2_TBOAllo_Seurat.rds")

######## Add NE Score for Fig. 2n ##########

### For NE score assignemnt

#### Converting to SCE to Add NE Score and create Violin plots
# Convert seurat object to SCE
library(SingleCellExperiment)
sce<-as.SingleCellExperiment(RPMvRPR2)
library(SummarizedExperiment)


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


###################### Violin plotting expression of NE score (Fig 2n) #####################
library(scater)

sub_cols<-c("turquoise3","maroon")
x<-scater::plotColData(sce, x = "Genotype", y = "NE_spearman", colour_by = "Genotype")+scale_discrete_manual(aesthetics = c("colour", "fill"),values=sub_cols)

#"RPM"="orange","RPMA"="firebrick2","SNL"="dodgerblue4")
# Add stats
fortest<-x+ geom_boxplot(fill=sub_cols, alpha=1/5, position = position_dodge(width = .2),size=0.2,color="black", notch=TRUE, notchwidth=0.3, outlier.shape = 2, outlier.colour=NA)
fortest


# Add p-values comparing groups
# Specify the comparisons you want
my_comparisons <- list( c("RPM", "RPR2"))
my_comparisons
fortest
library(ggplot2) 
library(ggpubr) 
library(tidyverse) 

fortest + stat_compare_means(comparisons = my_comparisons, method="t.test", paired=FALSE) # Add pairwise comparisons p-value
# Add global p-value

fortest + stat_compare_means(comparisons = my_comparisons, method="t.test", paired=FALSE,label = "p.signif") # Add pairwise comparisons p-value
# Add global p-value


# Add NE score to seurat object
ne_score<-sce$NE_spearman
saveRDS(RPMvRPR2,"111424_RPMvRPR2_TBO.rds")

RPMvRPR2@meta.data$NE_spearman<-ne_score
RPMvRPR2@meta.data$NE_spearman

# Plot for Fig 2n
VlnPlot(RPMvRPR2,features = c("NE_spearman"), same.y.lims=FALSE,group.by=c("Genotype"),cols=c("turquoise","turquoise4"),pt.size=.01, ncol=3)
FeaturePlot(RPMvRPR2, features = c("NE_spearman"), pt.size=0.2, reduction='umap',order=TRUE)+scale_color_viridis(option="rocket",direction=-1)& NoAxes()
table(RPMvRPR2@meta.data$leiden_scVI_1.3)

# END #







