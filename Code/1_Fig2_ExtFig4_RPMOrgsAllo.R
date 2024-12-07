## 08/23/24 RPM In vitro to in vivo TBO allograft for Fig 2 ##

library(Seurat)
library(SeuratObject)
library(zellkonverter)
library(devtools)
library(SummarizedExperiment)

setwd("~/Desktop/scRNAseq/")
getwd()

########### Reading in RPM TBO WT and Transformed data only (Extended Data Fig 4a,c) #########
# 09/24/24 AI

# Read in adata object as SCE 
adata<-readH5AD("061824_RPMTBO_OrgsOnly_Fig3c.h5ad")
table(adata$Cre)
adata
table(adata$leiden_scVI_1.2)
table(adata$UnID)
# RPM_Allo  
# 5627

dim(assay(adata,"counts"))
counts(adata)<-assay(adata,"counts")
counts(adata)
length(rownames(adata))
counts(adata)
assay(adata,"norm")
rowData(adata)
length(rowData(adata)$gene_ids)
adata

#Convert SCE to seurat
RPM_Orgs <- CreateSeuratObject(counts = counts(adata), meta.data = as.data.frame(colData(adata)))
RPM_Orgs
# An object of class Seurat 
# 55491 features across 2390 samples within 1 assay 
# Active assay: RNA (55491 features, 0 variable features)

DefaultAssay(RPM_Orgs)

#Create an assay of normalized gene expression
norm_assay <- CreateAssayObject(counts = assay(adata,"norm"))
norm_assay

RPM_Orgs
RPM_Orgs[["norm"]] <- norm_assay

# Set assay norm as default assay for seurat object
DefaultAssay(RPM_Orgs)<-'norm'

# Add embeddings of umap, and X_scVI_1.4 to assay norm
reducedDim(adata, "X_umap")
test<-reducedDim(adata, "X_umap")
colnames(test)<-c("UMAP_1","UMAP_2")
rownames(test)<-colnames(RPM_Orgs)
dim(test)
RPM_Orgs[['umap']] <- CreateDimReducObject(test, key="UMAP_", assay = "RNA")
RPM_Orgs[['umap']]@cell.embeddings


#colors<-c('#F39B7F99','brown1','#ff7f0e', 'turquoise', '#2ca02c', 'deepskyblue4', '#e7298a',  '#0aa6d8',  '#A0522D', '#a1c299', 'gold', '#984ea3', '#3C548899', '#377eb8',  '#d62728','navyblue', '#a1c299', 'purple','gray','maroon','gray20')
colors<-c('deepskyblue4', '#ff7f0e','#2ca02c', 'brown1', 'purple', '#A0522D', 'pink2',  '#a1c299',  '#0aa6d8', 'lightblue','salmon1','lightgreen')
DimPlot(RPM_Orgs,group.by='leiden_scVI_1.2',cols=colors, reduction='umap',label=FALSE,label.size=6)&NoAxes()

table(RPM_Orgs@meta.data$Batch)
# Org_Cre Org_No_Cre 
# 1230       1160 

DefaultAssay(RPM_Orgs)<-"RNA"
# Normalize and log counts for downstream signatures/visualization etc
RPM_Orgs<-NormalizeData(RPM_Orgs)

########################################################################
## Perform Cell Cycle Analysis for Ext Data Fig 4c and visualize ##
########################################################################
# lowercase phase is from python, upper will be from seurat
library(devtools)
library(hciR)
library(hciRdata)

s.genes<-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes

mouse94<-read.csv("Signatures/mouse94.csv")
s.genes.list<-subset(mouse94, mouse94$human_homolog %in% s.genes)
s.genes.mouse<-(s.genes.list$gene_name)
s.genes.mouse

g2m.genes.list<-subset(mouse94, mouse94$human_homolog %in% g2m.genes)
g2m.genes.mouse<-(g2m.genes.list$gene_name)
g2m.genes.mouse

RPM_Orgs <- CellCycleScoring(RPM_Orgs, s.features = s.genes.mouse, g2m.features = g2m.genes.mouse)
RPM_Orgs@meta.data$Phase<-factor(RPM_Orgs@meta.data$Phase,levels=c("G1","S","G2M"))

DimPlot(RPM_Orgs, reduction = "umap", group.by = "Phase",shuffle=TRUE,label=FALSE,pt.size=.4,cols=c("indianred3","green3","royalblue4"))+NoAxes()


############################# Read in scanpy object with organoids and allograft combined for Fig. 2d ################

# Read in adata object as SCE 
adata<-readH5AD("~/Desktop/090624_RPM_WT_Org_Allo_Fig2.h5ad")
table(adata$Cre)
table(adata$leiden_scVI_1.2)
table(adata$UnID)
# RPM_Allo  RPM_Org_Cre WT_Org_NoCre 
# 6594         1230         1160 

?readH5AD
dim(assay(adata,"counts"))
counts(adata)<-assay(adata,"counts")
counts(adata)
length(rownames(adata))
counts(adata)
assay(adata,"norm")
rowData(adata)
length(rowData(adata)$gene_ids)
adata

#Convert SCE to seurat
?CreateSeuratObject
RPM_TBO <- CreateSeuratObject(counts = counts(adata), meta.data = as.data.frame(colData(adata)))
RPM_TBO
# An object of class Seurat 
# 55491 features across 8984 samples within 1 assay 
# Active assay: RNA (55491 features, 0 variable features)

DefaultAssay(RPM_TBO)

#Create an assay of normalized gene expression
norm_assay <- CreateAssayObject(counts = assay(adata,"norm"))
norm_assay

RPM_TBO
RPM_TBO[["norm"]] <- norm_assay

# Set assay norm as default assay for seurat object
DefaultAssay(RPM_TBO)<-'norm'
table(RPM_TBO@meta.data$UnID)
# RPM_Allo  RPM_Org_Cre WT_Org_NoCre 
# 6594         1230         1160 

# Add embeddings of umap, and X_scVI_1.4 to assay norm
adata
library(SingleCellExperiment)

reducedDim(adata, "X_umap")
test<-reducedDim(adata, "X_umap")
colnames(test)<-c("UMAP_1","UMAP_2")
rownames(test)<-colnames(RPM_TBO)
dim(test)
RPM_TBO[['umap']] <- CreateDimReducObject(test, key="UMAP_", assay = "RNA")
RPM_TBO[['umap']]@cell.embeddings

colors<-c('deepskyblue4', '#ff7f0e','#2ca02c', 'brown1', 'purple', '#A0522D', 'pink2',  '#a1c299',  '#0aa6d8', 'lightblue','salmon1','lightgreen')

DimPlot(RPM_TBO,group.by='leiden_scVI_1.2',cols=colors, reduction='umap',label=FALSE,label.size=6)&NoAxes()


table(RPM_TBO@meta.data$UnID)
# RPM_Allo  RPM_Org_Cre WT_Org_NoCre 
# 6594         1230         1160 


RPM_TBO@meta.data$UnID<-factor(RPM_TBO@meta.data$UnID,levels=c("WT_Org_NoCre","RPM_Org_Cre","RPM_Allo"))
unid_cols<-c("orange1","darkorchid4", "turquoise")
DimPlot(RPM_TBO,group.by='UnID',cols=unid_cols,reduction='umap',shuffle=TRUE)

######## Before signatures, normalize and log count data############
DefaultAssay(RPM_TBO)<-'RNA'
RPM_TBO<-NormalizeData(RPM_TBO)
########################################################################
## Perform Cell Cycle Analysis and visualize for Ext Data Fig 4d ##
########################################################################

s.genes<-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes


mouse94<-read.csv("Desktop/Signatures/mouse94.csv")
s.genes.list<-subset(mouse94, mouse94$human_homolog %in% s.genes)
s.genes.mouse<-(s.genes.list$gene_name)
s.genes.mouse

g2m.genes.list<-subset(mouse94, mouse94$human_homolog %in% g2m.genes)
g2m.genes.mouse<-(g2m.genes.list$gene_name)
g2m.genes.mouse

RPM_TBO <- CellCycleScoring(RPM_TBO, s.features = s.genes.mouse, g2m.features = g2m.genes.mouse)
RPM_TBO@meta.data$Phase<-factor(RPM_TBO@meta.data$Phase,levels=c("G1","S","G2M"))

DimPlot(RPM_TBO, reduction = "umap", group.by = "Phase",shuffle=TRUE,label=FALSE,pt.size=.4,cols=c("indianred3","green3","royalblue4"))+NoAxes()

### Look at distribution of labeled tumors by cluster

##### What % of cells occupy what cluster for all types ####
Idents(RPM_TBO)<-'UnID'
Idents(RPM_TBO)

x<-table(Idents(RPM_TBO),RPM_TBO@meta.data$Phase)
x
proportions <- as.data.frame(100*prop.table(x, margin = 1))
proportions

colnames(proportions)<-c("Cluster", "Sample", "Frequency")


ggbarplot(proportions, x="Sample", y="Frequency", fill = "Sample", group = "Sample", ylab = "Frequency (percent)", xlab="Phase", palette =c("indianred3","green3","royalblue4"))+ theme_bw()+ facet_wrap(facets = "Cluster", scales="free_y", ncol =4)+ theme(axis.text.y = element_text(size=12)) +rotate_x_text(angle = 45)

# Stacked
p<-ggplot(proportions, aes(fill=Sample, y=Frequency, x=Cluster)) +
  geom_bar(position="stack", stat="identity")

p + scale_fill_manual(values=c("indianred3","green3","royalblue4"))+ theme_bw()+ theme(axis.text.y = element_text(size=20), axis.text.x=element_text(size=20), axis.title.x =element_text(size=20), axis.title.y = element_text(size=20), legend.text = element_text(size=20), legend.title = element_text(size=20))


########################################################################
######### Add Signature Data ##################
########################################################################
library(viridis)
library(ggplot2)
library(ggpubr)

######## Cell Type signatures ###########
# Using montoro consensus but minus Rpl and Rps genes
ct<-read.csv("Signatures/Montoro_Consensus.csv")
basal_noC<-ct$Basal_noC[1:41]
basal_noC
ne_noC<-ct$Neuroendocrine_noC[1:32]
ne_noC


RPM_TBO<-AddModuleScore(
  object = RPM_TBO,
  features = list(basal_noC),
  name = 'Basal_Consensus')

RPM_TBO<-AddModuleScore(
  object = RPM_TBO,
  features = list(ne_noC),
  name = 'NE_Consensus')


library(viridis)

FeaturePlot(RPM_TBO, features = c("NE_Consensus1"), pt.size=0.2, reduction='umap', order=TRUE)+ scale_color_gradientn(colors=rocket(10, direction=-1)) & NoAxes()
FeaturePlot(RPM_TBO, features = c("Basal_Consensus1"), pt.size=0.2, reduction='umap',order=TRUE)+ scale_color_gradientn(colors=rocket(10, direction=-1)) & NoAxes()


saveRDS(RPM_TBO,"090724_RPM_TBO_InVitrotoAllo_Sigs.rds")
RPM_TBO<-readRDS("090724_RPM_TBO_InVitrotoAllo_Sigs.rds")


##########################################################################################
############### Reading in and analyzing JUST the allograft cells, clustered in scvi/scanpy for Fig. 2e #####
#########################################################################################################
# 09/24/24 AI

# Read in adata object as SCE 
adata<-readH5AD("092824_RPM_Allo_Only_StringentQC_new.h5ad")
table(adata$Cre)
table(adata$leiden_scVI_1.4)
table(adata$UnID)
# RPM_Allo  
# 5627

dim(assay(adata,"counts"))
counts(adata)<-assay(adata,"counts")
counts(adata)
length(rownames(adata))
counts(adata)
assay(adata,"norm")
rowData(adata)
length(rowData(adata)$gene_ids)
adata

#Convert SCE to seurat
RPM_Allo <- CreateSeuratObject(counts = counts(adata), meta.data = as.data.frame(colData(adata)))
RPM_Allo
# An object of class Seurat 
# 55491 features across 6594 samples within 1 assay 
# Active assay: RNA (55491 features, 0 variable features)

DefaultAssay(RPM_Allo)

#Create an assay of normalized gene expression
norm_assay <- CreateAssayObject(counts = assay(adata,"norm"))
norm_assay

RPM_Allo
RPM_Allo[["norm"]] <- norm_assay

# Set assay norm as default assay for seurat object
DefaultAssay(RPM_Allo)<-'norm'

# Add embeddings of umap, and X_scVI_1.4 to assay norm
reducedDim(adata, "X_umap")
test<-reducedDim(adata, "X_umap")
colnames(test)<-c("UMAP_1","UMAP_2")
rownames(test)<-colnames(RPM_Allo)
dim(test)
RPM_Allo[['umap']] <- CreateDimReducObject(test, key="UMAP_", assay = "RNA")
RPM_Allo[['umap']]@cell.embeddings


colors<-c('deepskyblue4', '#ff7f0e','#2ca02c', 'brown1', 'purple', '#A0522D', 'pink2',  '#a1c299',  '#0aa6d8', 'lightblue','salmon1','lightgreen')
DimPlot(RPM_Allo,group.by='leiden_scVI_1.4',cols=colors, reduction='umap',label=FALSE,label.size=6)&NoAxes()

table(RPM_Allo@meta.data$Batch)
# RPM_Allo_New RPM_Allo_Old 
# 3476         2151 


### Normalize and log data before all signatures and for visualization #########
DefaultAssay(RPM_Allo)<-'RNA'
RPM_Allo<-NormalizeData(RPM_Allo)

########################################################################
## Perform Cell Cycle Analysis and visualize ##
########################################################################
# lowercase phase is from python, upper will be from seurat

s.genes<-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes

library(devtools)
library(hciR)
library(hciRdata)

mouse94<-read.csv("Desktop/Signatures/mouse94.csv")
s.genes.list<-subset(mouse94, mouse94$human_homolog %in% s.genes)
s.genes.mouse<-(s.genes.list$gene_name)
s.genes.mouse

g2m.genes.list<-subset(mouse94, mouse94$human_homolog %in% g2m.genes)
g2m.genes.mouse<-(g2m.genes.list$gene_name)
g2m.genes.mouse

RPM_Allo <- CellCycleScoring(RPM_Allo, s.features = s.genes.mouse, g2m.features = g2m.genes.mouse)
RPM_Allo@meta.data$Phase<-factor(RPM_Allo@meta.data$Phase,levels=c("G1","S","G2M"))

DimPlot(RPM_Allo, reduction = "umap", group.by = "Phase",shuffle=TRUE,label=FALSE,pt.size=.4,cols=c("indianred3","green3","royalblue4"))+NoAxes()

### Look at distribution of labeled tumors by cluster ###

##### What % of cells occupy what cluster for all types ####
Idents(RPM_Allo)<-'leiden_scVI_1.4'
Idents(RPM_Allo)

x<-table(Idents(RPM_Allo),RPM_Allo@meta.data$Phase)
x
proportions <- as.data.frame(100*prop.table(x, margin = 1))
proportions

colnames(proportions)<-c("Cluster", "Sample", "Frequency")


ggbarplot(proportions, x="Sample", y="Frequency", fill = "Sample", group = "Sample", ylab = "Frequency (percent)", xlab="Phase", palette =c("indianred3","green3","royalblue4"))+ theme_bw()+ facet_wrap(facets = "Cluster", scales="free_y", ncol =4)+ theme(axis.text.y = element_text(size=12)) +rotate_x_text(angle = 45)

# Stacked
p<-ggplot(proportions, aes(fill=Sample, y=Frequency, x=Cluster)) +
  geom_bar(position="stack", stat="identity")

p + scale_fill_manual(values=c("indianred3","green3","royalblue4"))+ theme_bw()+ theme(axis.text.y = element_text(size=20), axis.text.x=element_text(size=20), axis.title.x =element_text(size=20), axis.title.y = element_text(size=20), legend.text = element_text(size=20), legend.title = element_text(size=20))



####### Visualize A, N, P, P63 for Fig. 2e ##############

VlnPlot(RPM_Allo, features=c("Ascl1","Neurod1","Pou2f3","Trp63"), group.by=c("leiden_scVI_1.4"), cols=colors, alpha=0.7,ncol=4)

########################################################################
# Assign states for Fig. 2f
########################################################################

RPM_Allo@meta.data$Pheno<- ifelse(RPM_Allo@meta.data$leiden_scVI_1.4 %in% c("9"), "Tuft",
                                  ifelse(RPM_Allo@meta.data$leiden_scVI_1.4 %in% c("8"), "Basal",
                                         ifelse(RPM_Allo@meta.data$leiden_scVI_1.4 %in% c("1","3"), "NE/Neuronal",
                                                ifelse(RPM_Allo@meta.data$leiden_scVI_1.4 %in% c("0","2","6","7"), "NE",
                                                       ifelse(RPM_Allo@meta.data$leiden_scVI_1.4 %in% c("4","5"), "Neuronal", "Other")))))




table(RPM_Allo@meta.data$Pheno)
# NE NE/Neuronal    Neuronal        Tuft       Basal 
# 2938        1318        1034          12         325 


RPM_Allo@meta.data$Pheno<-factor(RPM_Allo@meta.data$Pheno, levels=c("NE","NE/Neuronal","Neuronal","Tuft","Basal"))
pheno_col<-c("brown2","darkorchid4","dodgerblue","orange","turquoise4")

DimPlot(RPM_Allo, group.by=c("Pheno"), cols=c("brown2","darkorchid4","dodgerblue","orange","turquoise4"), shuffle=TRUE, pt.size=0.6)+NoAxes()


######### Add Signature Data ##################
########################################################################
######### Look at other signatures from BioRxiv and cell types of OE
###########
## Look at ND1 vs ASCl1 vs POU2F3 ChIP targets as scores (Fig 2h)

chip<-read.csv("Signatures/ASCL1_NEUROD1_POU2F3_ChIP_Targets.csv")
achip<-chip$Borromeo_ASCL1_Targets
nchip<-chip$Borromeo_Oliver_NEUROD1_Targets
pchip<-chip$POU2F3_Targets
pchip<-subset(mouse94, mouse94$human_homolog %in% pchip)
pchip_mouse<-(pchip$gene_name)
pchip_mouse

RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(achip),
  name = 'ASCL1_Targets')

RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(nchip),
  name = 'NEUROD1_Targets')

RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(pchip_mouse),
  name = 'POU2F3_Targets')

FeaturePlot(RPM_Allo, features = c("ASCL1_Targets1"), pt.size=0.2, reduction='umap', order=TRUE)+ scale_color_gradientn(colors=rocket(10, direction=-1)) & NoAxes()
FeaturePlot(RPM_Allo, features = c("NEUROD1_Targets1"), pt.size=0.2, reduction='umap',order=TRUE)+ scale_color_gradientn(colors=rocket(10, direction=-1)) & NoAxes()
FeaturePlot(RPM_Allo, features = c("POU2F3_Targets1"), pt.size=0.2, reduction='umap',order=TRUE)+ scale_color_gradientn(colors=rocket(10, direction=-1)) & NoAxes()
VlnPlot(RPM_Allo,features = c("ASCL1_Targets1","NEUROD1_Targets1","POU2F3_Targets1"), same.y.lims=FALSE,group.by=c("Pheno"),cols=pheno_col,pt.size=.01,alpha=0.2)


######## Cell Type signatures (Fig 2i) ###########
# Using montoro consensus but minus Rpl and Rps genes
ct<-read.csv("Signatures/Montoro_Consensus.csv")
basal_noC<-ct$Basal_noC[1:41]
basal_noC
iono_noC<-ct$Ionocyte_noC[1:19]
iono_noC
ne_noC<-ct$Neuroendocrine_noC[1:32]
ne_noC
tuft_noC<-ct$Tuft_noC[1:60]
tuft_noC

RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(tuft_noC),
  name = 'Tuft_Consensus')

RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(basal_noC),
  name = 'Basal_Consensus')

RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(ne_noC),
  name = 'NE_Consensus')

RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(iono_noC),
  name = 'Ionocyte_Consensus')

VlnPlot(RPM_Allo,features = c("NE_Consensus1","Tuft_Consensus1","Basal_Consensus1"), same.y.lims=FALSE,group.by=c("Pheno"),cols=pheno_col,pt.size=.01,alpha=0.2)


## Look at archetype sigs as scores (Ext Data Fig 4f) ##

arch<-read.csv("Signatures/Archetype_Sigs_Maddox.csv")
a<-arch$SCLC.A
write.csv(a2,"Archetype_A2_signature_human.csv")
a2<-arch$SCLC.A2
a2
n<-arch$SCLC.N
p<-arch$SCLC.P
y<-arch$SCLC.Y

a<-a[1:987]
a_sc<-subset(mouse94, mouse94$human_homolog %in% a)
a_sc_mouse<-(a_sc$gene_name)
a_sc_mouse

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


RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(a_sc_mouse),
  name = 'A_Archetype')

RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(a2_sc_mouse),
  name = 'A2_Archetype')

RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(n_sc_mouse),
  name = 'N_Archetype')

RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(p_sc_mouse),
  name = 'P_Archetype')

RPM_Allo<-AddModuleScore(
  object = RPM_Allo,
  features = list(y_sc_mouse),
  name = 'Y_Archetype')

library(viridis)
RPM_Allo@meta.data$leiden_scVI_1.4
VlnPlot(RPM_Allo,features = c("A_Archetype1","N_Archetype1","P_Archetype1"),group.by=c("leiden_scVI_1.4"),cols=colors,pt.size=.01, alpha=0.1,ncol=4)


## Convert to generate NE score for Fig. 2g #
######## Add NE Score ##########

#### Converting to SCE to Add NE Score and create Violin plots
# Convert seurat object to SCE
library(SingleCellExperiment)
sce<-as.SingleCellExperiment(RPM_Allo)
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


###################### Violin plotting expression of interest #####################
# For example, NE score by phenotype #
library(scater)
colors
x<-scater::plotColData(sce, x = "Pheno", y = "NE_spearman", colour_by = "Pheno")+scale_discrete_manual(aesthetics = c("colour", "fill"),values=pheno_col)

# Add stats
fortest<-x+ geom_boxplot(fill=pheno_col, alpha=1/5, position = position_dodge(width = .2),size=0.2,color="black", notch=TRUE, notchwidth=0.3, outlier.shape = 2, outlier.colour=NA)
fortest


############# # For example, Basal Cell score by phenotype ##################
x<-scater::plotColData(sce, x = "Pheno", y = "Basal_Consensus1", colour_by = "Pheno")+scale_discrete_manual(aesthetics = c("colour", "fill"),values=pheno_col)
# Add stats
fortest<-x+ geom_boxplot(fill=pheno_col, alpha=1/5, position = position_dodge(width = .2),size=0.2,color="black", notch=TRUE, notchwidth=0.3, outlier.shape = 2, outlier.colour=NA)
fortest+ylim(-0.1, 1.2)
#+scale_y_break(c(3.5,8))


# Add p-values comparing groups
# Specify the comparisons you want
my_comparisons <- list( c("NE", "Tuft"),c("NE/Neuronal","Tuft"),c("Neuronal","Tuft"), c("Basal","Tuft"))
my_comparisons
fortest
library(ggplot2) 
library(ggpubr) 
library(tidyverse) 


fortest + stat_compare_means(comparisons = my_comparisons, method="t.test", paired=FALSE) # Add pairwise comparisons p-value
# Add global p-value

fortest + stat_compare_means(comparisons = my_comparisons, method="t.test", paired=FALSE,label = "p.signif") # Add pairwise comparisons p-value
# Add global p-value


# Add NE score to seurat object for Fig. 2g
ne_score<-sce$NE_spearman


RPM_Allo@meta.data$NE_spearman<-ne_score
RPM_Allo@meta.data$NE_spearman
RPM_Allo_Fig2e_meta<-RPM_Allo@meta.data

# Export meta data for statistical tests, etc. 
write.csv(RPM_Allo_Fig2e_meta,"RPM_Allo_Fig2e_meta.csv")

# Fig 2g plots #
FeaturePlot(RPM_Allo, features = c("NE_spearman"), pt.size=0.2, reduction='umap',order=TRUE)+scale_color_viridis(option="rocket",direction=-1)& NoAxes()
VlnPlot(RPM_Allo,features = c("NE_spearman"), same.y.lims=FALSE,group.by=c("Pheno"),cols=pheno_col,pt.size=.01)

# Save data #
# saveRDS(RPM_Allo,"110224_RPM_AlloOnly_Fig2_wFA.rds")
# RPM_Allo<-readRDS("092824_RPM_AlloOnly_Fig2_wFA.rds")

saveRDS(RPM_Allo,"092424_RPM_AlloOnly_Fig2.rds")
RPM_Allo<-readRDS("092424_RPM_AlloOnly_Fig2.rds")


table(RPM_Allo@meta.data$leiden_scVI_1.4)


## END ##
