# RPM + RPMA TBO Allograft +/- Cre in vitro to in vivo
# Abbie Ireland 05/31/2024
# Fig. 3-4, 6 and Extended Data Fig. 5-6, 9

#### Load packages##

#load necessary packages
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(SingleCellExperiment)
  library(ggplot2)
  library(cowplot)
  library(Seurat)
  library(SeuratObject)
  library(CellTagR)
})

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")

library("devtools")
devtools::install_github("morris-lab/CellTagR")

setwd("/Volumes/All_Staff/Abbie/scRNAseq/RPM_TBO_Allo2_CellTagAnalysis")

setwd("/Users/abbieireland/Library/Mobile Documents/com~apple~CloudDocs/")

########################################################################################################################
#### Let us test if there are shared clones in the In vitro minus Cre sample (not Cellplexed) and the In Vitro + Cre RPM sample AND the allograft
########################################################################################################################
setwd("/Users/abbieireland/Library/Mobile Documents/com~apple~CloudDocs/Documents/RPM_TBO_Allo2_CellTagAnalysis/0531_InVitroCretoInVivo")

#Set up the CellTag Object
bam.test.obj <- CellTagObject(object.name = "bam.cell.tag.obj", fastq.bam.directory = "./bams")
bam.test.obj
#Sample-1 is In Vitro + CMV Cre
#Sample-2 is RPM Allo In Vivo (Fastq 12/31)

# Extract the CellTag information
bam.test.obj <- CellTagExtraction(bam.test.obj, celltag.version = "v1")
bam.test.obj

# Check the bam file result
head(bam.test.obj@bam.parse.rslt[["v1"]])

# Generate the sparse count matrix
Barcode.Aggregate(list("./barcodes/0531_RPM_InVitro_Cre_barcodes.tsv", "./barcodes/0531_RPM_TBO_barcodes.tsv"), "./barcodes/all_barcodes.tsv")


# Generate the sparse count matrix
bam.test.obj <- CellTagMatrixCount(celltag.obj = bam.test.obj, barcodes.file ="./barcodes/all_barcodes.tsv")

bam.test.obj
# Object name:  bam.cell.tag.obj 
# Raw CellTag Counts =  1785 
# Raw Number of Cells with CellTag =  7418 
# Collapsed CellTag Counts =  0 
# Whitelisted CellTag Counts =  0 
# Whitelisted Number of Cells with CellTag =  0 

# Check the dimension of the raw count matrix
dim(bam.test.obj@raw.count)
tail(colnames(bam.test.obj@raw.count))
head(bam.test.obj@raw.count)
head(bam.test.obj@bam.parse.rslt)



##The collapsing steps are actually optional.. try with first ###
# Generating the collapsing file
bam.test.obj.collapsed <- CellTagDataForCollapsing(celltag.obj = bam.test.obj, output.file = "./collapsing.txt")


#Run starcode in terminal per sample
# starcode/starcode -s --print-clusters ~/Desktop/collapsing.txt > ~/Desktop/collapsing_result.txt

collapsed.rslt.dir <- "./collapsing_results/"

# Recount and generate collapsed matrix
bam.test.obj.collapsed <- CellTagDataPostCollapsing(celltag.obj = bam.test.obj.collapsed, collapsed.rslt.file = list.files(collapsed.rslt.dir, full.names = T))
bam.test.obj.collapsed
# Object name:  bam.cell.tag.obj 
# Raw CellTag Counts =  1785 
# Raw Number of Cells with CellTag =  7418 
# Collapsed CellTag Counts =  1447 
# Whitelisted CellTag Counts =  0 
# Whitelisted Number of Cells with CellTag =  0 

# Check the dimension of this collapsed count.
head(bam.test.obj@collapsed.count)
bam.test.obj.collapsed@collapsed.count
bam.test.obj.collapsed@raw.count


####PICK UP HERE W/ or W/O Collapse #######
# Calling binarization
bam.test.obj.collapsed <- SingleCellDataBinarization(bam.test.obj.collapsed, 2)


# Read the RDS file and get the object
dt.mtx.whitelist.path <- "./v1_whitelist_ours_082023_75percentile.csv"
bam.test.obj.collapsed <- SingleCellDataWhitelist(bam.test.obj.collapsed, dt.mtx.whitelist.path)
bam.test.obj.collapsed

# Object name:  bam.cell.tag.obj 
# Raw CellTag Counts =  1785 
# Raw Number of Cells with CellTag =  7418 
# Collapsed CellTag Counts =  1447 
# Whitelisted CellTag Counts =  1427 
# Whitelisted Number of Cells with CellTag =  4840 

MetricPlots(bam.test.obj.collapsed)
# Average:  1.542149 
# Frequency:  5.230554 

bam.test.obj.collapsed <- MetricBasedFiltering(bam.test.obj.collapsed, 20, comparison = "less")
bam.test.obj.collapsed <- MetricBasedFiltering(bam.test.obj.collapsed, 2, comparison = "greater")

MetricPlots(bam.test.obj.collapsed)

# Average:  2.705283 
# Frequency:  4.0911
bam.test.obj.collapsed
bam.test.obj.collapsed <- JaccardAnalysis(bam.test.obj.collapsed, fast = T)


# Call clones
bam.test.obj.collapsed <- CloneCalling(celltag.obj = bam.test.obj.collapsed, correlation.cutoff=0.8)

# Check them out!!
bam.test.obj.collapsed@clone.composition[["v1"]]
bam.test.obj.collapsed@clone.size.info[["v1"]]

comp<-bam.test.obj.collapsed@clone.composition[["v1"]]
size<-bam.test.obj.collapsed@clone.size.info[["v1"]]
tail(bam.test.obj.collapsed@whitelisted.count)

write.csv(comp,"./collapse_0.8_clonecomp.csv")
write.csv(size,"./collapse_0.8_clonesize.csv")


## Read in comp data to calculate quickly comp/clone 

comp_clone<-read.csv("collapse_0.8_clonecomp.csv")
head(comp_clone)
table(comp_clone$clone.id,comp_clone$Sample_ID)

table(comp_clone$Sample_ID)
# Allograft InVitro_Cre 
# 1488         516 

df<-table(comp_clone$clone.id, comp_clone$Sample_ID)
df
write.csv(df, "compperclonebySample.csv")


## Make an alluvial plot to track bottlenecks/changes in celltag populations over each sample
library(ggalluvial)
library(tidyverse)

# set.seed(0)
# data_bar <- data.frame(
#   stringsAsFactors = F,
#   Sample = rep(c("A", "B"), each = 10),
#   Percentage = runif(20),
#   Taxon = rep(1:10, by = 2)
# )
# data_bar
# data_bar %>%
#   group_by(Sample) %>%
#   mutate(perc = Percentage* 100/sum(Percentage)) %>%
#   ggplot(aes(y = perc, x = Sample, fill = as.character(Taxon))) +
#   geom_flow(aes(alluvium = Taxon), alpha= .5, color = "black",
#             curve_type = "linear", 
#             width = .1) +
#   geom_col(width = .5, color = "black") +
#   scale_y_continuous(NULL, expand = c(0,0)) +
#   cowplot::theme_minimal_hgrid() +
#   theme(legend.position="right")

# db<-read.csv("For_Connect_Bargraph_Plot_CloneTrack.csv")
# db<-data.frame(db)
# data_bar<-db
# data_bar
# table(data_bar$Sample)
# data_bar$Sample<-factor(data_bar$Sample, c("InVitro_NoCre","InVitro_Cre","Allograft"))
# 
# # So, sample should be... In Vitro No Cre, In Vitro Cre, Allograft
# # Percentage should be the percentage of the individual clone represented per sample... 
# # So if Clone 1 has 20% in Sample A, 2% in Sample B and 78% in Sample C they should add up to 100
# 
# # Taxons should be each clone.
# 
# ## first need to calcualte the actual percentage by group/Sample
# data_bar %>%
#   group_by(Sample) %>%
#   mutate(perc = Percentage* 100/sum(Percentage)) %>%
#   ggplot(aes(y = perc, x = Sample, fill = as.character(Clone.ID))) +
#   geom_flow(aes(alluvium = Clone.ID), alpha= .5, color = "black",
#             curve_type = "linear", 
#             width = .1) +
#   geom_col(width = .5, color = "black") +
#   scale_y_continuous(NULL, expand = c(0,0)) +
#   cowplot::theme_minimal_hgrid() +
#   theme(legend.position="none")
# 
# 
# data_bar %>%
#   ggplot(aes(y = Percentage, x = Sample, fill = as.character(Clone.ID))) +
#   geom_flow(aes(alluvium = Clone.ID), alpha= .5, color = "white",
#             curve_type = "linear", 
#             width = .5) +
#   geom_col(width = .5, color = "white") +
#   scale_y_continuous(NULL, expand = c(0,0)) +
#   cowplot::theme_minimal_hgrid() +
#   theme(panel.grid.major = element_blank(), 
#         axis.text.y = element_blank(), 
#         axis.ticks.y = element_blank(), legend.position="none")



### Showing All clones here 
db3<-read.csv("percent_clone_persample.csv")
db3<-data.frame(db3)
data_bar3<-db3
data_bar3
table(data_bar3$Sample)
data_bar3$Sample<-factor(data_bar3$Sample, c("InVitro_Cre","Allograft"))


## first need to calculate the actual percentage by group/Sample
data_bar3 %>%
  group_by(Sample) %>%
  ggplot(aes(y = Percentage, x = Sample, fill = as.character(Clone_ID))) +
  geom_flow(aes(alluvium = Clone_ID), alpha= .5, color = "black",
            curve_type = "sigmoid", 
            width = .2) +
  geom_col(width = .5, color = "black") +
  scale_y_continuous(NULL, expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() +
  theme(legend.position="right")

## Now plot only shared clones
db<-read.csv("clone_per_sample_sharedonly.csv")
db<-data.frame(db)
data_bar<-db
head(data_bar)
table(data_bar$Sample)
data_bar$Sample<-factor(data_bar$Sample, c("InVitro_Cre","Allograft"))

#"xspline", linear", "cubic", "quintic", "sine", "arctangent", and "sigmoid". 
## first need to calculate the actual percentage by group/Sample
data_bar %>%
  group_by(Sample) %>%
  ggplot(aes(y = Percentage, x = Sample, fill = as.character(Clone_ID))) +
  geom_flow(aes(alluvium = Clone_ID), alpha= .5, color = "black",
            curve_type = "sigmoid", 
            width = .5) +
  geom_col(width = .5, color = "black") +
  scale_y_continuous(NULL, expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() +
  theme(legend.position="right")



########################################################################################
########################################################################################
################# RPMA TBO CellTag analysis 032024 ########################
########################################################################################################################
#### Let us test if there are shared clones in the In vitro minus Cre sample (not Cellplexed) and the In Vitro + Cre RPM sample AND the allograft
########################################################################################################################
setwd("/Users/abbieireland/Library/Mobile Documents/com~apple~CloudDocs/Documents/RPMA_TBO_scRNAseq_031824/CellTagAnalysis/InVitroCretoAllo")

library("devtools")
devtools::install_github("morris-lab/CellTagR")
library("CellTagR")


#Set up the CellTag Object
bam.test.obj <- CellTagObject(object.name = "bam.cell.tag.obj", fastq.bam.directory = "./bams")
bam.test.obj
#Sample-1 is In Vitro + Cre 
#Sample-2 is RPMA Allo In Vivo

# Extract the CellTag information
bam.test.obj <- CellTagExtraction(bam.test.obj, celltag.version = "v1")
bam.test.obj

# Check the bam file result
head(bam.test.obj@bam.parse.rslt[["v1"]])

# Generate the sparse count matrix
Barcode.Aggregate(list("./barcodes/1-RPMA_TBO_CMV_barcodes.tsv", "./barcodes/2-RPMA_Allo_barcodes.tsv"), "./barcodes/all_barcodes.tsv")


# Generate the sparse count matrix
bam.test.obj <- CellTagMatrixCount(celltag.obj = bam.test.obj, barcodes.file ="./barcodes/all_barcodes.tsv")
bam.test.obj

# Object name:  bam.cell.tag.obj 
# Raw CellTag Counts =  4809 
# Raw Number of Cells with CellTag =  23689 
# Collapsed CellTag Counts =  0 
# Whitelisted CellTag Counts =  0 
# Whitelisted Number of Cells with CellTag =  0 

# Check the dimension of the raw count matrix
dim(bam.test.obj@raw.count)
tail(colnames(bam.test.obj@raw.count))
head(bam.test.obj@raw.count)
head(bam.test.obj@bam.parse.rslt)



##The collapsing steps are actually optional.. try with first ###
# Generating the collapsing file
bam.test.obj.collapsed <- CellTagDataForCollapsing(celltag.obj = bam.test.obj, output.file = "collapsing.txt")


#Run starcode in terminal per sample
# starcode/starcode -s --print-clusters ~/Desktop/collapsing.txt > ~/Desktop/collapsing_result.txt

collapsed.rslt.dir <- "collapsing_results"

# Recount and generate collapsed matrix
bam.test.obj.collapsed <- CellTagDataPostCollapsing(celltag.obj = bam.test.obj.collapsed, collapsed.rslt.file = list.files(collapsed.rslt.dir, full.names = T))
bam.test.obj.collapsed
# Object name:  bam.cell.tag.obj 
# Raw CellTag Counts =  4809 
# Raw Number of Cells with CellTag =  23689 
# Collapsed CellTag Counts =  4188 
# Whitelisted CellTag Counts =  0 
# Whitelisted Number of Cells with CellTag =  0 

# Check the dimension of this collapsed count.
head(bam.test.obj@collapsed.count)
bam.test.obj.collapsed@collapsed.count
bam.test.obj.collapsed@raw.count


####PICK UP HERE W/ or W/O Collapse #######
# Calling binarization
bam.test.obj.collapsed <- SingleCellDataBinarization(bam.test.obj.collapsed, 2)
bam.test.obj.collapsed

# Read the RDS file and get the object
dt.mtx.whitelist.path <- "./v1_whitelist_ours_082023_75percentile.csv"
bam.test.obj.collapsed <- SingleCellDataWhitelist(bam.test.obj.collapsed, dt.mtx.whitelist.path)
bam.test.obj.collapsed

# Object name:  bam.cell.tag.obj 
# Raw CellTag Counts =  4809 
# Raw Number of Cells with CellTag =  23689 
# Collapsed CellTag Counts =  4188 
# Whitelisted CellTag Counts =  4090 
# Whitelisted Number of Cells with CellTag =  15167 

MetricPlots(bam.test.obj.collapsed)
# Average:  1.786115 
# Frequency:  6.623472 

bam.test.obj.collapsed <- MetricBasedFiltering(bam.test.obj.collapsed, 20, comparison = "less")
bam.test.obj.collapsed <- MetricBasedFiltering(bam.test.obj.collapsed, 2, comparison = "greater")

MetricPlots(bam.test.obj.collapsed)

# W collapse 75 percentile of celltag whitelist
# Average:  3.480637 
# Frequency:  5.559658 

bam.test.obj.collapsed
bam.test.obj.collapsed <- JaccardAnalysis(bam.test.obj.collapsed, fast = T)


# Call clones
bam.test.obj.collapsed <- CloneCalling(celltag.obj = bam.test.obj.collapsed, correlation.cutoff=0.8)

# Check them out!!
bam.test.obj.collapsed@clone.composition[["v1"]]
bam.test.obj.collapsed@clone.size.info[["v1"]]

comp<-bam.test.obj.collapsed@clone.composition[["v1"]]
size<-bam.test.obj.collapsed@clone.size.info[["v1"]]
tail(bam.test.obj.collapsed@whitelisted.count)

write.csv(comp,"./collapse_0.8_clonecomp.csv")
write.csv(size,"./collapse_0.8_clonesize.csv")




### Read in comp data to calculate quickly comp/clone ###

getwd()
comp_clone<-read.csv("collapse_0.8_clonecomp.csv")

table(comp_clone$Sample)
# Allograft InVitro_Cre 
# 3916        2321 


df<-table(comp_clone$clone.id, comp_clone$Sample)
df
write.csv(df, "RPMA_compperclonebySample.csv")


###### Make alluvial plots to study dynamics over time ###########
db5<-read.csv("percent_clone_persample_all.csv")
db5<-data.frame(db5)
data_bar5<-db5
data_bar5
table(data_bar5$Sample)
data_bar5$Sample<-factor(data_bar5$Sample, c("InVitro_Cre","Allograft"))

## first need to calculate the actual percentage by group/Sample
data_bar5 %>%
  group_by(Sample) %>%
  ggplot(aes(y = Percentage, x = Sample, fill = as.character(Clone_ID))) +
  geom_flow(aes(alluvium = Clone_ID), alpha= .5, color = "black",
            curve_type = "sigmoid", 
            width = .1) +
  geom_col(width = .5, color = "black") +
  scale_y_continuous(NULL, expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() +
  theme(legend.position="none")


# Graph of only shared clones per sample
db5<-read.csv("percentclone_by_sample_shared.csv")
db5<-data.frame(db5)
data_bar5<-db5
data_bar5
table(data_bar5$Sample)
data_bar5$Sample<-factor(data_bar5$Sample, c("InVitro_Cre","Allograft"))

## first need to calculate the actual percentage by group/Sample
data_bar5 %>%
  group_by(Sample) %>%
  ggplot(aes(y = Percentage, x = Sample, fill = as.character(Clone_ID))) +
  geom_flow(aes(alluvium = Clone_ID), alpha= .5, color = "black",
            curve_type = "sigmoid", 
            width = .1) +
  geom_col(width = .5, color = "black") +
  scale_y_continuous(NULL, expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() +
  theme(legend.position="right")


######################################################################
#####################################################################
####################################################################
###################################################################


setwd("/Users/abbieireland/Library/Mobile Documents/com~apple~CloudDocs/Documents/RPM_TBO_Allo2_CellTagAnalysis")


######## MOVING INTO SEURAT NOW ############

#### Pull in adata object from Scanpy 
# scRNA seq RPMA TBO allograft plus RPM TBO allograft(s) # 
# Bringing in from SCANPY analysis.


# Load necessary packages
#load necessary packages
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(SingleCellExperiment)
  library(ggplot2)
  library(cowplot)
  library(Seurat)
  library(SeuratObject)
  library(viridis)
  library(ggpubr)
})

############## Seurat conversion and analysis using xenograft leiden clustering from scanpy ##########

########################################################################################################
########################################################################################################
######### Now, bring into Seurat, from my anndata object processed according to onb paper methods ########################################################################################################
########################################################################################################
################### START WITH ALLOGRAFT ONLY FROM ORIGINAL CLUSTERING ######################################
##############################################################################################################

library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(zellkonverter)


setwd("/Users/abbieireland/Library/Mobile Documents/com~apple~CloudDocs/Documents/RPM_TBO_Allo2_CellTagAnalysis")

# Read in adata object as SCE
adata<-readH5AD("040924_adata_RPM_RPMA_TBO_newandold_scVI_2.h5ad")
?readH5AD
table(adata$Cre)
table(adata$leiden_scVI_1.2)
table(adata$UnID)
# RPMA_Allo RPM_Allo_New RPM_Allo_Old 
# 10263         3660         2190 


table(adata$Cre)
table(adata$leiden_scVI_1.2)
table(adata$UnID)
# RPMA_Allo RPM_Allo_New RPM_Allo_Old 
# 11554         4414         2753 

adata


dim(assay(adata,"counts"))
counts(adata)<-assay(adata,"counts")
counts(adata)
length(rownames(adata))
counts(adata)
assay(adata,"norm")
rowData(adata)
length(rowData(adata)$gene_ids)

?CreateSeuratObject

#Convert SCE to seurat
TBO_seurat <- CreateSeuratObject(counts = counts(adata), meta.data = as.data.frame(colData(adata)))
TBO_seurat
# An object of class Seurat 
# 55800 features across 18721 samples within 1 assay 
# Active assay: RNA (55800 features, 0 variable features)
DefaultAssay(TBO_seurat)

#Create an assay of normalized gene expression
norm_assay <- CreateAssayObject(counts = assay(adata,"norm"))
norm_assay

TBO_seurat
TBO_seurat[["norm"]] <- norm_assay

# Set assay norm as default assay for seurat object
DefaultAssay(TBO_seurat)<-'norm'
table(TBO_seurat@meta.data$UnID)
# RPMA_Allo RPM_Allo_New RPM_Allo_Old 
# 11554         4414         2753 

# Add embeddings of umap, and X_scVI_1.4 to assay norm
adata
reducedDim(adata, "X_umap")
test<-reducedDim(adata, "X_umap")
colnames(test)<-c("UMAP_1","UMAP_2")
rownames(test)<-colnames(TBO_seurat)
dim(test)
TBO_seurat[['umap']] <- CreateDimReducObject(test, key="UMAP_", assay = "RNA")
TBO_seurat[['umap']]@cell.embeddings


# Fig. 3f #
#colors<-c('#F39B7F99','brown1','#ff7f0e', 'turquoise', '#2ca02c', 'deepskyblue4', '#e7298a',  '#0aa6d8',  '#A0522D', '#a1c299', 'gold', '#984ea3', '#3C548899', '#377eb8',  '#d62728','navyblue', '#a1c299', 'purple','gray','maroon','gray20')
colors<-c('deepskyblue4', '#ff7f0e','#2ca02c', 'brown1', 'purple', '#A0522D', 'pink2',  '#a1c299',  '#4ba4d3', '#b5d7e4','#ef9171','#a6eb9a')

DimPlot(TBO_seurat,group.by='leiden_scVI_1.2',cols=colors, reduction='umap',label=FALSE,label.size=6)&NoAxes()


# Begin assigning phenotypes/cell states

TBO_seurat@meta.data
TBO_seurat@meta.data$PhenoRPM<- ifelse(rownames(TBO_seurat@meta.data) %in% tuft, "Tuft",
                                    ifelse(rownames(TBO_seurat@meta.data) %in% ne, "NE",
                                           ifelse(rownames(TBO_seurat@meta.data) %in% n, "Neuronal",
                                                  ifelse(rownames(TBO_seurat@meta.data) %in% ne.n, "NE/Neuronal",
                                                         ifelse(rownames(TBO_seurat@meta.data) %in% basal, "Basal",
                                                                ifelse(rownames(TBO_seurat@meta.data) %in% TripleNeg, "Triple-Neg","RPMA"))))))
table(TBO_seurat@meta.data$PhenoRPM)
# NE NE/Neuronal    Neuronal        Tuft  Triple-Neg       Basal        RPMA 
# 1785        2468         985          15         168         307       10385 

TBO_seurat@meta.data$PhenoRPM<-factor(TBO_seurat@meta.data$PhenoRPM, levels=c("NE","NE/Neuronal","Neuronal","Tuft","Triple-Neg","Basal","RPMA"))
pheno_col<-c("brown2","darkorchid4","dodgerblue","orange","turquoise4","turquoise", rgb(.5,.5,.5,0.01 ))

DimPlot(TBO_seurat,group.by='PhenoRPM',cols=pheno_col, reduction='umap',label=FALSE,label.size=6, split.by="Genotype")&NoAxes()


# Assign phenotype for all tumors cells as in Fig. 3h #
TBO_seurat@meta.data$Pheno<- ifelse(TBO_seurat@meta.data$leiden_scVI_1.2 %in% c("11"), "Tuft",
                                       ifelse(TBO_seurat@meta.data$leiden_scVI_1.2 %in% c("3"), "NE",
                                              ifelse(TBO_seurat@meta.data$leiden_scVI_1.2 %in% c("4","9"), "Neuronal",
                                                     ifelse(TBO_seurat@meta.data$leiden_scVI_1.2 %in% c("2"), "NE/Neuronal",
                                                            ifelse(TBO_seurat@meta.data$leiden_scVI_1.2 %in% c("10"), "Basal",
                                                                   ifelse(TBO_seurat@meta.data$leiden_scVI_1.2 %in% c("0","1","6","5","7","8"), "Triple-Neg","NA"))))))




table(TBO_seurat@meta.data$Pheno)
TBO_seurat@meta.data$Pheno<-factor(TBO_seurat@meta.data$Pheno, levels=c("NE","NE/Neuronal","Neuronal","Tuft","Triple-Neg","Basal"))
pheno_col<-c("brown2","darkorchid4","dodgerblue","orange","turquoise4","turquoise")

DimPlot(TBO_seurat,group.by='Pheno',cols=pheno_col, reduction='umap',label=FALSE,label.size=6) & NoAxes()

#Basal- 10
#NE- 3
#N/NE- 2
#N- 9,4
#P - 11
#Triple-Neg/Mixed-0,1,6,5,7,8

DefaultAssay(TBO_seurat)<-'norm'

# Add CellTag annotations (starting by adding info for all clones with >10 cells total preQC)
write.csv(TBO_seurat@meta.data, "seurat_umap_metadata_RPM_RPMA_TBO_Allo.csv")


md<-read.csv("060124_seurat_umap_metadata_RPM_RPMA_TBO_Allo_wClones.csv")
md$CellTag_ID
md$CellTag_Bin
ct_clone<-md$CellTag_ID
ct_bin<-md$CellTag_Bin


TBO_seurat@meta.data$CellTag_Clone<-ct_clone
table(TBO_seurat@meta.data$CellTag_Clone)
# None   RPM_Clone_13   RPM_Clone_14   RPM_Clone_16   RPM_Clone_19    RPM_Clone_2   RPM_Clone_22 
# 11591             26              7              4              3             27              6 
# RPM_Clone_23   RPM_Clone_33   RPM_Clone_36   RPM_Clone_37   RPM_Clone_40    RPM_Clone_5   RPM_Clone_53 
# 829              6              9              9              5              3             11 
# RPM_Clone_55    RPM_Clone_6    RPM_Clone_9   RPMA_Clone_1 RPMA_Clone_100  RPMA_Clone_11  RPMA_Clone_12 
# 20             18              3              2             11            235             68 
# RPMA_Clone_13  RPMA_Clone_17   RPMA_Clone_2  RPMA_Clone_23  RPMA_Clone_24  RPMA_Clone_26  RPMA_Clone_28 
# 99             18             91             68             30             85            213 
# RPMA_Clone_3  RPMA_Clone_32  RPMA_Clone_36   RPMA_Clone_4  RPMA_Clone_42  RPMA_Clone_44   RPMA_Clone_5 
# 114            155             70              2             22             13              4 
# RPMA_Clone_53  RPMA_Clone_55  RPMA_Clone_58  RPMA_Clone_63  RPMA_Clone_64  RPMA_Clone_65  RPMA_Clone_66 
# 18             56            237            219            277             37            125 
# RPMA_Clone_67  RPMA_Clone_68  RPMA_Clone_69   RPMA_Clone_7  RPMA_Clone_72  RPMA_Clone_74  RPMA_Clone_75 
# 55            735            175              2              9             22             12 
# RPMA_Clone_76  RPMA_Clone_78   RPMA_Clone_8  RPMA_Clone_81  RPMA_Clone_82  RPMA_Clone_84  RPMA_Clone_86 
# 30             48             28             19             22             31             10 
# RPMA_Clone_87  RPMA_Clone_90  RPMA_Clone_93  RPMA_Clone_95  RPMA_Clone_96 
# 19             12              8             11             19 

TBO_seurat@meta.data$CellTag_Binary<-ct_bin
table(TBO_seurat@meta.data$CellTag_Binary, TBO_seurat@meta.data$UnID)

#                RPM_Allo_Old   RPM_Allo_New  RPMA_Allo
# Celltagged            0          986      3536
# None               2190         2674      6727


saveRDS(TBO_seurat,"060124_RPM_RPMA_TBO_Seurat_wCT_Clones_from_adata.rds")
TBO_seurat<-readRDS("060124_RPM_RPMA_TBO_Seurat_wCT_Clones_from_adata.rds")


########################################################################
# Before signature assessment, normalize and log count data 
########################################################################
DefaultAssay(TBO_seurat)<-'RNA'
TBO_seurat<-NormalizeData(TBO_seurat)

## Perform Cell Cycle Analysis and visualize ##
########################################################################
# lowercase phase is from python, upper will be from seurat

s.genes<-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes

s.genes.list<-subset(mouse94, mouse94$human_homolog %in% s.genes)
s.genes.mouse<-(s.genes.list$gene_name)
s.genes.mouse

g2m.genes.list<-subset(mouse94, mouse94$human_homolog %in% g2m.genes)
g2m.genes.mouse<-(g2m.genes.list$gene_name)
g2m.genes.mouse

TBO_seurat <- CellCycleScoring(TBO_seurat, s.features = s.genes.mouse, g2m.features = g2m.genes.mouse)
TBO_seurat@meta.data$Phase<-factor(TBO_seurat@meta.data$Phase,levels=c("G1","S","G2M"))

DimPlot(TBO_seurat, reduction = "umap", group.by = "Phase",shuffle=TRUE,label=FALSE,pt.size=.4,cols=c("indianred3","green3","royalblue4"))+NoAxes()

### Look at distribution of labeled tumors by cluster

##### What % of cells occupy what cluster for all types ####
Idents(TBO_seurat)<-'leiden_scVI_1.2'
Idents(TBO_seurat)

x<-table(Idents(TBO_seurat),TBO_seurat@meta.data$phase)
x
proportions <- as.data.frame(100*prop.table(x, margin = 1))
proportions

colnames(proportions)<-c("Cluster", "Sample", "Frequency")


ggbarplot(proportions, x="Sample", y="Frequency", fill = "Sample", group = "Sample", ylab = "Frequency (percent)", xlab="Phase", palette =c("indianred3","green3","royalblue4"))+ theme_bw()+ facet_wrap(facets = "Cluster", scales="free_y", ncol =4)+ theme(axis.text.y = element_text(size=12)) +rotate_x_text(angle = 45)

# Stacked
p<-ggplot(proportions, aes(fill=Sample, y=Frequency, x=Cluster)) +
  geom_bar(position="stack", stat="identity")

p + scale_fill_manual(values=c("indianred3","green3","royalblue4"))+ theme_bw()+ theme(axis.text.y = element_text(size=20), axis.text.x=element_text(size=20), axis.title.x =element_text(size=20), axis.title.y = element_text(size=20), legend.text = element_text(size=20), legend.title = element_text(size=20))


# Ext Data Fig. 5c,e)
VlnPlot(TBO_seurat, features = c("Ascl1","Neurod1","Pou2f3"),group.by="Genotype", cols=c("darkorchid4","orange"),pt.size=.01, ncol=3)
VlnPlot(TBO_seurat, features = c("Ascl1","Neurod1","Pou2f3"),group.by="leiden_scVI_1.2", cols=colors,pt.size=.01, ncol=3)

########################################################################
######### Add Signature Data ##################
########################################################################
library(viridis)
library(ggplot2)
library(ggpubr)

######### Look at other signatures from BioRxiv and cell types of OE
###########
## Look at ND1 vs ASCl1 ChIP targets as scores (Fig 3j)
getwd()
chip<-read.csv("Signatures/ASCL1_NEUROD1_ChIP_Targets.csv")
achip<-chip$Borromeo_ASCL1_Targets
nchip<-chip$Borromeo_Oliver_NEUROD1_Targets
pchip<-chip$POU2F3_Targets
achip<-chip$Borromeo_ASCL1_Targets
nchip<-chip$Borromeo_Oliver_NEUROD1_Targets

pchip<-subset(mouse94, mouse94$human_homolog %in% pchip)
pchip_mouse<-(pchip$gene_name)
pchip_mouse


TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(achip),
  name = 'ASCL1_Targets')

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(nchip),
  name = 'NEUROD1_Targets')

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(pchip_mouse),
  name = 'POU2F3_Targets')


library(viridis)

FeaturePlot(TBO_seurat, features = c("ASCL1_Targets1"), pt.size=0.2, reduction='umap')+scale_color_viridis(option="rocket",direction=-1)& NoAxes()
FeaturePlot(TBO_seurat, features = c("NEUROD1_Targets1"), pt.size=0.2, reduction='umap')+scale_color_viridis(option="rocket",direction=-1)& NoAxes()
FeaturePlot(TBO_seurat, features = c("POU2F3_Targets1"), pt.size=0.2, reduction='umap')+scale_color_viridis(option="rocket",direction=-1)& NoAxes()
VlnPlot(TBO_seurat, features = c("POU2F3_Targets1"),group.by="Genotype", cols=c("darkorchid4","orange"),pt.size=.01, ncol=1,alpha=0.05)


######## Cell Type signatures for Fig. 3k ###########
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

ct<-read.csv("Signatures/Montoro_Iono_Human.csv")
ct$Montoro_Iono_Human
ionoh<-subset(mouse94, mouse94$human_homolog %in% ct$Montoro_Iono_Human)
ionoh_hum<-(ionoh$gene_name)


ct<-read.csv("Signatures/Montoro_Iono_Extended.csv")
ct<-ct$Montoro_iono_Ext
ct<-subset(mouse94, mouse94$gene_name %in% ct)
ct<-ct$human_homolog
write.csv(ct,"Ionocyte_Mouse_HumanHomologs.csv")

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(tuft_noC),
  name = 'Tuft_Consensus')

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(basal_noC),
  name = 'Basal_Consensus')

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(ne_noC),
  name = 'NE_Consensus')

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(iono_noC),
  name = 'Ionocyte_Consensus')

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(ionoh_hum),
  name = 'Iono_Human')

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(ct),
  name = 'Iono_Mouse_Ext')


library(viridis)

FeaturePlot(TBO_seurat, features = c("NE_Consensus1"), pt.size=0.2, reduction='umap')+ scale_color_gradientn(colors=rocket(10, direction=-1)) & NoAxes()
FeaturePlot(TBO_seurat, features = c("Basal_Consensus1"), pt.size=0.2, reduction='umap')+ scale_color_gradientn(colors=rocket(10, direction=-1)) & NoAxes()
FeaturePlot(TBO_seurat, features = c("Tuft_Consensus1"), pt.size=0.2, reduction='umap')+ scale_color_gradientn(colors=rocket(10, direction=-1)) & NoAxes()
FeaturePlot(TBO_seurat, features = c("Ionocyte_Consensus1"), pt.size=0.2, reduction='umap')+ scale_color_gradientn(colors=rocket(10, direction=-1)) & NoAxes()

# Ext Data Fig 9d #

FeaturePlot(TBO_seurat, features = c("Iono_Human1"), pt.size=0.2, reduction='umap')+ scale_color_gradientn(colors=rocket(10, direction=-1)) & NoAxes()
FeaturePlot(TBO_seurat, features = c("Iono_Mouse_Ext1"), pt.size=0.2, reduction='umap')+ scale_color_gradientn(colors=rocket(10, direction=-1)) & NoAxes()


#########################################################
## Look at archetype sigs as scores (Fig. 3l)

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


TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(a_sc_mouse),
  name = 'A_Archetype')

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(a2_sc_mouse),
  name = 'A2_Archetype')

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(n_sc_mouse),
  name = 'N_Archetype')

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(p_sc_mouse),
  name = 'P_Archetype')

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(y_sc_mouse),
  name = 'Y_Archetype')


library(viridis)


a<-FeaturePlot(TBO_seurat, features = c("A_Archetype1"), pt.size=0.2, reduction='umap',order=TRUE)+scale_color_viridis(option="rocket",direction=-1)& NoAxes()
b<-FeaturePlot(TBO_seurat, features = c("A2_Archetype1"), pt.size=0.2, reduction='umap',order=TRUE)+scale_color_viridis(option="rocket",direction=-1)& NoAxes()
c<-FeaturePlot(TBO_seurat, features = c("N_Archetype1"), pt.size=0.2, reduction='umap',order=TRUE)+scale_color_viridis(option="rocket",direction=-1)& NoAxes()
d<-FeaturePlot(TBO_seurat, features = c("P_Archetype1"), pt.size=0.2, reduction='umap',order=TRUE)+scale_color_viridis(option="rocket",direction=-1)& NoAxes()
FeaturePlot(TBO_seurat, features = c("Y_Archetype1"), pt.size=0.2, reduction='umap',order=TRUE)+scale_color_viridis(option="rocket",direction=-1)& NoAxes()

grid.arrange(a,b,c,d,nrow=1)
VlnPlot(TBO_seurat, features = c("A_Archetype1","A2_Archetype1","N_Archetype1","P_Archetype1","Y_Archetype1"),group.by="Pheno", cols=pheno_col,pt.size=.01,alpha=0.01, ncol=5)
VlnPlot(TBO_seurat, features = c("A2_Archetype1"),group.by="leiden_scVI_1.2", cols=colors,pt.size=.01,alpha=0.01, ncol=1)

rpm_rpma_tbo<-TBO_seurat@meta.data
write.csv(rpm_rpma_tbo,"RPM_RPMA_TBO_Fig3e_metadata.csv")



########## Adding additional CARIS signatures (Fig. 6d) ###########
setwd("/Users/abbieireland/Desktop/scRNAseq")

caris<-read.csv("Signatures/Caris_Top100_HumanSCLC.csv")
a<-caris$Caris_SCLC.A
n<-caris$Caris_SCLC.N
p<-caris$Caris_SCLC.P
y<-caris$Caris_SCLC.Y
mixed<-caris$Caris_SCLC.Mixed
tn<-caris$Caris_SCLC.TN

# Convert to mouse
a<-subset(mouse94, mouse94$human_homolog %in% a)
a<-a$gene_name
n<-subset(mouse94, mouse94$human_homolog %in% n)
n<-n$gene_name
p<-subset(mouse94, mouse94$human_homolog %in% p)
p<-p$gene_name
y<-subset(mouse94, mouse94$human_homolog %in% y)
y<-y$gene_name
mixed<-subset(mouse94, mouse94$human_homolog %in% mixed)
mixed<-mixed$gene_name
tn<-subset(mouse94, mouse94$human_homolog %in% tn)
tn<-tn$gene_name


TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(a),
  name = 'Caris_A')
TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(n),
  name = 'Caris_N')
TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(p),
  name = 'Caris_P')
TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(y),
  name = 'Caris_Y')
TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(mixed),
  name = 'Caris_Mixed')
TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(tn),
  name = 'Caris_TN-')

FeaturePlot(TBO_seurat, features = c("Caris_TN-1"), pt.size=0.2, reduction='umap',)+scale_color_viridis(option="rocket",direction=-1)& NoAxes()

?FeaturePlot
table(TBO_seurat@meta.data$Pheno)
Ionocyte_Consensus1
VlnPlot(TBO_seurat,features = c("Caris_A1", "Caris_N1", "Caris_Mixed1", "Caris_Y1", "Caris_P1","Caris_TN-1"),group.by=c("Pheno"),cols=pheno_col,alpha=0, ncol=3)
table(TBO_seurat@meta.data$Pheno)

VlnPlot(TBO_seurat,features = c("Ionocyte_Consensus1", "Basal_Consensus1"),group.by=c("Pheno"),cols=pheno_col,alpha=0, ncol=3)


# Using liu consensus NMF lists, top 100 enriched (Fig 6f) #
liu<-read.csv("Signatures/Liu_NMF_Sigs.csv")
nmf1<-liu$NMF1
nmf2<-liu$NMF2
nmf3<-liu$NMF3
nmf4<-liu$NMF4


nmf1<-subset(mouse94, mouse94$human_homolog %in% nmf1)
nmf1_m<-(nmf1$gene_name)
nmf1_m

nmf2<-subset(mouse94, mouse94$human_homolog %in% nmf2)
nmf2_m<-(nmf2$gene_name)
nmf2_m

nmf3<-subset(mouse94, mouse94$human_homolog %in% nmf3)
nmf3_m<-(nmf3$gene_name)
nmf3_m

nmf4<-subset(mouse94, mouse94$human_homolog %in% nmf4)
nmf4_m<-(nmf4$gene_name)
nmf4_m

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(nmf1_m),
  name = 'NMF1-A')
TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(nmf2_m),
  name = 'NMF2-A/N')
TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(nmf3_m),
  name = 'NMF3-I')
TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(nmf4_m),
  name = 'NMF4-P')

VlnPlot(TBO_seurat,features = c("NMF1-A1","NMF2-A/N1","NMF3-I1","NMF4-P1"), group.by="Pheno",cols=pheno_col,pt.size=.01, ncol=2, alpha=0.02)
FeaturePlot(TBO_seurat,features = c("NMF1-A1","NMF2-A/N1","NMF3-I1","NMF4-P1"),pt.size=.01, ncol=2, alpha=0.2, order=TRUE)
saveRDS(TBO_seurat,"111424_TBO_seurat_RPM_RPMA.rds")



####### Add inflammatory signatures (Fig 6f) #######
inf<-read.csv("Signatures/Human_inflamed_Gay.csv")
t<-inf$T_Cell_Inflamed
mhc<-inf$MHC

t<-subset(mouse94, mouse94$human_homolog %in% t)
t_cellm<-t$gene_name
mhc<-subset(mouse94, mouse94$human_homolog %in% mhc)
mhc_m<-mhc$gene_name
mhc_m

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(mhc_m),
  name = 'MHC_Sig_Gay')

TBO_seurat<-AddModuleScore(
  object = TBO_seurat,
  features = list(t_cellm),
  name = 'T_Cell_Inflamed_Gay')

# Fig 6f feature plots

FeaturePlot(TBO_seurat, features = c("MHC_Sig_Gay1"), pt.size=0.2, reduction='umap',order=TRUE)+scale_color_viridis(option="rocket",direction=-1)& NoAxes()

FeaturePlot(TBO_seurat, features = c("NMF1-A1"), pt.size=0.2, reduction='umap',order=TRUE)+scale_color_viridis(option="rocket",direction=-1)& NoAxes()
FeaturePlot(TBO_seurat, features = c("NMF2-A/N1"), pt.size=0.2, reduction='umap',order=TRUE)+scale_color_viridis(option="rocket",direction=-1)& NoAxes()
FeaturePlot(TBO_seurat, features = c("NMF4-P1"), pt.size=0.2, reduction='umap',order=TRUE)+scale_color_viridis(option="rocket",direction=-1)& NoAxes()
FeaturePlot(TBO_seurat, features = c("NMF3-I1"), pt.size=0.2, reduction='umap',order=TRUE)+scale_color_viridis(option="rocket",direction=-1)& NoAxes()

FeaturePlot(TBO_seurat, features = c("T_Cell_Inflamed_Gay1"), pt.size=0.2, reduction='umap',order=TRUE)+scale_color_viridis(option="rocket",direction=-1)& NoAxes()

VlnPlot(TBO_seurat, features = c("MHC_Sig_Gay1"),group.by="Pheno", cols=pheno_col,pt.size=.01,alpha=0.01, ncol=1)
VlnPlot(TBO_seurat, features = c("T_Cell_Inflamed_Gay1"),group.by="Pheno", cols=pheno_col,pt.size=.01,alpha=0.01, ncol=1)

TBO_seurat@meta.data$`NMF2-A/N1`


#Clones by phenotype

Idents(TBO_seurat)<-'Pheno'
Idents(TBO_seurat)

x<-table(TBO_seurat@meta.data$CellTag_Clone, Idents(TBO_seurat))
x

proportions <- as.data.frame(100*prop.table(x, margin = 1))
proportions
colnames(proportions)<-c("Cluster", "Sample", "Frequency")

proportions 
# Stacked
p<-ggplot(proportions, aes(fill=Sample, y=Frequency, x=Cluster)) + 
  geom_bar(position="stack", stat="identity")

p + scale_fill_manual(values=pheno_col)+ 
  theme_bw()+ theme(axis.text.y = element_text(size=20), 
                    axis.text.x=element_text(size=14), axis.title.x =element_text(size=14), 
                    axis.title.y = element_text(size=18), legend.text = element_text(size=12), 
                    legend.title = element_text(size=18))+rotate_x_text(size=7,angle = 90)


write.csv(p$data,"comb_data_subtype_for_graph.csv")


#Genotype by leiden

Idents(TBO_seurat)<-'leiden_scVI_1.2'
Idents(TBO_seurat)

x<-table(TBO_seurat@meta.data$Genotype, Idents(TBO_seurat))
x

proportions <- as.data.frame(100*prop.table(x, margin = 1))
proportions
colnames(proportions)<-c("Cluster", "Sample", "Frequency")

proportions 
# Stacked
p<-ggplot(proportions, aes(fill=Sample, y=Frequency, x=Cluster)) + 
  geom_bar(position="stack", stat="identity")

p + scale_fill_manual(values=colors)+ 
  theme_bw()+ theme(axis.text.y = element_text(size=20), 
                    axis.text.x=element_text(size=14), axis.title.x =element_text(size=14), 
                    axis.title.y = element_text(size=18), legend.text = element_text(size=12), 
                    legend.title = element_text(size=18))+rotate_x_text(size=7,angle = 90)


# saveRDS(TBO_seurat,"060124_RPMRPMA_TBO_AlloOnly_Celltaganno_SCLCsubs.rds")
# setwd()
# 
# #Read in here 09 2024
# TBO_seurat<-readRDS("060124_RPMRPMA_TBO_AlloOnly_Celltaganno_SCLCsubs.rds")
# TBO_seurat

# NEWEST 
saveRDS(TBO_seurat,"092824_RPM_RPMA_TBO_Allos_CellTagAnno_Pheno.rds")
TBO_seurat<-readRDS("092824_RPM_RPMA_TBO_Allos_CellTagAnno_Pheno.rds")


######## Add NE Score (Fig. 3i) ##########

### For NE score assignment

#### Converting to SCE to Add NE Score and create Violin plots
?assay
# Convert seurat object to SCE
library(SingleCellExperiment)
sce<-as.SingleCellExperiment(TBO_seurat)
library(SeuratObject)
library(SeuratData)
library(SeuratDisk)
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
library(scater)

# Plot NE spearman, for example
x<-scater::plotColData(sce, x = "leiden_scVI_1.2", y = "NE_spearman", colour_by = "Pheno")+scale_discrete_manual(aesthetics = c("colour", "fill"),values=pheno_col)
x

# For example for POU2F3 targets
# Plot by genotype
x<-scater::plotColData(sce, x = "Genotype", y = "POU2F3_Targets1", colour_by = "Genotype")+scale_discrete_manual(aesthetics = c("colour", "fill"),values=c("darkorchid4","orange"))

#"RPM"="orange","RPMA"="firebrick2","SNL"="dodgerblue4")
# Add stats
fortest<-x+ geom_boxplot(fill=c("darkorchid4","orange"), alpha=1/5, position = position_dodge(width = .2),size=0.2,color="black", notch=TRUE, notchwidth=0.3, outlier.shape = 2, outlier.colour=NA)
fortest+ylim(-.75,1.25)
?geom_boxplot

# Add p-values comparing groups
# Specify the comparisons you want
my_comparisons <- list( c("RPM", "RPMA"))

library(ggplot2) 
library(ggpubr) 
library(tidyverse) 

fortest + stat_compare_means(comparisons = my_comparisons, method="t.test", paired=FALSE) # Add pairwise comparisons p-value
# Add global p-value

fortest + stat_compare_means(comparisons = my_comparisons, method="t.test", paired=FALSE,label = "p.signif") # Add pairwise comparisons p-value



## Back to seurat to add NE_spearman score (Fig. 3i)
TBO_seurat@meta.data$NE_spearman<-sce$NE_spearman

ne_spear<-TBO_seurat@meta.data$NE_spearman
ne_spear
write.csv(ne_spear,"ne_spearFig3.csv")



##############################################################################################################
######## Begin clonal analysis in FA space for Fig 4 #####################################################################
####### Add FA projections ######################################################################################
getwd()
# Read in adata object as SCE
adata2<-readH5AD("092824_RPM_RPMA_Allo_dpt_cellrank2.h5ad")
table(adata2$Cre)
table(adata2$leiden_scVI_1.2)
table(adata2$UnID)
# RPMA_Allo RPM_Allo_New RPM_Allo_Old 
# 10256         3653         2190 

dim(assay(adata2,"counts"))
counts(adata2)<-assay(adata2,"counts")
counts(adata2)
length(rownames(adata2))
counts(adata2)
assay(adata2,"norm")
rowData(adata2)
length(rowData(adata2)$gene_ids)

?CreateSeuratObject.Assay

#Convert SCE to seurat
TBO_seurat_fa <- CreateSeuratObject(counts = counts(adata2), meta.data = as.data.frame(colData(adata2)))
TBO_seurat_fa
# An object of class Seurat 
# 35841 features across 16099 samples within 1 assay 
# Active assay: RNA (35841 features, 0 variable features)
# 1 layer present: counts
DefaultAssay(TBO_seurat_fa)

#Create an assay of normalized gene expression
norm_assay <- CreateAssayObject(counts = assay(adata2,"norm"))
norm_assay

TBO_seurat_fa[["norm"]] <- norm_assay

# Set assay norm as default assay for seurat object
DefaultAssay(TBO_seurat_fa)<-'norm'


# Add embeddings of umap, and X_scVI_1.4 to assay norm
adata2
reducedDim(adata2, "X_umap")
test<-reducedDim(adata2, "X_umap")
colnames(test)<-c("UMAP_1","UMAP_2")
rownames(test)<-colnames(TBO_seurat_fa)
dim(test)
TBO_seurat_fa[['umap']] <- CreateDimReducObject(test, key="UMAP_", assay = "RNA")
TBO_seurat_fa[['umap']]@cell.embeddings


#colors<-c('#F39B7F99','brown1','#ff7f0e', 'turquoise', '#2ca02c', 'deepskyblue4', '#e7298a',  '#0aa6d8',  '#A0522D', '#a1c299', 'gold', '#984ea3', '#3C548899', '#377eb8',  '#d62728','navyblue', '#a1c299', 'purple','gray','maroon','gray20')
colors<-c('deepskyblue4', '#ff7f0e','#2ca02c', 'brown1', 'purple', '#A0522D', 'pink2',  '#a1c299',  '#4ba4d3', '#b5d7e4','#ef9171','#a6eb9a')

DimPlot(TBO_seurat_fa,group.by='leiden_scVI_1.2',cols=colors, reduction='umap',label=FALSE,label.size=6)&NoAxes()

# Add embeddings of umap, and X_scVI_1.4 to assay norm
# Add fa embeddings to RPM_Allo
reducedDim(adata2, "X_draw_graph_fa")
test<-reducedDim(adata2, "X_draw_graph_fa")
head(test)
colnames(test)
colnames(test)<-c("FA_1","FA_2")
rownames(test)<-colnames(TBO_seurat_fa)
dim(test)
head(colnames(TBO_seurat_fa))

TBO_seurat_fa[['fa']] <- CreateDimReducObject(test, key="FA_", assay = "RNA")
TBO_seurat_fa[['fa']]@cell.embeddings

# Fig. 4c
DimPlot(TBO_seurat_fa,group.by='leiden_scVI_1.2',cols=colors, reduction='fa',label=FALSE,label.size=6)&NoAxes()
write.csv(TBO_seurat_fa@meta.data,"TBO_fa_metadata_4clones.csv")

DimPlot(TBO_seurat_fa,group.by='Genotype',cols=c("darkorchid4","#ff7f0e"), reduction='fa',label=FALSE,label.size=6, shuffle=TRUE)&NoAxes()

######### Annotate CellTag Data #############

md<-read.csv("TBO_fa_metadata_4clones_anno.csv")
md$CellTag_ID
md$CellTag_Bin
ct_clone<-md$CellTag_ID
ct_bin<-md$CellTag_Bin
length(ct_clone)

TBO_seurat_fa@meta.data$CellTag_Clone<-ct_clone
table(TBO_seurat_fa@meta.data$CellTag_Clone)
# None   RPM_Clone_13   RPM_Clone_14   RPM_Clone_16   RPM_Clone_19    RPM_Clone_2   RPM_Clone_22   RPM_Clone_23 
# 11578             26              7              4              3             27              6            829 
# RPM_Clone_33   RPM_Clone_36   RPM_Clone_37   RPM_Clone_40    RPM_Clone_5   RPM_Clone_53   RPM_Clone_55    RPM_Clone_6 
# 6              9              9              5              3             11             20             18 
# RPM_Clone_9   RPMA_Clone_1 RPMA_Clone_100  RPMA_Clone_11  RPMA_Clone_12  RPMA_Clone_13  RPMA_Clone_17   RPMA_Clone_2 
# 3              2             11            235             68             99             18             91 
# RPMA_Clone_23  RPMA_Clone_24  RPMA_Clone_26  RPMA_Clone_28   RPMA_Clone_3  RPMA_Clone_32  RPMA_Clone_36   RPMA_Clone_4 
# 68             30             85            212            114            155             70              2 
# RPMA_Clone_42  RPMA_Clone_44   RPMA_Clone_5  RPMA_Clone_53  RPMA_Clone_55  RPMA_Clone_58  RPMA_Clone_63  RPMA_Clone_64 
# 22             13              4             18             56            237            219            277 
# RPMA_Clone_65  RPMA_Clone_66  RPMA_Clone_67  RPMA_Clone_68  RPMA_Clone_69   RPMA_Clone_7  RPMA_Clone_72  RPMA_Clone_74 
# 37            125             55            735            175              2              9             22 
# RPMA_Clone_75  RPMA_Clone_76  RPMA_Clone_78   RPMA_Clone_8  RPMA_Clone_81  RPMA_Clone_82  RPMA_Clone_84  RPMA_Clone_86 
# 12             30             48             28             19             22             31             10 
# RPMA_Clone_87  RPMA_Clone_90  RPMA_Clone_93  RPMA_Clone_95  RPMA_Clone_96 
# 19             12              8             11             19 

TBO_seurat_fa@meta.data$CellTag_Binary<-ct_bin
table(TBO_seurat_fa@meta.data$CellTag_Binary, TBO_seurat_fa@meta.data$UnID)
TBO_seurat_fa@meta.data$UnID<-factor(TBO_seurat_fa@meta.data$UnID,c("RPM_Allo_Old","RPM_Allo_New","RPMA_Allo"))
#                RPM_Allo_Old   RPM_Allo_New  RPMA_Allo
# Celltagged            0          986      3535
# None               2190         2667      6721


# Assign pheno/cell states as in Fig. 4c
TBO_seurat_fa@meta.data$Pheno<- ifelse(TBO_seurat_fa@meta.data$leiden_scVI_1.2 %in% c("11"), "Tuft",
                                    ifelse(TBO_seurat_fa@meta.data$leiden_scVI_1.2 %in% c("3"), "NE",
                                           ifelse(TBO_seurat_fa@meta.data$leiden_scVI_1.2 %in% c("4","9"), "Neuronal",
                                                  ifelse(TBO_seurat_fa@meta.data$leiden_scVI_1.2 %in% c("2"), "NE/Neuronal",
                                                         ifelse(TBO_seurat_fa@meta.data$leiden_scVI_1.2 %in% c("10"), "Basal",
                                                                ifelse(TBO_seurat_fa@meta.data$leiden_scVI_1.2 %in% c("0","1","6","5","7","8"), "Triple-Neg","NA"))))))




table(TBO_seurat_fa@meta.data$Pheno)
TBO_seurat_fa@meta.data$Pheno<-factor(TBO_seurat_fa@meta.data$Pheno, levels=c("NE","NE/Neuronal","Neuronal","Tuft","Triple-Neg","Basal","NA"))
pheno_col<-c("brown2","darkorchid4","dodgerblue","orange","turquoise4","turquoise")

DimPlot(TBO_seurat_fa,group.by='Pheno',cols=pheno_col, reduction='fa',label=FALSE,label.size=6) & NoAxes()
DimPlot(TBO_seurat_fa,group.by='Pheno',cols=pheno_col, reduction='umap',label=FALSE,label.size=6) & NoAxes()

saveRDS(TBO_seurat_fa, "11_2024_TBO_seurat_fa.rds")
########################################################################
######### Visualize CellTag Data ##################
########################################################################
cs<-table(TBO_seurat_fa@meta.data$CellTag_Clone, TBO_seurat_fa@meta.data$UnID)
cs
write.csv(cs,"092824_postQC_clonesize.csv")

########################################################################
######## Visualizing Celltag/clone data #########
########################################################################
Idents(TBO_seurat_fa)<-'CellTag_Binary'
# Robust defined as >5 cells per clone post-QC
TBO_seurat_fa@meta.data$Robust<-md$Robust
table(TBO_seurat_fa@meta.data$Robust, TBO_seurat_fa@meta.data$Genotype)
# RPM RPMA
# No     4870 6731
# Robust  973 3525


clones<-subset(TBO_seurat_fa,idents=c("Celltagged"))
table(clones@meta.data$Robust, clones@meta.data$Genotype)
#         RPM RPMA
# No       13   10
# Robust  973 3525


Idents(clones)<-'Robust'
Idents(clones)

clones<-subset(clones,idents=c("Robust"))
table(clones@meta.data$Robust, clones@meta.data$Genotype)
#       RPM RPMA
# Robust  973 3525

test<-table(clones@meta.data$CellTag_Clone)
test

#Clones by Leiden
Idents(clones)<-'leiden_scVI_1.2'
Idents(clones)

x<-table(clones@meta.data$CellTag_Clone,Idents(clones))
x
proportions <- as.data.frame(100*prop.table(x, margin = 1))

# proportions$Cluster
colnames(proportions)<-c("Cluster", "Sample", "Frequency")
proportions
# Stacked
p<-ggplot(proportions, aes(fill=Sample, y=Frequency, x=Cluster)) + 
  geom_bar(position="stack", stat="identity")


colors<-c('deepskyblue4', '#ff7f0e','#2ca02c', 'brown1', 'purple', '#A0522D', 'pink2',  '#a1c299',  '#0aa6d8', 'lightblue','salmon1','lightgreen')

p + scale_fill_manual(values=colors)+ 
  theme_bw()+ theme(axis.text.y = element_text(size=20), 
                    axis.text.x=element_text(size=14), axis.title.x =element_text(size=14), 
                    axis.title.y = element_text(size=18), legend.text = element_text(size=12), 
                    legend.title = element_text(size=18))+rotate_x_text(size=7,angle = 90)



# Cluster distribution in unbiased way...

# Uset this file to get x axis labels by removing dublicates
write.csv(p$data,"test_stacked_cluster_060124.csv")
df<-p$data

#Organize into frequency per clone
df
df <- df[order(df$Cluster),]
df

# Split every 12 into different column
# number of splits to do

# number of splits to do
n_cols <- nrow(df)/12
n_cols

# subsetting each split seperately, store in a list
df_l <- 
  lapply(1:n_cols,function(i){
    startrow <- (i-1) * 12 + 1
    endrow <- i*12
    
    df[startrow:endrow, c(1,2,3)]
  })


# bind list elements together
df_new <-
  do.call(cbind, df_l)
write.csv(df_new,"092824_cluster_percent_perclone_new.csv")


df2<-subset(df_new,colnames(df_new) %in% 'Frequency')
?subset

drop <- c("Sample","Cluster")
df2 = df_new[,!(names(df_new) %in% drop)]
df2
#Add colnames to this that match x axis labels of p stacked bargraph and good to go 
write.csv(df2,"092824_cluster_percent_perclone_onlyfreq.csv")
df2

# Perform unbiased hierarchical clustering and plot with pheatmap
library(readxl)
test<-read_excel("092824_cluster_percent_perclone_onlyfreq_labeled.xls", range=cell_cols("B:BA"))
test

library(pheatmap)
t<-pheatmap(test,cutree_rows = 1,cutree_cols = 5,cellwidth = 10, cellheight = 5,fontsize = 10, cluster_rows=FALSE,border_color=NA,color = colorRampPalette(c("darkturquoise","black","red2"))(30))
t$tree_col$labels

t<-pheatmap(test,cutree_rows = 1,cutree_cols = 1,cellwidth = 10, cellheight = 5,fontsize = 10, cluster_rows=FALSE,border_color=NA,color = colorRampPalette(c("darkturquoise","black","red2"))(30))


#Order of clones (obtain from copy paste of pdf export of t)
order<-read.csv("order_clones_092824.csv")
order<-order$order
order
proportions$Cluster

#Clones by Leiden (Fig. 4d)
Idents(clones)<-'leiden_scVI_1.2'

x<-table(clones@meta.data$CellTag_Clone,Idents(clones))
x
rownames(x)
proportions <- as.data.frame(100*prop.table(x, margin = 1))
proportions

proportions
proportions$Cluster
colnames(proportions)<-c("Cluster", "Sample", "Frequency")
table(proportions$Cluster)
proportions$Cluster<-factor(proportions$Cluster, order)
proportions$Sample<-factor(proportions$Sample, c("0","1","2","3","4","5","6","8","9","10","11","7"))
# Stacked
p<-ggplot(proportions, aes(fill=Sample, y=Frequency, x=Cluster)) + 
  geom_bar(position="stack", stat="identity")


colors<-c('deepskyblue4', '#ff7f0e','#2ca02c', 'brown1', 'purple', '#A0522D', 'pink2',  '#0aa6d8', 'lightblue','salmon1','lightgreen','#a1c299')

p + scale_fill_manual(values=colors)+ 
  theme_bw()+ theme(axis.text.y = element_text(size=20), 
                    axis.text.x=element_text(size=14), axis.title.x =element_text(size=14), 
                    axis.title.y = element_text(size=18), legend.text = element_text(size=12), 
                    legend.title = element_text(size=18))+rotate_x_text(size=7,angle = 90)


saveRDS(clones,"092824_RPM_RPMA_Allo_CelLTagClones_FA.rds")
clones<-readRDS("092824_RPM_RPMA_Allo_CelLTagClones_FA.rds")


##Annotate by group and look at dynamics##
write.csv(order,"order_clones_for_group_anno.csv")
write.csv(clones@meta.data, "new_meta_092824.csv")
# Add group data to metadata and read back in 
md<-read.csv("new_meta_092824_wGroups.csv")
table(md$Group)

clones@meta.data$Clone_Dynamics<-md$Group
Idents(clones)<-'Clone_Dynamics'
Idents(clones)
table(clones@meta.data$CellTag_Clone)


# Visualize different groups of clones in UMAP 
clones@meta.data$Pheno<- ifelse(clones@meta.data$leiden_scVI_1.2 %in% c("11"), "Tuft",
                                    ifelse(clones@meta.data$leiden_scVI_1.2 %in% c("3"), "NE",
                                           ifelse(clones@meta.data$leiden_scVI_1.2 %in% c("4","9"), "Neuronal",
                                                  ifelse(clones@meta.data$leiden_scVI_1.2 %in% c("2"), "NE/Neuronal",
                                                         ifelse(clones@meta.data$leiden_scVI_1.2 %in% c("10"), "Basal",
                                                                ifelse(clones@meta.data$leiden_scVI_1.2 %in% c("0","1","6","5","7","8"), "Triple-Neg","NA"))))))




table(clones@meta.data$Pheno)
clones@meta.data$Pheno<-factor(clones@meta.data$Pheno, levels=c("NE","NE/Neuronal","Neuronal","Tuft","Triple-Neg","Basal","NA"))
DimPlot(clones,group.by='Pheno',cols=pheno_col, reduction='fa',label=FALSE,label.size=6) & NoAxes()



# Visualize clonal patterns in FA 
# Patterns 1-4
g1<-subset(clones,idents=c('Group_1'))
g2<-subset(clones,idents=c('Group_2'))
g3<-subset(clones,idents=c('Group_3'))
g4<-subset(clones,idents=c('Group_4'))
g5<-subset(clones,idents=c('Group_5'))

table(g5$Pheno)

# For Fig. 4f overlap
# Export to overlap with plain gray fa map
clusterumap<-DimPlot(g5,group.by="Pheno",reduction='fa',order=TRUE, cols=pheno_col, pt.size=6, shuffle=TRUE)+ggtitle("Group 1") & NoAxes()&theme(legend.background = element_rect(fill = "transparent"),
                                                                                                           legend.box.background = element_rect(fill = "transparent"),
                                                                                                           panel.background = element_rect(fill = "transparent"),
                                                                                                           panel.grid.major = element_blank(),
                                                                                                           panel.grid.minor = element_blank(),
                                                                                                           plot.background = element_rect(fill = "transparent",
                                                                                                                                          color = NA))


ggsave("UMAPtest.png", plot= clusterumap, width=7, height=5, dpi=300, bg = "transparent")

# For Fig. 4h overlaps
table(g1$leiden_scVI_1.2)

clusterumap<-
  FeaturePlot(g5, features = c("dpt_pseudotime"), pt.size=4, reduction='fa')+scale_color_viridis(option="turbo",direction=-1)& NoAxes()&theme(legend.background = element_rect(fill = "transparent"),
                                                                                                                                                       legend.box.background = element_rect(fill = "transparent"),
                                                                                                                                                       panel.background = element_rect(fill = "transparent"),
                                                                                                                                                       panel.grid.major = element_blank(),
                                                                                                                                                       panel.grid.minor = element_blank(),
                                                                                                                                                       plot.background = element_rect(fill = "transparent",
                                                                                                                                                                                      color = NA))

ggsave("UMAPtest.png", plot= clusterumap, width=7, height=5, dpi=300, bg = "transparent")




table(g4@meta.data$leiden_scVI_1.2)
Idents(g4)<-'leiden_scVI_1.2'
c11<-subset(g4,ident=c("11") )
c11@meta.data$CellTag_Clone

g1<-rownames(g1@meta.data)
g2<-rownames(g2@meta.data)
g3<-rownames(g3@meta.data)
g4<-rownames(g4@meta.data)
g5<-rownames(g5@meta.data)

DimPlot(TBO_seurat_fa,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=g1,sizes.highlight=1, cols.highlight=c("red2"))+ggtitle("Group 1") & NoLegend() & NoAxes()
DimPlot(TBO_seurat_fa,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=g2,sizes.highlight=1, cols.highlight=c("orange"))+ggtitle("Group 2") & NoLegend() & NoAxes()
DimPlot(TBO_seurat_fa,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=g3,sizes.highlight=1, cols.highlight=c("green2"))+ggtitle("Group 3") & NoLegend() & NoAxes()
DimPlot(TBO_seurat_fa,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=g4,sizes.highlight=1, cols.highlight=c("royalblue2"))+ggtitle("Group 4") & NoLegend() & NoAxes()
DimPlot(TBO_seurat_fa,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=g1,sizes.highlight=.1, cols.highlight=c("gray"))+ggtitle("Group 1") & NoLegend() & NoAxes()

DimPlot(TBO_seurat,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=g5,sizes.highlight=1, cols.highlight=c("darkorchid3"))+ggtitle("Group 5") & NoLegend() & NoAxes()
?DimPlot

## Visualize individual clones for Ext Data Fig. 6b,c
test<-as.data.frame(clones$CellTag_Clone)
test$Barcodes<-rownames(test)
head(test)


test <- test %>% group_by(clones$CellTag_Clone)
test
test<-group_split(test)
test


DimPlot(TBO_seurat_fa,group.by="CellTag_Clone",reduction='umap',order=TRUE, cells.highlight=test[[52]]$Barcodes,sizes.highlight=2, cols.highlight=c("darkorchid3"))+ggtitle("Group 5") & NoLegend() & NoAxes()
DimPlot(TBO_seurat_fa,group.by="CellTag_Clone",reduction='umap',order=TRUE, cells.highlight=test[[1]]$Barcodes,sizes.highlight=2, cols.highlight=c("darkorchid3"))+ggtitle(paste0(as.data.frame(test[[52]][1])$`clones$CellTag_Clone`[2])) & NoLegend() & NoAxes()


plot_lst <- vector("list", length = 52)
for (i in 1:52) {
  g<-DimPlot(TBO_seurat_fa,group.by="CellTag_Clone",reduction='fa',order=TRUE, cells.highlight=test[[i]]$Barcodes,sizes.highlight=1, cols.highlight=c("#ff7f0e"))+ggtitle(paste0(as.data.frame(test[[i]][1])$`clones$CellTag_Clone`[2])) & NoLegend() & NoAxes()
  plot_lst[[i]] <- g
}

# Combine multiple plots for output, as desired
cowplot::plot_grid(plotlist = plot_lst[37:40], nrow=1)


# DEGs per cluster
Idents(TBO_seurat)<-'Pheno'
table(Idents(TBO_seurat))
TBO_seurat
ScaleData(TBO_seurat)
Cluster.Markers <- FindAllMarkers(TBO_seurat, only.pos = TRUE, logfc.threshold = 0.25, slot="scale.data")
?FindAllMarkers
test<-as.data.frame(Cluster.Markers)
table(test$cluster)

mouse_genes<-test$gene
length(mouse_genes)
mouse_genes<-subset(mouse94, mouse94$gene_name %in% mouse_genes)
human_genes<-(mouse_genes$human_homolog)
length(human_genes)
write.csv(test, "110724_RPM_RPMA_Allos_Pheno_Markers_Normalized.csv")

# DONE #
