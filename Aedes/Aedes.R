library(dplyr)
library(Seurat)
library(patchwork)
setwd('/home/a/asb55/asb55/aedes')

# Load the PBMC dataset
male.data <- Read10X(data.dir = "GSE160740_RAW/male/")
# Initialize the Seurat object with the raw (non-normalized data).
male <- CreateSeuratObject(counts = male.data, project = "male3k", min.cells = 3, min.features = 200)
male
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# Visualize QC metrics as a violin plot
VlnPlot(male, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
# Visualize QC metrics as a violin plot
VlnPlot(male, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

# plot1 <- FeatureScatter(male, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(male, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
plot2
male <- subset(male, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 )
male <- NormalizeData(male)







# Load the PBMC dataset
female.data <- Read10X(data.dir = "GSE160740_RAW/female/")
# Initialize the Seurat object with the raw (non-normalized data).
female <- CreateSeuratObject(counts = female.data, project = "female3k", min.cells = 3, min.features = 200)
female
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# Visualize QC metrics as a violin plot
VlnPlot(female, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
# Visualize QC metrics as a violin plot
VlnPlot(female, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

# plot1 <- FeatureScatter(female, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(female, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
plot2
female <- subset(female, subset = nFeature_RNA > 400 & nFeature_RNA < 2800 & nCount_RNA > 400)
female <- NormalizeData(female)



#merge two dataset
Aedes <- merge(male, y = c(female), add.cell.ids = c("male", "female"), project = "Aedes")





Aedes <- FindVariableFeatures(Aedes, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Aedes), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Aedes)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(Aedes)
Aedes <- ScaleData(Aedes, features = all.genes)
Aedes <- RunPCA(Aedes, features = VariableFeatures(object = Aedes))
# Examine and visualize PCA results a few different ways
print(Aedes[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Aedes, dims = 1:2, reduction = "pca")
DimPlot(Aedes, reduction = "pca")
DimHeatmap(Aedes, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(Aedes, dims = 1:15, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
Aedes <- JackStraw(Aedes, num.replicate = 100, dims = 50)
Aedes <- ScoreJackStraw(Aedes, dims = 1:50)
JackStrawPlot(Aedes, dims = 1:50)
ElbowPlot(Aedes, ndims = 50)
Aedes <- FindNeighbors(Aedes, dims = 1:35)
Aedes <- FindClusters(Aedes, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(Aedes), 5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
Aedes <- RunUMAP(Aedes, dims = 1:35)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(Aedes, reduction = "umap", group.by = "orig.ident")
DimPlot(Aedes, reduction = "umap")



# find markers for every cluster compared to all remaining cells, report only the positive
# ones
Aedes.markers <- FindAllMarkers(Aedes, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Aedes.markers
Aedes.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(Aedes.markers, file = "Aedes_MarkerGenes.tsv", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
Aedes31M <- which(Aedes.markers$cluster == 31)
Aedes31_Markers <- Aedes.markers[Aedes31M,]

write.csv(Aedes.markers, file = "/scratch/evassvis/asb55/aedes/Aedes_UNfiltered_MarkerGenes.tsv")

Aedesfilt <- subset(Aedes.markers, avg_log2FC > 1)
Aedesfilt
write.csv(Aedesfilt, file = "/scratch/evassvis/asb55/aedes/Aedes_afiltered_MarkerGenes.tsv")

save(Aedes, file = "/scratch/evassvis/asb55/aedes/Aedes_afterpca.Rdata")






load(file = "/scratch/evassvis/asb55/aedes/Aedes_afterpca.Rdata")

# GPCR Gene List
# read in Gene list or table
GPCR_OGs <- read.table("/scratch/evassvis/asb55/data/Aedes_aegypti.fasta.7TMDseqs.names.OGs.fix.sort", sep="\t")
# GPCR <- scran("/scratch/evassvis/mg478/SingleCellProteomes/Phobius/Aedes_aegypti.fasta.7TMDseqs.names", what="", sep="\n")

# get lis of names from table
GPCRs <- GPCR_OGs[,1]
GPCRs

AllCodes <- row.names(Aedes)
AllCodes
# GPCR_droso_irp = unique(grep(paste(GPCRs, collapse="|"), AllCodes, value=TRUE))
# GPCR_droso_irp
# GPCR_FiltList <- intersect(GPCR_droso_irp, GPCRs)
GPCR_FiltList <- intersect(GPCRs, AllCodes)
GPCR_FiltList

pdf("/scratch/evassvis/asb55/Dros_scRNA/figures/GPCR_heatmap.pdf", width = 50, height = 50)
DoHeatmap(Aedes, features = GPCR_FiltList,group.by = "ident", group.bar = TRUE, slot = "scale.data", label = TRUE)
dev.off()

dg=DotPlot(Aedes, features = GPCR_FiltList)
dg
gd=dg$data
gd
write.csv(gd, file = "/scratch/evassvis/asb55/Processed/Aedes/Aedes_GPcr_Dotplot.tsv")














Gpcr_Aedes_filt<- subset(gd,pct.exp> 25 & avg.exp.scaled >0.5)
write.csv(Gpcr_Aedes__filt, file = "/scratch/evassvis/asb55/Processed/Aedes/Aedes_GPcr_Filtered.tsv")


for (i in 1:nrow(Gpcr_Aedes_filt)) {
  match_index <- which(GPCR_OGs[[1]] == Gpcr_Aedes__filt[[i, 3]])
  if (length(match_index) > 0) {
    Gpcr_Aedes_filt[i, 6] <- as.character(GPCR_OGs[match_index, 2])
  }
}
GPCRs
Gpcr_Aedes_filt


write.csv(Gpcr_Aedes_filt, file = "/scratch/evassvis/asb55/Processed/Aedes/Aedes_GPcr_Filtered_with OG.tsv")

















Aedes_GPCR<- read.table("/scratch/evassvis/asb55/aedes/Aedes_GPcr_Filtered_with OG.tsv", sep=",")

Aedes_sero <- c("AAEL08360","AAEL017272","AAEL025125","AAEL027242","AAEL019804")




AAEL08360 <- which(Aedes_GPCR$V4=="AAEL08360")
AAEL017272 <- which(Aedes_GPCR$V4=="AAEL017272")
AAEL025125 <- which(Aedes_GPCR$V4=="AAEL025125")
AAEL027242 <- which(Aedes_GPCR$V4=="AAEL027242")
AAEL019804 <- which(Aedes_GPCR$V4=="AAEL019804")

HT<- c(AAEL08360, AAEL017272,AAEL025125,AAEL027242,AAEL019804)

Aedes_sero_Irp <- Aedes_GPCR[HT,]
Aedes_sero_Irp




write.csv(Aedes_sero_Irp, file = "/scratch/evassvis/asb55/Processed/Aedes_serotonin_.tsv")






dro_dop<- c("AAEL019437","AAEL005834","AAEL00266","AAEL019766")
AAEL019437 <- which(Aedes_GPCR$V4=="AAEL019437")
AAEL005834 <- which(Aedes_GPCR$V4=="AAEL005834")
AAEL00266 <- which(Aedes_GPCR$V4=="AAEL00266")
AAEL019766 <- which(Aedes_GPCR$V4=="AAEL019766")

dop <- c(AAEL019437,AAEL005834,AAEL00266,AAEL019766)

Aedes_dop <- Aedes_GPCR[dop,]
Aedes_dop
write.csv(Aedes_dop, file = "/scratch/evassvis/asb55/aedes/Aedes_dopamine_.tsv")
DotPlot(Aedes, features = Aedes_sero)
DotPlot(Aedes, features = dro_dop)


droso_dopdata=droso_dd$data
droso_dopdata





