library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
dorso.data <- Read10X(data.dir = "/scratch/evassvis/asb55/Dros_scRNA/57K_data")

# Initialize the Seurat object with the raw (non-normalized data).
dorso<- CreateSeuratObject(counts = dorso.data, project = "dorso3k", min.cells = 3, min.features = 200)
dorso
save(dorso, file = "/scratch/evassvis/asb55/Dros_scRNA/DORSOOBJJECT.RDATA")
# load("/scratch/evassvis/asb55/Dros_scRNA/DORSOOBJJECT.RDATA")
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats

pdf("/scratch/evassvis/asb55/Dros_scRNA/figures/vlnplot.pdf", width = 50, height = 50)



# Visualize QC metrics as a violin plot
VlnPlot(dorso, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
dev.off()
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

# plot1 <- FeatureScatter(dorso, feature1 = "nCount_RNA", feature2 = "percent.mt")

pdf("/scratch/evassvis/asb55/Dros_scRNA/figures/FeatureScatter.pdf", width = 50, height = 50)


plot2 <- FeatureScatter(dorso, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()
#plot1 + plot2
plot2
dorso <- subset(dorso, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 )
dorso <- NormalizeData(dorso)


dorso <- FindVariableFeatures(dorso, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(dorso), 10)

# plot variable features with and without labels

pdf("/scratch/evassvis/asb55/Dros_scRNA/figures/VariableFeaturePlot.pdf", width = 50, height = 50)



plot1 <- VariableFeaturePlot(dorso)
dev.off()

pdf("/scratch/evassvis/asb55/Dros_scRNA/figures/Labelpoints.pdf", width = 50, height = 50)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
dev.off()

pdf("/scratch/evassvis/asb55/Dros_scRNA/figures/combined.pdf", width = 50, height = 50)
plot1 + plot2
dev.off()

all.genes <- rownames(dorso)
dorso <- ScaleData(dorso, features = all.genes)
dorso <- RunPCA(dorso, features = VariableFeatures(object = dorso), npcs = 100)
# Examine and visualize PCA results a few different ways
print(dorso[["pca"]], dims = 1:5, nfeatures = 5)



pdf("/scratch/evassvis/asb55/Dros_scRNA/figures/VizDimLoadings.pdf", width = 50, height = 50)

VizDimLoadings(dorso, dims = 1:2, reduction = "pca")
dev.off()
pdf("/scratch/evassvis/asb55/Dros_scRNA/figures/DimPlot.pdf", width = 50, height = 50)

DimPlot(dorso, reduction = "pca")
dev.off()

pdf("/scratch/evassvis/asb55/Dros_scRNA/figures/DimHeatmap.pdf", width = 50, height = 50)
DimHeatmap(dorso, dims = 1, cells = 500, balanced = TRUE)
dev.off()

pdf("/scratch/evassvis/asb55/Dros_scRNA/figures/DimHeatMap.pdf", width = 50, height = 50)
DimHeatmap(dorso, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
dorso <- JackStraw(dorso, num.replicate = 100, dims = 100)
dorso <- ScoreJackStraw(dorso, dims = 1:100)

pdf("/scratch/evassvis/asb55/Dros_scRNA/figures/JackStrawplot", width = 50, height = 50)


JackStrawPlot(dorso, dims = 1:100)
dev.off()

pdf("/scratch/evassvis/asb55/Dros_scRNA/figures/ElbowPlot", width = 50, height = 50)



ElbowPlot(dorso, ndims = 100)
dev.off()

dorso <- FindNeighbors(dorso, dims = 1:100)
dorso <- FindClusters(dorso, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(dorso), 5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pdf("/scratch/evassvis/asb55/Dros_scRNA/figures/Umap", width = 50, height = 50)
dorso <- RunUMAP(dorso, dims = 1:100)
dev.off()
save(dorso, file = "/scratch/evassvis/asb55/Dros_scRNA/DORSOOBJJECTprocessed.RDATA")
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
pdf("/scratch/evassvis/asb55/Dros_scRNA/figures/Dimplot", width = 50, height = 50)
DimPlot(dorso, reduction = "umap")
dev.off()


save(dorso, file = "/scratch/evassvis/asb55/Dros_scRNA/DORSOOBJJECTprocessed.RDATA")
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
dorso.markers <- FindAllMarkers(dorso, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
dorso.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
cluster0.markers <- FindMarkers(dorso, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

dorsofilt <- subset(dorso.markers, avg_log2FC > 1)
dorsofilt
write.csv(dorsofilt, file = "/scratch/evassvis/asb55/Dros_scRNA/Dorso_filtered_MarkerGenes.tsv")

