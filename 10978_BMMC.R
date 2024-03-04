library(dplyr)
library(Seurat)
library(patchwork)

# 1. Open Data and create Seurat Object
# Download data, decompress with Zip7 and store in your desired folder

# Load the PBMC dataset
pbmc.data <- ReadMtx(mtx = "C:/Users/Admin/Desktop/R-script/10978_BMMC/Raw_Data/matrix.mtx",
                     features = "C:/Users/Admin/Desktop/R-script/10978_BMMC/Raw_Data/features.tsv",
                     cells =  "C:/Users/Admin/Desktop/R-script/10978_BMMC/Raw_Data/barcodes.tsv")


# To check, if the cell count matrix is correct, look at the first 10 columns and rows
# columns = barcodes ; rows = features/genes
pbmc.data[1:10,1:10]

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "10978_PBMC", min.cells = 10, min.features = 200)
# An object of class Seurat
# 18680 features across 17228 samples within 1 assay
# Active assay: RNA (18680 features, 0 variable features)
# 1 layer present: counts
pbmc
str(pbmc)

# Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

# 5064020568 bytes
dense.size <- object.size(as.matrix(pbmc.data))
dense.size

# 316769016 bytes
sparse.size <- object.size(pbmc.data)
sparse.size

# 16 bytes
dense.size/sparse.size

# To check how your meta data looks right now
# matric, with columns = number of counts and number features; and rows = barcodes/individual cells
View(pbmc@meta.data)

# 2. Pre-Processing the Data

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# So first add the percentage of mitochondrial genes to our metadata
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot to check quality of cells before filtering
# Export as 10978_PBMC_PreQC
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Idents(pbmc) <- "seurat_clusters"
# pbmc <- subset(pbmc, idents = "6", invert = T)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
# nFeature to nCount should follow a straight line, so that a good number of genes is sequenced to a satisfying amount
# Export as 10978_PBMC_PreQC2
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# After this preliminary analysis, we use these data to filter cells
# We filter cells that have unique feature (=gene) counts over 5000 (or less than 200) (more seems like duplicates, less like empty droplets)
# We filter cells that have >10% mitochondrial counts (dying cells have abberant high mitochondrial gene expression)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# visualize subset now
# An object of class Seurat 
# 18680 features across 17025 samples within 1 assay 
# Active assay: RNA (18680 features, 0 variable features)
# 1 layer present: counts
pbmc

# 3. Normalization of Data (feature expression measurements for each cell by the total expression)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#These settings are the default settings, which can also be achieved by this command
pbmc <- NormalizeData(pbmc)

# To check out which commands you have already performed on your Seurat Object
str(pbmc)

# Identify highly variable features (genes, which have a high cell to cell variance)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# 4. Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# See the names of the 10 most variable genes
top10

# plot variable features with and without labels
# Export as 10978_PBMC_HighlyVariableFeatures
plot1 <- VariableFeaturePlot(pbmc)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1

# 5. Scaling the Data
# Shifts the expression of each gene, so that the mean expression across cells is 0
# Scales the expression of each gene, so that the variance across cells is 1
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
# The results of this are stored in pbmc[["RNA"]]$scale.data
# By default, only variable features are scaled.
# You can specify the features argument to scale additional features
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# 6. Perform linear dimensional reduction (creates a list with genes with the most positive and negative loadings)
# By default it just covers the 2000 variable features which we determined before
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA (Principle Component Analysis) results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca") + NoLegend()
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the dimensionality of the dataset
ElbowPlot(pbmc)

# 7. Clustering the cells (dims = number of PC to include in downstream analysis)
pbmc <- FindNeighbors(pbmc, dims = 1:15)

# Understanding the resolution (how many cluster at which resolution)
pbmc <- FindClusters(pbmc, resolution = c(0.1,0.3,0.5,0.7,1))

# Find clusters can also be performed in a different resolution (0.1 - 1), the lower the number the fewer clusters
# To look at the metadata and see the different resolutions in columns to choose which one works best for you
View(pbmc@meta.data)

# Visualize data reduction plot at different resolutions to identify which one fits best
# Avoid overlapping clusters in higher resolutions
DimPlot(pbmc, group.by = "RNA_snn_res.0.1", label = TRUE)

# Plot a ClusTree to visualize the clustering of your cells
ClusTree <- clustree(pbmc)
ClusTree

# To set identity of clusters (just do this, if you want to change resolution)
Idents(pbmc) <- "RNA_snn_res.0.1"

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

# 8. Non-dimensional reduction (tSNE/UMAP)
set.seed(0)
pbmc <- RunUMAP(pbmc, dims = 1:20)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE, label.size = 6)

# Save this object from here, that it can easily be opened from here
saveRDS(pbmc, file = "C:/Users/Admin/Desktop/R-script/10978_PBMC/Output/10978_PBMC_UMAP.rds")

# 9. Finding differentially expressed features (cluster biomarkers)

# Activate presto package for faster DE (differential expression) analysis
library(presto)

# The default assay needs to be RNA, check this by the following code
DefaultAssay(pbmc)

# Find all markers (that takes long and is not necessary)
FindAllMarkers(pbmc, logfc.threshold = 0, min.pct = 0.1, only.pos = TRUE, test.use = 'roc', slot = 'counts')

# find all markers of a specific cluster compared to all other cluster
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
# visualize the expression of the top 2 markers of this cluster throughout all cluster
VlnPlot(pbmc, features = c(row.names(cluster2.markers)[1], row.names(cluster2.markers)[2]))

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

#Visulaize the top10 markers of each cluster
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n=10, wt = avg_log2FC)
top10
top10_cluster_2 <- top10 %>% filter(cluster == 2)

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.0, test.use = "roc", only.pos = TRUE)

# Plotting a DotPlot of the gene expression of the features of interest throughout the different cluster
# First, Define the genes you want to plot
genes_of_interest <- c("CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "NCAM1", "KLRB1", "NKG7", "CD19", "MS4A1", "CD38", "CD14", "CD68", "FCGR3A", "CD1C", "LILRA4", "PPBP", "HBB", "CD34")

# Create dot plot for the selected genes across cell clusters
dotplot <- DotPlot(object = pbmc, features = genes_of_interest)
dotplot + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot expression probability distributions of individual gene markers across clusters
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

# visualizes feature expression on a tSNE or PCA plot
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

# DoHeatmap() generates an expression heatmap for given cells and features.
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# 10. Assigning cell type identity to clusters
# Known markers for specific cell types
#   Cluster ID 	Markers 	    Cell Type
#   0 	        IL7R, CCR7 	  Naive CD4+ T
#   1 	        CD14, LYZ 	  CD14+ Mono
#   5 	        IL7R, S100A4 	Memory CD4+
#   3 	        MS4A1 	      B
#   4 	        CD8A 	CD8+    T
#   7 	        FCGR3A, MS4A7 FCGR3A+ Mono
#   2 	        GNLY, NKG7 	  NK
#   8 	        FCER1A, CST3 	DC
#   6 	        PPBP 	        Platelet

#renaming clusters just individually because you already know the names (Resolution 0.1)
new.cluster.ids <- c("CD4+ T", "CD14+ Mono", "NK", "Mix", "B", "CD8+ T", "CD16+ Mono", "Megakaryocytes", "DC", "Low Quality")
#renaming clusters just infividually because you alread know the names (Resolution 0.7)
new.cluster.ids <- c("Naive CD4+ T", "CD14+ Mono1", "CD14+ Mono2", "Effector Memory CD4+ T", "Mixed", "Central Memory CD4+ T", "CD56 high NK", "Naive B", "CD16+ Monocytes", "CD16 high NK", "Memory CD8+ T", "Naive CD8+ T", "Memory B", "Megakaryocytes", "Plasmacytoid DC", "Myeloid DC", "?", "??", "???", "Proliferating T", "Mixed2", "HPSCs", "Mixed3")

names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE, label.size = 6, pt.size = 0.5)

#Plotting and saving plot as image
library(ggplot2)
plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "../10978_PBMC_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)

#Saving final document
saveRDS(pbmc, file = "C:/Users/Admin/Desktop/R-script/Seurat Tutorial JB/Output/pbmc3k_final.rds")