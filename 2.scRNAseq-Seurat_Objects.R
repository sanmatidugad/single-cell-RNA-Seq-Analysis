# 4. Create Seurat object ----
datadir <-  'counts_filtered'
list.files(datadir)

expression_matrix <- Read10X(datadir, gene.column = 2, cell.column = 1, unique.features = TRUE, strip.suffix = FALSE)

# hc.1k.seurat <- CreateSeuratObject(counts = expression_matrix, min.cells = 3)  %>% 
#   NormalizeData(verbose = FALSE) %>% 
#   FindVariableFeatures(verbose = FALSE)

hc.1k.seurat <- CreateSeuratObject(counts = expression_matrix, 
                                   min.cells = 3, min.features = 300,
                                   project = "1K_HCM")
# 1K Heart Cells in Mouse

# Calculating percent of mitochondrial reads
# NOTE: change 'mt' to 'MT' in human for mito genes, change R(pl|n)[0-9] to RP[S|L] in human for Ribo genes

hc.1k.seurat[["percent.mt"]] <- PercentageFeatureSet(object = hc.1k.seurat, pattern = "^mt-") 
hc.1k.seurat[["percent.ribo"]] <- PercentageFeatureSet(object = hc.1k.seurat, pattern = "^R(pl|n)[0-9]") 
#some other quality control features are, ERCC genes, housekeeping genes.

# in the violin plot, features = genes detected, while counts = total molecules detected
# nCount_RNA - number of UMIs per cell, each point is a cell, It is basically the Column sum for each cell.
# nCount_RNA is just Summing the counts of all RNA molecules for each cell across all genes. 
# nFeature_RNA - number of genes, each point is a cell - nrows()
# percent.mito - % reads in each cell that mapped to mitochondria. (mitochondrial genes) 

# Make violin plot
vln_plot = VlnPlot(hc.1k.seurat, c("nCount_RNA", "nFeature_RNA", "percent.mt","percent.ribo"),
                   pt.size = 0.1, ncol = 2) + NoLegend()

vln_plot

# Filter your data
hc.1k.seurat <- subset(hc.1k.seurat, subset = 
                         nCount_RNA < 60000  & 
                         nCount_RNA > 500 & 
                         nFeature_RNA > 500 & 
                         percent.mt < 20 &
                         percent.ribo <  10)


# Observe Data After Filtering.
VlnPlot(hc.1k.seurat, c("nCount_RNA", "nFeature_RNA", "percent.mt","percent.ribo"), pt.size = 0.1, ncol = 2)

# NOTE: you need to be careful when setting cut-offs that you're not losing unique cell populations


## 5. Expression Normalization and Find Highly Variable Genes.
##    After removing unwanted gene cells from the data, next step is to normalize the data.
##    Log Normalize normalizes each gene expression for each cell by the total expression and multiplies this by a scale factor of 10,000 by default/
##    This data is then log transformed. LogNormalize also stabilizes the variance.
##    Other methods of normalization are CLR - That applies a Central Log ratio transformation,
##            RC - Relative Counts - Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor.
##                                   No log-transformation is applied. For counts per million (CPM) set scale.factor = 1e6


hc.1k.seurat = NormalizeData(hc.1k.seurat, normalization.method = "LogNormalize", scale.factor = 1e4) %>%
  FindVariableFeatures(selection.method = "vst", verbose = TRUE)

top10 = head(VariableFeatures(hc.1k.seurat), 10)

plot1 = VariableFeaturePlot(hc.1k.seurat)
plot2 = LabelPoints(plot = plot1, points = top10, repel=TRUE, xnudge = 0, ynudge = 0)

plot1+ plot2

# another QA plot
ggplot(hc.1k.seurat@meta.data, aes(nCount_RNA, nFeature_RNA)) +
  geom_point(alpha = 0.7, size = 0.5) +
  labs(x = "Total UMI counts per cell", y = "Number of genes detected")

# Potential things to look for in the type of QA plot produced above:
# 1. Data points in the bottom LEFT hand quadrant = low genes and UMIs per cell. May represent poor quality cells.
# 2. Data points in the bottom RIGHT hand quadrant = low genes but high UMIs per cell. These could be dying cells, but also could represent a population of a low complexity celltype (i.e red blood cells).

# Visualizing Independent Component Analysis in different ways. (ICA) - Linear Dimension Reduction
# Plot UMAP ----
# it is standard practice to apply a linear transformation ('scaling') before PCA. For single cell data this includes:
# 1. Shifting the expression of each gene, so that the mean expression across cells is 0
# 2. Scaling the expression of each gene, so that the variance across cells is 1
# This gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
hc.1k.seurat <- ScaleData(hc.1k.seurat, verbose = FALSE)
hc.1k.seurat <- RunPCA(hc.1k.seurat, npcs = 40, verbose = FALSE)
hc.1k.seurat <- RunUMAP(hc.1k.seurat, reduction = "pca", dims = 1:40)
hc.1k.seurat <- RunTSNE(hc.1k.seurat, reduction.use = "pca", dims.use = 1:50, perplexity=20)

hc.1k.seurat <- FindNeighbors(hc.1k.seurat, reduction = "pca", dims = 1:40)
hc.1k.seurat <- FindClusters(hc.1k.seurat, resolution = 0.5)

DimPlot(hc.1k.seurat, reduction = "umap", split.by = "orig.ident", label = TRUE)
#Plot Heatmap
DimHeatmap(hc.1k.seurat, dims = 1:9, cells = 100, reduction = "pca", balanced = TRUE)

#  Jackstraw Plot - The process involves plotting the significance of each gene's association with each principal component (using QQ-plots)
# and then calculating a p-value to determine if each PC explains a significant amount of variation in the dataset.
# This helps identify which PCs are most relevant for understanding the underlying patterns in the data.

hc.1k.seurat = JackStraw(hc.1k.seurat, reduction = "pca")
hc.1k.seurat = ScoreJackStraw(hc.1k.seurat, dims = 1:20)

JackStrawPlot(hc.1k.seurat, dims = 1:40)
PCASigGenes(hc.1k.seurat, pcs.use = 1, pval.cut = 0.001)[1:50]

## Also known as Cluster bend plot - Used to determine the optimal number of clusters in a clustering analysis
ElbowPlot(hc.1k.seurat, ndims = 30, reduction = "pca")
