library(tidyverse)
library(Seurat)
library(patchwork)
library(tximport)
library(EnsDb.Hsapiens.v79)
library(biomaRt)
library(SeqGSEA)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)

#import alevin data, report unfiltered cell and gene counts
txi <- tximport('/projectnb/bf528/users/vangogh2022/project_4/salmon.output/alevin/quants_mat.gz', type="alevin")

#create Seurat object
seurat_obj <- CreateSeuratObject(counts = txi$counts, min.cells=3, min.features = 200)
unfiltered_num_genes <- seurat_obj@assays$RNA@counts@Dim[2]
unfiltered_num_cells <- seurat_obj@assays$RNA@counts@Dim[1]

#add mitocondrial genes and log10genesperUMI
seurat_obj[['percent.mt']] <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")
seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)

#visualize filtering parameters
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave('violin.png', device='png', bg='transparent')
plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(seurat_obj, feature1 = 'nFeature_RNA', feature2 = 'percent.mt')
plot4 <- seurat_obj@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
figure <- ggarrange(plot1, plot2, plot3, labels=c('A', 'B', 'C'), ncol = 3)
ggsave('featurescatter.png', plot=figure, device='png', bg='transparent', width = 18, height = 8)
plot4

#Filter cells
seurat_obj <- subset(seurat_obj, 
                     subset = (nFeature_RNA > 200) & (nFeature_RNA < 4000) & (nCount_RNA > 500) & (percent.mt < 20) & (log10GenesPerUMI > 0.85))

#filter genes
# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = seurat_obj, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
seurat_obj <- CreateSeuratObject(filtered_counts, meta.data = seurat_obj@meta.data)

filtered_num_cells <- seurat_obj@assays$RNA@counts@Dim[2]
filtered_num_genes <- seurat_obj@assays$RNA@counts@Dim[1]

#Normalize
seurat_obj <- NormalizeData(seurat_obj)

#Find variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(seurat_obj)

#Scale the data
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

#dimension reduction
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")
DimPlot(seurat_obj, reduction = "pca")
DimHeatmap(seurat_obj, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat_obj, dims = 1:4, cells = 500, balanced = TRUE)

#determine 'dimensionality' of the dataset
seurat_obj <- JackStraw(seurat_obj, num.replicate = 100)
seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:20)
JackStrawPlot(seurat_obj, dims = 1:20)
elbow <- ElbowPlot(seurat_obj)
ggsave('elbow.png', plot=elbow, device='png', bg='transparent')

#clustering and pie chart
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:8)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

saveRDS(seurat_obj, file = "seurat.rds")


seurat_obj <- readRDS('seurat.rds')

clusters <- as.data.frame(table(Idents(seurat_obj)))
clusters$percent <- clusters$Freq / sum(clusters$Freq) * 100 
clusters$percent <- formatC(clusters$percent, digits=2, format='f') 
clusters$percent <- paste(clusters$percent, '%', sep='')

colnames(clusters)[1] <- 'Cluster'

myPalette <- brewer.pal(12, "Paired") 

png('pierchart.png', width=600, height = 400)
pie_chart <- pie(clusters$Freq, labels=clusters$percent, col = myPalette, border = 'white') +
  legend(1.5, .7, clusters$Cluster, cex = .7, fill = myPalette, title='Cluster') +
  par(mar=c(5,4,4,7)+0.3)
dev.off()


#save filtered and unfiltered gene/cell counts in dataframe

df <- data.frame('rows' = c('Before Filtering', 'After Filtering'),
                 'Cells'=c(unfiltered_num_cells, filtered_num_cells), 
                 'Genes'=c(unfiltered_num_genes, filtered_num_genes),
                 row.names = 'rows')
write.csv(df, 'cell_gene_counts.csv', row.names = TRUE)
