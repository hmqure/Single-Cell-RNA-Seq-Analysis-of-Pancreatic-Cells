# Load Packages
library(tidyverse)
library(Seurat)
library(dplyr)

setwd("/projectnb/bf528/users/vangogh2022/project_4/analyst/")

# Load Data
changed <- readRDS("/projectnb/bf528/users/vangogh2022/project_4/analyst/seurat.rds")

# Identify mks
changed_mks <- Findsummks(changed, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_mks <- changed_mks %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
view(changed_mks)
view(top_mks)

# Marker genes used in original analysis
marker_genes <- read_csv("Marker_Genes.csv")
sumgenes <- unlist(strsplit(marker_genes$Genes, ","))

#UMAP
cell <- RunUMAP(changed, dims = 1:10)
DimPlot(cell, reduction = "umap", label = TRUE)

# Heatmap
png("heatmap_clusters.png", height = 800, width = 1000)
DoHeatmap(changed, features = top_mks$gene) + NoLegend()
dev.off()

# Clusters
# clusters (0,2,3), (1,4), (5,9), and (7,9) overl
panglao <- read_tsv("https://panglaodb.se/mks/PanglaoDB_mks_27_Mar_2020.tsv.gz") %>%  filter(str_detect(species,"Hs")) %>% select(c("official gene symbol", "cell type", "organ"))

view(panglao)

# Violin Plots 
pdf("violins.pdf")
violin(changed, features = c("SST"))
violin(changed, features = c("GCG")) 
violin(changed, features = c("CPA1"))
violin(changed, features = c("KRT19")) 
violin(changed, features = c("PDGFRB")) 
violin(changed, features = c("INS"))
violin(changed, features = c("PPY")) 
violin(changed, features = c("VWF", "PECAM1","CD34")) 
violin(changed, features = c("SDS", "CD163")) 
dev.off()

# Assigning cell types
new_cluster_ids <- c("Delta", "Alpha", "Alpha", "Acinar", "Ductal", "Stellate", "Beta", "Gamma", "Beta", "Acinar", "Vascular", "Macrophage")
names(new_cluster_ids) <- levels(changed)
changed <- RenameIdents(changed, new_cluster_ids)
DimPlot(changed, reduction = "umap", label = TRUE, pt.size = 0.5)

# Clustered UMAP
cell <- RunUMAP(changed, dims = 1:10)
DimPlot(cell, reduction = "umap", label = TRUE)

png("features.png", width = 1000, height = 1000)
FeaturePlot(cell, features = sumgenes)
dev.off()

# Clustered Heatmap
png("heatmap_celltype.png", height = 900, width = 1500)
DoHeatmap(changed, features = top_mks$gene) + NoLegend()
dev.off()

# Save Seurat Object
saveRDS(changed, file = "analyst.rds")

# Save Marker Genes
write_csv(changed_mks, "marker_genes_final.csv")
deseq_changed <- Findsummks(changed, test.use = "DESeq2", max.changed.per.ident = 50)
