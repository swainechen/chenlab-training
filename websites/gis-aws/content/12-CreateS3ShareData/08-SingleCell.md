+++
title = "Add'l-4: Single-cell analysis"
weight = 210
tags = ["S3", "Single-cell", "RCAv2"]
+++

PAGE UNDER UPDATION....

In this tutorial we will be processing and analysing a 10x Genomics PBMC single cell gene expression dataset using the Seurat package to - identify marker genes and major cell types, visualise quality control metrics of single cell data. Additionally, we will also see an example of how an R package, RCAv2, can be used to identify major cell types by projecting single cell gene expression profiles against reference transcriptomes.

####	1. Processing a 10x Genomics PBMC single cell gene expression dataset by Seurat ([Adapted from Satija Lab Seurat Guided Clustering](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html))

1.	Download, extract and load the 10x Genomics PBMC single cell gene expression dataset
```bash
cd \
mkdir single-cell
cd single-cell

wget https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar zxvf pbmc3k_filtered_gene_bc_matrices.tar.gz
```

```bash
R

library(dplyr)
library(Seurat)
library(RCAv2)

pbmc.data <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")
```

2.	Create a Seurat object containing the raw data
```bash
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
```

3.	Add percent of mitochondrial reads per cell as a quality control (QC) filter
```bash
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
```

4.	Visualise QC metrics using violin plots
```bash
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

5.	Visualise QC metrics using scatter plots
```bash
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```

6.	Subset data based on nFeature_RNA (number of detected genes, NODG) and pMito metrics
```bash
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```

7.	Log-normalise data
```bash
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
```

8.	Perform feature selection (2,000 highly variable genes)
```bash
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
```

9.	Identify the 10 most highly variable genes - note how these correspond to marker genes
```bash
top10 <- head(VariableFeatures(pbmc), 10)
```

10.	Plot variable features
```bash
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
```

11.	Scale data using all genes
```bash
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
```

12.	Perform PCA on scaled data using variable features
```bash
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
```

13.	Visualise PCA results
```bash
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
```

14.	Determine number of PCs to use in following analyses
```bash
ElbowPlot(pbmc)
```

15.	Cluster cells
```bash
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
```

16.	Perform UMAP (non-linear dimensional reduction) analysis and visualisation
```bash
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
dev.off()
q()
```
```bash
mv Rplots.pdf Seurat.pdf
```
You can now inspect the results by transferring this file to your local machine. Refer to [Local machine to EC2 instance and vice-versa - Step4](http://slchen-lab-training.s3-website-ap-southeast-1.amazonaws.com/12-creates3sharedata/07-datatransfer.html)

####	2. Example of [RCAv2 package](https://github.com/prabhakarlab/RCAv2) for identifying major cell types by projecting single cell gene expression profiles against reference transcriptomes.

Now let us download an RDS that contains a pre-processed Seurat object with cell cluster and RCA annotation information. 

```bash
wget https://gisxaws-training.s3.ap-southeast-1.amazonaws.com/pbmc3k_final_withRCA2.rds
``` 

1.	 Identify marker genes for every cluster compared to all other cells
```bash
library(dplyr)
library(Seurat)

pbmc <- readRDS(file = "pbmc3k_final_withRCA2.rds")
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
```

2.	 Plot raw counts of marker genes of interest, and visualise marker gene expression
```bash
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
FeaturePlot(pbmc, features = c("MS4A1", "CD3E", "CD4", "CD8A",
                               "GNLY", "CD14", "FCGR3A", "LYZ", 
                               "FCER1A", "PPBP"))
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
```

3.	Perform major cell type annotation using RCA2
```bash
RCA_from_Seurat<-RCAv2::createRCAObject(pbmc@assays$RNA@counts, pbmc@assays$RNA@data)
RCA_from_Seurat<-RCAv2::dataProject(rca.obj = RCA_from_Seurat)
RCA_from_Seurat<-RCAv2::dataClust(RCA_from_Seurat)
RCA_from_Seurat<-RCAv2::estimateCellTypeFromProjection(RCA_from_Seurat)
```

4.	Add annotations from RCA2 to Seurat object
```bash
df_RCA_cell <- data.frame(barcode = colnames(RCA_from_Seurat$projection.data),
                          cell = unlist(RCA_from_Seurat$cell.Type.Estimate))
rownames(df_RCA_cell) <- df_RCA_cell$barcode
pbmc$RCA.proj.cell <- df_RCA_cell[colnames(pbmc),"cell"]
```

5.	Aggregate annotations from RCA2 per Seurat cluster
```bash
df_RCA_results <- FetchData(pbmc, vars = c("seurat_clusters","RCA.proj.cell"))
df_RCA_results$count <- 1
df_RCA_results <- aggregate(count ~ ., df_RCA_results, FUN = sum)
df_RCA_results <- df_RCA_results %>% group_by(seurat_clusters) %>% top_n(n = 5, wt = count)
df_RCA_results <- df_RCA_results[order(df_RCA_results$seurat_clusters),]

vec_RCA_cluster_annotations <- (df_RCA_results %>% filter(count == max(count)))$RCA.proj.cell

for (i in 0:(length(vec_RCA_cluster_annotations)-1)){
  print(df_RCA_results %>% filter(seurat_clusters == i))
}
```

6.	Annotate each cell with its Seurat cluster's RCA2 annotation
```bash
df_umap_cluster_RCA2 <- FetchData(pbmc, vars = c("UMAP_1","UMAP_2","seurat_clusters"))
df_umap_cluster_RCA2$cell <- NA

for (i in 0:(length(vec_RCA_cluster_annotations)-1)){
  df_umap_cluster_RCA2[which(df_umap_cluster_RCA2$seurat_clusters==i),"cell"] <- vec_RCA_cluster_annotations[i+1]
}

pbmc$RCA.annotation <- df_umap_cluster_RCA2[colnames(pbmc),"cell"]
```

7.	Visualise UMAP with RCA2 annotations
```bash
DimPlot(pbmc, reduction = "umap", group.by = "RCA.annotation", label = TRUE, label.size = 5)+
  NoLegend()
q()
```
You can once again inspect the results by transferring this file to your local machine.
