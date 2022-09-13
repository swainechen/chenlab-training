+++
title = "Add'l-4: Single-cell analysis"
weight = 210
tags = ["S3", "Single-cell", "RCAv2"]
+++

In this tutorial we will be processing and analysing a 10x Genomics PBMC single cell gene expression dataset using the Seurat package to - identify marker genes and major cell types, visualise quality control metrics of single cell data. Additionally, we will also see an example of how an R package, RCAv2, can be used to identify major cell types by projecting single cell gene expression profiles against reference transcriptomes.  

####	1. A 10x Genomics PBMC single cell gene expression dataset can be processed by the [Satija Lab Seurat Guided Clustering](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)
Today we will use a pre-processed Seurat object to illustrate the difference between using the Seurat package along with RCAv2. 

```bash
cd \
mkdir single-cell
cd single-cell
wget https://gisxaws-training.s3.ap-southeast-1.amazonaws.com/pbmc3k_final_withoutRCA2.rds

R

library(dplyr)
library(Seurat)
pbmc_withoutRCA <- readRDS(file = "pbmc3k_final_withoutRCA2.rds")
pbmc_withoutRCA <- RunUMAP(pbmc_withoutRCA, dims = 1:10)
Seurat <- DimPlot(pbmc_withoutRCA, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 5)+
  NoLegend()
 
head(pbmc_withoutRCA)
```
Observe that the pbmc\_withoutRCA file is without annotations prior to running RCAv2
![Seurat\_pbmc](/images/rcav2/pbmc_withoutRCA.png)


####	2. Example of [RCAv2 package](https://github.com/prabhakarlab/RCAv2) for identifying major cell types by projecting single cell gene expression profiles against reference transcriptomes.  
Now let us process the above Seurat object with cell cluster and RCA annotation information.  

```bash
library(RCAv2)
pbmc <- pbmc_withoutRCA
```

1.	Perform major cell type annotation using RCAv2 package on the Seurat object without RCA annotation
```bash
RCA_from_Seurat<-RCAv2::createRCAObject(pbmc@assays$RNA@counts, pbmc@assays$RNA@data)
RCA_from_Seurat<-RCAv2::dataProject(rca.obj = RCA_from_Seurat)
RCA_from_Seurat<-RCAv2::dataClust(RCA_from_Seurat)
RCA_from_Seurat<-RCAv2::estimateCellTypeFromProjection(RCA_from_Seurat)
```  

2.	Add annotations from RCAv2 to Seurat object
```bash
df_RCA_cell <- data.frame(barcode = colnames(RCA_from_Seurat$projection.data),
                          cell = unlist(RCA_from_Seurat$cell.Type.Estimate))
rownames(df_RCA_cell) <- df_RCA_cell$barcode
pbmc$RCA.proj.cell <- df_RCA_cell[colnames(pbmc),"cell"]
```  

3.	Aggregate annotations from RCA2 per Seurat cluster
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

4.	Annotate each cell with its Seurat cluster's RCA annotation
```bash
df_umap_cluster_RCA2 <- FetchData(pbmc, vars = c("UMAP_1","UMAP_2","seurat_clusters"))
df_umap_cluster_RCA2$cell <- NA

for (i in 0:(length(vec_RCA_cluster_annotations)-1)){
  df_umap_cluster_RCA2[which(df_umap_cluster_RCA2$seurat_clusters==i),"cell"] <- vec_RCA_cluster_annotations[i+1]
}

pbmc$RCA.annotation <- df_umap_cluster_RCA2[colnames(pbmc),"cell"]

head(pbmc)
```
Observe the pbmc\_withRCA file which now has RCA annotations
![SeuratRCA\_pbmc](/images/rcav2/pbmc_withRCA.png)

5.	Visualise UMAP with RCA annotations
```bash
RCA <- DimPlot(pbmc, reduction = "umap", group.by = "RCA.annotation", label = TRUE, label.size = 5)+
  NoLegend()
gridExtra::grid.arrange(Seurat,RCA,ncol=1)
q()
# type in 'n' in response to "Save workspace image?"

mv Rplots.pdf SeuratVSRCA.pdf
```
![Seurat\_pbmc](/images/rcav2/SeuratVSRCA.png)

You can now inspect the plots in SeuratVSRCA.pdf by transferring this file to your local machine. Refer to [Local machine to EC2 instance and vice-versa - Step4](http://slchen-lab-training.s3-website-ap-southeast-1.amazonaws.com/12-creates3sharedata/07-datatransfer.html)


>_**Optional**: Now let us review single cell datasets and try few more visualizations using a pre-processed Seurat object with RCA annotation information   
	1) How to identify marker genes in single cell datasets   
	2) Visualisation of marker gene expression in single cell datasets   
	3) Visualisation of major cell type annotations made using RCA2_   

```bash
wget https://gisxaws-training.s3.ap-southeast-1.amazonaws.com/pbmc3k_final_withRCA2.rds

R

library(dplyr)
library(Seurat)
pbmc <- readRDS(file = "pbmc3k_final_withRCA2.rds")

### Identify marker genes for every cluster compared to all other cells
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

### Plot raw counts of marker genes of interest, and visualise marker gene expression, and major cell type annotations made using RCA2
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
FeaturePlot(pbmc, features = c("MS4A1", "CD3E", "CD4", "CD8A",
                               "GNLY", "CD14", "FCGR3A", "LYZ",
                               "FCER1A", "PPBP"))
dev.off()
q()
```   

