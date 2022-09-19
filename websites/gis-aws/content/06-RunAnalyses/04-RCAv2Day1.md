+++
title = "c. Single Cell RNA-Seq Data Summary"
date = 2022-09-19T09:05:54Z
weight = 50
tags = ["RCAv2", "RNA-seq", "single-cell"]
+++

Recent advances in sequencing technology has allowed for single-cell transcriptomics, in which the RNA from a single cell can be quantified and measured. This allows for the identification of cell type simply from the levels of marker gene expression. Day 2 will involve a [more in-depth analysis](http://slchen-lab-training.s3-website-ap-southeast-1.amazonaws.com/12-creates3sharedata/08-singlecell.html) of single-cell data; for this exercise, we will just get used to using R in the AWS EC2 terminal. 

First, let's make a directory for this exercise and download some data:

```bash 
cd ~
mkdir explore 
cd explore

aws s3 cp s3://gisxaws-training/pbmc3k_final_withoutRCA2.rds .
```

The analysis is run with built-in R, a different language from what you have been using in the rest of the tutorial (bash) and those commands will not work. This session will require a AWS machine of 8GB RAM.

Next, start the R environment:
```bash
R
```

Then, we will load our libraries and data file:
```R
library(dplyr)
library(Seurat)

data <- readRDS(file = "pbmc3k_final_withoutRCA2.rds")
```

We will study this dataset in-depth later, but for now familiarise yourself with navigating R. R has many strengths, including statistical analyses and data visualisation. There is a graphic user interface known as Rstudio, but for this workshop we will be using our terminals only.

First, familiarise yourself with this data- what does the data structure look like? How large is it?:

```R
# Let's take a look:
head(data)

nrow(data)

colnames(data)

class(data)

typeof(data)
```

In this file, the columns represent the follows:

**orig.ident**          -> The name of the original dataset  
**nCount\_RNA**         -> The number of UMI reads per cell  
**nFeature\_RNA**       -> The number of detected genes per same cell  
**percent.mt**          -> Percent of Features that are mitochondial  
**RNA\_snn\_res.0.5**   -> Cluster ID based on SNN, with a resolution of 0.5  
**seurat\_clusters**    -> Seurat cluster ID   

Each row represents a single feature- a single molecule of genetic information.  
We can also see that this is a _Seurat Object_, and is S4 class. 

Let's take a look at the number of detected genes per same cell, and see how the different clusters compare:

```R

```

We will continue our analysis of this file in Day 2. If you have remaining time, please feel free to explore the [Seurat tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) hosted by the Satija Lab. Here, you can download the raw data used to generate this file and explore the different QC measures used to process and cluster the 10x sequencing data.
