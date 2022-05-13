***Bambu*** is a R package for multi-sample transcript discovery and quantification using long read RNA-Seq data. You can use ***Bambu*** after read alignment to obtain expression estimates for known and novel transcripts and genes. The output from ***Bambu*** can directly be used for visualisation and downstream analysis such as differential gene expression or transcript usage.

Goals:
- Count the number of different full length transcripts int he sample
- Detect novel transcripts and genes in the samples
- Identify differentially expressed transcripts
- Detect alternate isoform usage of novel genes

Software Requirements:
 - R packages: [***Bambu***](http://bioconductor.org/packages/bambu/), [***DeSeq2***](http://bioconductor.org/packages/DESeq2/), [***DEXseq***](http://bioconductor.org/packages/DEXSeq/)
 - Other softwares : [***Minimap2***](https://github.com/lh3/minimap2), [***Samtools***](http://www.htslib.org/), [***Gffcompare***](https://github.com/gpertea/gffcompare)

In terms of machine specification, for this analysis, please spin off an AWS machine of at least 8 Gb RAM. 

### Installation
The AMI you were provided should have all the requirements installed. But in the case if you want to install locally after the tutorial, we have also included below on how ***Bambu*** can be installed.

***Bambu*** is available through GitHub and Bioconductor

Downloading through Bioconductor:
```rscript
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("bambu")
```

Downloading directly from GitHub:
First download the repository ***Bambu*** from github
```bash
git clone https://github.com/GoekeLab/bambu.git
```
Within R, load the package (be sure to put the correct path to ***Bambu*** below)
```rscript
devtools::load_all("/path/to/bambu")
library(bambu)
```

### Downloading reference genome and annotations

The default mode to run ***Bambu*** is using a set of aligned reads (bam files), reference genome annotations (gtf file, TxDb object, or bambuAnnotation object), and reference genome sequence (fasta file or BSgenome). In this tutorial we will be using data from human cancer cell-lines from the [SG-NEx project](https://github.com/GoekeLab/sg-nex-data). 

```bash
# create work directory
mkdir bambu_training
cd bambu_training

# download genome fasta file 
aws s3 cp --no-sign-request s3://sg-nex-data/data/annotations/genome_fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa ./

# download gtf file
aws s3 cp --no-sign-request s3://sg-nex-data/data/annotations/gtf_file/Homo_sapiens.GRCh38.91.sorted.gtf ./
```

### Downloading and aligning reads

The SG-NEx project provides aligned reads that can be used directly with Bambu, however sometimes, raw fastq reads might also be provided instead of aligned reads. In that case, alignment can be done as follows:

Note that we are aligning the reads to the genome fasta and not the transcriptome fasta as ***Bambu*** uses intron junctions to distingush novel transcripts.

```bash
# download fastq from s3 bucket
aws s3 cp --no-sign-request s3://sg-nex-data/data/sequencing_data/fastq/SGNex_H9_directRNA_replicate1_run1/SGNex_H9_directRNA_replicate1_run1.fastq.gz ./

# align using Minimap2 
minimap2 -ax splice -uf -k14 Homo_sapiens.GRCh38.dna.primary_assembly.fa SGNex_H9_directRNA_replicate1_run1.fastq.gz > SGNex_H9_directRNA_replicate1_run1.sam  #needs more than 4gb ram
samtools view -Sb SGNex_H9_directRNA_replicate1_run1.sam | samtools sort -o SGNex_H9_directRNA_replicate1_run1.bam
samtools index SGNex_H9_directRNA_replicate1_run1.bam
```
Depending on your data source, you may need to run minimap2 with different paramters. See below taken from the ***Minimap2*** [documentation](https://github.com/lh3/minimap2)

```bash
minimap2 -ax splice ref.fa rna-reads.fa > aln.sam       # spliced long reads (strand unknown)  
minimap2 -ax splice ref.fa nanopore-cdna.fa > aln.sam   # Nanopore 2D cDNA-seq              
minimap2 -ax splice -uf -k14 ref.fa reads.fa > aln.sam  # noisy Nanopore Direct RNA-seq            
minimap2 -ax splice:hq -uf ref.fa query.fa > aln.sam    # Final PacBio Iso-seq or traditional cDNA        
```


To save time we will download the aligned bam files for the remaining samples from the S3 bucket.

```bash
# download aligned bam files for H9 samples 
aws s3 cp s3://sg-nex-data/data/sequencing_data/bam/GIS_H9_directRNA_Rep1_Run2_R1.bam ./
aws s3 cp s3://sg-nex-data/data/sequencing_data/bam/GIS_H9_directRNA_Rep1_Run2_R1.bam.bai ./
aws s3 cp s3://sg-nex-data/data/sequencing_data/bam/GIS_H9_directRNA_Rep2_Run2_R1.bam ./
aws s3 cp s3://sg-nex-data/data/sequencing_data/bam/GIS_H9_directRNA_Rep2_Run2_R1.bam.bai ./
aws s3 cp s3://sg-nex-data/data/sequencing_data/bam/GIS_H9_directRNA_Rep3_Run2_R1.bam ./
aws s3 cp s3://sg-nex-data/data/sequencing_data/bam/GIS_H9_directRNA_Rep3_Run2_R1.bam.bai ./

# download aligned bam files for HEYA8 samples 
aws s3 cp s3://sg-nex-data/data/sequencing_data/bam/GIS_HEYA8_directRNA_Rep1_Run2_R1.bam ./
aws s3 cp s3://sg-nex-data/data/sequencing_data/bam/GIS_HEYA8_directRNA_Rep1_Run2_R1.bam.bai ./
aws s3 cp s3://sg-nex-data/data/sequencing_data/bam/GIS_HEYA8_directRNA_Rep2_Run2_R1.bam ./
aws s3 cp s3://sg-nex-data/data/sequencing_data/bam/GIS_HEYA8_directRNA_Rep2_Run2_R1.bam.bai ./
aws s3 cp s3://sg-nex-data/data/sequencing_data/bam/GIS_HEYA8_directRNA_Rep3_Run1_R1.bam ./
aws s3 cp s3://sg-nex-data/data/sequencing_data/bam/GIS_HEYA8_directRNA_Rep3_Run1_R1.bam.bai ./

```
### Quality Control

[***Samtools stats***](http://www.htslib.org/doc/samtools-stats.html) can be used to quickly check the read accuracy and read mappability in the bam files, in the example below, we checked the aligned bam file for the **GIS_H9_directRNA_Rep1_Run2** sample.

```bash
# quickly check the number of reads mapped in the bam file, and also read accuracy 
samtools stats GIS_H9_directRNA_Rep1_Run2_R1.bam
```


### Initializing R 

To summarise the data. We have 6 direct RNA Nanopore long-read samples, 3 replicates each from H9 and HEYA8 human cancer cell line. 
H9 are embryonic stem stells and HEYA8 are ovary cells extracted from adenocarcinoma tissue. Now that we have all our data, let's start up R

```bash
R
```



### Preparing data for ***Bambu***
```rscript
# set work directory 
setwd("bambu_training")


library(bambu)
fa.file <- "./Homo_sapiens.GRCh38.dna.primary_assembly.fa"
gtf.file <- "Homo_sapiens.GRCh38.91.sorted.gtf"
annotations <- prepareAnnotations(gtf.file) # this will prepare annotation object for Bambu, gtf file can be loaded directly, but for best practice, you can create annotation object first, and then you can use it repeatedly with Bambu without the need to prepare it again.


samplePath  = "./"
samples.H9 = c("GIS_H9_directRNA_Rep1_Run2_R1.sorted.bam", "GIS_H9_directRNA_Rep2_Run2_R1.sorted.bam", "GIS_H9_directRNA_Rep3_Run2_R1.sorted.bam")
samples.HEYA8 = c("GIS_HEYA8_directRNA_Rep1_Run2_R1.sorted.bam", "GIS_HEYA8_directRNA_Rep2_Run2_R1.sorted.bam", "GIS_HEYA8_directRNA_Rep3_Run1_R1.sorted.bam")
```



### Running ***Bambu***

The default mode to run ***Bambu*** is using a set of aligned reads (bam files), reference genome annotations (gtf file, TxDb object, or bambuAnnotation object), and reference genome sequence (fasta file or BSgenome). ***Bambu*** will return a ***summarizedExperiment*** object with the genomic coordinates for annotated and new transcripts and transcript expression estimates. We highly recommend to use the same annotations that were used for genome alignment. 

```rscript
se.H9 = bambu(reads = samples.H9, annotations = annotations, genome = fa.file, rcOutDir = "./rc")
se.HEYA8 = bambu(reads = samples.HEYA8, annotations = annotations, genome = fa.file, rcOutDir = "./rc")
```
If you have a gtf file and fasta file you can simply replace `annotations` with `gtf.file` in the above codes.


#### Output
***Bambu*** returns a ***SummarizedExperiment*** object which can be accessed as follows:

```bash
assays(se.H9) #returns the transcript abundance estimates as counts or CPM
rowRanges(se.HEYA8) #returns a GRangesList with all annotated and newly discovered transcripts
rowData(se.H9) #returns additional information about each transcript such as the gene name and the class of newly discovered transcript
```

#### Save output

As we used the `rcOutDir` argument with ***Bambu*** the read class files were automatically saved, allowing us to restart transcript discovery or quantification. To save the final results we can save the output ***SummarizedExperiment*** object. 

```rscript
saveRDS(se.H9, file = "bambu.H9.rds")
saveRDS(se.HEYA8, file = "bambu.HEYA8.rds")
```

To reload them
```rscript
se.H9 = loadRDS("bambu.H9.rds")
se.HEYA8 = loadRDS("bambu.HEYA8.rds")
```

### Output new annotations as a gtf file

The `.gtf` format is used by bioinformatic tools to describe the location of genomic features, such as genes. Below is how output the newly discovered annotations from ***Bambu*** as a `.gtf` for use in downstream tools. Keep in mind that the novel annotations are added to the reference annotations, so all original annotations will still be present in the output. To remove the known annotations, filtering needs to be done and will be shown in a later step.

```rscript
writeToGTF(rowRanges(se.H9), "./H9_novel_annotations.gtf")
writeToGTF(rowRanges(se.HEYA8), "./HEYA8_novel_annotations.gtf")
```

### Identify novel transcripts in one sample and not another
Now let's find some novel transcripts that occur uniquely in H9 but not in HEYA8 and visualise them. First we will filter the annotations to remove the reference annotations. Then we can sort them by the `txNDR`, which ranks them by how likely they are a real transcript. 
```rscript

se.H9.filtered <- se.H9[mcols(rowRanges(se.H9))$newTxClass != "annotation",]
se.HEYA8.filtered <- se.HEYA8[mcols(rowRanges(se.HEYA8))$newTxClass != "annotation",]

#investigate high confidence novel isoforms
head(mcols(rowRanges(se.filtered))[order(mcols(rowRanges(se.filtered))$txNDR),])
plotBambu(se, type = "annotation", gene_id = "gene.101")
#no novel genes
se.filtered.noNovelGene = se.filtered[rowData(se.filtered)$newTxClass != "newGene-spliced",]
mcols(rowRanges(se.filtered.noNovelGene))[order(mcols(rowRanges(se.filtered.noNovelGene))$txNDR),]
plotBambu(se, type = "annotation", gene_id = "ENSG00000141429")

#find which ones are unique
se.unique <- !bambu:::isReadClassCompatible(rowRanges(se.filtered), rowRanges(se2.filtered))$equal
se2.unique <- !bambu:::isReadClassCompatible(rowRanges(se2.filtered), rowRanges(se.filtered))$equal

```

![gene.101](/images/bambu/gene.101.png)

![ENSG00000141429](/images/bambu/ENSG00000141429.png)

### gffCompare 

?????

### UCSC genome browser

We can visualise these novel transcripts using a genome browser such as IGV or the UCSC genome browser. For this tutorial we will use the UCSC genome browser. 

First we need to produce the bed files from the new transcripts:
```rscript
annotations.H9.UCSC = rowRanges(se.H9.filtered)
seqlevelsStyle(annotations.H9.UCSC) <- "UCSC" #this reformats the chromosome names 
annotations.H9.UCSC = keepStandardChromosomes(annotations.H9.UCSC, pruning.mode="coarse") #removes chromosomes UCSC doesn't have
writeToGTF(annotations.H9.UCSC, "./H9_novel_annotations.UCSC.gtf")
```

Now we need to create a bucket to store this data so it can be accessed by the UCSC genome browser. We need to upload this so that it is publically accessible. Note: Anyone can access this file now so do not use this method for any files that need to remain private or that you do not want to be accessed by anyone. 
```bash
aws s3 cp H9_novel_annotations.UCSC.gtf s3://bucket/path/ --acl public-read
```
Now go to https://genome.ucsc.edu/ in your browser                 
My data > Custom Tracks > add custom tracks                   
In the box labeled "Paste URLs or data" copy in path of the file you copied onto your bucket                             
Remember to replace "bucket" and "path" with the real names and path                    
https://"bucket".s3.ap-southeast-1.amazonaws.com/"path"/H9_novel_annotations.gtf

![UCSC brower](/images/bambu/UCSC.png)

### Bioconductor archive files (rcFiles)

`rcFiles` represent a ***SummerizedExperiment*** object holding the read classes calculated from an individual file. They can be used to skip the read class construction part of ***Bambu***. `rcFiles` can be generated in two ways, either as a direct output of the `bambu()` function when `quant` and `discovery` are FALSE, or as written outputs when a path is provided to the `rcOutdir` argument. When `rcFiles` are output using `rcOutdir` this is done using ***BiocFileCache***. For more details on how to access, use, and identify these files see (https://bioconductor.org/packages/release/bioc/html/BiocFileCache.html). A short example is shown below.

```rscript
library(BiocFileCache)
bfc = BiocFileCache("path/to/rcOutput/", ask = FALSE)
info = bfcinfo(bfc)
```
The info object is a ***tibble*** object which associates the filename (fpath) with the sample (rname) to help you identify which `.rds` file you need.

### Combine quantification

For many downstream tools different samples need to use the same annotations so expression values can be compared. ***Bambu*** can also be run with multiple sample types. Here we will use the `rcFiles` we generated earlier to speed up the run time. This may take some time to run. 

```rscript
se.combined = bambu(rcFile = info$rpath, annotations = annotations, genome = fa.file)
```

### Visualising quantification
***Bambu*** allows us to quickly visualise the quantification of the different samples in heatmap and pca form. In the heatmap the correlation between the transcript expression of each of the samples is compared. You should see that samples from the same cell line are clustered together. This should also be the case in the PCA plot.
```rscript
library(ggplot2)
plotBambu(se, type = "heatmap")

plotBambu(se, type = "pca")
```

![heatmap](/images/bambu/heatmap.png)

![PCA](/images/bambu/PCA.png)


### Identifying differentially expressed genes 

One of the most common tasks when analysing RNA-Seq data is the analysis of differential gene expression across a condition of intertest. Here we use ***DESeq2*** to find the differentially expressed genes between H9 and HEYA8 cell lines. As ***DESeq2*** requires integer estimate input, so we will first round the estimates from ***Bambu***, similar to how results from ***Salmon*** are used, which is the most famous transcript quantification tool for short read RNA-Seq data.

```rscript
library(DESeq2)

# add a metadata column to define the differential condition, the cell line here
colData(se.combined)$condition = c("H9", "H9", "H9", "HEYA8", "HEYA8", "HEYA8")

# convert transcript expression to gene expression
se.combined.gene = transcriptToGeneExpression(se.combined)

# prepare counts to a usable object for `DESeq` function
dds <- DESeqDataSetFromMatrix(round(assays(se.combined.gene)$counts), colData = colData(se.combined.gene),
    design = ~condition)
    
# apply DESeq function
dds.deseq <- DESeq(dds)

# extract the results tables and visualise 
deGeneRes <- DESeq2::results(dds.deseq, independentFiltering = FALSE)
head(deGeneRes[order(deGeneRes$padj), ])
summary(deGeneRes)

# Log fold change shrinkage for visualization and ranking
library(apeglm)
resLFC <- lfcShrink(dds.deseq, coef = "condition_HEYA8_vs_H9", type = "apeglm")
plotMA(resLFC, ylim = c(-3, 3))
```
![MA_plot](/images/bambu/DESeq2_MA_plot.png)

For more details on how ***DESeq2*** can be used, please visit the ***DESeq2*** [tutorial](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

### Detecting alternatively used isoforms

We can also used ***DEXseq*** to look at differential isoform expression

```rscript
library(DEXSeq)
dxd <- DEXSeqDataSet(countData = round(assays(se)$counts), sampleData = as.data.frame(colData(se)),
    design = ~sample + exon + condition:exon, featureID = rowData(se)$TXNAME, groupID = rowData(se)$GENEID)
dxr <- DEXSeq(dxd)
head(dxr)

plotMA(dxr, cex = 0.8)
```

![MA_plot](/images/bambu/DEXSeq_MA_plot.png)

For more details on how ***DEXSeq*** can be used, please visit the ***DEXSeq*** [tutorial](https://bioconductor.org/packages/devel/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html).
