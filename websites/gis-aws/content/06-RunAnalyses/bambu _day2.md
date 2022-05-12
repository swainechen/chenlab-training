***bambu*** is a R package for multi-sample transcript discovery and quantification using long read RNA-Seq data. You can use ***bambu*** after read alignment to obtain expression estimates for known and novel transcripts and genes. The output from ***bambu*** can directly be used for visualisation and downstream analysis such as differential gene expression or transcript usage.

### Downloading reference genome and annotations

The default mode to run ***bambu** is using a set of aligned reads (bam files), reference genome annotations (gtf file, TxDb object, or bambuAnnotation object), and reference genome sequence (fasta file or BSgenome). In this tutorial we will be using data from human cancer cell-lines from the SG-NEx project. 

Goals:
Detect novel transcripts and genes in the samples
Identify differentially expressed transcripts
Detect alternate isoform usage of novel genes

```bash

aws s3 cp genome/path ./
aws s3 cp annotation/path ./

```

### Downloading and aligning reads

The SG-NEx project provides aligned reads that can be used directly with Bambu, however many datasets will only provide the raw reads (fastq). Therefore we will download the raw reads from SG-NEx and align them ourselves.

Note that we are aligning the reads to the genome and not the transcriptome as bambu uses intron junctions to distingush novel transcripts.

```bash

aws s3 cp fastq/path

#./minimap2 -ax splice ref.fa rna-reads.fa > aln.sam       # spliced long reads (strand unknown)
#./minimap2 -ax splice ref.fa nanopore-cdna.fa > aln.sam # Nanopore 2D cDNA-seq
./minimap2 -ax splice -uf -k14 ref.fa reads.fa > aln.sam  # noisy Nanopore Direct RNA-seq
#./minimap2 -ax splice:hq -uf ref.fa query.fa > aln.sam    # Final PacBio Iso-seq or traditional cDNA

samtools view -Sb XXX.sam | samtools sort -o XXX.bam
samtools index XXX.bam $>> XXX.bam.bai

```

To download the aligned reads from SG-NEx (skipping the above step)

```bash

aws s3 cp bam/path ./
aws s3 cp bam.bai/path ./

```

Now that we have all our data, let's start up R
```bash
R
```

### Preparing data for Bambu
```rscript
library(bambu)
gtf.file <- "PATH"
annotation <- prepareAnnotations(gtf.file)
fa.file <- "PATH"
sample <- "PATH"
```

### Running Bambu

The default mode to run bambu is using a set of aligned reads (bam files), reference genome annotations (gtf file, TxDb object, or bambuAnnotation object), and reference genome sequence (fasta file or BSgenome). bambu will return a summarizedExperiment object with the genomic coordinates for annotated and new transcripts and transcript expression estimates. We highly recommend to use the same annotations that were used for genome alignment. If you have a gtf file and fasta file you can run bambu with the following options:

```rscript

se = bambu(reads = sample, annotations = annotations, genome = fa.file)
```

### Output
bambu returns a SummarizedExperiment object which can be accessed as follows:

assays(se) returns the transcript abundance estimates as counts or CPM
rowRanges(se) returns a GRangesList with all annotated and newly discovered transcripts
rowData(se) returns additional information about each transcript such as the gene name and the class of newly discovered transcript

### Performing transcript discovery

Some research questions focus on only indentifying novel transcripts and genes and quantification is not important. We can stop the quantification module of bambu as follows:
```rscript

newAnnotations <- bambu(reads = sample, annotations = annotations, genome = fa.file, quant = FALSE)
```
### Visualize novel transcripts

```rscript
plotBambu(se, type = "annotation", gene_id = "ENSG00000099968")
```

### Identify novel transcripts in one sample and not another

### UCSC genome browser

### Bioconductor archive files (rcFiles)

rcFiles represent a SummerizedExperiment object holding the read classes calculated from an individual file. They can be used to skip the read class construction part of Bambu. rcFiles can be generated in two ways, either as a direct output of the bambu() function when quant and discovery are FALSE, or as written outputs when a path is provided to the rcOutdir argument. When rcFiles are output using rcOutdir this is done using BiocFileCache. For more details on how to access, use, and identify these files see (https://bioconductor.org/packages/release/bioc/html/BiocFileCache.html). A short example is shown below.

```rscript
library(BiocFileCache)
bfc = BiocFileCache("path/to/rcOutput/", ask = FALSE)
info = bfcinfo(bfc)
```
The info object is a tibble which associates the filename (fpath) with the sample (rname) to help you identify which .rds file you need.

### Visualising quantification

```rscript
library(ggplot2)
plotBambu(se, type = "heatmap")

plotBambu(se, type = "pca")
```

### Identifying differentially expressed genes

One of the most common tasks when analysing RNA-Seq data is the analysis of differential gene expression across a condition of intertest. Here we use DESeq2 to find the differentially expressed genes between MCF7 and K562 cell lines. Similar to using results from Salmon, estimates from bambu will first be rounded.

```rscript
library(DESeq2)
dds <- DESeqDataSetFromMatrix(round(assays(seGene)$counts), colData = colData(se),
    design = ~condition)
dds.deseq <- DESeq(dds)
deGeneRes <- DESeq2::results(dds.deseq, independentFiltering = FALSE)
head(deGeneRes[order(deGeneRes$padj), ])
```

### Detecting alternatively used isoforms

We used DEXSeq to detect alternative used isoforms.

```rscript
library(DEXSeq)
dxd <- DEXSeqDataSet(countData = round(assays(se)$counts), sampleData = as.data.frame(colData(se)),
    design = ~sample + exon + condition:exon, featureID = rowData(se)$TXNAME, groupID = rowData(se)$GENEID)
dxr <- DEXSeq(dxd)
head(dxr)
```
