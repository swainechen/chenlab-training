+++
title = "b. Long-Read RNA-Seq data"
date = 2019-01-24T09:05:54Z
weight = 40
tags = ["bambu", "RNA-seq", "transcript"]
+++

Long read RNA sequencing, Nanopore sequencing, has been widely adopted in different areas. The ability to produce full-length reads of RNAs has been a great improvement to the previous sequencing technologies. This is shown to have the potential to uncover the full genome picture, same as for transcriptomics, where varying isoforms can be discovered in full, fusion isoforms even. 

Here we use this long read RNA sequencing analysis tool, Bambu, a R package that can easily process multiple samples in one command that does both transript discovery and quantification. Comparing to other softwares, this software is comparably more efficient while being more sensitive and precise in transcript discovery, which is very important in transcript quantification as well. 

```bash 
mkdir -p /tmp/RNASeq
cd /tmp/RNASeq
```

The analysis is run with built-in R, a different language from what you have been using in the rest of the tutorial (bash) and those commands will not work. This session will require a AWS machine of 8GB RAM.

To start an R session
```bash
R
```

A workspace image keeps track of all the variables you have stored and can be reopened at a later point. For this part of the tutorial it is not needed

To load the bambu packages, which is already pre-installed:
```rscript
library(bambu)
```

Now we will load in a small test dataset including aligned reads, a genome, and annotations. This will be explained in further detail on day 2.

```rscript
test.bam <- system.file("extdata", "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam", package = "bambu")
  
fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa", package = "bambu")

gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")
```

Lets try and detect if we can find any novel transcripts in this long-read data and quantify their expression

```
bambuAnnotations <- prepareAnnotations(gtf.file)

se <- bambu(reads = test.bam, annotations = bambuAnnotations, genome = fa.file)
```

You can explore the data with a few commands but we will go into more detail with this in a later tutorial
```rscript
se
head(rowData(se))
head(rowRanges(se))
head(assays(se))
head(assays(se)$counts)
```

The next command will produce a plot showing the different transcripts found for the gene ENSG00000107104 and their expression. We will show you later in the tutorial on how to download this plot to your local machine so you can view it
```rscript
pdf(file = "ENSG00000107104_CPM.pdf", width = 6, height = 8)
plotBambu(se, type = "annotation", gene_id = "ENSG00000107104")
dev.off()
```

***NOTE:*** To leave the R session at any point once you have opened it. Do not type this in until you have completed the rest of the steps or else you will need to reopen R. 
```rscript
quit()
# type in 'n' in response to "Save workspace image?
```

Once you have left R you an see if the file was made
```bash
ls
```
You should see a file called ENSG00000107104_CPM.pdf.

Congratulations, you now know how to run transcript discovery and quantification using R on an AWS EC2 instance. Next, we'll look at how to make an AMI from your own machine. 
