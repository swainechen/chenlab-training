+++
title = "b. RNA-Seq data"
date = 2019-01-24T09:05:54Z
weight = 40
tags = ["bambu", "RNA-seq", "transcript"]
+++


We will now run some sample analysis using R

First lets test that everything is installed correctly by running bambu using a small test dataset and see that it runs through without error.

We will open an R session. This will change the terminal to R which is a different language from what you have been using in the rest of the tutorial (bash) and those commands will not work. 

This will show us where R is installed
```bash
which R
```

This will start the R session
```bash
R
```

A workspace image keeps track of all the variables you have stored and can be reopened at a later point. For this part of the tutorial it is not needed

The bambu package has already been installed on your machine and all you need to do is load it into R
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

You can explore the data with a few commands but we will go into more detail with this later
```rscript
se
```
```rscript
rowData(se)
```
```rscript
rowRanges(se)
```
```rscript
assays(se)
```
```rscript
assays(se)$counts
```
```rscript
pdf(file = "ENSG00000107104_CPM.pdf", width = 6, height = 8)
plotBambu(se, type = "annotation", gene_id = "ENSG00000107104")
dev.off()
```
The last command will produce a plot showing the different transcripts found for the gene ENSG00000107104 and their expression. We will show you later in the tutorial on how to download this plot to your local machine so you can view it

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

Congratulations, you now know how to run transcript discovery and quantification using R on an AWS EC2 instance. Next we will...
