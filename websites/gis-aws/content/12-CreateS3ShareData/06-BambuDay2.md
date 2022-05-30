+++
title = "e. RNA-Seq data-Bambu"
weight = 150
tags = ["S3", "RNA-Seq", "Bambu"]
+++


## Usage examples of Bambu

In this tutorial we will be using data from human cancer cell-lines from the [SG-NEx project](https://github.com/GoekeLab/sg-nex-data). We have 6 direct RNA Nanopore long-read samples, 3 replicates each from A549 and HepG2 human cancer cell line. A549 cell line are extracted from lung tissues from a patient with lung cancer and HepG2 are extracted from hepatocellular carcinoma from a patient with liver cancer. We will use Bambu to identify and quantify novel isoforms in these two cell lines.  

### Downloading reference genome and annotations

As mentioned in #Bambu-Day 1, the default mode to run ***Bambu*** is using a set of aligned reads (bam files), reference genome annotations (gtf file, TxDb object, or bambuAnnotation object), and reference genome sequence (fasta file or BSgenome). Here we use the ensembl Grch 38 genome sequence and annotation files, which we have downloaded and stored in the SG-NEx data open S3 bucket. 

```bash
# create work directory
mkdir bambu_training
cd bambu_training

# download genome fasta file 
aws s3 cp --no-sign-request s3://sg-nex-data/data/annotations/genome_fasta/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa ./
aws s3 cp --no-sign-request s3://sg-nex-data/data/annotations/genome_fasta/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai ./

# download gtf file
aws s3 cp --no-sign-request s3://sg-nex-data/data/annotations/gtf_file/Homo_sapiens.GRCh38.91.gtf ./
```

### Downloading and aligning reads 

The SG-NEx project provides aligned reads that can be used directly with Bambu, however sometimes, raw fastq reads might also be provided instead of aligned reads. In that case, alignment can be done as follows:

Note that we are aligning the reads to the genome fasta and not the transcriptome fasta as ***Bambu*** uses intron junctions to distingush novel transcripts.

```bash
# download fastq from s3 bucket
aws s3 cp --no-sign-request s3://sg-nex-data/data/bambu_training/fastq/A549_directRNA_sample1.fastq.gz ./

# align using Minimap2 
minimap2 -ax splice -uf -k14 Homo_sapiens.GRCh38.dna.primary_assembly.fa A549_directRNA_sample1.fastq.gz > A549_directRNA_sample1.sam  #needs more than 4gb ram
samtools view -Sb A549_directRNA_sample1.sam | samtools sort -o A549_directRNA_sample1.bam
samtools index A549_directRNA_sample1.bam
```


For this analysis, you can directly download the bam files needed from the open S3 bucket

```bash
# download aligned bam files for A549 samples 
aws s3 sync --no-sign-request s3://sg-nex-data/data/bambu_training/bam/ ./ 
```
### Quality Control

[***Samtools stats***](http://www.htslib.org/doc/samtools-stats.html) can be used to quickly check the read accuracy and read mappability in the bam files, in the example below, we checked the aligned bam file for the **SGNex_A549_directRNA_replicate1_run1** sample.

```bash
# quickly check the number of reads mapped in the bam file, and also read accuracy 
samtools stats A549_directRNA_sample1.bam | head -46
```
![stats_output](/images/bambu/samtools_stats_examples.png)

This bam file has a read mappability of `150602/184107` by taking the ratio between reads mapped and raw total sequences and a read accuracy of `1-0.178599` by minusing 1 with the error rate. Usually we consider a sample with a mappability of 80% above as relatively good, and the read accuracy is also common for Nanopore reads. 



### Preparing data for ***Bambu***
```rscript
# set work directory 
setwd("bambu_training")


library(bambu)
fa.file <- "Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
gtf.file <- "Homo_sapiens.GRCh38.91.gtf"
annotations <- prepareAnnotations(gtf.file) # this will prepare annotation object for Bambu, gtf file can be loaded directly, but for best practice, you can create annotation object first, and then you can use it repeatedly with Bambu without the need to prepare it again.


samples.A549 <- list.files(".", pattern = "A549.*.bam$", full.names = TRUE)
samples.HepG2 <- list.files(".", pattern = "HepG2.*.bam$", full.names = TRUE)
```



### Running ***Bambu***

To run ***Bambu*** is using, you will need 

- a set of aligned reads (bam files), 
- reference genome annotations (gtf file, TxDb object, or bambuAnnotation object), 
- and reference genome sequence (fasta file or BSgenome). 

***Bambu*** will return a ***summarizedExperiment*** object with 

- the genomic coordinates for annotated and new transcripts 
- and transcript expression estimates. 

We highly recommend to use the same annotations that were used for genome alignment. 

```rscript
rcFolder <-"./rc"
se.A549 <- bambu(reads = samples.A549, annotations = annotations, genome = fa.file, rcOutDir = rcFolder) # suggest to run with multiple cores 
se.HepG2 <- bambu(reads = samples.HepG2, annotations = annotations, genome = fa.file, rcOutDir = rcFolder) # suggest to run with multiple cores 
```
If you have a gtf file and fasta file you can simply replace `annotations` with `gtf.file` in the above codes.


#### Output
***Bambu*** returns a ***SummarizedExperiment*** object which can be accessed as follows:

```bash
assays(se.A549) #returns the transcript abundance estimates as counts or CPM
rowRanges(se.A549) #returns a GRangesList with all annotated and newly discovered transcripts
rowData(se.A549) #returns additional information about each transcript such as the gene name and the class of newly discovered transcript
```

#### Save output

As we used the `rcOutDir` argument with ***Bambu*** the read class files were automatically saved, allowing us to restart transcript discovery or quantification. To save the final results we can save the output ***SummarizedExperiment*** object. 

```rscript
saveRDS(se.A549, file = "bambu.A549.rds")
saveRDS(se.HepG2, file = "bambu.HepG2.rds")
```

To reload them
```rscript
se.A549 <- loadRDS("bambu.A549.rds")
se.HepG2 <- loadRDS("bambu.HepG2.rds")
```

### Output new annotations as a gtf file

The `.gtf` format is used by bioinformatic tools to describe the location of genomic features, such as genes. Below is how to output the newly discovered annotations from ***Bambu*** as a `.gtf` for use in downstream tools. Keep in mind that the novel annotations are added to the reference annotations, so all original annotations will still be present in the output. To remove the known annotations, filtering needs to be done and will be shown in a later step.

```rscript
writeToGTF(rowRanges(se.A549), "./A549_novel_annotations.gtf")
writeToGTF(rowRanges(se.HepG2), "./HepG2_novel_annotations.gtf")
```

### Identify novel transcripts in one sample and not another
Now let's find some novel transcripts that occur uniquely in A549 but not in HepG2 and visualise them. First we will filter the annotations to remove the reference annotations. Then we can sort them by the `txNDR`, which ranks them by how likely they are a real transcript. 


```rscript

se.A549.filtered <- se.A549[mcols(rowRanges(se.A549))$newTxClass != "annotation",]
se.HepG2.filtered <- se.HepG2[mcols(rowRanges(se.HepG2))$newTxClass != "annotation",]

#investigate high confidence novel isoforms
head(mcols(rowRanges(se.A549.filtered))[order(mcols(rowRanges(se.A549.filtered))$txNDR),])
plotBambu(se, type = "annotation", gene_id = "gene.101")

# visualize their expression in both 

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



### UCSC genome browser

We can visualise these novel transcripts using a genome browser such as IGV or the UCSC genome browser. For this tutorial we will use the UCSC genome browser. 

First we need to produce the bed files from the new transcripts:
```rscript
annotations.A549.UCSC = rowRanges(se.A549.filtered)
seqlevelsStyle(annotations.A549.UCSC) <- "UCSC" #this reformats the chromosome names 
annotations.A549.UCSC = keepStandardChromosomes(annotations.A549.UCSC, pruning.mode="coarse") #removes chromosomes UCSC doesn't have
writeToGTF(annotations.A549.UCSC, "./A549_novel_annotations.UCSC.gtf")
```

Now we can upload the annotations created to the public accessible bucket that we have created before (links).  
```bash
aws s3 cp A549_novel_annotations.UCSC.gtf s3://bucket/path/ --acl public-read
```
Now go to https://genome.ucsc.edu/ in your browser                 
My data > Custom Tracks > add custom tracks                   
In the box labeled "Paste URLs or data" copy in path of the file you copied onto your bucket                             
Remember to replace "bucket" and "path" with the real names and path                    
https://"bucket".s3.ap-southeast-1.amazonaws.com/"path"/A549_novel_annotations.gtf

![UCSC brower](/images/bambu/UCSC.png)


