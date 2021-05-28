---
title: "Run your analyses"
date: 2019-01-24T09:05:54Z
weight: 300
pre: "<b>VII ⁃ </b>"
tags: ["HPC", "Overview", "Batch"]
---
We will now run some sample code.

**First, let's check our tools:**


```bash
which bwa
```

Output shows where bwa is installed.

```bash
which samtools
```
Output shows where samtools is installed.



**Basic Analyses**


1. Get reference sequence

```bash
mkdir -p /tmp/outbreaks/SG-M1
```
```bash
cd /tmp/outbreaks/SG-M1
```
```bash
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/275/545/GCF_001275545.2_ASM127554v2/GCF_001275545.2_ASM127554v2_genomic.fna.gz
```
```bash
gunzip GCF_001275545.2_ASM127554v2_genomic.fna.gz
```
```bash
mv GCF_001275545.2_ASM127554v2_genomic.fna SG-M1.fna
```

2. Map and call SNPs

```bash
bwa index SG-M1.fna
```
```bash
bwa mem SG-M1.fna /tmp/fastq/SRR6327950/SRR6327950_1.fastq.gz /tmp/fastq/SRR6327950/SRR6327950_2.fastq.gz | samtools view -bS - > SRR6327950.bam
```
```bash
samtools sort SRR6327950.bam -o SRR6327950-sort.bam
```
```bash
samtools index SRR6327950-sort.bam
```
```bash
lofreq faidx SG-M1.fna
```
```bash
lofreq call -f SG-M1.fna -r NZ_CP012419.2:400000-500000 SRR6327950-sort.bam > SRR6327950-400k.vcf
```
Mapping takes ~5 min on a t2.medium. Sorting takes ~2 min, lofreq ~1 min. If need to downsample:

```bash
seqtk sample /tmp/fastq/SRR6327950/SRR6327950_1.fastq.gz 0.25 | gzip -c > SRR6327950-reduced_1.fastq.gz
```
```bash
seqtk sample /tmp/fastq/SRR6327950/SRR6327950_2.fastq.gz 0.25 | gzip -c > SRR6327950-reduced_2.fastq.gz
```

3. Assembly (takes ~4 min then will run out of RAM):
```bash
spades.py -t 2 -1 /tmp/fastq/SRR6327950/SRR6327950_1.fastq.gz -2 /tmp/fastq/SRR6327950/SRR6327950_2.fastq.gz -o SRR6327950_spades
```

**NOTE**: This assembly above will complete on a t3a.large, takes about 5 hours.
