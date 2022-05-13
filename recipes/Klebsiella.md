# Klebsiella assembly and annotation

There's a great tool that unifies analysis of *Klebsiella pneumoniae* called [Kleborate](https://github.com/katholt/Kleborate).
This will do all the basics - species (somewhat complicated for what we used to call K. pneumoniae), resistance, MLST, key virulence factors, and serotype.
All it requires is an assembly.

## Download data

We'll take an example dataset from NCBI:
* Run `SRR6505299` (WGS of resistant E. coli and K. pneumoniae: GASRECK48)
* BioSample `SRS2866648` (aka `SAMN08382736`)
* Study `SRP131304` (3rd generation cephalosporin-resistant E. coli and K. pneumonia (GASREC) Genome sequencing and assembly)
* BioProject `PRJNA431029` (3rd generation cephalosporin-resistant E. coli and K. pneumonia (GASREC) Genome sequencing and assembly)

Downloading this data set requires ~4GB of storage due to the intermediate files, but the final fastq files are only ~360MB.
```
mkdir -p /home/ubuntu/fastq/SRR6505299
aws s3 --no-sign-request sync s3://sra-pub-run-odp/sra/SRR6505299/ /home/ubuntu/fastq/SRR6505299
# the files are stored in .sra format, so we need to convert to fastq
cd /home/ubuntu/fastq/SRR6505299
fasterq-dump SRR6505299
for i in *.fastq; do gzip $i; done
# clean up to save some space
rm /home/ubuntu/fastq/SRR6505299/SRR6505299
```

## Assembly
```
# run the assembly with velvet
SLC-wgs.pl -q1 SRR6505299_1.fastq.gz -q2 SRR6505299_2.fastq.gz -name SRR6505299 -output_dir /home/ubuntu/ -assembler velvetoptimizer
# sample commands for SPAdes or SKESA
#SLC-wgs.pl -q1 SRR6505299_1.fastq.gz -q2 SRR6505299_2.fastq.gz -name SRR6505299 -output_dir /home/ubuntu/ -assembler spades
#SLC-wgs.pl -q1 SRR6505299_1.fastq.gz -q2 SRR6505299_2.fastq.gz -name SRR6505299 -output_dir /home/ubuntu/ -assembler skesa
```

## Annotation with Kleborate
```
# we'll use the .assembly file, which has the contigs (the .fna file joins them with 100 N's between contigs)
cd /home/ubuntu
tar xvzf SRR6505299.tgz
cd SRR6505299
kleborate -a SRR6505299.assembly --all -o SRR6505299.kleborate.txt
```
This is an ST307 *K. pneumoniae* strain with 12 predicted resistance genes, including OXA-1 and CTX-M-15, fluoroquinolone resistance mutations in gyrA and parC, a zero (0) virulence score.
