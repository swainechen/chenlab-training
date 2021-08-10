# Lacer
## FASTQ quality score recalibration
From the earliest 

Get some sample data
```
# MG1655 genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz
mv GCF_000005845.2_ASM584v2_genomic.fna MG1655.fna

# UTI89 genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/265/GCF_000013265.1_ASM1326v1/GCF_000013265.1_ASM1326v1_genomic.fna.gz
gunzip GCF_000013265.1_ASM1326v1_genomic.fna.gz
mv GCF_000013265.1_ASM1326v1_genomic.fna UTI89.fna

# sample fastq data - UTI89 Illumina NextSeq 500, ~260MB total
kingfisher get --run-identifier SRR7349974 -m ena-ascp
```

Map the data and run lacer
```
bwa index MG1655.fna
bwa mem MG1655.fna SRR7349974_1.fastq.gz SRR7349974_2.fastq.gz | samtools view -bS > SRR7349974-MG1655.bam
samtools sort SRR7349974-MG1655.bam -o SRR7349974-MG1655-sort.bam
samtools index SRR7349974-MG1655-sort.bam
lacer.pl -bam SRR7349974-MG1655-sort.bam -ref MG1655.fna -stopbases 3000002 -outfile SRR7349974-MG1655.lacer.txt
```

There are several options for applying the recalibration. One of the easier ones is to recalibrate the fastq files directly:
```
# or you can recalibrate the fastq files directly
lacepr --fastq SRR7349974_1.fastq.gz --recal SRR7349974-MG1655.lacer.txt --pair 1 --out SRR7349974-recal_1.fastq.gz
lacepr --fastq SRR7349974_2.fastq.gz --recal SRR7349974-MG1655.lacer.txt --pair 2 --out SRR7349974-recal_2.fastq.gz
```

The "traditional" way is to recalibrate the bam file. This can be done with the include lacepr program or with GATK:
```
# recalibrate the bam file using lacepr
lacepr --bam SRR7349974-MG1655-sort.bam --recal SRR7349974-MG1655.lacer.txt --out SRR7349974-MG1655-sort-recal.bam

# recalibrate the bam file using GATK - requires sorting out read groups
gatk AddOrReplaceReadGroups -I SRR7349974-MG1655-sort.bam -O SRR7349974-MG1655-RG.bam -RGID 1 -RGLB lib -RGPL ILLUMINA -RGPU Unit -RGSM Sample
gatk CreateSequenceDictionary -R MG1655.fna -O MG1655.dict
gatk ApplyBQSR --input SRR7349974-MG1655-RG.bam --output SRR7349974-MG1655-sort-recal.bam --reference MG1655.fna --bqsr-recal-file SRR7349974-MG1655.lacer.txt
```
