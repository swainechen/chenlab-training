# Lacer
## FASTQ quality score recalibration
From the earliest 

Get some sample data. The Illumina data is for UTI89, which is about 2% divergent from MG1655. Lacer will recalibrate the qualities nearly identically regardless of whether you use the UTI89 or MG1655 references. Because GATK requires further assumptions for recalibration, the recalibrated qualities will vary significantly if you do the same experiment with GATK, and the ones obtained using the "perfect" reference (UTI89) will match those given by Lacer.
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
samtools faidx MG1655.fna
lacer.pl -bam SRR7349974-MG1655-sort.bam -ref MG1655.fna -stopbases 3000000 -outfile SRR7349974-MG1655.lacer.txt
```

There are several options for applying the recalibration. One of the easier ones is to recalibrate the fastq files directly:
```
# or you can recalibrate the fastq files directly
lacepr --fastq SRR7349974_1.fastq.gz --recal SRR7349974-MG1655.lacer.txt --pair 1 --out SRR7349974-recal_1.fastq.gz
lacepr --fastq SRR7349974_2.fastq.gz --recal SRR7349974-MG1655.lacer.txt --pair 2 --out SRR7349974-recal_2.fastq.gz
```

The "traditional" way is to recalibrate the bam file. This can be done with the include lacepr program or with GATK. For most microbial genomes (which GATK really can't recalibrate, due to the requirement for a perfect reference or known SNPs), read groups aren't commonly needed or used, so it's a little easier to use lacepr:
```
# recalibrate the bam file using lacepr - force the read group
lacepr --bam SRR7349974-MG1655-sort.bam --recal SRR7349974-MG1655.lacer.txt --out SRR7349974-MG1655-sort-recal.bam --rg NULL
```

The equivalent using GATK requires adding read groups. NULL is the read group string used by lacer when no read group is in the bam file during recalibration, so you have to add that (or add the read groups before recalibration):
```
gatk AddOrReplaceReadGroups -I SRR7349974-MG1655-sort.bam -O SRR7349974-MG1655-RG.bam -RGID 1 -RGLB lib -RGPL ILLUMINA -RGPU NULL -RGSM Sample
gatk CreateSequenceDictionary -R MG1655.fna -O MG1655.dict
gatk ApplyBQSR --input SRR7349974-MG1655-RG.bam --output SRR7349974-MG1655-sort-recal.bam --reference MG1655.fna --bqsr-recal-file SRR7349974-MG1655.lacer.txt
```
