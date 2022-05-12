# Assembly and analysis of a Salmonella isolate

We'll use a publicly available data set as an example.
The organization of data in the main public databases associated quite a few pieces of information with each data set.
Many of these have their own accession number.
For this, we'll use (from most specific to least specific):
* Run `SRR12151671` (WGS of Salmonella: SGEHI2016-PSU-BS-067SL)
* BioSample `SRS6951036` (aka `SAMN15447839`)
* Study `SRP270279` (Whole Genome Sequencing of Salmonella isolates from Thailand)
* BioProject `PRJNA644105` (Whole Genome Sequencing of Salmonella isolates from Thailand)

The "raw" sequencing data is accessed by the "Run" accession.
Note that a single BioSample may have multiple Runs associated with it (and, in turn, a single Study or BioProject can have multiple BioSamples).
Fortunately, AWS now has the entire SRA on S3, organized by Run accessions - this can be accessed without credentials (`--no-sign-request` option on the aws CLI).

It's always good to keep your data organized.
For this example, we'll just use our home directory, but you may want to consider using a dedicated location on additionally mounted storage, for example, to store fastq files.

This data set requires ~450MB of storage for the final gzip'd fastq files.
However, the total storage required for the following commands (due to needing to temporarily store the initial sra file and intermediate fastq files) is ~4GB.
```
mkdir -p /home/ubuntu/fastq/SRR12151671
aws s3 --no-sign-request sync s3://sra-pub-run-odp/sra/SRR12151671/ /home/ubuntu/fastq/SRR12151671
# the files are stored in .sra format, so we need to convert to fastq
cd /home/ubuntu/fastq/SRR12151671
fasterq-dump SRR12151671
for i in *.fastq; do gzip $i; done
# clean up to save some space
rm /home/ubuntu/fastq/SRR12151671/SRR12151671
```

## Assembly and genome annotation
I have a wrapper script in my SLC Closet repository that helps do a basic assembly and annotation.
It uses velvet by default (doing a full optimization of the kmer parameter), but also can use SPAdes or skesa currently.
It outputs a .tgz file which has the assembly, annotation, and log files.
```
# check the help
SLC-wgs.pl -help

# run the assembly with velvet
SLC-wgs.pl -q1 SRR12151671_1.fastq.gz -q2 SRR12151671_2.fastq.gz -name SRR12151671 -output_dir /home/ubuntu/ -assembler velvetoptimizer
# sample commands for SPAdes or SKESA
#SLC-wgs.pl -q1 SRR12151671_1.fastq.gz -q2 SRR12151671_2.fastq.gz -name SRR12151671 -output_dir /home/ubuntu/ -assembler spades
#SLC-wgs.pl -q1 SRR12151671_1.fastq.gz -q2 SRR12151671_2.fastq.gz -name SRR12151671 -output_dir /home/ubuntu/ -assembler skesa
```

## Calling MLSTs
One can use both raw short reads and assembled sequences to call MLSTs.
The databases from PubMLST for Salmonella are already installed on this machine.
You'll likely just do one of these, but both are shown here:
```
# using short reads and SRST2
srst2 --input_pe /home/ubuntu/fastq/SRR12151671/SRR12151671_1.fastq.gz /home/ubuntu/fastq/SRR12151671/SRR12151671_2.fastq.gz --output SRR12151671-MLST --log --mlst_db /usr/local/lib/SRST2/MLST/Salmonella_enterica.fasta --mlst_definitions /usr/local/lib/SRST2/MLST/senterica.txt --mlst_delimiter _

# using the assembly we just did
# the following script is from the SLC Closet repository
# it looks for the SRST2 databases, and also will detect and automatically process assembly .tgz files from the SLC-wgs.pl script
blast-mlst.pl -species Senterica /home/ubuntu/SRR12151671.tgz > SRR12151671-assembly.mlst
```
This is an ST1549 strain.

## Calling resistance genes and virulence factors
These again leverage SRST2 for short reads and a blast-based script (against the same databases) for assemblies.
Both are shown here but you will likely just use one depending on your own preferences.
These require some database setup that is done on the `CHENLAB-PUBLIC` AMI already and described in the [SRST2](/sysadmin/bioinformatics.md) section under [Bioinformatics software setup](/sysadmin/bioinformatics.md).
```
# using short reads and SRST2
# find the correct file
ls /usr/local/lib/SRST2/0.2.0
ls /usr/local/lib/SRST2/Senterica
# there is a combined database set up for Senterica
srst2 --input_pe /home/ubuntu/fastq/SRR12151671/SRR12151671_1.fastq.gz /home/ubuntu/fastq/SRR12151671/SRR12151671_2.fastq.gz --output SRR12151671-VF --log --gene_db /usr/local/lib/SRST2/Senterica/Senterica-combined-2021-06-02.fasta

# using the assembly we just did
# the following script is from the SLC Closet repository
# it looks for the SRST2 databases, and also will detect and automatically process assembly .tgz files from the SLC-wgs.pl script
blast-vf.pl -species Senterica /home/ubuntu/SRR12151671.tgz > SRR12151671-assembly.genes
# resistances
grep 'DB:ARGannot' SRR12151671-assembly.genes
# others - like virulence genes
grep -v 'DB:ARGannot' SRR12151671-assembly.genes

# another alternative for assemblies - using abricate
# the .assembly file has all the contigs
# the .fna file has all the contigs joined (in arbitrary order) by 100 N's
# either can be used
tar xvzf SRR12151671.tgz
abricate SRR12151671/SRR12151671.fna
abricate SRR12151671/SRR12151671.assembly
```
`SRST2` finds one predicted resistance gene (Aac6-Iaa, encoding resistance to aminoglycosides; this is common on the Salmonella chromosome but often doesn't confer phenotypic resistance), and quite a few genes classified as virulence factors.
`Abricate` doesn't find any resistance genes.

## Calling serotype
There is a dedicated program called `SeqSero` for calling Salmonella serotypes, this can use either short reads or assemblies.
This strain ends up being 
```
# call using short reads
SeqSero.py -m 2 -i /home/ubuntu/fastq/SRR12151671/SRR12151671_1.fastq.gz /home/ubuntu/fastq/SRR12151671/SRR12151671_2.fastq.gz -d SRR12151671-fastq
# results are in this file
cat SRR12151671-fastq/Seqsero_result.txt

# pull out the assembly fasta file and call serotypes from there
tar xvzf SRR12151671.tgz
# check the files available
ls -l SRR12151671/
# SeqSero requires an input assembly to have an extension of .fasta, .fna, .fa, or .fsa
# we can use either the .assembly file (after linking or renaming it) or the .fna file here
SeqSero.py -m 4 -i SRR12151671/SRR12151671.fna -d SRR12151671-assembly
# results are in this file
cat SRR12151671-assembly/Seqsero_result.txt
```
This strain is predicted to be be an Altona serotype strain.
