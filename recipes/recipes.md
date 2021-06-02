# Recipes
Here are some examples using the resources described in this GitHub repository (and associated [training website]()).

## Basic assembly and analysis of a Salmonella isolate

We'll use a publicly available data set as an example.
The organization of data in the main public databases associated quite a few pieces of information with each data set.
Many of these have their own accession number.
For this, we'll use (from most specific to least specific):
* Run `SRR12151671` (WGS of Salmonella: SGEHI2016-PSU-BS-067SL)
* BioSample `SRS6951036` (aka `SAMN15447839`)
* Study `SRP270279` (Whole Genome Sequencing of Salmonella isolates from Thailand)
* BioProject `PRJNA644105`

The "raw" sequencing data is accessed by the "Run" accession.
Note that a single BioSample may have multiple Runs associated with it (and, in turn, a single Study or BioProject can have multiple BioSamples).

It's always good to keep your data organized.
For this example, we'll just use our home directory, but you may want to consider using a dedicated location on additionally mounted storage, for example, to store fastq files.
```
mkdir -p /home/ubuntu/fastq/SRR12151671
cd /home/ubuntu/fastq/SRR12151671
kingfisher --run-identifier SRR12151671 -m ena-ascp
```

### Assembly and genome annotation
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

### Calling MLSTs
One can use both raw short reads and assembled sequences to call MLSTs.
The databases from PubMLST for Salmonella are already installed on this machine.
You'll likely just do one of these, but both are shown here:
```
# using short reads and SRST2
srst2 --input_pe /home/ubuntu/fastq/SRR12151671/SRR12151671_1.fastq.gz /home/ubuntu/fastq/SRR12151671/SRR12151671_2.fastq.gz --output SRR12151671-MLST --log --mlst_db /usr/local/lib/SRST2/MLST/Salmonella_enterica.fasta --mlst_definitions /usr/local/lib/SRST2/MLST/senterica.txt --mlst_delimiter _

# using the assembly we just did
# the following script is from the SLC Closet repository
# it looks for the SRST2 databases, and also will detect and automatically process assembly .tgz files from the SLC-wgs.pl script
blast-mlst.pl -species Senterica /home/ubuntu/SRR12151671.tgz
```

### Calling resistance genes and virulence factors
These again leverage SRST2 for short reads and a blast-based script (against the same databases) for assemblies.
Both are shown here but you will likely just use one depending on your own preferences.
```
# using short reads and SRST2
# find the correct file
ls /usr/local/lib/SRST2/0.2.0
ls /usr/local/lib/SRST2/Senterica
# there is a combined database set up for Senterica
srst2 --input_pe /home/ubuntu/fastq/SRR12151671/SRR12151671_1.fastq.gz /home/ubuntu/fastq/SRR12151671/SRR12151671_2.fastq.gz --output SRR12151671-VF --log --gene_db /usr/local/lib/SRST2/Senterica/<.fasta>

# using the assembly we just did
# the following script is from the SLC Closet repository
# it looks for the SRST2 databases, and also will detect and automatically process assembly .tgz files from the SLC-wgs.pl script
blast-vf.pl -species Senterica /home/ubuntu/SRR12151671.tgz
```
