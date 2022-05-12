# Metagenomic classification using Kraken 2

Classifying short reads from a metagenomic sequencing sample is a very common application. There are a few ways to do this, we'll be using Kraken 2 here.

The AMI already has everything installed, including the "minikraken" databases (which are still kind of big). Overall, this should make it easy to just transfer in your data and run an analysis!

For details on the installation of Kraken 2, see [this section of the bioinformatics sysadmin page](https://github.com/swainechen/chenlab-training/blob/main/sysadmin/bioinformatics.md#Kraken-2). For details on transferring files in to the AMI, see [this section of the AWS-GIS training website](http://slchen-lab-training.s3-website-ap-southeast-1.amazonaws.com/12-creates3sharedata/05-datatransfer-local2ec2.html). And of course you should be familiar with the [Kraken 2](https://ccb.jhu.edu/software/kraken2/) website.

## Get some sample data and run Kraken 2
The following will grab a relatively small data set from the public databases as an example. This data set is run [`SRR5978047`](https://www.ncbi.nlm.nih.gov/sra/?term=SRR5978047), sample `SRS2467593` aka `SAMN078560824`, from project `SRP112732` aka `PRJNA394877`. The experiment is `SRX3134635: human gut metagenome in setting of bone marrow transplant`.
```
mkdir -p /home/ubuntu/fastq/SRR5978047
aws s3 --no-sign-request sync s3://sra-pub-run-odp/sra/SRR5978047/ /home/ubuntu/fastq/SRR5978047
# the files are stored in .sra format, so we need to convert to fastq
cd /home/ubuntu/fastq/SRR5978047
fasterq-dump SRR5978047
for i in *.fastq; do gzip $i; done
# clean up to save some space
rm /home/ubuntu/fastq/SRR5978047/SRR5978047

kraken2 -db /usr/local/lib/Kraken2/k2_standard_8gb_20210517 SRR5978047.fastq.gz --report SRR5978047.kraken2.report > SRR5978047.kraken2.raw.txt
```

The main report will come out in `SRR5978047.kraken2.report`, and all the raw classification data in the `SRR5978047.kraken2.raw.txt` file.
