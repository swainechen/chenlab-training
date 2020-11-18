+++
title = "b. Download data from S3 bucket"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 50
tags = ["tutorial", "aws console", "ec2"]
+++

Now that you have access to shared S3 buckets we are now ready to download data from S3 buckets to your machine.

1.	Run the command to list all the S3 buckets

```bash
aws --profile training s3 ls --region ap-southeast-1
```
you will find the bucket named **slchen-lab-transfer** in the list, this bucket contains multiple folders. 

**NOTE**: Only a specific few folders may have been shared,

We will now copy files from a specific location within this folder to your machine

Before that, lets create a folder on your machine to where the data will be copied to.

2.	Run the following commands to create a directory to hold the data and further cd to that directory.

```bash
mkdir /tmp/fastq/SRR6327950
```

```bash
cd /tmp/fastq/SRR6327950
```

3.	Run the AWS CLI command to copy files from the S3 bucket to the directory just created

```bash
aws --profile training s3 cp s3://slchen-lab-transfer/GIS-training/SRR6327950/SRR6327950_1.fastq.gz . --region ap-southeast-1
```

```bash
aws --profile training s3 cp s3://slchen-lab-transfer/GIS-training/SRR6327950/SRR6327950_2.fastq.gz . --region ap-southeast-1
```

4.	Verify the MD5 checksum as you would normally do.

5.	Download additional data using regular utilities

```bash
ena-fast-download.py SRR6327950
```	

We are now ready to run the analyses.
