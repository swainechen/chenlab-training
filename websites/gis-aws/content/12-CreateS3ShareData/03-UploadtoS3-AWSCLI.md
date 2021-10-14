+++
title = "Using AWS CLI"
weight = 80
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++

1.	Run the AWS CLI command to list all the S3 buckets associated with the account in the region configured with **aws configure**.

```bash
aws s3 ls
```

Output will be a list of buckets including the new bucket which was created above.

2.	Run AWS CLI commands.   

**List the contents of the bucket**

```bash
aws s3 ls s3://{Bucket Name}
```

Verify that the output includes your uploaded files.


**Download files from the S3 bucket**

To download files we will use the AWS CLI **s3 cp** command.

```bash
aws s3 cp s3://{Bucket Name}/{object name} ./
```

Verify the download by listing the directory and confirm the object is present.

For more detailed steps on AWS CLI for S3 refer to Section VI (**Download data shared using Amazon S3**) of the workshop.


**Upload files to the S3 bucket**

 To upload files we will use the AWS CLI **s3 cp** command.

```bash
aws s3 cp ./{your-file} s3://{Bucket Name}/
```

3.	Verify the contents of the S3 bucket after the copy is completed via command line (below) or by using the AWS Console.

```bash
aws s3 ls s3://{Bucket Name}
```
