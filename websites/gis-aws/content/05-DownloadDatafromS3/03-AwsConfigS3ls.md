+++
title = "a. AWS CLI commands for S3"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 50
tags = ["tutorial", "aws console", "ec2"]
+++

Now that you have access to the AWS CLI, you can use it to download files from an S3 bucket. For more information on Amazon S3, see Amazon Simple Storage Service Documentation. 

1.	Run an AWS CLI command to list all the S3 buckets associated with the account.

```bash
aws s3 ls
```

The output will list all the S3 buckets (if any) associated with the account.

For the purpose of this workshop we will be using a **named profile** with AWS CLI which is a collection of settings and credentials that you can apply to a AWS CLI command. When you specify a profile to run a command, the specific settings and credentials associated with the profile are used to run that command. This is typically followed when credentials are shared to access resources. Use the AWS CLI with the following named profile to access the S3 bucket folders shared with us for this workshop.  
Access Key ID: `AKIAWTEVRD7DDRXIXUVV`  
Secret Access Key: `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`  
Default region: `ap-southeast-1`  

2.	Enter the AWS CLI command to configure AWS Credentials with a profile named "training":

```bash
aws configure  --profile training
```

3.	Enter the **aws access key id** provided for the training.

```bash
AWS Access Key ID [None]: [Access Key ID]
```

4.	Next Enter the **aws secret key** provided for the training.

```bash
AWS Secret Access Key [None]: [Secret Access Key]

```

5.	Enter a **Default region name**.

```bash
Default region name [None]: ap-southeast-1
```

Hit enter to accept the defaults for **output format**


6. Rerun the AWS CLI S3 command to list the folders within an S3 bucket, but this time with the profile. This will show you all folders of the bucket that share the account credentials.

```bash
aws s3 ls --profile training s3://slchen-lab-transfer/GIS-training/
```

The output is a list of folders within the S3 bucket (slchen-lab-tranfer) that is shared.


Expected output:
```bash
PRE SRR6327875/
PRE SRR6327950/
```
