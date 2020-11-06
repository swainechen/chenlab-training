+++
title = "a. AWS CLI commands for S3"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 50
tags = ["tutorial", "aws console", "ec2"]
+++

Now that you have access to the AWS CLI, you can use it to download files from an S3 bucket. Fore more information on Amazon S3, see Amazon Simple Storage Service Documentation. 

1.	Run an aws cli command to list all the S3 buckets associated with the account.

```bash
aws s3 ls
```

The output will list all the S3 buckets associated with the account.


For the prupose of this workshop we will be using a **named profile** with AWS CLI which is a collection of settings and credentials that you can apply to a AWS CLI command. When you specify a profile to run a command, the specific settings and credentials associated with the profile are used to run that command.This is typically followed when credentials are shared to access resources. We will be using a named profile to use AWS CLI to access the S3 bucket folders shared with us for this workshop.


2.	Enter the AWS CLI command to configure AWS Credentials with a profile named "training" :

```bash
aws configure  --profile training
```

3.	Enter the **aws access key id** provided for the trianing.

```bash
aws_access_key_id=[Access Key ID]
```

4.	Next Enter the **aws secret key** provided for the training.

```bash
aws_secret_access_key=[Secret Access Key]

```

5.	Enter a **Default region name**

```bash
Default region name=us-east-1
```

Hit enter to accept the defaults for **output format**


6. Rerun the aws cli S3 command run earlier but this time with the profile to list all the S3 buckets associated with the shared account credentials.

```bash
aws --profile training s3 ls
```

The output is a list of S3 buckets that are shared.


E.g. output:
```bash
2013-07-11 17:08:50 mybucket
2013-07-24 14:55:44 mybucket2
```
