+++
title = "AWS Configure - General Use"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 60
tags = ["tutorial", "Prerequisite", "ec2"]
+++

1.	Run the **aws configure** command

```bash
aws configure
```

2.	Enter your account's **aws access key id**

```bash
aws_access_key_id=[Access Key ID]
```

3.	Enter your account's **aws secret access key**

```bash
aws_secret_access_key=[Secret Access Key]
```

4.	Enter a **default region name**. For the purposes of this workshop, use _ap-southeast-1_ (the region for Singapore).

```bash
Default region name=ap-southeast-1
```

Hit enter to accept the defaults for **output format**


5.	Run the ec2 **describe-instances** command and check the output.

```bash
aws ec2 describe-instances
```

This gives a description of all the EC2 instances in the account for the specified region.

6.	In addition, take a look at the key-pairs we have for the selected region.

```bash
aws ec2 describe-key-pairs
```

**Note:** This is a very brief introduction to AWS CLI. With great power comes great responsibility, so familiarise yourself and practice before using it to automate instances. Your wallet will thank you.  

Next, you will use the AWS CLI to interact with Amazon S3.
