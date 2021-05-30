+++
title = "a. Introduction to AWS CLI"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 40
tags = ["tutorial", "Prerequisite", "ec2"]
+++

The [AWS CLI](https://aws.amazon.com/cli/) allows you to manage services using the command line and control services through scripts. Many users choose to conduct some level of automation using the AWS CLI.

Connect to an EC2 instance and paste the following command. AWS CLI is already installed on your machine, but this will initiate an update. 

```bash
sudo apt update && sudo apt upgrade; dpkg --list | grep awscli
```

**Check Existing Amazon EC2 Instances**

Use the following commands to explore AWS CLI (use **q** to exit windows):

- general AWS CLI help
- help related to Amazon EC2 commands
- the list of your existing instances with their key characteristics
- the list of your registered SSH key-pairs

```bash
aws help
```
```bash
aws ec2 help
```
```bash
aws ec2 describe-instances
```
```bash
aws ec2 describe-key-pairs
```

The above commands may fail when a region is not specified, or there may be a permissions error as it requires your AWS credentials to access other AWS resources such as EC2, S3, and so on. You can enter this information manually. If you are using Event Engine for this workshop, you can copy the information on the Console page and paste it instead.

For general use, the **aws configure** command is the fastest way to set up your AWS CLI installation.

1.	run the **aws configure** command

```bash
aws configure
```

2.	Enter your account **aws access key id**

```bash
aws_access_key_id=[Access Key ID]
```

3.	Enter your account **aws secret access key**

```bash
aws_secret_access_key=[Secret Access Key]

```

4.	Enter a **Default region name**. For the purposes of this lab, use ap-southeast-1.

```bash
Default region name=ap-southeast-1
```

Hit enter to accept the defaults for  **output format**


5.	Let us rerun the ec2 **describe-instances** command again and check the output.

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
