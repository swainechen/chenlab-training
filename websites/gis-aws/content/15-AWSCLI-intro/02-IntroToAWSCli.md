+++
title = "a. Introduction to AWS CLI"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 40
tags = ["tutorial", "Prerequisite", "ec2"]
+++

The [AWS CLI](https://aws.amazon.com/cli/) allows you to manage services using the command line and control services through scripts. Many users choose to conduct some level of automation using the AWS CLI.

Open a Terminal window and paste the following command. The AWS CLI updates the version, if necessary.

AWS CLI is already installed on the image, however we will update the package.

```bash
sudo apt update && sudo apt upgrade; dpkg --list | grep awscli
```

**Check Existing Amazon EC2 Instances**

Use the following commands to display:

- the general AWS CLI help,
- the help related to Amazon EC2 commands,
- the list of your existing instances with their key characteristics and
- the list of your registered SSH key-pairs.

```bash
aws help
```
```bash
aws ec2 help
```
Let us run a command to have a look at the EC2 instances in our account.

```bash
aws ec2 describe-instances
```
The command fails with a region not specified or a permissions error as it requires your AWS credentials to access other AWS resources such as EC2, S3 etc..

For general use, the aws configure command is the fastest way to set up your AWS CLI installation.

1.	run the **aws configure** command

```bash
aws configure
```

2.	Enter your accounts **aws access key id**

```bash
aws_access_key_id=[Access Key ID]
```

3.	Enter your account **aws secret key**

```bash
aws_secret_access_key=[Secret Access Key]

```

4.	Enter a **Default region name**

```bash
Default region name=us-east-1
```

Hit enter to accept the defaults for  **output format**


5.	Let us rerun the ec2 describe-instances command again and check the output.

```bash
aws ec2 describe-instances
```

The command now runs and shows as output a description of all the EC2 instances in the account for the specified region.

6.	Let us run another command to have a look at the key-pairs we have for the region selected.

```bash
aws ec2 describe-key-pairs
```
**NOTE**:Type **q** to exit the help pages.

Next, you use the AWS CLI to interact with Amazon S3.
