+++
title = "a. Introduction to AWS CLI"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 40
tags = ["tutorial", "Prerequisite", "ec2"]
+++

The [AWS CLI](https://aws.amazon.com/cli/) allows you to manage services using the command line and control services through scripts. Many users choose to conduct some level of automation using the AWS CLI.

Connect to an EC2 instance and paste the following command. AWS CLI is already installed on your machine, but this will initiate an update. There are a few ways to check this, depending on how it was installed

```bash
sudo apt update && sudo apt upgrade; dpkg --list | grep awscli
pip3 list | grep awscli
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

The above `describe-instances` and `describe-key-pairs` commands will likely fail because we did not set our AWS account's credentials: the access keys and the region. This results in a permissions error as credentials are required to access AWS resources such as EC2, S3, and so on (without them AWS CLI does not know which account you are referring to). You can enter those manually using the `aws configure` command.
  
As we have not yet configured the instance with the AWS account's credentials, the **~/.aws/** folder will not be available. Now let us go ahead with **aws configure**.
  
_Important: if you are using Event Engine for this workshop, you can copy the code snippet with the temporary account's credentials from the Console page and paste it in your command line instead. This will be analogous to going through the manual setting up described in points 1-4 below (you can skip to point 5 directly)._

![EE console](/images/hpc-aws-parallelcluster-workshop/EE_console_login.png)    

For general use, the **aws configure** command is the fastest manual way to set up your AWS CLI and connect it to your account.

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
