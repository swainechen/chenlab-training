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
dpkg --list | grep awscli
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

[View Point 2 of the visual for difference between "Event Engine" and "Direct sign-in"]({{< ref "01-hpc-overview/03-awsaccess.html" >}}).
