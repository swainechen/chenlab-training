---
title: " Start a Linux machine from a shared AMI"
date: 2019-01-24T09:05:54Z
weight: 40 
pre: "<b>III ⁃ </b>"
tags: ["HPC", "Introduction", "EC2", "Optional"]
---

This workshop section walks you through setting up your own EC2 Linux instance from a premade AMI, which you will use to run analyses in later sections of the workshop.

We'll be using an image built by the Chen Lab at NUS and GIS. This has a lot of common tools useful for bacterial (and general) genomics installed. We have documented [what is on this AMI and how it was set up](https://github.com/swainechen/chenlab-training/tree/main/sysadmin) on our [accompanying GitHub repository](https://github.com/swainechen/chenlab-training). Those instructions may be useful if you want to set up the same software on another system and may help you with some hints when installing other software as well for your own work.

Specifically, you learn how to:

**a.	[Sign in to the AWS Management Console and explore it.]({{< ref "03-ec2fromsharedami/04-aws-ec2-dashboard.html" >}})**

**b.	[Create an Amazon EC2 instance from an AMI.]({{< ref "03-ec2fromsharedami/06-launchec2instance.html" >}})**

**c.	[SSH into an EC2 instance and run Linux commands.]({{< ref "03-ec2fromsharedami/07-sshec2instance.html" >}})**

<!-- Total cost to run this lab is expected to be less than **two dollars should you run it on your own or company account** if the guide is followed as discussed, in one sitting, and including the deletion of all resources outlined through the guide. -->
