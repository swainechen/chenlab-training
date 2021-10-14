+++
title = "g. FAQ"
date = 2021-05-23T16:24:30-04:00
draft = false 
weight = 80
tags = ["help", "FAQ"]
+++

## FAQ

**Which amount of CPU and RAM should I choose for my machine? Which EC2 instance do I choose?**  
Of course, that depends on your application but just for a simple use case such as installing packages and moving small files even a small instance from _T2_ family of 2-8 GB RAM and 1-2 vCPUs might be enough. If you are looking to run memory-intensive applications (e.g. processing big amounts of data), look into memory-optimized _R5_ or _R5a_ instances, where 16-64 GB RAM and 2-8 vCPUs are reasonable to start with. If you are running applications that are designed to utilize multiple cores, compute-optimized _C5_ instances offer a combination of higher number vCPUs vs RAM, with options like 16 vCPUs per 32 GB of memory.

**What is an SSH client? Which one should I use to access my machine?**  
SSH client is a tool that allows you to connect to the remote machine (such as an initiated AWS EC2 instance - a computer in the cloud) from your local computer. On Mac and Linux users, ssh clients should be already preinstalled as part of core utilities and accessible via the terminal with the command `ssh`. On Windows, [Putty](https://www.putty.org/) is a popular choice, however, if you want the full functionality of  a UNIX-like command line, newer Windows versions come with the ability to install [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10), which emulates a Linux machine using your Windows resources.

**What is the username I should use for the different AMIs?**
This depends on the AMI you're using. A good summary of the main ones can be found at [this link](https://alestic.com/2014/01/ec2-ssh-username/)

**Which of these services are availability zone or region limited?**  
Essentially, all of them - EC2 instances availability, EBS volumes and S3 buckets pricing might differ between geographic regions. Most of the regions, however, would have a very similar list of available services and you should be able to run your core setup of EC2-EBS-S3 anywhere in the world where AWS is available.

**What is the current pricing list?**  
AWS provide a convenient tool for exploring your costs (you will need to sign into your AWS account to use it): [AWS cost explorer](https://aws.amazon.com/aws-cost-management/aws-cost-explorer/?c=4&pt=1). However, probably the most expensive and most widely used part of AWS is the EC2 instances in so-called "on-demand mode", where you select a machine according to your preferred specifications pay for it for as long as it is running. Here is the [on-demand pricing page for EC2 instances](https://aws.amazon.com/ec2/pricing/on-demand/) where you can familiarize yourself with the costs of different machines. Do not forget to set the appropriate region and be sure to check out various instance types (e.g. compute- or memory- optimized).  

**I liked the AMI we used during the course. Can I use it in my own AWS projects later?**  
Absolutely! Right now we are polishing our AMI and making it public, and we are also determined to keep supporting it in the future. Feel free to use is as a base for your own experiments and to advertise it further. You will always be able to find the latest version by searching for "CHENLAB-PUBLIC" when you're launching an instance.

