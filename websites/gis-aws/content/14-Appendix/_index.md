---
title: " Appendix"
date: 2020-04-24T7:05:54Z
weight: 650
pre: "<b>XVI ⁃ </b>"
tags: ["HPC", "EFA", "Elastic Fabric Adapter", "Network", "MPI"]
---

## Glossary

**AMI** is short for Amazon Machine Image - an encrypted template of your cloud machine, functionally similar to the snapshot of a computer's root drive. The image includes an operating system, additional software and layers of your application.

**AWS CLI** is AWS Command Line Interface (CLI) - a unified tool to manage your AWS services. With just one tool to download and configure, you can control multiple AWS services from the command line and automate them through scripts. Click here to know more - [AWS CLI](https://aws.amazon.com/cli/) 

**Bucket** is, in Amazon terms, is a container for stored objects. 

**Cloud** is a combination of a data center and a system of computational resources available to many users over the internet. Click here to know more - [Cloud](https://en.wikipedia.org/wiki/Cloud_computing)

**CPU**, also called a processor, is a central processing unit of your machine, performing calculations and executing instructions from RAM. Modern processors almost always have multiple _cores_, which divide the workload of CPU and execute tasks in parallel. Some processors also employ _virtual_ cores, which is essentially executing several streams of instructions concurrently on a single physical core to utilize the available resources more efficiently.

**EC2** is short for Amazon Elastic Compute Cloud (Amazon EC2) - a web service that provides secure, resizable compute capacity in the cloud. Click here to know more - [EC2](https://aws.amazon.com/ec2/?ec2-whats-new.sort-by=item.additionalFields.postDateTime&ec2-whats-new.sort-order=desc)

**EC2 instance** is a copy of AMI running as a virtual server in the AWS cloud.

**EC2 key-pair**, consisting of a private key and a public key, is a set of security credentials that you use to prove your identity when connecting to an instance. Click here to know more - [EC2 Key-Pairs](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-key-pairs.html)

**Private & Public keys** are a part of a widespread encryption system, in which each user controls a pair of keys - _public_ (which may be made available to others) and _private_ (which should never be known by anyone but the owner). In this system, a message is encrypted using the recipient's public key, but can only be decrypted and read using the recipient's private key. Click here to know more - [Public-key cryptography](https://en.wikipedia.org/wiki/Public-key_cryptography)

**RAM** is short for Random Access Memory, the working memory of your computer (as opposed to _storage_ memory), which has almost instant access to the CPU and is used to stored working data and machine code. In other words, the amount of RAM directly defines the maximum size of workload the computer can handle at any given moment.

**Regions & Availability Zones**: Amazon EC2 is hosted in multiple locations worldwide. These locations are composed of Regions, Availability Zones, Local Zones, and Wavelength Zones. Each Region is a separate geographic area. Click here to know more - [Regions & Availability Zones](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/using-regions-availability-zones.html) 

**S3** is short for Amazon Simple Storage Service, a storage service for the internet. Used to deposit and retrieve any amount of data from anywhere on the web.

**SSH**, or a Secure Socket Shell, is a network protocol that is used to provide secure access to a computer over an unsecured network.

## FAQ

**Which amount of CPU and RAM should I choose for my machine? Which EC2 instance do I choose?**  
Of course, that depends on your application but just for a simple use case such as installing packages and moving small files even a small instance from _T2_ family of 2-8 GB RAM and 1-2 vCPUs might be enough. If you are looking to run memory-intensive applications (e.g. processing big amounts of data), look into memory-optimized _R5_ or _R5a_ instances, where 16-64 GB RAM and 2-8 vCPUs are reasonable to start with. If you are running applications that are designed to utilize multiple cores, compute-optimized _C5_ instances offer a combination of higher number vCPUs vs RAM, with options like 16 vCPUs per 32 GB of memory. 

**What is an SSH client? Which one should I use to access my machine?** SSH client is a tool that allows you to connect to the remote machine (such as an initiated AWS EC2 instance - a computer in the cloud) from your local computer. On Mac and Linux users, ssh clients should be already preinstalled as part of core utilities and accessible via the terminal with the command `ssh`. On Windows, [Putty](https://www.putty.org/) is a popular choice, however, if you want the full functionality of  a UNIX-like command line, newer Windows versions come with the ability to install [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10), which emulates a Linux machine using your Windows resources. 

