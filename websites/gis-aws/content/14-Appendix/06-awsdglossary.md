+++
title = "f. AWS Glossary"
date = 2021-05-23T16:24:30-04:00
draft = false 
weight = 70
tags = ["help", "glossary", "terms"]
+++

## Glossary

This is a short list of the most widely used cloud computing and AWS-related terms in this workshop that you might want to familiarize yourself with. It should be enough for the purposes of the current course, however, here is [the complete curated list of AWS references](https://docs.aws.amazon.com/general/latest/gr/glos-chap.html) by Amazon team.

**Access key** is, in AWS terms, the combination of an access key ID (for example, `AKIAIOSFODNN7EXAMPLE`) and a secret access key (for example, `wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY`). You use access keys to sign API requests that you make to AWS. It is a unique pair of cryptographic identifiers used as a signature of sender for server requests, a way of protecting the request from being altered.

**Access key ID** is a unique identifier that is associated with a secret access key.

**AMI** is short for Amazon Machine Image - an encrypted template of your cloud machine, functionally similar to the snapshot of a computer's root drive. The image includes an operating system, additional software and layers of your application.

**AWS CLI** is AWS Command Line Interface (CLI) - a unified tool to manage your AWS services. With just one tool to download and configure, you can control multiple AWS services from the command line and automate them through scripts. Click here to know more - [AWS CLI](https://aws.amazon.com/cli/)

**Bucket** is, in Amazon terms, is a container for stored objects.

**Cloud** is a combination of a data center and a system of computational resources available to many users over the internet. Click here to know more - [Cloud](https://en.wikipedia.org/wiki/Cloud_computing)

**CPU**, also called a processor, is a central processing unit of your machine, performing calculations and executing instructions from RAM. Modern processors almost always have multiple _cores_, which divide the workload of CPU and execute tasks in parallel. Some processors also employ _virtual_ cores, which is essentially executing several streams of instructions concurrently on a single physical core to utilize the available resources more efficiently.

**EBS volume** is a storage unit of Elastic Block Store - a block storage system employed by AWS, typically used in conjunction with EC2 instances. Functionally, this an analogy of the hard drive of your typical desktop machine.

**EC2** is short for Amazon Elastic Compute Cloud (Amazon EC2) - a web service that provides secure, resizable compute capacity in the cloud. Click here to know more - [EC2](https://aws.amazon.com/ec2/?ec2-whats-new.sort-by=item.additionalFields.postDateTime&ec2-whats-new.sort-order=desc)

**EC2 instance** is a copy of AMI running as a virtual server in the AWS cloud.

**EC2 key-pair**, consisting of a private key and a public key, is a set of security credentials that you use to prove your identity when connecting to an instance. Click here to know more - [EC2 Key-Pairs](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-key-pairs.html)

**Private & Public keys** are a part of a widespread encryption system, in which each user controls a pair of keys - _public_ (which may be made available to others) and _private_ (which should never be known by anyone but the owner). In this system, a message is encrypted using the recipient's public key, but can only be decrypted and read using the recipient's private key. Click here to know more - [Public-key cryptography](https://en.wikipedia.org/wiki/Public-key_cryptography)

**RAM** is short for Random Access Memory, the working memory of your computer (as opposed to _storage_ memory), which has almost instant access to the CPU and is used to stored working data and machine code. In other words, the amount of RAM directly defines the maximum size of workload the computer can handle at any given moment.

**Regions & Availability Zones**: Amazon EC2 is hosted in multiple locations worldwide. These locations are composed of Regions, Availability Zones, Local Zones, and Wavelength Zones. Each Region is a separate geographic area. Click here to know more - [Regions & Availability Zones](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/using-regions-availability-zones.html)

**S3** is short for Amazon Simple Storage Service, a storage service for the internet. Used to deposit and retrieve any amount of data from anywhere on the web.

**Secret access key** is a key that's used in conjunction with the access key ID to cryptographically sign programmatic AWS requests.

**Session token** is a credential for temporary, limited-privilege API access during a particular session for [AWS Identity and Access Management (IAM)](https://docs.aws.amazon.com/general/latest/gr/glos-chap.html#IAM) users or for users that you authenticate. Used in conjunction with an access key and a secret access key.

**SSH**, or a Secure Socket Shell, is a network protocol that is used to provide secure access to a computer over an unsecured network.

