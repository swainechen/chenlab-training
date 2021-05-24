+++
title = "c. Tagging 101"
date = 2021-05-23T16:24:30-04:00
draft = false
weight = 55 
tags = ["help", "ami", "tagging"]
+++

Tagging is an often overlooked but incredibly crucial tool to using AWS. When managing large groups or multiple projects, this will allow you to track your resources in multiple dimensions.  

**Tagging Basics**  
Tags are composed of keys and values. The key the overall category that you want to track- say "grant", or "name", or "user". The value is the associated category value, such at "0037592", "PopgenMachine", or "Fred". Tags are not added automatically by AWS, and are instead dependent on the user. Within this tutorial, all resources can be tagged (EC2, AMIs, Snapshots, EBS Volumes). You can use up to 50 tags on each resource, so go wild. AWS has a handy [tag editor function](https://docs.aws.amazon.com/ARG/latest/userguide/tag-editor.html) that can help you search for resources and apply tags.  

**Tagging Examples**  
Tagging resources by hand  
- giving an instance tags (Name as default tag, add grant tag, add user tag, add project tag)  
- giving an EBS volume and AMI a tag (link back to lesson for how to make/attach these guys)  
- pulling all resources with a tag in cost explorer  

Tagging resources using the tag editor function  
- finding the tag editor  
- searching for resources  
- adding tags  
- pulling all resources in cost explorer  

**Further Reading**  
As with many things, AWS has a page for that. [Here](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/Using_Tags.html), you can read about tagging in more detail.  
