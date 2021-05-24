+++
title = "b. Best Practices for AWS"
date = 2021-05-23T10:46:30-04:00
draft = false 
weight = 50
tags = ["help", "best practices"]
+++

There's a lot to track when running an AWS machine. This page has a concise list of best practices that you can use as a guide. The key topics covered here include:

- Security

- Tagging

- Budgeting

- Storage and instance upkeep


**Security**
    This one is pretty straight forward.

- Don't share your account access information.

- Make sure your machine's security permissions are locked down first- then you can open them up as needed.

- Avoid making resources public whenever possible.

- Make S3 buckets private by default. When you want to share that information, be careful not to make resources public.


**Tagging**
    This one is quite important. Tagging is a handy tool that AWS has developed as a method for tracking resource usage. The most practical application of tags includes the tagging of resources to a specific user, or to a specific grant. For a more in-depth tutorial on tagging, please see this link. In short:

- Make sure your tags are consistent.

- Make sure that all uses of an account udnerstand the tagging conventions.


**Budgeting**
    Briefly, it's best to use the in-built AWS budgeting tools and alerts. These will enable you to set an overall account budget and make sure you don't surpass it. 
    For more information, check out this link: https://aws.amazon.com/aws-cost-management/


**Storage and instance upkeep**
    There's a saying about AWS- you pay most for the resources you forget about. Cost is best managed when you make sure that you're using precisely what you need- and you toss what you don't. In short, make sure of the following:

- When you're done using your instance, stop it.

- Make sure the instance shoe fits. Don't use an overpowered machine for a small job.`

- Only use the storage you need. Don't keep 10GiB of data on a 100GiB volume.

- When saving large datasets that are rarely accessed, keep on S3 and allow to go to glacier.
