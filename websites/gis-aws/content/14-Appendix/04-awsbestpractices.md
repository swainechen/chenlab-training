+++
title = "d. Best Practices for AWS"
date = 2021-05-23T10:46:30-04:00
draft = false 
weight = 50
tags = ["help", "best practices"]
+++

There is a lot to track when running an AWS machine. This page has a concise list of best practices that you can use as a guide:

**Security**  
    This one is pretty straightforward.

- Do not share your account access information.
- Handle your private security keys with caution. Do not share your machine's security keys unnecessarily. Remember, if someone knows the private key associated with your machine, they can access it.
- Make sure your machine's security permissions are locked down first - then you can open them up as needed.
- Avoid making resources public whenever possible. Do not store sensitive data in publicly available storage.
- Make S3 buckets private by default. When you want to share the information there, be careful not to make resources public.

**Tagging**  
    This one is quite important. Tagging is a handy tool that AWS has developed as a method for tracking resource usage. The most practical application of tags includes the tagging of resources to a specific user, or to a specific grant. For a more in-depth tutorial on tagging, please see [this official guide](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/Using_Tags.html) and check out the Tagging section in this appendix. Most importantly:

- Make sure your tags are consistent.
- Make sure that all users of an account understand the tagging conventions.

**Budgeting**  
    Briefly, it is best to use the in-built AWS budgeting tools and alerts. These will enable you to set an overall account budget and make sure you don't surpass it. 
    For more information, check out this link: https://aws.amazon.com/aws-cost-management/

**Storage and instance upkeep**  
    There is a saying about AWS - you pay most for the resources you forget about. Cost is best managed when you make sure that you're using precisely what you need - and you toss what you do not. In short, make sure of the following:

- When you are done using your instance, stop it.
- Make sure the instance shoe fits. Don't use an overpowered machine for a small job.
- Only use the storage you need. Don't keep 10GiB of data on a 100GiB volume. Remember that you are paying for all of the storage attached to your machine, not only the part with data. 
- When saving large datasets that are rarely accessed, keep it in S3 buckets and allow to go to Glacier storage class.
