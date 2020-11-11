---
title: "Cleaning up"
date: 2019-01-24T09:05:54Z
weight: 600 
pre: "<b>XV ⁃ </b>"
tags: ["HPC", "Overview", "Batch"]
---

In this final segment of the workshop, you will learn how to clean up after yourself. In AWS, it's easy to spin up instances, make storage volumes, and security groups. The only problem? They can drain your budget if you forget to delete the resources you don't use. In this section, we will:

-   Terminate EC2 instances
-	Delete S3 buckets
-	Make room in the budget for pizza

**Terminating EC2 Instances**
We've covered this before, but just to make sure:

1.  open the menu for instance, choose instance state, and terminate.

2.  When prompted for verification, choose "Yes, Terminate".

It'll still be visible for a bit, but then the entry will be deleted.

**Deleting S3 buckets**

If you delete a bucket, it and its contents are gone forever. If you enabled versioning, those versions are gone, too. As with termination, be *absolutely sure* you're ready to delete these resources.

-  show how to delete bucket

**Deleting Snapshots and Images**

- show how to delete snapshots and image

**Are we done yet?**

Revisit billing and costs- are your resources gone?

Website will stay publically available.

Contact us for future info.

And now, we have arrived at the end of the workshop!
