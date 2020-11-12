+++
title = "a. Delete all of your extra resources"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 50
tags = ["tutorial", "aws console", "ec2"]
+++

**Terminating EC2 Instances**
We've covered this before, but just to make sure:

1.  Select your EC2 instance from the available list. Click **Actions**, **Instance State**, and then **Terminate**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2InstanceTerminate-2.png)

2.  When prompted for verification, select choose "Yes, Terminate".

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2InstanceTerminateConfirm-2.png)

3.  The instance state will now show "Terminated".

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
