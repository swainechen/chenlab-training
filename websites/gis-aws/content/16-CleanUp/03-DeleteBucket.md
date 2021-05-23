+++
title = "b. Delete S3 Buckets"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 50
tags = ["tutorial", "aws console", "ec2"]
+++

**Deleting S3 buckets**

If you delete a bucket, it and its contents are gone forever. If you enabled versioning, those versions are gone, too. As with termination, be *absolutely sure* you're ready to delete these resources.

1.  Open the S3 dashboard and select your bucket. In the upper right-hand corner, click **Delete**.

![S3 Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3DeleteButton.png)

2.  In the dialogue box, type the name of your bucket in the space provided. Note the warnings on the page- when the bucket is deleted, it's gone forever, and other people have the opportunity to use the bucket name.

![S3 Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3delete.png)

3.  Your bucket will no longer show on the S3 menu. Success!
