+++
title = "c. Make bucket public"
weight = 90
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++

You will now learn how to share an S3 Bucket with your collaborator and make the S3 Objects public.

1.	In the AWS Management Console search bar, type **S3**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3Service.png)

2.	Click on **S3** to open the **Amazon S3** Console.

3.	Scroll through the list of buckets and click on the bucket you created.
**NOTE**: Make sure you click on the bucket name link and not just select the bucket using the checkbox.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3SelectBucketLink.png)

4.	Click on **Permissions** to view the bucket access configurations.
Observe under the **Block Public Access (bucket settings)** the **Block all public access** setting is **on**. Note this only allows objects in this bucket to be public; but you still have to make the individual objects public as well.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissions.png)

5.	Click on **Edit**.

6.	Uncheck the **Block all public access** check box and click on **Save changes**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic1.png)

7.	Type in "confirm" in the text box and click on **Confirm**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic2.png)

You will now be shown a message indicating that public access settings to the bucket have been updated successfully.

Next we will update the Access Control List to provide read and write permissions for individual objects.

8.	Now click on **Edit** against **Object Ownership**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic2-2.png)

9.	Select **ACLs enabled**, tick the **I acknowledge that ACLs will be restored** and **Save changes**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic2-3.png)

10.	Scroll down to the **Access control list (ACL)** section and click on **Edit**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic3.png)

11.	Next under **Grantee** - **Everyone (public access)** check **List** (under Objects).

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic4.png)

12.	Check the **I understand the effects of these changes on my objects and buckets** box and then click on **Save changes**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic5.png)

At this point you have effectively shared the S3 bucket. Being able to list an object in your S3 bucket is different from being able to download it. This needs object permissions via the ACL (Access control list).

_We will now run a larger analyses using the [long-read RNA-Seq data]({{< ref "12-creates3sharedata/07-bambuday2.html" >}}) and [Single-cell data]({{< ref "12-creates3sharedata/08-singlecell.html" >}}) to explore the following features of S3
-Downloading data from a publicly shared S3 bucket
-S3 access control lists (ACLs) to enable public read access
-Accessing an object in the S3 bucket through its URL_
