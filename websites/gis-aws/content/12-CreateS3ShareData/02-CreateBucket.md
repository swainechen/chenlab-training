+++
title = "a. Create Amazon S3 bucket"
weight = 60
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++

1.	In the AWS Management Console search bar, type **S3**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3Service.png)


2.	Choose and click on **S3** to open the **Amazon S3** Console.

3.	Click on **Create bucket** 

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3CreateBucket.png)

4.	In the **Create bucket** dialog, fill in your **Bucket name**. 
**NOTE**: The bucket name provided must be globally unique as the namespace is shared by all AWS accounts. This means that after a bucket is created, the name of that bucket cannot be used by another AWS account in any AWS Region until the bucket is deleted.   

If you get an error, change the bucket name until a unique name is accepted.

Next, select the region from the **AWS Region** field.
**NOTE**: S3 is a regional service and not associated with any specific availability zone.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3CreateBucket1.png)

Scroll down.

Read through the **Block Public Access settings for this bucket** to get an idea on bucket access policies.
Leave **Block all public access** checked (this is the default configuration).

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3CreateBucket2.png)


Scroll down.

Under the **Tags(0)-optional** section provide a Key-Value for the Tag, e.g. Key="Bucket" Value=[firstname+initials-S3], leave the rest of the configuration options with default values, and scroll down.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3CreateBucket3.png)

Next click **Create bucket**

The Bucket will now be created.

5.	On the list of buckets, scroll and verify that your new bucket has been created. Selecting the bucket will also display the bucket properties, permissions & management configurations. 

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3SelectBucket.png)

6.      Click on the bucket to access and view the bucket properties, permissions & management configurations.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3SelectBucket1.png)
