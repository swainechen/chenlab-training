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

Next select the region form the **Region** field.
**NOTE**: S3 is a regional service and not associated with any specific availability zone.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3CreateBucket1.png)

Click on **Next**

On the **Configuration options** section provide a Key-Value for the Tag, e.g. Key="Bucket" Value=[firstname+initials-S3], leave the rest of the configuration options with default values and click on **Next**

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3CreateBucket2.png)


On the **Set permissions** section read through the **Block public access (bucket settings)** to get an idea on bucket access policies.
Leave the **Block all public access** as checked which is the default configuration.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3CreateBucket3.png)

Click on **Next**

On the **Review** section do a quick review of the configurations selected.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3CreateBucket4.png)

Next click **Create bucket**

The Bucket will now be created.

5.	On the list of buckets scroll through the list of buckets and verify that your new bucket has been created. Selecting the bucket will also display the buckets Properties, Permissions & Managment configurations. 

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3SelectBucket.png)
