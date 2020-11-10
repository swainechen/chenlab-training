+++
title = "c. Share S3 bucket with another Researcher"
weight = 60
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++


We will now share files with your collaborators via an S3 bucket.

1.	In the AWS Management Console search bar, type **S3**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3Service.png)

2.	Choose and click on **S3** to open the **Amazon S3** Console.

3.	Scroll through the list of buckets and choose and click on the bucket name created. **NOTE**: Make sure you click on the bucket name link and not just select the bucket using the checkbox.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3SelectBucketLink.png)

4.	Click on **Permissions** to view the bucket access configurations.

Observe the first button **Bucket Public Access** is selected by default which shows the **Block all public access** setting is **on**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissions.png)


5.	Click on **Edit**

6.	Uncheck the **Block all public access** check box and click on **Save**

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic1.png)

7.	Type in "confirm" in the text box and click on **Confirm**

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic2.png)

You will now be shown a message indicating that Public access settings to the bucket have been updated successfully.

**Next we will update the Access control list to provide permisions such as allowing to list and write objects**


8.	Next click on the **Access Control List** button.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic4.png)


9.	Next select **Everyone** Group under **Public access** section and further check **List objects** on the dialog titled **Everyone** and click on **Save**

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic5.png)

Verify that the value under **List Objects** column for the **Everyone** Group under the **Public access** section shows **Yes**


10.	Select the file object from the list you wish to make public and share, for others to be able to download it. 

A dialog box appears with the Overview, Properties and Permissions of the object.

Click on **Permissions**

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsObjectPublic1.png)

11.	On the Permissions page of the object select **Everyone** under the **Public access** section and further check the **Read object** check box in the dialog box labelled **Everyone** and click on **Save**

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsObjectPublic2.png)


Verify that the **Read Object** value is set to **Yes** for the **Everyone** Group under **Public access**

At this point you have effectively shared the file object. You can check with your collaborators to check if they can list the objects in your shared bucket. Try and download the files you shared by running the AWS CLI commands.
