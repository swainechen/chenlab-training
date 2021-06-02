+++
title = "c. Share S3 bucket with another Researcher"
weight = 60
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++


We will now share files with your collaborators via an S3 bucket.

1.	In the AWS Management Console search bar, type **S3**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3Service.png)

2.	Click on **S3** to open the **Amazon S3** Console.

3.	Scroll through the list of buckets and click on the bucket you created.  
  **NOTE**: Make sure you click on the bucket name link and not just select the bucket using the checkbox.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3SelectBucketLink.png)

4.	Click on **Permissions** to view the bucket access configurations.

Observe the first button **Bucket Public Access** is selected by default and shows the **Block all public access** setting is **on**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissions.png)


5.	Click on **Edit**.

6.	Uncheck the **Block all public access** check box and click on **Save changes**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic1.png)

7.	Type in "confirm" in the text box and click on **Confirm**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic2.png)

You will now be shown a message indicating that public access settings to the bucket have been updated successfully.

**Next we will update the Access Control List to provide read and write permissions for objects**


8.	Scroll down to the **Access control list (ACL)** section and click on **Edit**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic3.png)


9.	Next under **Grantee** - **Everyone (public access)** check **List** (under Objects).

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic4.png)


10.	Check the **I understand the effects of these changes on my objects and buckets** box and then click on **Save changes**. 

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic5.png)

At this point you have effectively shared the S3 bucket. Check with your collaborators to see if they can list the objects in your shared bucket.   


11.	Now let us make a file object downloadable. Select the file you would like to share. Click on **Actions** and **Make public**.  

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic6.png)  


12.	Click on **Make public**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic7.png)

Only the files you made public should now be downloadable by running the AWS CLI commands.
