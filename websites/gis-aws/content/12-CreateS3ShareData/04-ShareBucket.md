+++
title = "c. Share S3 bucket with another Researcher"
weight = 60
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++


We will now share files with another fellow team member/researcher via an S3 bucket.

1.	In the AWS Management Console search bar, type **S3**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3Service.png)

2.	Choose and click on **S3** to open the **Amazon S3** Console.

3.	Scroll through the list of buckets and choose and click on the bucket name created.**NOTE**: make sure you click on the bucket name link and not just selecting the bucket using the checkbox.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3SelectBucketLink.png)

4.	Click on **Permissions** to view the bucket access configurations.

Observe the the first button **Bucket Public Access** is selected by default which shows the **Block all public access** setting is **on**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissions.png)


5.	Click on **Edit**

6.	Uncheck the **Bloc all public acces** check box and click on **Save**

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic1.png)

7.	Type in "confirm" in the text box and click on **Confirm**

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic2.png)

You will now be shown a message indicating that Public access settings to the bucket have been updated successfully.

**Next we will update the Access control list to provide permisions such as allowing to list object, write obejcts**


9.	Next click on the **Access Control List** button.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic4.png)


10.	Next select **Everyone** Group under **Public access** section and further check **List objects** on the dialog titled **Everyone** and click on **Save**

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic5.png)

Verify that the value under **List Objects** column for the **Everyone** Group under the **Public access** section shows **Yes**


11.	Select the file object from the list you wish make public and share for others to be able to download it. 

You will get slide in Dialog showing the Overview, Properties and Permissions of the object.

Click on **Permissions**

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsObjectPublic1.png)

12.	On the Permissions page of the object select **Everyone** under the **Public access** section and further check the **Read object** check box on the Dialog labled **Everyone** and click on **Save**

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsObjectPublic2.png)


Verify that the **Read Object** valye is set to **Yes** for the **Everyone** Group under **Public access**

At this point you have effectivly shared the file object, you can check with your team/fellow Researcher to check if they can list the objects in your shared bucket. and try and dwonload the files you have shared by running AWS CLI commands.
