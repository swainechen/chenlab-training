+++
title = "c. Make objects public"
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
Observe the first button **Bucket Public Access** is selected by default and shows the **Block all public access** setting is **on**. Note this only allows objects in this bucket to be public; but you still have to make the individual objects public as well.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissions.png)

5.	Click on **Edit**.

6.	Uncheck the **Block all public access** check box and click on **Save changes**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic1.png)

7.	Type in "confirm" in the text box and click on **Confirm**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic2.png)

You will now be shown a message indicating that public access settings to the bucket have been updated successfully.

Next we will update the Access Control List to provide read and write permissions for individual objects.

8.	Scroll down to the **Access control list (ACL)** section and click on **Edit**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic3.png)  

9.	Next under **Grantee** - **Everyone (public access)** check **List** (under Objects).

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic4.png)  

10.	Check the **I understand the effects of these changes on my objects and buckets** box and then click on **Save changes**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic5.png)  

At this point you have effectively shared the S3 bucket. Check with your collaborators to see if they can list the objects in your shared bucket using the following AWS CLI command.

```bash
aws s3 ls s3://{Collaborator-Bucket-Name}
```

11.	Being able to list an object in your S3 bucket is different from being able to download it. Now let us make a file object downloadable. Select the file you would like to share. Click on **Actions** and **Make public**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic6.png)

12.	Click on **Make public**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic7.png)

13.	Now let’s try inspecting the image you made public, using a web browser. First, **click on the image name** in the S3 bucket you created and then copy the **Object URL** into a web browser.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsObjectPrivate1.png)

14.	Next, the files you made public should now be downloadable by your collaborator by running the AWS CLI commands.

```bash
aws s3 cp s3://{Collaborator-Bucket-Name}/{object name} ./
```

15.	Now let us grant access to a specific AWS account i.e. of your collaborator via **Access for other AWS accounts**. For this we will need to update the Access Control List to provide read and write permissions for objects only to a specific collaborator. **Follow Step 3 above** and click on **Permissions**. Scroll down to the **Access control list (ACL)** section and click on **Edit**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsObjectPrivate2.png)

16.	Next under Grantee - Everyone (public access) uncheck List (under Objects) and click on **Add grantee** under **Access for other AWS accounts**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsObjectPrivate3.png)

17.	In order to provide S3 Bucket access to a specific collaborator, ask your collaborator to provide you with an AWS Canonical ID by using the following command

```bash
aws s3api list-buckets --query Owner.ID --output text
```

18.	Enter your collaborator’s AWS Canonical ID to grant access only to your collaborator. You may provide/ revoke access to multiple AWS accounts by clicking on Add/ Remove grantee. Check the **List** and **Write** boxes under Objects and click on **Save changes**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsObjectPrivate4.png)

You have now successfully shared your S3 Bucket privately with a specific collaborator.

**Transferring a file securely to your collaborator’s bucket**

Now transfer any specific file to your collaborator’s S3 Bucket which you have been granted access to by using the following command

```bash
aws s3 cp path/to/filetocopy s3://{Collaborator-Bucket-Name}/ --acl public-read
```

