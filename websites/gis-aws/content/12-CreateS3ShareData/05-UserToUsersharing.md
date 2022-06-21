+++
title = "d. S3 exercises: Add'l-1: Sharing data with a collaborator"
weight = 120
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++


In this section we will use the AWS CLI (Command Line Interface) to share data between users:

1.	Check with your collaborators to see if they can list the objects in your shared bucket using the following AWS CLI command.

```bash
aws s3 ls s3://{Collaborator-Bucket-Name}
```

2.	Now let us make a file object in your bucket downloadable. Select the file you would like to share. Click on **Actions** and **Make public via ACL**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic6.png)

3.	Click on **Make public**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic7.png)

4.	Let’s try inspecting the file you made public, using a web browser. First, **click on the image name** in the S3 bucket you created and then copy the link under **Object URL** into a web browser.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsObjectPrivate1.png)

5.	Next, the files you made public should now be downloadable by your collaborator by running the AWS CLI commands.

```bash
aws s3 cp s3://{Collaborator-Bucket-Name}/{object name} ./
```

6.	Now let us grant access to a specific AWS account i.e. of your collaborator via **Access for other AWS accounts**. For this we will need to update the Access Control List to provide read and write permissions for objects only to a specific collaborator. **Follow Step 3 above** and click on **Permissions**. Scroll down to the **Access control list (ACL)** section and click on **Edit**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsObjectPrivate2.png)

7.	Next under Grantee - Everyone (public access) uncheck List (under Objects) and click on **Add grantee** under **Access for other AWS accounts**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsObjectPrivate3.png)

8.	In order to provide S3 Bucket access to a specific collaborator, ask your collaborator to provide you with an AWS Canonical ID by using the following command

```bash
aws s3api list-buckets --query Owner.ID --output text
```

9.	Enter your collaborator’s AWS Canonical ID to grant access only to your collaborator. You may provide/ revoke access to multiple AWS accounts by clicking on Add/ Remove grantee. Check the **List** and **Write** boxes under Objects and click on **Save changes**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsObjectPrivate4.png)

You have now successfully shared your S3 Bucket privately with a specific collaborator.

**Transferring a file securely to your collaborator’s bucket**

Now transfer any specific file to your collaborator’s S3 Bucket which you have been granted access to by using the following command

```bash
aws s3 cp path/to/filetocopy s3://{Collaborator-Bucket-Name}/ --acl public-read
```

