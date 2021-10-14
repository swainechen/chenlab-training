+++
title = "Access to a specific account"
weight = 110
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++

Earlier we made the S3 bucket available to **Everyone (public access)**. Now let us grant access to a specific AWS account i.e. of your collaborator via **Access for other AWS accounts**

1.	Once again click on the bucket created earlier and then on Permissions to view the bucket access configurations. 

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsObjectPrivate1.png)

2.	Next, we will update the Access Control List to provide read and write permissions for objects only to a specific collaborator. Scroll down to the Access control list (ACL) section and click on Edit.  

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsObjectPrivate2.png)

3.	Next under Grantee - Everyone (public access) uncheck List (under Objects) and click on **Add grantee** under **Access for other AWS accounts**.  

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsObjectPrivate3.png)

4.	In order to provide S3 Bucket access to a specific collaborator, ask your collaborator to provide you with an AWS Canonical ID by using the following command  

```bash
aws s3api list-buckets --query Owner.ID --output text
```

5.	Enter your collaborator’s AWS Canonical ID to grant access only to your collaborator. You may provide/ revoke access to multiple AWS accounts by clicking on Add/ Remove grantee. Check the **List** and **Write** boxes under Objects and click on **Save changes**.  

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsObjectPrivate4.png)

You have now successfully shared your S3 Bucket privately with a specific collaborator.

    
   
**Transferring a file securely to your collaborator’s bucket**

Now transfer any specific file to your collaborator’s S3 Bucket which you have been granted access to by using the following command

```bash
aws s3 cp path/to/filetocopy s3://{Collaborator-Bucket-Name}/
```

>_**Note**: Though you will be able to transfer a file securely to your collaborator's bucket, you will not be able to **copy, download or view** the file transferred to your bucket (by your collaborator), which needs special **bucket permissions**._
