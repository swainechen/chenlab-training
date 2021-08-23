+++
title = "c. Share S3 bucket with another Researcher"
weight = 60
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++


A. Before we share an S3 bucket with a collaborator, let us first  
-Copy data from your local machine to your AWS instance and vice versa  
-Copy data from one instance to another instance. 

**Copy data from your local machine to your AWS instance and vice versa**
1.	On your local computer, open a new terminal and create a sample text file. 

```bash
echo “Hello, world” > test.txt
```
  
2.	Use scp to copy this sample test file into the home directory of your remote machine/ instance, by running the following command on your local machine. 
```bash
scp -i path/to/key.pem path/to/test.txt ubuntu@<ip-address>:/home/ubuntu/
```

If you are using windows command prompt and have the PuTTY suite of utilities installed, run the following command   
```bash
pscp -P 22 -i path/to/key.ppk path/	to/test.txt ubuntu@<ip-address>:/home/ubuntu/
```

Alternative: use rsync - the options here preserve timestamps and permissions (-a), and only transfer if it would result in an updated version (-u). This can help with tracking what the latest file is; using scp as above will just overwrite the destination even if the file exists at the destination and is newer.
```bash
rsync -e 'ssh -i /path/to/key.pem' -avxu /path/to/test.txt ubuntu@<ip-address>:/home/ubuntu/
```

3.	Open the test.txt file using a text editor on your remote machine and edit it - add random text and save it.

4.	Copy the test.txt file back to your local machine by using the scp/pscp command on your local terminal.

```bash
scp -i path/to/key.pem ubuntu@<ip-address>:/home/ubuntu/test.txt path/to/folder/localmachine
```
If you are using windows command prompt run the following command 
```bash
pscp -P 22 -i path/to/key.ppk ubuntu@<ip-address>:/home/ubuntu/test.txt path/to/folder/localmachine 
```

Alternative: use rsync:
```bash
rsync -e ‘ssh -i /path/to/key.pem’ -avux ubuntu@<ip-address>:/home/ubuntu/test.txt /path/to/folder/localmachine
```

5.	Open the copied file to review your edits.

**Copy data from one instance to another instance**  
Now let us copy the sample test.txt file from your current instance (CI) to another instance (AI) you have earlier started. 

1.	Open a new terminal, and use it to ssh into the instance you started at the beginning of this workshop (AI); if it is unavailable, quickly spin up another instance. The specs do not matter.

-Make sure to copy your .pem key (for AI) to your instance (CI). This is the private key file that you are using to log in to the instances, i.e. the key pair you generated in step xx. Be sure to change the permissions as well (chmod 0400 keyname.pem). You can use scp to copy it to your new CI instance. (You will need the key file in order to communicate (SSH or SCP) with another instance with the same key pair)  
-SSH into your AWS instance (AI) from your local machine where you have the text.txt file.  
-Run the following command to copy the file from your current instance to the other (AI)  
```bash
scp -i path/to/key.pem path/to/test.txt ec2-user@<ip-address>:/home/ec2-user/ 
```
NOTE: Username to be used here are **ec2-user** or **ubuntu** if your instance was spun up from an Amazon Linux OS or Ubuntu Linux OS images, respectively.

2.	SSH into the instance you copied the file into from your current instance using the following command

```bash
ssh -i path/to/key.pem ec2-user@<ip-address>
```

3.	Review the file you copied to this instance. Type “exit” to simply return to your earlier instance.

---
B. You will now learn to share an S3 Bucket with your collaborator i.e.  
-Transferring data to your private bucket, but giving a collaborator access to one file in that bucket.  
-Giving a collaborator access to a public bucket with their files (and only them)  
-Transferring a file securely to your collaborator’s bucket.

**Transferring data to your private bucket, but giving a collaborator access to one file in that bucket**

We will now share files with your collaborators via an S3 bucket.  

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

**Next we will update the Access Control List to provide read and write permissions for objects**

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

11.	Now let us make a file object downloadable. Select the file you would like to share. Click on **Actions** and **Make public**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic6.png)

12.     Click on **Make public**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3BucketPermissionsBucketPublic7.png)

Only the files you made public should now be downloadable by running the AWS CLI commands.

```bash
aws s3 cp s3://{Collaborator-Bucket-Name}/{object name} ./
```

**Giving a collaborator access to a public bucket with their files (and only them)**

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
