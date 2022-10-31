+++
title = "b. S3 data transfer"
weight = 60
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++

Let us now upload and download data from S3, using:

1.	AWS Console   

2.	AWS CLI  

---
**Using AWS Console (via GUI)**

1.	In the AWS Management Console search bar, type **S3**.  

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3Service.png)

2.	Click on **S3** to open the **Amazon S3** Console.

3.	Scroll through the list of buckets and click on your new bucket. **NOTE**: make sure you click on the bucket name link, and not just select the bucket using the checkbox.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3SelectBucketLink.png)

4.	Click the **Upload** button.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3UploadClickUpload.png)

5.	Either drag and drop files or click  **Add Files** to browse and select the files to be uploaded. For the sake of this workshop, add a small text file or any other dummy file / image that does **NOT** contain sensitive information.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3UploadAddFiles.png)

6.	Click on **Upload**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3UploadAddFilesUpload.png)

The files will be uploaded to the S3 folder and the file object names will appear in the list of the objects.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3UploadAddFilesUploaded.png)


---
**Using AWS CLI (from the EC2 instance)** 
>_**IMPORTANT**: Make sure you have configured the AWS CLI credentials (Refer to [Section V.a]({{< ref "15-awscli-intro/02-introtoawscli.html" >}}) for AWS CLI configure commands)._

1.	Run the AWS CLI command to list all the S3 buckets associated with the account in the region configured with **aws configure**.

```bash
aws s3 ls
```

Output will be a list of buckets including the new bucket which was created above.

2.	Run AWS CLI commands.

**List the contents of the bucket**

```bash
aws s3 ls s3://{Bucket Name}
```

Verify that the output includes your uploaded files.


**Download files from the S3 bucket**

To download files we will use the AWS CLI **s3 cp** command.

```bash
aws s3 cp s3://{Bucket Name}/{object name} ./
```

Verify the download by listing the directory and confirm the object is present.

For more detailed steps on AWS CLI for S3 refer to Section VI ([**Download data shared using Amazon S3**]({{< ref "05-downloaddatafroms3.html" >}})) of the workshop.


**Upload files to the S3 bucket**

 To upload files we will use the AWS CLI **s3 cp** command.

```bash
aws s3 cp ./{your-file} s3://{Bucket Name}/
```

3.	Verify the contents of the S3 bucket after the copy is completed via command line (below) or by using the AWS Console.

```bash
aws s3 ls s3://{Bucket Name}
```
