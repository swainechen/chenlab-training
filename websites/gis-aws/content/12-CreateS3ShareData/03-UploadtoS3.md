+++
title = "b. Upload data to S3"
weight = 60
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++

**We will first use the Console to upload files:**

1.	In the AWS Management Console search bar, type **S3**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3Service.png)

2.	Choose and click on **S3** to open the **Amazon S3** Console.

3.	Scroll through the list of buckets and choose and click on the bucket name created. **NOTE**: make sure you click on the bucket name link and not just select the bucket using the checkbox.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3SelectBucketLink.png)

4.	Click on the **Upload** button.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3UploadClickUpload.png)

5.	Either drag and drop files or click on **Add Files** to browse and select the files to be uploaded.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3UploadAddFiles.png)

6.	Click on **Upload**

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3UploadAddFilesUpload.png)

The files will be uploaded to the S3 folder and the file object names would appear in the list box of the objects.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3UploadAddFilesUploaded.png)


**We will now use AWS CLI to list objects, upload and download files from an EC2 instance:**

**NOTE**: Make sure you have configured the AWS CLI credentials (Refer to section V.a for AWS CLI configure commands).


7.	Run the AWS CLI command to list all the S3 buckets associated with the account in the region configured with **aws configure** commands

```bash
aws s3 ls
```

Output will be a list of buckets including the new bucket which was created above.

8.	Run AWS CLI command to list the contents of the bucket

```bash
aws s3 ls s3://{Bucket Name}
```

Verify that the output includes files uploaded in the previous step via the console.

**We will now download files from the S3 bucket:**

To download files we will use the AWS CLI **s3 cp** command.

```bash
aws s3 cp s3://{Bucket Name}/{object name} ./
```

Verify the download by listing the directory and confirm the object present in the list-box of objects for the bucket.

For more detailed steps on AWS CLI for S3 refer to section VI (**Download data shared using Amazon S3**) of the workshop.

**We will now upload files to the S3 bucket**

 To upload files we will use the AWS CLI **s3 cp** command.

```bash
aws s3 cp ./{your-file} s3://{Bucket Name}/
```

9.	After the copy is completed. Verify the contents of the S3 bucket.

```bash
aws s3 ls s3://{Bucket Name}
```

Verify the output includes the uploaded file. You could optionally use the AWS Console to view the bucket contents.

