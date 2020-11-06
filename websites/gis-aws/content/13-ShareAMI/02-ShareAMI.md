+++
title = "a. Share an AMI and EBS Volume snapshot"
weight = 60
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++

We will now share the AMI we earlier created with another team memeber/fellow Researcher. 
**Note**: Ensure you have the individuals or groups AWS Account Number for this purpose. 

1.	In the AWS Management Console search bar, type **EC2**.

2.	Choose and click on **EC2** to open the **EC2 dashboard** Console.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/aws-ec2.png)

3.	On the Left navigation pane click on **AMIs** under the **Images** section 

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2AMIs.png)

4.	From the list of AMIs select the AMI you created in the previous section of the workshop.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2AMI.png)

5.	Next click on **Actions** and choose and click on **Modify Image Permissions**

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2AMIPermissions.png)

6.	In the **Modify Image Permissions** dialog type in the AWS Account Number in the **AWS Account Number** field and click on **Add Permissions**. 

Observe that the default selection of the **This image is currently** attribute is set to **Private** indicating that you are sharing this with a specific AWS account and not publishing it the public domain.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2AMIPermissionsDialog1.png)

7.	click on the **Save** button.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2AMIPermissionsDialog2.png)

Your AMI has now been shared with a specific AWS Account you provided the AWS Account Number for. you could ask your team-memeber/Fellow Researcher to create an EC2 instance from a shared AMI, to have a look at those steps see the earlier sections in the workshop that desctibe those.


**It is also possible to share an EBS Volume snapshot in a similar manner**

8.	On the Left navigation pane of the **EC2 dashboard** click on **Snapshots** under the **Elastic Block Store** section.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2Snapshots.png)

9.	Next select the EBS Snapshot you would like to share and further click on **Actions** and choose and select **Modify Permissions**. Follow the same steps as sharing the AMI i.e. 6 and 7 above.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2SnapshotsPermissions.png)

You would now have also shared the snapshot with a  AWS Account you provided the AWS Account Number for. you could ask your team-memeber/Fellow Researcher to create a Volume out of the shared snapshot.


