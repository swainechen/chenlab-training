+++
title = "a. Share an AMI and EBS Volume snapshot"
weight = 60
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++

We will now share the AMI we earlier created with another collaborator. 
**Note**: Ensure you have the individual's or group's AWS Account Number for this purpose. 

1.	In the AWS Management Console search bar, type **EC2**.

2.	Choose and click on **EC2** to open the **EC2 dashboard** Console.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/aws-ec2.png)

3.	On the Left navigation pane click on **AMIs** under the **Images** section 

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2AMIs.png)

4.	From the list of AMIs select the AMI you created in the previous section of the workshop.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2AMI.png)

5.	Next click on **Actions** and choose and click on **Modify Image Permissions**

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2AMIPermissions.png)

6.	In the **Modify Image Permissions** dialog type in the AWS Account Number in the **AWS Account Number** field and click on **Add Permissions**. Make sure to clarify the correct AWS Account Number with your collaborator.

Observe that the default selection of the **This image is currently** attribute is set to **Private** indicating that you are sharing this with a specific AWS account and not publishing it the public domain.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2AMIPermissionsDialog1.png)

7.	Click on the **Save** button.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2AMIPermissionsDialog2.png)

Your AMI has now been shared with a specific AWS Account you provided the AWS Account Number for. You could now invite your collaborator to create an EC2 instance from your shared AMI (Refer to section X about setting up a linux machine from a previously created AMI). 


**It is also possible to share an EBS Volume snapshot in a similar manner**

8.	On the left navigation pane of the **EC2 dashboard**, under the **Elastic Block Store** section, click on **Snapshots**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2Snapshots.png)

9.	Next select the EBS Snapshot you would like to share and click on **Actions**. Select **Modify Permissions** from the drop-down list. Follow the same steps on sharing the AMI from steps 6 and 7 above.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2SnapshotsPermissions.png)

You would now have also shared the snapshot with an AWS Account you provided the AWS Account Number for. You could ask your collaborator to create a Volume from the shared snapshot.


