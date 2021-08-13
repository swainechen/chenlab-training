+++
title = "a. Share an AMI and EBS Volume snapshot"
weight = 60
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++

We will now share the AMI we earlier created with another collaborator. 

**Note**: Ensure you have the individual's or group's AWS Account Number for this purpose. 

1.	In the AWS Management Console search bar, type **EC2**.

2.	Click on **EC2** to open the **EC2 dashboard** Console.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/aws-ec2.png)

3.	On the Left navigation pane click on **AMIs** under the **Images** section.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2AMIs.png)

4.	From the list of AMIs, select the AMI you created in the previous section of the workshop.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2AMI.png)

5.	Next click on **Actions** and click on **Modify Image Permissions**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2AMIPermissions.png)

6.	In the **Modify Image Permissions** dialog type in the AWS Account Number in the **AWS Account Number** field and click on **Add Permissions**. Be absolutely sure to clarify the correct AWS Account Number with your collaborator.  

The AWS Account Number can be found here

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2AMIAWSAccNum.png)

Observe that by default, the image is set to Private. Keep that option checked, and enter the AWS Account Number in the space provided below. This will ensure that your image is still private (not available to public access), but the specific account mentioned will have access.  

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2AMIPermissionsDialog1.png)

7.	Click on the **Save** button.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2AMIPermissionsDialog2.png)

Your AMI has now been shared with this specific AWS Account Number provided. Invite your collaborator to click on **Launch Instance** on the **EC2 Dashboard** and search for the AMI you recently shared under **My AMIs** > **Shared with me**.  

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2AMIshared.png)

Your collaborator will now be able to create an EC2 instance from your shared AMI (Refer to section X about setting up a linux machine from a previously created AMI). 


**It is also possible to share an EBS Volume snapshot in a similar manner**.

8.	On the left navigation pane of the **EC2 dashboard**, under the **Elastic Block Store** section, click on **Snapshots**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2Snapshots.png)

9.	Next select the EBS Snapshot you would like to share and click **Actions**. Select **Modify Permissions** from the drop-down list. Follow the same steps on sharing the AMI from steps 6 and 7 above.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2SnapshotsPermissions.png)

Now, you have successfully shared the snapshot with the provided AWS Account. Invite your collaborator to view the shared snapshot by first clicking on **Snapshots** under the **Elastic Block Store** on the left navigation pane and then on the drop-down list that says **Owned By Me** and select **Private Snapshots**.  

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2Snapshotshared.png)

You can ask your collaborator to create a Volume from the shared snapshot.
