+++
title = "a. Create EBS Volume from snapshot"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 40
tags = ["tutorial", "Prerequisite", "ec2"]
+++

1.	In the AWS Management Console search bar, type **EC2** to navigate to the EC2 Service.


2.	On the Left Navigation Bar click on **Snapshots** under **Elastic Block Store**

![Copy Paste](/images/hpc-aws-parallelcluster-workshop/EC2Snapshots.png)


3.	Select the Snapshot you created in the section IX of the workshop. **NOTE**: You can identify the Snapshot by the Snapshot Name Tag you had provided earlier.

![Copy Paste](/images/hpc-aws-parallelcluster-workshop/EC2SelectSnapshot.png)

4.	Next click on **Actions** and choose **Create Volume** from the drop-down list.

![Copy Paste](/images/hpc-aws-parallelcluster-workshop/EC2CreateVolume.png)

5.	After selecting **Create Volume**, make note that the **Size** field is prepopulated but you will have the option to increase the size (however, you cannot decrease the size). **Important**: Make sure that the **Availability Zone** is the same as the EC2 Instance you are going to attach the Volume to.

6.	Click on **Add Tag** to uniquely tag the resource. Provide "Name" for the **Key** and "[your initials]-EBS" for the **Value**.

![Copy Paste](/images/hpc-aws-parallelcluster-workshop/EC2VolumeTag.png)

7.	**Create Volume** will create a new EBS Volume. If the Volume is created you will see the following screen:

![Copy Paste](/images/hpc-aws-parallelcluster-workshop/EC2VolumeCreated.png)
