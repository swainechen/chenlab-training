+++
title = "b. Attach EBS Volume to Instance"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 50
tags = ["tutorial", "aws console", "ec2"]
+++

1.	On the left navigation bar click on **Volumes** under **Elastic Block Store** to view all the Volumes. Select your newly created Volume.  

2.	Click on **Actions** and select **Attach Volume** from the drop-down list.

3.	In the **Attach Volume** dialog, click on the **Instance** field and select the EC2 Instance either by looking for the Instance ID or Name Tag. Click **Attach** to attach the Volume. In the list of displayed Volumes you should see your new Volume is **In Use**.

![Copy Paste](/images/hpc-aws-parallelcluster-workshop/EC2AttachVolume-2.png)

4.	On the left navigation pane click **EC2 Dashboard** at the very top. 

5.	Select **Running Instances**, then click your EC2 Instance with the newly attached volume.  

![Copy Paste](/images/hpc-aws-parallelcluster-workshop/EC2RunningInstances.png)  
 

6.	In the **Description** tab below scroll down to **Block Devices** to view the new block device added (e.g. /dev/sdh). **Important**: This device name could be different from what is shown below. Keep track of the device name and size; you'll need that information for mounting the volume later. 

![Copy Paste](/images/hpc-aws-parallelcluster-workshop/EC2AttachedVolume-2.png)

