+++
title = "b. Attach EBS volume to instance"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 50
tags = ["tutorial", "aws console", "ec2"]
+++

1.	On the left navigation bar click on **Volumes** under **Elastic Block Store** to view all the volumes. Locate and select your newly created volume searching the list for the unique Name Tag provided in the previous step.

2.	Click on **Actions** and select **Attach Volume** from the drop-down list.

3.	In the **Attach Volume** dialog, click on the **Instance** field and select the EC2 instance either by looking for the Instance ID or Name Tag. Click **Attach** to attach the volume. In the list of displayed volumes you should see your new volume displaying an  **In Use** status.

![Copy Paste](/images/hpc-aws-parallelcluster-workshop/EC2AttachVolume.png)

4.	On the left navigation pane click on  **EC2 Dashboard** on the very top of the list of items. 

5.	Select **Running instances** and further select the EC2 instance from the list of EC2 instances to which you attached the volume.

![Copy Paste](/images/hpc-aws-parallelcluster-workshop/EC2RunningInstances.png)
 

6.	In the **Description** tab below scroll down to section **Block Devices** to view the new block device added (e.g. /dev/sdh). **Important**: This device name could be different from what is shown below. Keep track of the device name and size; you'll need that information for mounting the volume later. 

![Copy Paste](/images/hpc-aws-parallelcluster-workshop/Ec2AttachedVolume.png)

