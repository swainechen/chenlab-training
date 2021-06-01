+++
title = "b. Attach EBS Volume to Instance"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 50
tags = ["tutorial", "aws console", "ec2"]
+++

1.	On the Left Navigation Bar click on **Volumes** under **Elastic Block Store** to view all the volumes. Select your newly created volume searching the list for the unique Name tag provided in the previous step.

2.	Click on **Actions** and further click on **Attach Volume**

3.	In the **Attach Volume** dialog, click on the **Instance** field and select the EC2 instance by looking for the Instance ID or Name Tag in the list. Click **Attach** to attach the volume.

![Copy Paste](/images/hpc-aws-parallelcluster-workshop/EC2AttachVolume.png)

If attached successfully - In the list of displayed volumes you should see your new volume having a status (under the State Column) indicating **In Use**.

4.	On the left Navigation pane click **EC2 Dashboard** on the very top of the list of items. 

5.	Select **Running instances** and then select your EC2 Instance from the list of instances to which you attached the Volume.

![Copy Paste](/images/hpc-aws-parallelcluster-workshop/EC2RunningInstances.png)
 

6.	In the EC2 Details section select **Storage** to see the **Block devices** attached to the instance. Notice the new volume just created attached as an additional Block Device with a Volume ID, Device Name and Volume Size.

![Copy Paste](/images/hpc-aws-parallelcluster-workshop/Ec2AttachedVolume.png)

**Important**: Note down the specific device name for the next step of mounting the volume. The drive name may differ from what's shown. **NOTE**: Newer Linux kernels may rename your devices to **/dev/xvdf through /dev/xvdp** internally, even when the device name entered here (and shown in the details) is /dev/sdf through /dev/sdp.
