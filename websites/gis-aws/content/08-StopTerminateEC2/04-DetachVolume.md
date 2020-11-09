+++
title = "b. Detach the additional EBS Volume"
weight = 60
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++

1.	In the AWS Management Console search bar, type **EC2**.

2.	In the left navigation pane click on **Instances** under the **Instances** section.

3.	Select your EC2 Instance and In the Instance Details pane below select **Storage** to view the **Block Devices** and further click on **Volume ID** of the additional EBS Volume device we attached earlier (e.g. /dev/sdh **NOTE**: the device name may be different in your case).

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2SelectVolume.png)

This will open the Volumes page showing the specific Volumes details.

4.	Click on the **Actions** button and further click on **Detach Volume** to detach the volume from the EC2 instance.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2VolumeDetach.png)

5.	On the **Detach Volume** Dialog click **Yes, Detach**, the volume will now be detached and in a few moments the **Status** of the volume will change to **available** indicating that it could be attached to another instance if required.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2VolumeDetachDialog.png)





