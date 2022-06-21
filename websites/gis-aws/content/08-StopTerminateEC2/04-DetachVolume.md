+++
title = "b. Detach the additional EBS Volume"
weight = 60
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++

1.	In the AWS Management Console search bar, type **EC2**.

2.	In the left Navigation pane click on **Instances** under the **Instances** section.

3.	Select your EC2 Instance and in the Instance Details pane below, select **Storage** to view the **Block Devices** and then click on the **Volume ID** of your additional EBS Volume attached previously (e.g. /dev/sdh; the device name may be different in your case).

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2SelectVolume.png)

This will open the Volumes page showing the specific Volume details.

4.	Click on the **Actions** button and further click on **Detach Volume** to detach the volume from the EC2 Instance.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2VolumeDetach.png)

5.	In the **Detach Volume** dialog box, click **Detach**; the volume will now be detached and in a few moments the **Status** of the volume will change to **available**, indicating that it could be attached to another instance if required.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2VolumeDetachDialog.png)

