+++
title = "c. Delete EBS Volumes"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 50
tags = ["tutorial", "aws console", "ec2"]
+++

**Deleting EBS Volumes**

Time to send our extra storage to the bin.

1.  Using the EC2 dashboard, select **Elastic Block Store -> Volumes**. Select your volumes to be deleted. Before proceeding, make sure that the **Status** column shows "Available"; if the status shows "In Use", go to **Actions -> Detach** and wait for the volumes to change status.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2VolumeWindow.png)

2.  Go to **Actions -> Delete Volumes**. 

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2VolumeActionDelete.png)

3.  In the prompt window, click **Yes, Delete**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2VolumeDeleteWindow.png)

These volume status should now show "Deleting". They will be gone shortly. 
