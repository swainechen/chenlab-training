+++
title = "d. Delete the detached Volume"
weight = 60
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++


Now that a Snapshot has been taken of the volume, we will delete the volume itself.

1.	In the left Navigation pane click on **Volumes** under the **Elastic Block Store** section to view all the volumes.	 

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2VolumeDelete.png)

2.	Click on **Actions** and select **Delete Volume**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2VolumeDeleteDialogConfirm.png)

3.	On the **Delete Volume** dialog click on **Delete** to confirm the delete. The EBS Volume will now be deleted. This volume no longer exists on AWS, and you will no longer be charged for it. 

