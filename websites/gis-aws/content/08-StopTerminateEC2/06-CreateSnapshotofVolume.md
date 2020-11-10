+++
title = "c. Create snapshot of detached Volume"
weight = 60
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++

1.	Click on **Actions** and choose and select **Create Snapshot**. 

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2VolumeCreateSnapshot.png)


2.	In the **create snapshot** screen provide a **Description** for the Volume, this will be identifier to search for the snapshot of the volume.


![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2VolumeCreateSnapshotScreen.png)


3.	Click on **Add Tag** to provide a unique name tag for better identificaiton of the snapshot. Enter "[your initials]-training" in the **Key** field and a unique value in the **Value** fild.

**NOTE**: for the purpose of this training, use your previous tag and provide "[your initials]-EBS" as the value.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EBSVolumeSnapshot.png)

4.	Click on **Create Snapshot**


5.	A Snapshot creation task will be started and a message indicating the same will be displayed. Click on the **Close** button.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2VolumeCreateSnapshotSuccess.png)

6.	Click on  **Snapshots** under the **Elastic Block Store** section in the left navigation pane to view all the snapshots created.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2VolumeSnapshotCreated.png)
