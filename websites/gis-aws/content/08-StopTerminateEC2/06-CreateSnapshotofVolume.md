+++
title = "c. Create snapshot of detached Volume"
weight = 60
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++

1.	Click on **Actions**, then select **Create Snapshot**. 

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2VolumeCreateSnapshot.png)


2.	In the **Create Snapshot** screen, provide a **Description** for the Volume; this will be an identifier to search for the Snapshot of the Volume.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2VolumeCreateSnapshotScreen.png)


3.	Click on **Add Tag** to provide a unique name tag for better identificaiton of the snapshot. Enter "Name" in the **Key** field, and "[your initials]-snapshot" for the **Value** field. Add the "User" and "Grant" tag you've previously used as well.


![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EBSVolumeSnapshot.png)

4.	Click on **Create Snapshot**.


5.	A Snapshot creation task will be started, and a message indicating the same will be displayed in a green bar on the top of the screen. 

6.	Click on **Snapshots** under the **Elastic Block Store** section in the left Navigation pane to view all the snapshots created.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2VolumeSnapshotCreated.png)
