+++
title = "a. Share an EBS snapshot"
weight = 60
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++

We will now EBS Volume snapshot in a similar manner.

8.	On the left navigation pane of the **EC2 dashboard**, under the **Elastic Block Store** section, click on **Snapshots**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2Snapshots.png)

9.	Next select the EBS Snapshot you would like to share and click **Actions**. Select **Modify Permissions** from the drop-down list. Follow the same steps on sharing the AMI from steps 6 and 7 above.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2SnapshotsPermissions.png)

Now, you have successfully shared the snapshot with the provided AWS Account. Invite your collaborator to view the shared snapshot by first clicking on **Snapshots** under the **Elastic Block Store** on the left navigation pane and then on the drop-down list that says **Owned By Me** and select **Private Snapshots**.  

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2Snapshotshared.png)

You can ask your collaborator to now create a Volume from the shared snapshot.
