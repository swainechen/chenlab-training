+++
title = "a. Delete all of your extra resources"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 50
tags = ["tutorial", "aws console", "ec2"]
+++

**Terminating EC2 Instances**
We've covered this before, but just to make sure:

1.  Select your EC2 instance from the available list. Click **Actions**, **Instance State**, and then **Terminate**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2InstanceTerminate-2.png)

2.  When prompted for verification, select **Terminate**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2InstanceTerminateConfirm-2.png)

3.  The instance state will now show **Terminated**.

It'll still be visible for a bit, but then the entry will be deleted.


**Deleting S3 buckets**
If you delete a bucket, it and its contents are gone forever. If you enabled versioning, those versions are gone, too. As with termination, be *absolutely sure* you're ready to delete these resources.

1.  Open the S3 dashboard and select your bucket. In the upper right-hand corner, click **Delete**.

![S3 Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3DeleteButton.png)

2.  In the dialogue box, type the name of your bucket in the space provided. Note the warnings on the page- when the bucket is deleted, it's gone forever, and other people have the opportunity to use the bucket name.

![S3 Management Console](/images/hpc-aws-parallelcluster-workshop/S3/S3delete.png)

3.  Your bucket will no longer show on the S3 menu. Success!

**Deleting EBS Volumes**
Time to send our extra storage to the bin.

1.  Using the EC2 dashboard, select **Elastic Block Store -> Volumes**. Select your volumes to be deleted. Before proceeding, make sure that the **Status** column shows "Available"; if the status shows "In Use", go to **Actions -> Detach** and wait for the volumes to change status.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2VolumeWindow.png)

2.  Go to **Actions -> Delete Volumes**. 

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2VolumeActionDelete.png)

3.  In the prompt window, click **Yes, Delete**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2VolumeDeleteWindow.png)

These volume status should now show "Deleting". They will be gone shortly. 


**Deregistering AMIs**  
Next, we'll deregister our AMI. As with all previous removals, this is *permanent*. Be careful, and make sure you're deregistering the correct AMI. 

1.  Using the EC2 dashboard, select **Images -> AMI**. Select your AMI to be de-registered. 

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2FindAMI.png)

2.  Go to **Actions -> Deregister**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2FindAMIActions.png)

3.  Select **Deregister**. In the prompt window, select **Continue**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2DeregisterAMI.png)

4.  Your AMI should be gone. 

**Deleting Snapshots**
1.  Using the EC2 dashboard, select **Elastic Block Store -> Snapshots**. Select your snapshot to be deleted.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2FindSnapshot.png)

2.  Go to **Actions -> Delete** 

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2SnapshotActions.png)

3.  Select **Delete**. In the prompt window, select **Yes, Delete**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2SnapshotDelete.png)

4. Your Snapshot should no longer exist! 


**Conclusion**
At this point, all of your resources should be gone. You can revisit the Billing and Costs segment and, using your tags, verify that all chargeable resources have been removed. 

We hope this AWS workshop was constructive. This website will remain up and will be improved over time as a resource and a reference for your use in the future. Please contact your instructors with any questions! We're always happy to help.
