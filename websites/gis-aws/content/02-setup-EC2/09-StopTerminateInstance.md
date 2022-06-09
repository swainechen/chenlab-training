+++
title = "h. Stop and Terminate instance"
weight = 90
tags = ["tutorial", "cloud9", "s3"]
+++

Now let us stop and terminate the instance - **before following these steps, listen to your instructor**. For educational purposes, they may ask you to review the steps, but not actually terminate the instance.

1.	Return to the AWS console. In the AWS Management Console search bar, type "EC2", or click the **EC2 Dashboard** item in the left-hand menu.

2.	Choose **EC2** to open the **EC2 Dashboard**.

3.	Click on **Instances** under the **Instances** section in the left navigation pane to view all the instances that are Running, Stopped and Terminated.

4.	Select the instance created and click on **Instance State**, then **Stop instance**.


![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2StopInstance.png)

5.	Click on **Stop** to confirm stopping the instance.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2StopConfirm.png)

The instance will now attempt to stop, wait for a few moments for the instance to stop, until **Instance State** column shows **Stopped**.

We will now terminate the EC2 instance.

6.	Select the EC2 instance from the list of instances. Click on **Instance State**, and then click on **Terminate instance**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2InstanceTerminate.png)


7.	In the Terminate Instance dialog, observe the note **On an EBS-backed instance, the default action is for the root EBS volume to be deleted when the instance in terminated. Storage on any local drives will be lost**. To proceed, select **Terminate**.


![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2InstanceTerminateConfirm.png)

The instance will now be terminated.

8.	After a minute, the **Instance State** of the instance which will now show **Terminated**.

