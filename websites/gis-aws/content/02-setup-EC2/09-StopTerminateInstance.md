+++
title = "g. Stop and Terminate instance"
weight = 90
tags = ["tutorial", "cloud9", "s3"]
+++

Now let us stop and terminate the instance

1.	Login to the AWS console, in the AWS Managament Console search bar, type "EC2".

2.	Choose **EC2** to open the **EC2 Dashboard**.

3.	Click on **Instances** under the **Instances** section in the left Navigation pane to view all the instances that are Running, Stopped and Terminated.

4.	Select the instance created and click on **Actions** button and further choose and click on **Instance State** and futher choose and click on **Stop instance**


![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2StopInstance.png)

5.	Click on **Stop** to confirm stopping the instance.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2StopConfirm.png)

The instance will now attempt to stop, wait for a few moments for the instance to stop, until **Instance State** column shows **Stopped**

We will now terminate the EC2 instance.

6.	Select the EC2 instance from the list of instances. Click on **Actions** and choose and click on **Instance State** and further choose and click on **Terminate instance**

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2InstanceTerminate.png)


7.	On the Terminate Instance dialog, observe the note **Associated resources may incur costs after these instances are terminated** and further clickon the **Delete EBS volume** to read the additional note. To proceed  select **Terminate**.


![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2InstanceTerminateConfirm.png)

The instance will now be terminated.

8.	After a minute, view the **Instance State** of the instance which will show **Terminated**

