+++
title = "e. Terminate EC2 Instance"
weight = 60
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++

We will now terminate the EC2 instance.

**NOTE**:that we have already created an Image of the instance in the previous section of the workshop.

1.	Select the EC2 instance from the list of instances. Click on **Actions** and choose and click on **Instance State** and further choose and click on **Terminate**

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2InstanceTerminate.png)


2.	On the Terminate Instance dialog, observe the note **Associated resources may incur costs after these instances are terminsated** and further clickon the **Delete EBS volume** to read the additional note.To porceed  select **Yes, terminate**.


![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2InstanceTerminateConfirm.png)

The instance will now be terminated.

3.	After a minute view the **Instance State** of the instance, it will indicate as **Terminated**
