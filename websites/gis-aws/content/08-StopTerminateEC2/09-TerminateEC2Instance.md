+++
title = "e. Terminate EC2 Instance"
weight = 60
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++

We will now Terminate the EC2 Instance.

**NOTE**: Ensure your previous EC2 Image was successful; we'll be using it in the next portion of the workshop. In this segment, you will be terminating your EC2 Instance. This is permanent and irreversible. However, using this EC2 Image, it is possible to spin up an identical machine. In general, be very cautious when terminating machines.

1.	Select the EC2 Instance from the list of instances. Click on **Instance State** > **Terminate Instance**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2InstanceTerminate-2.png)


2.	On the Terminate Instance dialog, observe the note "Associated resources may incur costs after these instances are terminated" and click  **Delete EBS Volume**. Read the additional note. To proceed, select **Yes, terminate**. This is an important safeguard, and never click this automatically!  


![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2InstanceTerminateConfirm-2.png)

The instance will now be terminated.

3.	After a minute view the **Instance State** of the instance, it will indicate **Terminated**.

  
