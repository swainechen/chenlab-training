+++
title = "a. Terminate your EC2 Instance"
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

