+++
title = "a. Exploring EC2 Instance Lifecycle"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 50
tags = ["tutorial", "aws console", "ec2"]
+++

An Amazon EC2 instance transitions through different states throughout its lifecycle.

To view the current instance state of an EC2 Instance:

1.	In the AWS Management Console search bar, type **EC2**.

2.	Choose EC2 to open the **EC2 Dashboard**.

3.	In the left navigation pane select **Instances** under the **Instances** section.

4.	To see the current instance state, select your EC2 Instance and observe the value under the **Instance state** column. The status will indicate Running, Stopped or Terminated.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2InstanceState.png)


Below is a diagram of an EC2 instance lifecycle. Note that the instance state affects the costs while using EC2. For more information see [Instance Lifecycle](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-instance-lifecycle.html)


![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2Lifecycle.png)  
  
Very briefly, an EC2 instance will accumulate cost while it is running. When you stop the EC2 instance, you no longer accumulate costs associated with it.  
  
However, you will always be charged for the storage on your machine. Whether it's the storage on the root drive or EBS volumes attached to your root, these costs do not stop when your machine is stopped. They will only cease when you permanently delete these volumes. Furthermore, you pay for the entire storage present on the disk, not just what you use. Make sure to select resources appropriate for your project. 


