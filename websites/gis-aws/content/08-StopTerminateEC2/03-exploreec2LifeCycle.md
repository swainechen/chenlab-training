+++
title = "a. Exploring EC2 Instance Lifecycle"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 50
tags = ["tutorial", "aws console", "ec2"]
+++

An Amazon EC2 instance transitions through different states from the moment you launch it through to its termination.

To view the current instance state of an EC2 Instance:

1.	In the AWS Management Console search bar, type **EC2**.

2.	Choose EC2 to open the **EC2 Dashboard**.

3.	In the left navigation pane select **Intances** under the **Instances** section.

4.	To see the current instance state select your earlier created  EC2 Instance and observe the value under the **Instance state** column. It should be either indicating - Running, Stopped or Terminated.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2InstanceState.png)


Here is the state diagram of an EC2 instance. Also notice that the state affects the costs while using EC2. For more information see [Instance Lifecycle](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-instance-lifecycle.html)


![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2Lifecycle.png)




