+++
title = "d. Deregister your AMI"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 50
tags = ["tutorial", "aws console", "ec2"]
+++

**Deregistering AMIs**  

Next, we'll deregister our AMI. As with all previous removals, this is *permanent*. Be careful, and make sure you're deregistering the correct AMI. 

1.  Using the EC2 dashboard, select **Images -> AMI**. Select your AMI to be de-registered. 

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2FindAMI.png)

2.  Go to **Actions -> Deregister**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2FindAMIActions.png)

3.  Select **Deregister**. In the prompt window, select **Continue**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2DeregisterAMI.png)

4.  Your AMI should be gone. 
