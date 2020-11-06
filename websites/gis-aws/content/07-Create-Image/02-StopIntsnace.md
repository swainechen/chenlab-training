+++
title = "a. Stop EC2 Instance"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 50
tags = ["tutorial", "aws console", "ec2"]
+++

We will now stop the running instance. Before we proceed it is a good idea to gracefully stop any running applications and ensure that data is saved. Stopping an instance is the equivalent of turning off your computer/server.

**NOTE**: It is best practices but not mandatory to stop the instance before creating the image. However while creating the image the instance will udnergo a reboot, hence make sure there are no applications running.

1.	Before stopping let us copy all the required data from /tmp to a more persistant storage as stopping the instance will clear out /tmp. Copy all the required files to the EBS volume we mounted at **/mnt/volume1** before proceeding.

```bash
cp -ax /tmp/fastq /mnt/volume1/
```

```bash
cp -ax /tmp/outbreaks /mnt/volume1/
```

2.	Login to the AWS console, in the AWS Management Console search bar, type **EC2**.

3.	Choose **EC2** to open the **EC2 Dashboard**.

4.	In the left Navigation pane, Click on **Instances** under the **Instances** section to view all the Intsnces that are Running, Stopped and Terminated.

5.	Select your Instance. 

**NOTE**: The **Instance State** of the EC2 instance shows **running** state.

6.	Click on **Actions** button and further choose and click on  **Instance State** and further choose and click on **Stop**

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2StopInstance.png)

7.	On the Dialog box click on **Stop** to confirm stopping the instance.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2StopConfirm.png)

8.	Note that the instance will now attemped to be stopped, and after a short duration the **Instance State** will change to **Stopped**
