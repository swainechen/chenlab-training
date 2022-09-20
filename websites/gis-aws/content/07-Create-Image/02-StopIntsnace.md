+++
title = "a. Stop EC2 Instance"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 50
tags = ["tutorial", "aws console", "ec2"]
+++

We will now stop the running instance. Before we proceed, it is a good idea to stop any running applications and ensure that data is saved. Stopping an instance is equivalent to turning off your computer/server.

**NOTE**: It is best practices, but not mandatory, to stop the instance before creating the image. However while creating the image the instance will undergo a reboot, so make sure there are no applications running.

1.	Before stopping, copy all the required data from /tmp to a more persistent storage as stopping the instance will clear out /tmp. Copy all the required files to the EBS volume we mounted previously, at **/mnt/volume1**, before proceeding.

Make sure to set the permissions to the directory prior to copying.  

```bash
sudo chmod -R 777 /mnt/volume1/
```

```bash
cp -ax /tmp/fastq /tmp/outbreaks /tmp/RNASeq /tmp/SingleCell /mnt/volume1/
```
2.	Login to the AWS Console, and in the AWS Management Console search bar, type **EC2**.

3.	Choose **EC2** to open the **EC2 Dashboard**.

4.	In the left Navigation pane, Click on **Instances** to view all the instances that are Running, Stopped and Terminated.

5.	Select your instance. 

**NOTE**: The **Instance State** of the EC2 instance shows the machine is **Running**.

6.	Click on the **Actions** button and then click on **Instance State**, and finally click on **Stop**

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2StopInstance-2.png)

7.	On the dialog box click on **Stop** to confirm stopping the instance.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2StopConfirm-2.png)

8.	Note that the instance will now attempt to stop, and after a short duration the **Instance State** will change to **Stopped**. You may have to refresh the window.
