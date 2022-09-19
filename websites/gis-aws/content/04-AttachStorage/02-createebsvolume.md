+++
title = "a. Create an EBS Volume"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 40
tags = ["tutorial", "createebsvolume", "ec2"]
+++

1.	In the AWS Management Console search bar, type EC2 to navigate to the EC2 Service.


2.	On the Left Navigation Bar click on **Volumes** under **Elastic Block Store**

![Copy Paste](/images/hpc-aws-parallelcluster-workshop/EC2Volumes.png)


3.	Click on **Create Volume** (on the top right corner) to create a new volume

4.	In the **Create Volume** page fill out the required size of the volume in GB, e.g. 10 GB or more depending on the data set size.

**NOTE**: For the purposes of this lab 10 GB will suffice.   

**Important**: Make sure the **Availability Zone** is the same as the EC2 instance you are going to attach the volume to.   

![Copy Paste](/images/hpc-aws-parallelcluster-workshop/CreateVolumeAZ.png)

![Copy Paste](/images/hpc-aws-parallelcluster-workshop/CreateVolume.png)


5.	Click on **Add Tag** to uniquely tag the resource. Provide "Name" for the **Key** and "[your initials]-EBS" for the value.

![Copy Paste](/images/hpc-aws-parallelcluster-workshop/EC2VolumeTag.png)

6.	Click on **Create Volume** to Create the volume. 
