+++
title = "a. Create EBS Volume from snapshot"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 40
tags = ["tutorial", "Prerequisite", "ec2"]
+++

1.	In the AWS Management Console search bar, type **EC2** to navigate to the EC2 Service.


2.	On the Left Navigation Bar click on **Snapshots** under **Elastic Block Store**

![Copy Paste](/images/hpc-aws-parallelcluster-workshop/EC2Snapshots.png)


3.	From the list of snapshots select the snapshot you created in the previous section of the workhop, **NOTE**: you can identify the snapshot by the Snapshot Name Tag you had provided earlier.

![Copy Paste](/images/hpc-aws-parallelcluster-workshop/EC2SelectSnapshot.png)

4.	Next click on **Actions** button and choose and click on **Create Volume**

![Copy Paste](/images/hpc-aws-parallelcluster-workshop/EC2CreateVolume.png)

4.	In the **Create Volume** notice that the **Size** field is already populated but you could increase the size (however, you cannot decrease the size). **Important**:Make sure the **Availability Zone** is the same as the EC2 instance you are going to attach the volume to.

5.	Click on **Add Tag** to uniquely tag the resource. Provide a value "Name" for the **Key** and a unique name for the **Value**.
**NOTE**:For the purposes of the lab please use your [Name+Lastname-EBS] as the **Value**  

![Copy Paste](/images/hpc-aws-parallelcluster-workshop/EC2VolumeTag.png)

6.	Click on **Create Volume** to Create the volume. If the volume is created you would get the following screen:

![Copy Paste](/images/hpc-aws-parallelcluster-workshop/EC2VolumeCreated.png)
