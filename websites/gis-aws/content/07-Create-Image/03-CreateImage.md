+++
title = "b. Create EC2 Image"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 50
tags = ["tutorial", "aws console", "ec2"]
+++

We will now create an image of the EC2 Instance.

1.	Select the EC2 instance we just stopped,click on **Actions** button, next click on **Image** and choose and click **Create Image** 


![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2CreateImage.png)


2.	On the **Create Image** page provide an **Image Name** and **Image description**. Ensure you remove the extra EBS volume we had attached to the EC2 Instance by clicking on the "x" button (see image below) . Click on **Create Image** to start the Image creation task.


![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2CreateImageDialog.png)

An image creation task is started in the backgroud and you would see a message indicating "Successfully created ..." on the top of the EC2 Dashboard page.

**NOTE**:Incase you had not stopped the instance and were still logged into the instance, the Instance would have rebooted and your ssh session would need to be restarted to connect back to the instance.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2CreateImageSuccess.png)

3.	Pause for a minute to let the image creation taks to progress. To check if the image has been created choose and click on **AMIs** under the **Images** section in the left Navigation pane to view all the AMI's created earlier as well as the new one being created.

You will see that the new AMI is either being created or already created and ready for use with the name supplied in the **Create Image** page.


![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2AMI.png)
