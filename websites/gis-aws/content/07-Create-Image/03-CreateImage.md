+++
title = "b. Create EC2 Image"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 50
tags = ["tutorial", "aws console", "ec2"]
+++

We will now create an image of the EC2 instance.

1.	Select the EC2 instance we stopped, click on the **Actions** button, click on **Image**, and click **Create Image**.


![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2CreateImage.png)


2.	On the **Create Image** page, provide an **Image Name** and **Image description**. Ensure you remove the extra EBS volume we had attached to the EC2 instance by clicking on the "x" button (see image below) . Click on **Create Image** to start the Image creation task.  


![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2CreateImageDialog.png)

An image creation task is started in the backgroud and you will see a message stating "Successfully created ..." on the top of the EC2 Dashboard page when it's ready to go.

**NOTE**: In case you had not stopped the instance and were still logged in, the instance would have rebooted and your ssh session would have ended. You will need to connect again (verify the IP address).

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2CreateImageSuccess.png)

3.	Pause for a minute to let the Image creation task progress. To check if the image has been created click on **AMIs** under the **Images** section in the left Navigation pane to view all the AMI's created earlier as well as the new one in progress.

You will see that the new AMI is either being created or already created and ready for use with the name you've given.


![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2AMI.png)

4.	Select your new AMI instance, click on the **Actions** button, then **Manage Tags**. Make sure to tag your image with a Name, User, and Grant tag, as before. For this instance, you can use the Name value "[your-initials]-AMI".  

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2AMI2.png)

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/EC2AMIAddTags.png)

