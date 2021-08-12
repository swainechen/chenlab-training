+++
title = "e. Launch an EC2 Instance"
weight = 80
tags = ["tutorial", "cloud9", "aws cli", "s3"]
+++

We will now launch an EC2 Linux based instance.


1.	Return to the AWS Management Console and open the Amazon EC2 Dashboard.

2.	In the upper-right corner of the AWS Management Console, confirm you are in the desired AWS region **(e.g., Singapore)**.

![EC2LaunchInstance](/images/hpc-aws-parallelcluster-workshop/ConsoleRegion.png)


3.	Click **Launch instance**, then click **Launch instance** again from the drop down menu.

![EC2LaunchInstance](/images/hpc-aws-parallelcluster-workshop/EC2LaunchInstance.png)

4.	Next in the **Step 1** page, select the first **Amazon Linux 2 AMI(HVM), SSD Volume Type** and check **64-bit (x86)** under the **Select** button and click on the **Select**. 

**NOTE**: The ami-xxxxxxxxx label and specific versions of the installed package may be different than in the image below.

![EC2SelectAMI](/images/hpc-aws-parallelcluster-workshop/EC2SelectAMI.png)

5.	In the **Step 2** page, choose an Instance type, select the **c5.2xlarge** instance from the list and click **Next: Configure Instance Details**.  
**Note:** This is a large instance - we're doing this on purpose. This will not be a default instance you should use until you've considered the resources you need for your project.

![EC2ChooseInstancetype](/images/hpc-aws-parallelcluster-workshop/EC2ChooseInstanceType-1.png)

6.	On **Step 3** page, Configure Instance Details - leave the default settings. Note that the **Subnet** field can be configured to launch the instance in a specific **Availability Zone**; while we are keeping the default for this workshop, this gives you control over the location of your machine. Click the **Next: Add Storage** button in the bottom right corner. [Learn More](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/configuring-instance-metadata-service.html)

![EC2ChooseInstancetype](/images/hpc-aws-parallelcluster-workshop/EC2StepConfigureInstance.png)


7.	On the **Step 4** page, you have the ability to modify or add storage and disk drives to the instance. For this lab, we will simply accept the storage defaults. Take note that the default setting for **Delete on Termination** is affirmative. This indicates that if the machine is terminated, the root volume associated with the instance will be deleted. You need to uncheck this if you plan to store data on the root volume which you would want to access later. Next, click **Next: Add Tags**.


![EC2ChooseInstancetype](/images/hpc-aws-parallelcluster-workshop/EC2AddStorage.png)


8.	On **Step 5** page, you can choose a friendly name for your instance by clicking "Add Tag", and entering "Name" for the key, and "my-first-EC2" for the value. This Name key, more correctly known as a tag, will appear in the console once the instance launches. It makes it easy to keep track of running machines in a complex environment. Create additional tags similar for the ones you used for the key pair previously - give a user and grant key to this machine, and enter the same values. When ready, click **Next: Configure Security Group**.

![EC2AddTags](/images/hpc-aws-parallelcluster-workshop/EC2AddTags-3.png)

9.	On **Step 6** page, you will be prompted to create a new security group, which will be your firewall rules. Provide a name for your new security group.

**NOTE**: For the purpose of this lab use your "[Your initials]-sec-group".

![EC2ConfigSecGroups](/images/hpc-aws-parallelcluster-workshop/EC2ConfigSecGroups.png)

10.	Confirm an existing SSH rule exists which allows TCP port 22. To accept connections from Anywhere select the drop-down box under the **Source** column and select **Anywhere** which will correspond to **0.0.0.0/0, ::/0**. 


![EC2ConfigSecGroups](/images/hpc-aws-parallelcluster-workshop/EC2ConfigureSecGroupInbound.png)

11.	Click the **Review and Launch** button.

12.	Review your configuration and choices, and then click Launch.

13.	Select the key pair that you created in the beginning of this lab from the drop-down and check the **I acknowledge** checkbox. Then click the Launch Instances button.

![EC2ConfigSecGroups](/images/hpc-aws-parallelcluster-workshop/EC2SelectKeypair.png)


Your instance will now start, which may take a moment. You will be shown the **Launch Status** page with the message that your **instances are now launching**


14.	On the lower right of the page click on **View Instances** to view the list of EC2 instances. Click on your instance. It will go through an initialization process. Once your instance has launched, you will see your Linux server as well as the Availability Zone the instance is in, and the publicly routable DNS name.

![EC2ViewInstances](/images/hpc-aws-parallelcluster-workshop/EC2ViewInstances.png)
