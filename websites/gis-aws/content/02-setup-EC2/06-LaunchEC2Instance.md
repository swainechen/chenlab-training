+++
title = "e. Launch an EC2 Instance"
weight = 80
tags = ["tutorial", "cloud9", "aws cli", "s3"]
+++

We will now launch an EC2 Linux based instance.


1.	Sign into the AWS Management Console and open the Amazon EC2 console at https://console.aws.amazon.com/ec2

2.	In the upper-right corner of the AWS Management Console, confirm you are in the desired AWS region **(e.g., N. Virginia)**.

![EC2LaunchInstance](/images/hpc-aws-parallelcluster-workshop/ConsoleRegion.png)


3.	Click on  **Launch instance** and further choose and click **Launch instance**

![EC2LaunchInstance](/images/hpc-aws-parallelcluster-workshop/EC2LaunchInstance.png)

4.	Next in the **Step 1** page, select the first **Amazon Linux 2 AMI(HVM), SSD Volume Type** and further check **64-bit (x86)** under the **Select** button and click on the **Select**. 

**NOTE**: that the ami-xxxxxxxxx label and specific versions of the installed package may be different than in the image below.

![EC2SelectAMI](/images/hpc-aws-parallelcluster-workshop/EC2SelectAMI.png)

5.	In the **Step 2** page, choose an Instance type, select the **t3a.medium** instance from the list and click **Next: Configure Instance Details**.

![EC2ChooseInstancetype](/images/hpc-aws-parallelcluster-workshop/EC2ChooseInstanceType.png)

6.	On **Step 3** page, Configure Instance Details - leave the default settings, One aspect to note here is the **Subnet** field that can be configured to launch the instance in a specific **Availability Zone**. Click on **Next: Add Storage** button at the bottom right corner . [Learn More](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/configuring-instance-metadata-service.html)

![EC2ChooseInstancetype](/images/hpc-aws-parallelcluster-workshop/EC2StepConfigureInstance.png)


7.	On **Step 4** page you have the ability to modify or add storage and disk drives to the instance. For this lab, we will simply accept the storage defaults and click **Next: Add Tags**.

**NOTE**: Take notice of the default selection of the **Delete on Termination** box checked - this indicates that if the machine is terminated the Root Volume associated with the Instance will be deleted, you need to uncheck this if you plan to store data on the root volume which you would want to access later.

![EC2ChooseInstancetype](/images/hpc-aws-parallelcluster-workshop/EC2AddStorage.png)


8.	On **Step 5** page, you can choose a friendly name for your instance by clicking "Add Tag", and entering "[your initials]-training" for the Key part and "[Your initials]-linux-machine" for the Value part. This Name key, more correctly known as a tag, will appear in the console once the instance launches. It makes it easy to keep track of running machines in a complex environment. Click **Next: Configure Security Group**

![EC2AddTags](/images/hpc-aws-parallelcluster-workshop/EC2AddTags-2.png)

9.	On **Step 6** page, you will be prompted to create a new security group, which will be your firewall rules. Provide a name for your new security group.

**NOTE** for the purpose of this lab use your "[Your initials]-sec-group" for uniqueness and easy searchability.

![EC2ConfigSecGroups](/images/hpc-aws-parallelcluster-workshop/EC2ConfigSecGroups.png)

10.	Confirm an existing SSH rule exists which allows TCP port 22. To accept connections from Anywhere select the drop-down box under the **Source** column and select **Anywhere** which will corrospond to **0.0.0.0/0, ::/0**. 


![EC2ConfigSecGroups](/images/hpc-aws-parallelcluster-workshop/EC2ConfigureSecGroupInbound.png)

11.	Click the **Review and Launch** button.

12.	Review your configuration and choices, and then click Launch.

13.	Select the key pair that you created in the beginning of this lab from the drop-down and check the **I acknowledge** checkbox. Then click the Launch Instances button.

![EC2ConfigSecGroups](/images/hpc-aws-parallelcluster-workshop/EC2SelectKeypair.png)


Your instance will now be starting, which may take a moment. You will be shown the **Launch Status** page with the message that your **instances are now launching**


14.	On the lower right of the page click on **View Instances** to view the list of EC2 instances.  The instance will go through an initialization process. Once your instance has launched, you will see your Linux server as well as the Availability Zone the instance is in, and the publicly routable DNS name.

![EC2ViewInstances](/images/hpc-aws-parallelcluster-workshop/EC2ViewInstances.png)
