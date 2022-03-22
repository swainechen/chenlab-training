+++
title = "b. Launch an EC2 Instance from AMI"
weight = 80
tags = ["tutorial", "cloud9", "aws cli", "s3"]
+++

We will now launch an EC2 Linux-based Instance.

1.	In the upper-right corner of the AWS Management Console, confirm you are in the desired AWS region (e.g., Singapore).

![EC2LaunchInstance](/images/hpc-aws-parallelcluster-workshop/ConsoleRegion.png)

2.	In the left navigation pane click on **AMIs** under the **Images** section.

![EC2LaunchInstance](/images/hpc-aws-parallelcluster-workshop/EC2AMIs.png)

3.	Select the AMI created in the previous section of the workshop then click the **Launch** button.

![EC2SelectAMI](/images/hpc-aws-parallelcluster-workshop/EC2LaunchAMI.png)

4.	Click on the **Select** Button to select the AMI.

5.	In the **Step 2** page, choose an Instance type, select the **t3a.medium** Instance from the list and click **Next: Configure Instance Details**.

![EC2ChooseInstancetype](/images/hpc-aws-parallelcluster-workshop/EC2ChooseInstanceType.png)

6.	On **Step 3** page, Configure Instance Details - leave the default settings. Note that the **Subnet** field can be configured to launch the instance in a specific **Availability Zone**. Make sure this availability zone matches your previous resources. Click on **Next: Add Storage** button at the bottom right corner . [Learn More](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/configuring-instance-metadata-service.html)

![EC2ChooseInstancetype](/images/hpc-aws-parallelcluster-workshop/EC2StepConfigureInstance.png)


7.	On the **Step 4** page, you have the ability to modify or add storage and disk drives to the instance. For this lab, we will simply accept the storage defaults and click **Next: Add Tags**.

**NOTE**: Take note that **Delete on Termination** box is checked by default- this indicates that if the machine is terminated, the Root Volume associated with the Instance will be deleted. You need to uncheck this if you plan to store data on the Root Volume which you want to preserve after termination. 

![EC2ChooseInstancetype](/images/hpc-aws-parallelcluster-workshop/EC2AddStorage1.png)


8.	On the **Step 5** page, you can choose a friendly name for your instance by clicking "Add Tag", and entering "Name" for the Key, and "[your initials]-new-server" for the Value. This Name key, more correctly known as a tag, will appear in the console once the instance launches. It makes it easy to keep track of running machines in a complex environment. Click **Next: Configure Security Group**.

![EC2AddTags](/images/hpc-aws-parallelcluster-workshop/EC2AddTags.png)

9.	On the **Step 6** page, you will be prompted to create a new security group, which will set your firewall rules. Select your previous security group from the drop-down list.   

![EC2ConfigSecGroups](/images/hpc-aws-parallelcluster-workshop/EC2ConfigSecGroups-2.png)

10.	Confirm an existing SSH rule exists which allows TCP port 22. To accept connections from anywhere select the drop-down box under the **Source** column and select **Anywhere** which will correspond to **0.0.0.0/0, ::/0**.


![EC2ConfigSecGroups](/images/hpc-aws-parallelcluster-workshop/EC2ConfigureSecGroupInbound.png)

11.	Click the **Review and Launch** button.

12.	Review your configuration and choices, and then click **Launch**.

13.	Select the key pair that you created in the beginning of this lab from the drop-down list and check the **I acknowledge** checkbox. Then click the **Launch Instances** button.

![EC2ConfigSecGroups](/images/hpc-aws-parallelcluster-workshop/EC2SelectKeypair-2.png)

Your Instance will now start, which may take a moment. The **Launch Status** page will show the message that your **instances are now launching**.


14.	On the lower right of the page click on **View Instances** to view the list of EC2 Instances.  The Instance will go through an initialization process. Once your Instance has launched, you will see your Linux server as well as the Availability Zone the Instance is in, and the publicly routable DNS name. The page may not automatically refresh, so check periodically. 

![EC2ViewInstances](/images/hpc-aws-parallelcluster-workshop/EC2ViewInstances-3.png)
