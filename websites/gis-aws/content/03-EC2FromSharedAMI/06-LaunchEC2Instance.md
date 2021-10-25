+++
title = "b. Launch an EC2 Instance"
weight = 80
tags = ["tutorial", "cloud9", "aws cli", "s3"]
+++

We will now launch an EC2 Linux-based instance.

1.	Click on **Launch Instance**.

![EC2LaunchInstance](/images/hpc-aws-parallelcluster-workshop/EC2LaunchInstance.png)

2.	In **Step1: Choose an Amazon Machine Image(AMI)**, enter the AMI ID shared for the training into the **Search for an AMI...** search box. The **AMI** may differ from the one pictured; refer to instructor's directions.  
For this training, search for `CHENLAB-PUBLIC` and choose the AMI which has ID `ami-07b3bbf5c844d85c9`. You can find out more about this AMI in Appendix: About the AMI. 

![EC2SelectAMI](/images/hpc-aws-parallelcluster-workshop/EC2SearchAMI.png)

The results of the search appear  **in Community AMIs** and **in AWS Marketplace**.

![EC2SelectAMI](/images/hpc-aws-parallelcluster-workshop/EC2SearchAMIResult.png)


3.	Click on the results **in Community AMIs**.

The AMI ID will now appear.

![EC2SelectAMI](/images/hpc-aws-parallelcluster-workshop/EC2SearchAMISharedWithMe.png)

4.	Click on the **Select** button to select the AMI.

5.	In the **Step 2** page, choose an instance type. Select the **t3a.medium** instance from the list and click **Next: Configure Instance Details**.

![EC2ChooseInstancetype](/images/hpc-aws-parallelcluster-workshop/EC2ChooseInstanceType.png)

6.	On the **Step 3** page, Configure Instance Details - leave the default settings. Note the **Subnet** field can be configured to launch the instance in a specific **Availability Zone**. Click on **Next: Add Storage** button at the bottom right corner. [Learn More](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/configuring-instance-metadata-service.html)

![EC2ChooseInstancetype](/images/hpc-aws-parallelcluster-workshop/EC2StepConfigureInstance.png)


7.	On **Step 4** page you have the ability to modify or add storage and disk drives to the instance. For this lab, we will simply accept the storage defaults and click **Next: Add Tags**.

**NOTE**: Notice the default selection for the **Delete on Termination** box is checked - this indicates that if the machine is terminated the root volume associated with the instance will be deleted, you need to uncheck this if you plan to store data on the root volume which you would want to access later.

![EC2ChooseInstancetype](/images/hpc-aws-parallelcluster-workshop/EC2AddStorage1.png)


8.	On the **Step 5** page, you can choose a friendly name for your instance by clicking "Add Tag", and entering "Name" for the **Key** and "[your initials]-linux-server" as the **value**. This Name key, more correctly known as a **tag**, will appear in the console once the instance launches. It makes it easy to keep track of running machines in a complex environment. Add the same tags that we used before- one tag for user, and one for grant. Keep the values consistent. Click **Next: Configure Security Group**.

![EC2AddTags](/images/hpc-aws-parallelcluster-workshop/EC2AddTags.png)

9.	On **Step 6** page, you will be prompted to create a new security group, which will define your firewall rules. If you have not yet made a security group, provide a name for your new security group. Otherwise, click **Select an existing security group**, select your previous key, and skip ahead to **Step 11**.

**NOTE**: For the purpose of this lab use your "[your-initials]-sec-group" for uniqueness and easy searchability.

![EC2ConfigSecGroups](/images/hpc-aws-parallelcluster-workshop/EC2ConfigSecGroups-2.png)

10.	Confirm an existing SSH rule exists which allows TCP port 22. To accept connections from Anywhere select the drop-down box under the **Source** column and select **Anywhere** which will corrospond to **0.0.0.0/0, ::/0**.


![EC2ConfigSecGroups](/images/hpc-aws-parallelcluster-workshop/EC2ConfigureSecGroupInbound.png)

11.	Click the **Review and Launch** button.

12.	Review your configuration and choices, and then click **Launch**.

13.	Select the key pair that you created in the beginning of this lab from the drop-down and check the **I acknowledge** checkbox. Then click the **Launch Instances** button.

![EC2ConfigSecGroups](/images/hpc-aws-parallelcluster-workshop/EC2SelectKeypair-2.png)


Your instance will now start, which may take a moment. You will be shown the **Launch Status** page with the message that your **instances are now launching**.


14.	On the lower right of the page click on **View Instances** to view the list of EC2 instances.  The instance will go through and initialization process. Once your instance has launched, you will see your Linux server as well as the Availability Zone the instance is in, and the publicly routable DNS name.

![EC2ViewInstances](/images/hpc-aws-parallelcluster-workshop/EC2ViewInstances-2.png)
