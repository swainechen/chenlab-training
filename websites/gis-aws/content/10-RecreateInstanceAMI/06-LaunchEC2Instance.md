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

3.	Select the AMI created in the previous section of the workshop then click the **Launch instance from Image** button.

![EC2SelectAMI](/images/hpc-aws-parallelcluster-workshop/EC2LaunchAMI.png)

4.	On the **Launch an instance** page, you can choose a friendly name for your instance.

![EC2LaunchInstance](/images/hpc-aws-parallelcluster-workshop/EC2InstanceName-2.png)

5.	Click on "Add additional tags" and  "Add Tag". You will find the "Name" you provided for your instance. Now enter the key and value. These keys, more correctly known as a tags, will appear in the console once the instance launches. It makes it easy to keep track of running machines in a complex environment. Create additional tags similar for the ones you used for the key pair previously - give a user and grant key to this machine, and enter the same values. When ready, select **Instances, Volumes, Network interfaces** under the **Resource types**.

![EC2AddTags](/images/hpc-aws-parallelcluster-workshop/EC2AddTags-2.png)

6.	Under **Instance type**, click on the drop-down arrow and type **t2.2xlarge** into the search bar. 
**Note:** This is a large instance - we're doing this on purpose. This will be required for analysis in Section XII. This is not to be considered as a default instance until you've considered the resources based on your project needs.

![EC2ChooseInstancetype](/images/hpc-aws-parallelcluster-workshop/EC2ChooseInstanceType-3.png)

7.	Select the key pair that you created in the beginning of this lab from the drop-down list under **Key pair (login)**.

![EC2ConfigSecGroups](/images/hpc-aws-parallelcluster-workshop/EC2SelectKeypair-2.png)

8.	Next click on **Edit** against **Network settings**. You will be prompted to provide **Subnet** and **Security group** details. The security groups will be your firewall rules.
&nbsp;&nbsp;&nbsp; a.   Note that the Subnet field can be configured to launch the instance in a specific Availability Zone; while we are keeping the default for this workshop, this gives you control over the location of your machine.
&nbsp;&nbsp;&nbsp; b.   Choose **Select existing security group**

![EC2ConfigSecGroups](/images/hpc-aws-parallelcluster-workshop/EC2ConfigSecGroups-3.png)

9.	Under **Configure storage**, you have the ability to modify or add storage and disk drives to the instance. For this lab, we will simply accept the storage defaults.

![EC2ChooseInstancetype](/images/hpc-aws-parallelcluster-workshop/EC2AddStorage-2.png)

10.	Review your configuration under **Summary** and click **Launch Instance**.

![EC2ChooseInstancetype](/images/hpc-aws-parallelcluster-workshop/EC2StepConfigureInstance-3.png)

Your instance will now start, which may take a moment. You will be shown the **Launch Status** page with the message **Successfully initiated launch of instance**.

11.	On the lower right of the page click on **View all Instances** to view the list of EC2 instances. Click on your instance. It will go through an initialization process. Once your instance has launched, you will see your Linux server as well as the Availability Zone the instance is in, and the publicly routable DNS name.

![EC2ViewInstances](/images/hpc-aws-parallelcluster-workshop/EC2ViewInstances-3.png)

