+++
title = "b. Launch an EC2 Instance"
weight = 80
tags = ["tutorial", "cloud9", "aws cli", "s3"]
+++

We will now launch an EC2 Linux-based instance.

1.	Click on **Launch Instance**, then click **Launch instance** again from the drop down menu.

![EC2LaunchInstance](/images/hpc-aws-parallelcluster-workshop/EC2LaunchInstance.png)

2.	On the **Launch an instance** page, you can choose a friendly name for your instance.

![EC2LaunchInstance](/images/hpc-aws-parallelcluster-workshop/EC2InstanceName-2.png)

3.	Click on **Add additional tags** and  "Add Tag". You will find the "Name" you provided for your instance. Now enter the key and value. These keys, more correctly known as a tags, will appear in the console once the instance launches. It makes it easy to keep track of running machines in a complex environment. Create additional tags similar for the ones you used for the key pair previously - give a user and grant key to this machine, and enter the same values. When ready, select **Instances, Volumes, Network interfaces** under the **Resource types**.

![EC2AddTags](/images/hpc-aws-parallelcluster-workshop/EC2AddTags-2.png)

4.	Under **Application and OS images (Amazon Machine Image)** enter the AMI ID shared for the training into the search box. The **AMI** may differ from the one pictured; refer to instructor's directions.  

**Note:**
> For this training, search for `CHENLAB-PUBLIC` and choose the most recent AMI version; as of June, 2022 this is v2.10, which has ID `ami-0e688975feec36afb`. You can find out more about CHENLAB-PUBLIC AMIs in [Appendix: About the AMI](http://slchen-lab-training.s3-website-ap-southeast-1.amazonaws.com/14-appendix/06-aboutyourami.html).

![EC2SelectAMI](/images/hpc-aws-parallelcluster-workshop/EC2SearchAMI.png)

The results of the search appear  **in Community AMIs** and **in AWS Marketplace**.

![EC2SelectAMI](/images/hpc-aws-parallelcluster-workshop/EC2SearchAMIResult.png)


4.	Click on the results **in Community AMIs**.

The AMI ID will now appear.

![EC2SelectAMI](/images/hpc-aws-parallelcluster-workshop/EC2SearchAMISharedWithMe.png)

5.	Click on the **Select** button to select the AMI. 

6.	Under **Instance type**, click on the drop-down arrow and type **t3a.medium** into the search bar.
**Note:** We're doing this on purpose. This will not be a default instance you should use until you've considered the resources you need for your project.

![EC2ChooseInstancetype](/images/hpc-aws-parallelcluster-workshop/EC2ChooseInstanceType-2.png)

7.	Select the key pair that you created in the beginning of this lab from the drop-down list under **Key pair (login)**.

![EC2ConfigSecGroups](/images/hpc-aws-parallelcluster-workshop/EC2SelectKeypair-2.png)

8.	Next click on **Edit** against **Network settings**. You will be prompted to provide **Subnet** and **Security group** details. The security groups will be your firewall rules.  
&nbsp;&nbsp;&nbsp; a.   Note that the Subnet field can be configured to launch the instance in a specific Availability Zone; while we are keeping the default for this workshop, this gives you control over the location of your machine.  
&nbsp;&nbsp;&nbsp; b.   Provide a name for your new security group. **NOTE**: For the purpose of this lab use your "[Your initials]-sec-group".  
&nbsp;&nbsp;&nbsp; c.   Confirm an existing SSH rule exists which allows TCP port 22.  
&nbsp;&nbsp;&nbsp; d.   To accept connections from anywhere, select the drop-down box under the **Source type** column and select **Anywhere** which will correspond to **0.0.0.0/0**.  

![EC2ConfigSecGroups](/images/hpc-aws-parallelcluster-workshop/EC2ConfigSecGroups-2.png)

9.	Under **Configure storage**, you have the ability to modify or add storage and disk drives to the instance. For this lab, we will simply accept the storage defaults.

![EC2ChooseInstancetype](/images/hpc-aws-parallelcluster-workshop/EC2AddStorage-2.png)

10.	Review your configuration under **Summary** and click **Launch Instance**.

![EC2ChooseInstancetype](/images/hpc-aws-parallelcluster-workshop/EC2StepConfigureInstance-2.png)

Your instance will now start, which may take a moment. You will be shown the **Launch Status** page with the message **Successfully initiated launch of instance**.

11.	On the lower right of the page click on **View all Instances** to view the list of EC2 instances. Click on your instance. It will go through an initialization process. Once your instance has launched, you will see your Linux server as well as the Availability Zone the instance is in, and the publicly routable DNS name.

![EC2ViewInstances](/images/hpc-aws-parallelcluster-workshop/EC2ViewInstances-2.png)

