+++
title = "c. Create a new Key Pair"
weight = 70
tags = ["tutorial", "cloud9", "ParallelCluster"]
+++

In this section of the lab you will create a Linux insnatce and login to it to run a few linux commands. To manage the instance, you need to be able to connect to it via SSH. The following steps outline how to create a unique SSH key pair for this purpose.

1.	Sign into the AWS Management Console and open the Amazon EC2 console at 
[https://console.aws.amazon.com/ec2](https://console.aws.amazon.com/ec2)

2.	In the upper-right corner of the AWS Management Console, confirm you are in the desired AWS region **(e.g., N. Virginia)**.

![ConsoleEC2KeyPair](/images/hpc-aws-parallelcluster-workshop/ConsoleRegion.png)

3.	In the left navigation pane click on **Key Pairs** under the **Network & Security** section.  This will display a page to manage your SSH key pairs. 

![ConsoleEC2KeyPair](/images/hpc-aws-parallelcluster-workshop/EC2KeyPair.png)


4.	On the **Key Paris** page click the **Create Key Pair** button at the top of the browser window.

![ConsoleEC2CreateKeyPair](/images/hpc-aws-parallelcluster-workshop/EC2CreateKeyPair.png)

5.	In the resulting pop up window, provide a key pair name of your choosing and select **.ppk** or **.pem** depending on which ssh client you will use to connect to the EC2 Instance (for Putty select **.ppk** .For  Windows Susbsystem Linux, MacOS-X ssh Client, Linux ssh client,  Windows ssh client select **.pem**) and click on **Create key pair**  


**NOTE**: For the purpose of this workshop please provide a unique key pair name, so that it might not clash with key pair names provided by other participants using the same AWS account. You could use your [first name + last name] in the key pair for uniqueness easy identification.


![ConsoleEC2CreateKeyPairDialog](/images/hpc-aws-parallelcluster-workshop/EC2CreateKeyPairDialog.png)

6.	The key pair you created should automatically download to your system.  Follow any browser instructions to save the file to the default download location.  Verify that the file has downloaded and you can locate it. Remember the full path to this .ppk/pem file you just downloaded. This file contains your private key for future SSH connections. 

7.	You will see a message appear at the top of the screen that says **Successfully created key pair**.  You will see the key pair you created listed.
