+++
title = "c. SSH into an EC2 instance"
weight = 90
tags = ["tutorial", "cloud9", "s3"]
+++

The operating system of your local computer determines the options that you have to connect from your local computer to your Linux instance. Here are the directions for **Putty**. Click on the links to select the directions for other ssh clients:
[Windows Susbsystem for Linux](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/WSL.html),
[MacOS-X ssh Client](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/AccessingInstancesLinux.html),
[Linux ssh client](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/AccessingInstancesLinux.html)


**Directions for PuTTY**

We will now use the key pair with the PuTTy SSH client to open an SSH session to the EC2 instance.

1.	Start the PuTTy application.

2.	Enter in the IP address of the EC2 instance in the **Host Name (or IP address)** field.

3.	To save the session (so we dont have to repeat the configuration everytime) provide a name in the **Saved Sessions** field.

4.	In the left navigation pane under **Connections** click on **SSH** and next click on **Auth**.

5.	In the **Private key file for authentication** field click **Browse** and browse in the **.ppk** file downloaded when creating the key pair earlier.

6.	Next in the left navigation pane click on **Session** (**scroll up if you need to**) and and further click on **Save** to save the session.

7.	Next click on **Open** to start the SSH session.

8.	PuTTy will now attempt to start an SSH Session.

**NOTE**: You will get a security alert **The Servers host key is not cached in the registry**, click on **Yes** to continue.

9.	You will now be presented with the **login as** prompt.

10.	At the **login as** prompt type **ec2-user** and hit enter.

You would have now logged into the EC2 instance on AWS.
