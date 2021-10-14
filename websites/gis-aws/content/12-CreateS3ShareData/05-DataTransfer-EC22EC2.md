+++
title = "EC2 to EC2 instance"
draft = false
weight = 130
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++

Now let us copy the sample test.txt file from your current instance (CI) to another instance (AI) you have earlier started.

1.	Open a new terminal, and use it to ssh into the instance you started at the beginning of this workshop (AI); if it is unavailable, quickly spin up another instance. The specs do not matter.

- Make sure to copy your .pem key (for AI) to your CI. This is the private key file that you are using to log in to the instances, i.e. the key pair you have been specifying with the `-i` option in SSH. Be sure to change the permissions as well (chmod 0400 keyname.pem). You can use scp to copy it to your new CI. (You will need the key file in order to communicate (SSH or SCP) with another instance with the same key pair)
- SSH into your AWS instance (AI) from your local machine where you have the text.txt file.
- Run the following command to copy the file from your current instance to the other instance (AI)
```bash
scp -i path/to/key.pem path/to/test.txt ec2-user@<ip-address>:/home/ec2-user/
```
NOTE: Username to be used here are **ec2-user** or **ubuntu** if your instance was spun up from an Amazon Linux OS or Ubuntu Linux OS images, respectively.

2.	SSH into the instance you copied the file into from your current instance using the following command

```bash
ssh -i path/to/key.pem ec2-user@<ip-address>
```

3.	Review the file you copied to this instance. Type “exit” to simply return to your earlier instance.

