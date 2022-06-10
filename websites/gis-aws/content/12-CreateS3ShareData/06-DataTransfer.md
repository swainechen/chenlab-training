+++
title = "Add'l-2: Remote data transfer"
weight = 150
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++


In this section we will use the AWS CLI (Command Line Interface) to copy data from:

1.	Your local machine to your AWS EC2 instance and vice versa.
2.	One EC2 instance to another EC2 instance.

---
**Local machine to EC2 instance and vice-versa**
1.	On your local computer, open a new terminal and create a sample text file.

```bash
echo “Hello, world” > test.txt
```

2.	Use scp to copy this sample test file into the home directory of your remote machine/ instance, by running the following command on your local machine.
```bash
scp -i path/to/key.pem path/to/test.txt ubuntu@<ip-address>:/home/ubuntu/
```

If you are using windows command prompt and have the PuTTY suite of utilities installed, run the following command
```bash
pscp -P 22 -i path/to/key.ppk path/to/test.txt ubuntu@<ip-address>:/home/ubuntu/
```

Alternative: use rsync - the options here preserve timestamps and permissions (-a), and only transfer if it would result in an updated version (-u). This can help with tracking what the latest file is; using scp as above will just overwrite the destination even if the file exists at the destination and is newer.
```bash
rsync -e 'ssh -i /path/to/key.pem' -avxu /path/to/test.txt ubuntu@<ip-address>:/home/ubuntu/
```

3.	Open the test.txt file using a text editor on your remote machine and edit it - add random text and save it.

4.	Copy the test.txt file back to your local machine by using the scp/pscp command on your local terminal.

```bash
scp -i path/to/key.pem ubuntu@<ip-address>:/home/ubuntu/test.txt path/to/folder/localmachine
```
If you are using windows command prompt run the following command
```bash
pscp -P 22 -i path/to/key.ppk ubuntu@<ip-address>:/home/ubuntu/test.txt path/to/folder/localmachine
```

Alternative: use rsync:
```bash
rsync -e 'ssh -i /path/to/key.pem' -avux ubuntu@<ip-address>:/home/ubuntu/test.txt /path/to/folder/localmachine
```

5.	Open the copied file to review your edits.

---
**EC2 to EC2 instance**   
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

