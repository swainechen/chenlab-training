+++
title = "AWS Configure - Event Engine"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 50
tags = ["tutorial", "Prerequisite", "ec2"]
+++
   
_Important: As you are using Event Engine for this workshop, you will be using the code snippet with the temporary account's credentials from the Console page. This will be analogous to going through the manual setting up described next in points 1-4 under "AWS Configure - General Use"._   

1.	Copy the **Credentials/CLI Snippet** onto the command line.  
![EE console](/images/hpc-aws-parallelcluster-workshop/EE_console_login.png)    

2.	Now rerun the ec2 **describe-instances** command again and check the output.

```bash
aws ec2 describe-instances
```

This gives a description of all the EC2 instances in the account for the specified region.

3.	 Let's take a look at the key-pairs we have for the selected region.

```bash
aws ec2 describe-key-pairs
```
