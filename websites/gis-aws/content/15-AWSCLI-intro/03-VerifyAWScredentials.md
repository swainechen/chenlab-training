+++
title = "b. Verifying AWS credentials"
date = 2021-06-22T10:46:30-04:00
draft = false
weight = 70
tags = ["tutorial", "Prerequisite", "ec2"]
+++

Before we interact with the Amazon S3 let us take a look at the importance of the AWS credentials. 

AWS security credentials are used to verify
1.	Who you are
2.	 Your permission to access the resources that you are requesting  
AWS uses these security credentials to authenticate and authorize your requests.  
  
We configured **who you are** in the earlier section. Let us inspect the credentials and config files in the **~/.aws/** folder now.
**_(Not applicable for AWS configure with Event Engine code snippet)_**  
  
```bash
cat ~/.aws/credentials
```

```bash
cat ~/.aws/config
```  
  
  
  
>_IMPORTANT: If you are using Event Engine for this workshop, make sure to have the **AWS ACCESS KEY ID**, **AWS SECRET ACCESS KEY**, **AWS DEFAULT REGION** match the temporary account's credentials from the Console page on Event Engine. The **AWS SESSION TOKEN** is unique for each Event Engine session (this session token is NOT APPLICABLE for an original AWS account)._   
   
   
    
We now successfully configured **who you are**. Next, we will configure **your permission to access the resources that you are requesting (i.e. with a "named profile")** to interact with Amazon S3.
