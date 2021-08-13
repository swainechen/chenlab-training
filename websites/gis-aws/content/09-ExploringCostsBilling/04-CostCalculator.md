+++
title = "c. AWS Pricing Calculator"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 50
tags = ["tutorial", "aws console", "ec2"]
+++
  
AWS provides a pricing calculator that helps estimate your monthly AWS bill more efficiently. Using this tool, you can add, modify and remove services from an estimated bill and it will recalculate estimated monthly charges automatically.  
You can click [here](https://calculator.aws/#/) to link to the pricing calculator. Once you’ve navigated to the site, click **Create estimate**.   

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/AWSPricingCalculator.png)  
   
**Generating an EC2 and EBS Estimate**  
Now we’ll make a cost estimate based on our own machines. First, we’ll run an estimate using the specs for the AMI that we’ve currently got running. As a reminder, these are the resources attached:  
-  EC2 instance type: t3a.medium (2 vCPU, 4 GiB)  
-  EBS volume is 100 GB  
 
1.  On the **Select service** page, enter “EC2” into the search bar. Click **Configure** for the first option that just says “Amazon EC2”. 

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/AWSPricingCalculator_EC2.png)
 
2.  On the **Configure Amazon EC2** page, enter a useful description for our first estimate. For this first one, we can use “AMI estimate”.  

3.  Select **Quick estimate**. This will serve our purposes, but later you can play with the Advanced estimate option later.  
  
4.  Under the **EC2 instance specifications section**, we’ll search for our instance type. Here you can either search for instances by system requirements or by name. Since we know we want **t3a.medium**, we’ll search by name. For quantity, keep the value 1. Utilization, by default, is given by % Utilization / Month. For this initial test, we’ll keep it at 100%- this means that your instance will be running full time for the entire period.  

5.  For the pricing strategy, we’ll select **On-Demand Instances**. This is the more expensive option, and it also allows for a good amount of flexibility when running your analyses. You’ll use this when you have an infrequent or inconsistent workload. At this point, you can click the drop down arrow and show the calculation for projected cost, but we’ll continue through to the EBS prediction to cover all our bases.  

6. Under **Amazon Elastic Block Storage (EBS)**, keep the default General Purpose SSD (gp2) but change the storage to reflect our actual EBS volume’s size, 100gb.  

7.  Now we can check the total pricing. Both EBS and EC2 are calculated individually (can click “show calculations” in the respective section to view the specific breakdown), and a whole estimate is provided at the bottom of the page.  

8.  Now use the calculator to explore the following scenarios:  
-  Instead of leaving your instance on for the entire month, you only use it during work hours (8 hours a day).  
-  Remember that c5.2xlarge instance we started in the very beginning? What would happen if we left that running for the entire month?  

-  How would different pricing plans affect the cost of leaving the c5.2xlarge instance running?  
  
  
**Pricing S3 Costs**  
Next, we’ll look at the costs of storing data on AWS’s Simple Storage Service (S3). Let’s hypothetically say you’ve got a dataset of 100,000 bacterial sequences, with each sequence file being around 5.25 megabases, and these are all saved as fasta files. This will be approximately 160 GB of data. We’ll use the S3 pricing calculator to see how storage costs of this database change depending on storage tier. 
  
1.  Return to the **Pricing Calculator** main page. Enter “S3” into the search bar. Click **Configure** for the first option that just says “Amazon Simple Storage Service (S3)”. 

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/AWSPricingCalculator_S3.png)  

2.  On the **Configure Amazon Simple Storage Service (S3)** page, enter a useful description for our first estimate. For this first one, we can use “S3 estimate”.  

3.  Make sure your region is **Asia Pacific (Singapore)**. 
 
4.  Under **Select S3 Storage classes and other features**, make sure the following are selected: S3 Standard, S3 Glacier, and S3 Glacier Deep Archive.  

5.  For all three sections, put 160 GB as our **storage**. As you scroll down, you can see how the storage class affects the cost for the same amount of data.  
  
6.  Feel free to explore with the pricing calculator further- how much would it cost to store 160 GB of sequences on an EBS? How does this compare to S3?  
 
For more details on how to use the calculator, navigate to the [AWS user guide](https://docs.aws.amazon.com/pricing-calculator/latest/userguide/getting-started.html).  


