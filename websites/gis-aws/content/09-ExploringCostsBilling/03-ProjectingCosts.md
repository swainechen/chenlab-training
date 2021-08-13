+++
title = "b. Service Pricing Pages"
date = 2019-09-18T10:46:30-04:00
draft = false
weight = 50
tags = ["tutorial", "aws console", "ec2"]
+++

Projecting costs before utilizing resources is not only required for budgeting, but can also help optimise your usage of AWS resources. There are several online tools that can help project costs. Here is a quick overview on how to view pricing for specific AWS services.

**Exploring the specific AWS Service pages**  
Every AWS service has its own web page with a **pricing** section to provide more information on costs associated with the service. Here are the pricing web page links for EC2, S3, and EBS. 


- [EC2 Pricing page](https://aws.amazon.com/ec2/pricing/)  

- [EBS Pricing page](https://aws.amazon.com/ebs/pricing/)  

- [S3 Pricing page](https://aws.amazon.com/s3/pricing/)  

**Exploring EC2 Pricing**  

1.  On the **EC2 Pricing page** scroll down and click **On-Demand Pricing**  

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/OnDemandPricing.png)  

2.  On the **Amazon EC2 On-Demand Pricing** page you can explore the pricing for EC2 instance per hour for vatious operating systems. Notice that the pricing is specific on the regional level.  

3.  To view the pricing for a specific region click on the **Region** drop-down menu and select **Asia-Pacific (Singapore)**; the pricing list will dynamically change to reflect the pricing for that region.  

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/OnDemandPricingRegion.png)  
  
  
**Exploring EBS Pricing**  

1.  At the top of the **EBS Pricing page**, change the **Region** to **Asia-Pacific (Singapore)**.

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/OnDemandPricingRegion_EBS.png)  
 
2.  Scroll down to explore the division of resources; here, it’s divided into Amazon EBS Volumes, EBS Snapshots, restoring costs, and practical examples. At the bottom of the page is a pricing calculator; you’ll read more about this in our next section.   
  
  
**Exploring S3 Pricing** 
 
1.  The first thing to note for the **AWS S3 Pricing page** is the series of tabs at the top- these let you explore the cost not only for storage, but for moving, retrieving, and replicating your data.  

![AWS Management Console](/images/hpc-aws-parallelcluster-workshop/OnDemandPricingRegion_S3.png)

2.  On the **Storage** tab, change your **Region** to **Asia-Pacific (Singapore)** (as shown above). Scroll down. 
 
3.  Here you can see the pricing tiers for different levels of storage; it’s important to note how S3’s intelligent tiering can move your data to different classes, and save you money.  
 
   
Next, we’ll synthesize all of this information by using the AWS pricing calculator to estimate the cost of running a few jobs on AWS. 
