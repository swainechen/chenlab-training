---
title: "Attach additional storage"
date: 2019-01-24T09:05:54Z
weight: 60 
pre: "<b>IV ⁃ </b>"
tags: ["HPC", "Introduction", "EC2", "Optional"]
---

The created EC2 instance has a root volume with enough storage capacity to support the operating system, as well as a few additional applications required for supporting the analyses.
We will now need additional storage to accommodate the large datasets we are about to download (e.g. Raw reads, aligned sequences etc.).
Hence we need to attach additional storage.

In this section we will:

**a.	[Create an Elastic Block Store (EBS) volume.]({{< ref "04-attachstorage/02-createebsvolume.html" >}})**

**b.	[Attach the new volume to an EC2 instance.]({{< ref "04-attachstorage/03-attachvolume.html" >}})**

**c. 	[Format and mount the new volume.]({{< ref "04-attachstorage/04-mountvolume.html" >}})**

&nbsp;&nbsp;&nbsp;	1.	Format the new volume with a filesystem.

&nbsp;&nbsp;&nbsp;	2.	Mount the volume and copy some data into it.
