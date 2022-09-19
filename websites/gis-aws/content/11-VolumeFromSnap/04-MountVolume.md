+++
title = "c. Mount the Volume and check data"
weight = 60
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++

1.	Login to your EC2 Instance.

2.	List the available disks using the following command:

```bash
lsblk
```
	
The output will list the disks you attached to your instance.

**NOTE**: Depending on the Linux version and the machine type, the device names may differ. The EC2 Console will generally show **/dev/sdX**, where X is a lower-case letter, but you may see **/dev/xvdX** or **/dev/nvmeYn1**. The following table may help with translating. Another way to help track is to pick different sizes for your EBS volumes (such as 151, 152, 153 GB for different volumes).
| Device name (Console) | Alternate 1 | Alternate 2 |
| --- | --- | --- |
| /dev/sda | /dev/xvda | /dev/nvme0n1 |
| /dev/sdb | /dev/xvdb | /dev/nvme1n1 |
| /dev/sdc | /dev/xvdc | /dev/nvme2n1 |
| /dev/sdd | /dev/xvdd | /dev/nvme3n1 |
| /dev/sde | /dev/xvde | /dev/nvme4n1 |
| /dev/sdf | /dev/xvdf | /dev/nvme5n1 |

3.  For example, **nvme1n1** is our disk, shown below. Not only does it have the proper size, but there is no preexisting mount point.

```bash
NAME        MAJ:MIN RM  SIZE RO TYPE MOUNTPOINT
nvme0n1     259:0    0  150G  0 disk
└─nvme0n1p1 259:1    0  150G  0 part /
nvme1n1     259:2    0   10G  0 disk
```

4.	Check if the volume has any data using the following command, using the disk identified in **step 3**:

```bash
sudo file -s /dev/nvme1n1
```

If the above command output shows "/dev/nvme1n1: data", it means your Volume is empty.

5.	Create a directory of your choice to mount our new ext4 volume. Let's use the name "volume1".

```bash
sudo mkdir /mnt/volume1
```

6.	Mount the volume to “volume1” directory using the following command:

```bash
sudo mount /dev/nvme1n1 /mnt/volume1
```

7.	cd into volume1 directory and check the disk space to confirm the volume mount.

```bash
cd /mnt/volume1
df -h .
```

8.	The above command will show the free space in the volume1 directory.

9.	Ensure your copied data is present on the volume.

10.	You will need to set ownership of this drive from root to user with **chown**. Replace **user** with the username for your current OS (likely ubuntu or ec2-user).  

```bash
sudo chown -R <user> /mnt/volume1
```
