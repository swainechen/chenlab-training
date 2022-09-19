+++
title = "c. Mount the Volume"
weight = 60
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++

1.	Log in to your EC2 Instance.

2.	List the available disks using the following command:

```bash
lsblk
```
	
The output will list the disks attached to your instance.

**NOTE**: Depending on the Linux version and the machine type, the device names may differ. The EC2 Console will generally show **/dev/sdX**, where X is a lower-case letter, but you may see **/dev/xvdX** or **/dev/nvmeYn1**. The following table may help with translating. Another way to help track is to pick different sizes for your EBS volumes (such as 151, 152, 153 GB for different volumes).

| Device name (Console) | Alternate 1 | Alternate 2 |
| --- | --- | --- |
| /dev/sda | /dev/xvda | /dev/nvme0n1 |
| /dev/sdb | /dev/xvdb | /dev/nvme1n1 |
| /dev/sdc | /dev/xvdc | /dev/nvme2n1 |
| /dev/sdd | /dev/xvdd | /dev/nvme3n1 |
| /dev/sde | /dev/xvde | /dev/nvme4n1 |
| /dev/sdf | /dev/xvdf | /dev/nvme5n1 |

3.	Check for an unmounted filesystem of size 10GB. For example: "nvme1n1" as shown below:

```bash
NAME        MAJ:MIN RM  SIZE RO TYPE MOUNTPOINT
nvme0n1     259:0    0  150G  0 disk
└─nvme0n1p1 259:1    0  150G  0 part /
nvme1n1     259:2    0   10G  0 disk
```

4.	Check if the volume has any data using the following command:

```bash
sudo file -s /dev/nvme1n1
```

Where "nvme1n1" is the device you noted from the previous section after attaching the device to the EC2 Instance.

	If the above command output shows "/dev/nvme1n1: data", it means your volume is empty.

5.	Format the volume to the ext4 filesystem using the following command.

```bash
sudo mkfs -t ext4 /dev/nvme1n1
```

**NOTE**: This file-system formatting step is only for a new device, **DO NOT** run this step while mounting an existing volume as it will wipe out all data on the device.

6.	Create a directory of your choice to mount our new ext4 volume. Let's use the name "volume1".

```bash	
sudo mkdir /mnt/volume1
```


7.	Mount the volume to “volume1” directory using the following command.

```bash
sudo mount /dev/nvme1n1 /mnt/volume1
```


8.	cd into the volume1 directory and check the disk space for confirming the volume mount.

```bash
cd /mnt/volume1
df -h .
```
The above command would show the free space in the volume1 directory.


9.  At this point, the drive is owned by root and not user. We will want to change ownership of the drive, so that you can change the contents of the drive (Add/remove files, etc).  

```bash  
sudo chown -R ubuntu /mnt/volume1
```  

10.	For your own information, it is possible to later remove this device. Practice unmounting, then remounting it. To unmount the volume, you have to use the following command. Make sure to be outside the directory to unmount the volume.

```bash
sudo umount /dev/nvme1n1
```  
But we'll need this device for later, so remember to re-mount it.  

```bash
sudo mount /dev/nvme1n1 /mnt/volume1
```  
