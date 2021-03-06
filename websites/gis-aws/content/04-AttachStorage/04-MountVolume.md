+++
title = "c. Mount the Volume"
weight = 60
tags = ["tutorial", "dashboard", "ParallelCluster"]
+++

1.	Login to your EC2 Instance.

2.	list the available disks using the following command

```bash
lsblk
```
	
The output will list the disk you attached to your instance.

**NOTE**: Newer Linux kernels may rename your devices to **/dev/xvdf** through **/dev/xvdp** internally, even when the device name entered here (and shown in the details) is /dev/sdf through /dev/sdp.

3.	Check for an unmounted filesystem of size 10GB. For example: "nvem1n1" as shown below.

```bash
NAME        MAJ:MIN RM  SIZE RO TYPE MOUNTPOINT
nvme0n1     259:0    0  150G  0 disk
└─nvme0n1p1 259:1    0  150G  0 part /
nvme1n1     259:2    0   10G  0 disk
```

4.	Check if the volume has any data using the following command.

```bash
sudo file -s /dev/nvme1n1
```

Where "nvem1n1" is the device you noted from the previous section after attaching the device to the EC2 Instance.

	If the above command output shows "/dev/nvem1n1: data", it means your volume is empty.

5.	Format the volume to the ext4 filesystem using the following command.

```bash
sudo mkfs -t ext4 /dev/nvem1n1
```

**NOTE**: This file-system formatting step is only for a new device, **DO NOT** run this step while mounting an existing volume as it will wipe out all data on the device.

6.	Create a directory of your choice to mount our new ext4 volume. Let's use the name "volume1".

```bash	
sudo mkdir /mnt/volume1
```


7.	Mount the volume to “volume1” directory using the following command.

```bash
sudo mount /dev/nvem1n1 /mnt/volume1
```


8.	cd into volume1 directory and check the disk space for confirming the volume mount.

```bash
cd /mnt/volume1
df -h .
```
The above command would show the free space in the volume1 directory.

**NOTE**: At this point, the drive is owned by root and not user. If you would like to edit files in this directory you'll need to change ownership with **chown**. We'll go over this later. 

9.	To unmount the volume, you have to use the following command. Make sure to be outside the directory to unmount the volume.

```bash
sudo umount /dev/nvem1n1
```
