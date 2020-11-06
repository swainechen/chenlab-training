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

3.	Check if the volume has any data using the following command.

```bash
sudo file -s /dev/xvdf
```

Where "xvdf" is the device you noted from the previous section after attaching the device to the EC2 Instance.

	If the above command output shows "/dev/xvdf: data", it means your volume is empty.

4.	Format the volume to the ext4 filesystem using the following command.

```bash
sudo mkfs -t ext4 /dev/xvdf
```

**NOTE**: This file-system formatting step is only for a new device, **DO NOT** run this step while mounting an existing volume as it will wipe out all data on the device.

5.	Create a directory of your choice to mount our new ext4 volume. I am using the name "newvolume"

```bash	
sudo mkdir /mnt/volume1
```


6.	Mount the volume to “newvolume” directory using the following command.

```bash
sudo mount /dev/xvdf /mnt/volume1
```


7.	cd into newvolume directory and check the disk space for confirming the volume mount.

```bash
cd /volume1
df -h .
```

8.	The above command would show the free space in the newvolume directory.

9.	To unmount the volume, you have to use the following command.

```bash
umount /dev/xvdf
```
