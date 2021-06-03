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
**NOTE**: Newer Linux kernels may rename your devices to **/dev/xvdf** through **/dev/xvdp** internally, even when the device name entered here (and shown in the details) is /dev/sdf through /dev/sdp. To identify your disk for mounting, you can also check the size and whether there is an existing mount point. 

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

10.	You will need to set ownership of this drive from root to user with **chown**.

```bash
sudo chown -R user /mnt/volume1
```
