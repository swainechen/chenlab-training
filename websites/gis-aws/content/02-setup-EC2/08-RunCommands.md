+++
title = "g. Run Linux commands "
weight = 90
tags = ["tutorial", "cloud9", "s3"]
+++

Let us run a few Linux commands to install a few software packages. Enter the following commands, and press **y** and **enter** when necessary. The instance that we have started is completely blank and does not contain any downloaded programs or software. We will need to build it up from scratch.  

From yum **(nmap)**: just to install a network tool

```bash
sudo yum install nmap
```

From source **(seqtk)**:
```bash
sudo su -
```
```bash
yum install git
```
```bash
cd /usr/local/src
```
```bash
git clone https://github.com/lh3/seqtk.git
```
```bash
cd seqtk/
```
```bash
make
```
Whoops! That command failed. We're still missing something. Now install gcc zlib-devel: 
```bash
yum install gcc zlib-devel
```
```bash
make
```
Now the command should be recognized! There we go. Try seqtk again:
```bash
./seqtk
```

As you can see, the default "base" images are exactly that - clean initial installations of just the operating system. These are perfect for customizing so that you have exactly the software you might need - but some setup will be required. Another option is to use an **AMI** (Amazon Machine Image) that already has some software installed. We've made our lab's AMI, which has quite a bit of genomics sofware pre-installed, available for all to use. We'll do this in the next section.

