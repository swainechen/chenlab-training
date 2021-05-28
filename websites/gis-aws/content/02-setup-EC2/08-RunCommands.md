+++
title = "g. Run Linux commands "
weight = 90
tags = ["tutorial", "cloud9", "s3"]
+++

Let us run a few Linux commands to install a few software packages. Enter the following commands, and press **y** and **enter** when necessary. The instance that we've started is completely blank and does not contain any downloaded programs or software. We'll need to build it up from scratch.  

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
  
As you can see, working from an empty instance can be quite troublesome. You'll need to download absolutely everything you need. AWS offers AMIs- Amazon Machine Instances- that can act as templates for pre-built machines. Let's get ready to try that.

