+++
title = "f. Run Linux commands "
weight = 90
tags = ["tutorial", "cloud9", "s3"]
+++

Let us run a few Linux commands to install a few software packages

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
If it failed, good! It was supposed to. Now install gcc zlib-devel: 
```bash
yum install gcc zlib-devel
```
```bash
make
```
Now the command should be recognized! Try seqtk again:
```bash
./seqtk
```


