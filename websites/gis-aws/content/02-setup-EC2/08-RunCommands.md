+++
title = "f. Run Linux commands "
weight = 90
tags = ["tutorial", "cloud9", "s3"]
+++

Let us run a few Linux commands to install a few software packages

From apt **(nmap)**: just to install a network tool

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
```bash
yum install gcc zlib-devel
```
```bash
make
```
```bash
./seqtk
```


