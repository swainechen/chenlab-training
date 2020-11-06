+++
title = "f. Run Linux commands "
weight = 90
tags = ["tutorial", "cloud9", "s3"]
+++

Let us run a few Linux commands to install a few software packages

From apt **(kraken)**:

```bash
sudo yum install kraken
```

From source **(abricate)**:
```bash
sudo su -
```
```bash
cd /usr/local/src
```
```bash
git clone https://github.com/tseemann/abricate.git
```
```bash
cd abricate/
```
```bash
ln -s /usr/local/src/abricate/bin/* /usr/local/bin
```
```bash
abricate --check
```
```bash
abricate --setupdb
```


