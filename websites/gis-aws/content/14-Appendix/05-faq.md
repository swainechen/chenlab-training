+++
title = "f. FAQ"
date = 2021-05-23T16:24:30-04:00
draft = false 
weight = 80
tags = ["help", "FAQ"]
+++

## FAQ

**Which amount of CPU and RAM should I choose for my machine? Which EC2 instance do I choose?**
Of course, that depends on your application but just for a simple use case such as installing packages and moving small files even a small instance from _T2_ family of 2-8 GB RAM and 1-2 vCPUs might be enough. If you are looking to run memory-intensive applications (e.g. processing big amounts of data), look into memory-optimized _R5_ or _R5a_ instances, where 16-64 GB RAM and 2-8 vCPUs are reasonable to start with. If you are running applications that are designed to utilize multiple cores, compute-optimized _C5_ instances offer a combination of higher number vCPUs vs RAM, with options like 16 vCPUs per 32 GB of memory.

**What is an SSH client? Which one should I use to access my machine?** SSH client is a tool that allows you to connect to the remote machine (such as an initiated AWS EC2 instance - a computer in the cloud) from your local computer. On Mac and Linux users, ssh clients should be already preinstalled as part of core utilities and accessible via the terminal with the command `ssh`. On Windows, [Putty](https://www.putty.org/) is a popular choice, however, if you want the full functionality of  a UNIX-like command line, newer Windows versions come with the ability to install [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10), which emulates a Linux machine using your Windows resources.

**Which of these services are availability zone or region limited?**  

**What is the current pricing list?**  
Keep it updated, but date it.  


