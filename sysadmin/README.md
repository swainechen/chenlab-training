# Basic Sysadmin
This contains the basics of setting up a "clean" install that you might launch on AWS. The base AMI images (for RedHat or Ubuntu, for example) come quite bare - which is a good thing! Doing updates and getting the basic software installed can be daunting when you start out, but the ability to customize your workstation can open up lots of new possibilities. It will also help you upgrade and troubleshoot issues when (not if!) they come up.

## Ubuntu
All v2.0x versions of the CHENLAB-PUBLIC are built on:
`Ubuntu Server 20.04 LTS (HVM), SSD Volume Type - ami-0d058fe428540cd89 (64-bit x86)`

All v2.1x versions of the CHENLAB-PUBLIC AMI are built on:
`Ubuntu Server 22.04 LTS (HVM), SSD Volume Type - ami-04d9e855d716f9c99 (64-bit (x86))`

The default username is `ubuntu`. Many things can be run with just 1-2 CPUs and 2-4 GB of RAM (i.e. a `t3a.small` or `t3a.medium`).
However, if you want to do assembly or metagenomic classification (such as with Kraken2), you need more RAM (16-32 GB), and I often use an `r5a.large` or `r5a.xlarge`.

Setting up a system is split into the following steps (that progressively get more and more specific / bespoke:

* [General system updates](system.md)
* [R and Bioconductor](R-Bioconductor.md)
* [Bioinformatics software](bioinformatics.md)

## Cleaning an AMI for release
Some tasks need to be done to clean up an AMI for public use.
This isn't so useful for most users but is here just in case.
* [AMI Cleaning](AMI-cleaning.md)
