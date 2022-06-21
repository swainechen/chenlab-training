# Bioinformatics software
Bioinformatics is moving very fast. It's hard to keep up with version changes, and because of some details of academic funding and training, there is a lot of variation in organization and design.

From a setup point of view, there are two categories of software out there. Some software that can be installed in a relatively "standard" Unix-like manner, while others require some custom updates or custom configuration (sometimes for compilation, nonstandard paths, or for needing special libraries / data to be downloaded).

Generally the strategy is to keep all the downloads and sources in `/usr/local/src` and then link into `/usr/local/bin` or other appropriate directories in `/usr/local`.
Also, software tends to go through a lifecycle where there are lots of updates, then eventually it stabilizes and/or becomes less updated for whatever reason.
The standard Ubuntu LTS (Focal, 20.04 in this case) repositories don't keep up with the latest software - they are made to support system stability, and the latest versions are incorporated into the interim point releases (like 20.10, 21.04, etc).
So for some of these software packages, you have a choice between installing from the Ubuntu repositories and installing from a binary or source release.

Some considerations for installing from Ubuntu - you'll get automatic updates with your system, including security updates, the setup is about as pain-free as it can be, uninstallation is clean and easy, and it's easier to track versions for reproducibility. However, you may get stuck at an old version. For some software, that isn't a big deal, because it's not changing much; for others, it can mean that you simply aren't able to do things or some bugs won't have been fixed.

Installing from a binary or source package has the opposite characteristics. You can ensure you're using the latest version, including the latest updates to the "live" codebase, but you may have to do adaptation of the installation to your system, compiling can sometimes require troubleshooting, uninstallation can be difficult to do cleanly, and you may not always know how to mark which version you're using (especially if you are routinely updating to the current git source tree instead of using a release version).

I've made some recommendations for the software below based on these considerations. There are no hard and fast rules, and for most I provide instructions for multiple ways to install these packages. One thing is I've tried to avoid most of the python managers like conda, which can lead to lots of complex dependencies as more software gets installed - this should therefore help keep it relatively "clean" and simple for you to then install other software through that route!

### Another reminder about using `sudo` and the root account:
The commands below are written to do most things as root - this is not the best practice. The [standard recommendation](https://tldp.org/HOWTO/Software-Building-HOWTO-3.html) is to only use sudo or the root account when it's absolutely needed. There are a couple ways to do this:
- download into `/home/ubuntu/src`, build there, then sudo install
- make a `staff` or `admin` or `src` group, add the `ubuntu` user to that group, `chown -R root:staff /usr/local/src` (or even the whole `/usr/local` tree), then do all the builds as ubuntu +/- the install as root

Here, though, we're doing the simple version to just run around as root. This takes away a layer of potential problems and if you're on a fresh cloud machine, it's not a huge deal to mess something up since there really shouldn't be other users you're concerned about and you can just make a new instance.


Note that the version numbers are included in the commands below - you'll have to update those to whatever the current version is at the time.

## Overview
### General
* [ASCP](#ASCP)
* [BLAST+](#BLAST)
* [genome-tools](#genome-tools)
* [HYPHY](#HYPHY)
* [Kingfisher](#Kingfisher)
* [meme](#meme)
* [MUMmer](#MUMmer)
* [pbbam](#pbbam)
* [samtools](#samtools-including-htslib)
* [sra-tools](#sra-tools)
* [SLC Closet](#SLC-Closet)

### Read processing
* [bowtie](#bowtie)
* [bwa](#bwa)
* [deML](#deML)
* [fastp](#fastp)
* [lacer](#lacer)
* [minimap2](#minimap2)
* [Porechop](#Porechop)
* [poretools](#poretools)
* [seqmagick](#seqmagick)
* [seqtk](#seqtk)
* [STAR](#STAR)
* [Trimmomatic](#Trimmomatic)

### Assembly
* [A5 assembler](#A5-assembler)
* [canu](#canu)
* [Flye](#Flye)
* [miniasm](#miniasm)
* [raven](#raven)
* [redbean](#redbean)
* [SGA](#SGA)
* [skesa](#skesa)
* [SPAdes](#SPAdes)
* [velvet](#velvet)
* [Trycycler](#Trycycler)
* [Unicycler](#Unicycler)
* [FinIS](#FinIS)
* [OPERA-LG](#OPERA-LG)

### Post-processing, variant calling
* [BLASR](#BLASR)
* [breseq](#breseq)
* [GATK](#GATK)
* [GraphMap2](#GraphMap2)
* [lofreq](#lofreq)
* [mash](#mash)
* [medaka](#medaka)
* [nanopolish](#nanopolish)
* [pilon](#pilon)
* [racon](#racon)
* [snippy](#snippy)

### Annotation and classification
* [abricate](#abricate)
* [EzClermont](#EzClermont)
* [GBS-SBG](#GBS-SBG)
* [Kleborate](#Kleborate)
* [Kraken 2](#Kraken-2)
* [prokka](#prokka)
* [referenceseeker](#referenceseeker)
* [Roary](#Roary)
* [SeqSero](#SeqSero)
* [SRST2](#SRST2)

### Visualization
* [bandage](#bandage)
* [BRIG](#BRIG)
* [Circos](#Circos)
* [EasyFig](#EasyFig)
* [SeqFindr](#SeqFindr)
* [slcview](#slcview)

## Detailed instructions

### General
* [ASCP](#ASCP)
* [BLAST+](#BLAST)
* [genome-tools](#genome-tools)
* [HYPHY](#HYPHY)
* [Kingfisher](#Kingfisher)
* [meme](#meme)
* [MUMmer](#MUMmer)
* [pbbam](#pbbam)
* [samtools](#samtools-including-htslib)
* [sra-tools](#sra-tools)
* [SLC Closet](#SLC-Closet)

#### [ASCP](http://downloads.asperasoft.com/connect2/)
This is useful for fast downloads, such as from ENA or Genbank. This seems to be best (and easiest) to install as the user (ubuntu).
The download page is at [https://www.ibm.com/aspera/connect/](https://www.ibm.com/aspera/connect/), use this to update the link and version below.
```
# make sure you're in a user (ubuntu) shell and not root
cd /home/ubuntu
wget https://d3gcli72yxqn2z.cloudfront.net/downloads/connect/latest/bin/ibm-aspera-connect_4.2.0.42_linux.tar.gz
tar xvzf ibm-aspera-connect_4.2.0.42_linux.tar.gz
./ibm-aspera-connect_4.2.0.42_linux.sh
# update your path to include this (for ex. in .bashrc)
echo '# path for ascp' >> /home/ubuntu/.bashrc
echo 'export PATH=$PATH:/home/ubuntu/.aspera/connect/bin' >> /home/ubuntu/.bashrc
```
Note that, if needed, the standard key required is at `/home/ubuntu/.aspera/connect/etc/asperaweb_id_dsa.openssh`.

#### [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
Ubuntu LTS version: 2.12.0<br/>
Online version: 2.13.0-1

The Ubuntu repositories have this, but it's at version 2.9.0-2. The current version right now is at 2.13.0-1. Generally for most users, BLAST has been pretty stable, but here's how to update it if you need the latest version (you'll have to repeat these for each update as well after you do this manual install).

The latest version can be found at [https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) - you'll have to update the links and version numbers appropriately below.

First check on what we already have and where it is
```
# check the current version
blastn -h
# this should be 2.9.0+ if it's from the Ubuntu repositories
# this will likely be 2.13.0+ or higher if you manually updated it as below

# check where it's installed
which blastn
# this will be /usr/bin/blastn for the version from the Ubuntu repositories
# this should be /usr/local/bin/blastn for a manually installed/updated version

# if you need to, make sure you know the order of directories in your path
env | grep PATH
# if you haven't changed this manually, /usr/local/bin should come before /usr/bin, which is why we can do a manual update, link it into /usr/local/bin, and that will take precedence over the version from the Ubuntu repositories
```

Install the latest version:
```
sudo su -
cd /usr/local/src
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz
tar xvzf ncbi-blast-2.13.0+-x64-linux.tar.gz
cd ncbi-blast-2.13.0+
# always good to look at the README, though not much in this one
less README
cd /usr/local/src/ncbi-blast-2.13.0+/bin
for i in *; do ln -s /usr/local/src/ncbi-blast-2.13.0+/bin/$i /usr/local/bin; done
```

You can verify that the default version is what you expect with the commands in the first block above. You may also want to remove the older Ubuntu version to avoid any potential confusion: `sudo apt purge ncbi-blast+ ncbi-blast+-legacy blast2`.

#### [genome-tools](https://github.com/swainechen/genome-tools)
This is the first open source project I wrote - it was a set of scripts to help manage the data files for multiple genomes.
Some of the file formats have changed but I still find this useful, particularly if you have just a few individual genomes you use frequently (for me, there are just a few E. coli strains I use heavily).
There are also a few useful scripts included.
The main setup is pretty standard, but then you need to add data.
```
sudo su -
cd /usr/local/src
wget https://github.com/swainechen/genome-tools/archive/refs/tags/1.3.tar.gz
tar xvzf 1.3.tar.gz
cd genome-tools-1.3
# check the docs
less INSTALL
./setup.pl
# accept all defaults, and copy Orgmap.pm to /usr/local/lib/site_perl

# download some genomes and make them available
```

#### [HYPHY](https://github.com/veg/hyphy)
Ubuntu LTS version: 2.5.36<br/>
Latest GitHub version: 2.5.39

I've generally avoided conda when the setup is otherwise simple, so that it can be used for other software later without much trouble. This setup is intended to be used for a single person or to drive a specific pipeline. Therefore, we'll compile HYPHY from source using the latest release tarball.
```
sudo su -
cd /usr/local/src
wget https://github.com/veg/hyphy/archive/refs/tags/2.5.39.tar.gz
tar xvzf 2.5.39.tar.gz
cd hyphy-2.5.39
# check the instructions
less README.md

# then do the installation. We have openMPI installed from the initial system setup
cmake .
make HYPHYMPI
# this defaults to install to /usr/local
make install
```

#### [Kingfisher](https://github.com/wwood/kingfisher-download)
**Note: this seems to be changing quite a bit in terms of install and syntax, as is the organization of the public repositories. Other options include [ENA's own tool](https://github.com/enasequence/ena-ftp-downloader) or prefetch + fasterq-dump from NCBI's [SRA Tools](#sra-tools). A good place to check for other options is https://www.biostars.org/p/9493498/**

This is a convenient tool for downloading public data sets, such as from GenBank or ENA.
Again, the setup is pretty simple so I've done it without conda.
Fortunately, this only requires one additional library to install (extern).

```
sudo su -
cd /usr/local/src
git clone https://github.com/wwood/kingfisher-download
cd kingfisher-download
# check the docs
less README.md
pip3 install bird_tool_utils'>='0.2.17
ln -s /usr/local/src/kingfisher-download/bin/kingfisher /usr/local/bin
ln -s /usr/local/src/kingfisher-download/ena-fast-download.py /usr/local/bin
```
If you've installed [ASCP](#ASCP) as described in this guide (or you're using the CHENLAB-PUBLIC AMI), then the `-m ena-ascp` method will just work, as will the more standard methods.

#### [meme](http://meme-suite.org/meme/)
Current version is 5.4.1.
```
sudo su -
cd /usr/local/src
wget https://meme-suite.org/meme/meme-software/5.4.1/meme-5.4.1.tar.gz
tar xvzf meme-5.4.1.tar.gz
cd meme-5.4.1

# always good to check the docs
less README
less INSTALL

# this one we have to explicitly set the prefix
./configure --prefix=/usr/local --enable-build-libxml2 --enable-build-libxslt
# check everything looks ok for the output from the configure script, then:
make
# this has a test suite which is good to do
make test
# some meme?? tests fail if you run as root
# they also fail if you don't have enough slots for the MPI version, i.e. CPUs
make install

# then, as per instructions, add things to your path (for the user, not necessarily root)
echo '# path for meme' >> /home/ubuntu/.bashrc
echo 'export PATH=$PATH:/usr/local/libexec/meme-5.4.1' >> /home/ubuntu/.bashrc
```

#### [MUMmer](https://github.com/mummer4/mummer)
Ubuntu LTS version: 3.23<br/>
Latest GitHub version: 4.0.0 release candidate 1

MUMmer is a venerable suite of programs for aligning DNA and protein sequences. It has recently seen a new incarnation with version 4. We'll therefore install the latest GitHub release, since the Ubuntu LTS repositories have the previous version 3 generation.

```
sudo su -
cd /usr/local/src
wget https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz
tar xvzf mummer-4.0.0rc1.tar.gz
cd mummer-4.0.0rc1
# check the docs - these are on the website, particularly INSTALL.md
./configure --prefix=/usr/local
make
make install

# test it - this needs LD_LIBRARY_PATH set, which should be fine if you run this as the ubuntu user
nucmer -h
```

#### [pbbam](https://github.com/PacificBiosciences/pbbam)
Ubuntu LTS version: 2.0.0<br/>
GitHub version: 2.1.0

This is a set of tools for handling the specific bam files that PacBio used to use.

The lastest version from the Github repository has some issues with compiling with new gcc versions, so we'll install from the repositories.

Install from the Ubuntu repositories (**Recommended**):
```
sudo apt install pbbamtools
```

Install from the GitHub release (May 2022 - doesn't work with :
```
sudo su -
cd /usr/local/src
wget https://github.com/PacificBiosciences/pbbam/archive/refs/tags/v2.1.0.tar.gz
tar xvzf v2.1.0.tar.gz
cd pbbam-2.1.0
# check the docs
less README.md
less INSTALL.md

# build and test - use meson and ninja
mkdir build
cd build
# there's a minor change in the option for tests (it's -Dtests, not -Denable-tests)
meson --prefix /usr/local/ -Dtests=true ..
ninja
ninja test
# there are some git errors because we used a release tarball, but everything else works
ninja install
# this installs some shared libraries in /usr/local/lib (among other things in /usr/local), so we need to update LD_LIBRARY_PATH
echo "# For pbbam" >> /home/ubuntu/.bashrc
echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/x86_64-linux-gnu" >> /home/ubuntu/.bashrc

# test it
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/x86_64-linux-gnu
pbmerge
pbindex
```

#### [samtools](http://www.htslib.org/) (including htslib)
Ubuntu LTS version: 1.13-4<br/>
GitHub version: 1.15.1

The HTSlib / samtools suite provides core tools and libraries that are also used by many other bioinformatics software tools. Again there is a version in the Ubuntu repositories (1.10), but the current version as of this writing is 1.15. This may or may not matter for what you're doing. If it does, here's how to update it.

Figuring out what version you have currently:
```
dpkg --list | grep samtools
which samtools
# this generally will be /usr/bin/samtools if it's from the Ubuntu packages
# this generally will be /usr/local/bin/samtools if you manually installed it
samtools
```

Installing the latest version (check the version numbers throughout). We'll install all the packages including the development htslib libraries that other software might need. First, the regular samtools utilities, which will pull in the matching htslib:
```
sudo su -
cd /usr/local/src
wget https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2
tar xvjf samtools-1.15.1.tar.bz2
cd samtools-1.15.1
# always check the basic docs
less README
less INSTALL
# do the main configuration and installation
./configure --enable-configure-htslib
make
# the default prefix is to install to /usr/local already
make install

# htslib is configured - go make and install
cd /usr/local/src/samtools-1.15.1/htslib-1.15.1
make
make install
```

Then bcftools:
```
sudo su -
cd /usr/local/src
wget https://github.com/samtools/bcftools/releases/download/1.15.1/bcftools-1.15.1.tar.bz2
tar xvjf bcftools-1.15.1.tar.bz2
cd bcftools-1.15.1
# always check the basic docs
less README
less INSTALL
# do the main configuration and installation
./configure
make
# the default prefix is to install to /usr/local already
make install
```

Again, you can go back above to verify what version you're now using by default (you may need to do a `hash -r` to have bash update program locations). Also, you can remove the version from the Ubuntu repositories if you want to avoid potential confusion: `sudo apt purge samtools libhts3 bcftools`. But these may be pulled in by other packages leading to some complex dependencies, so this step isn't really necessary.

Note a couple useful pieces of information for the hts libraries:
```
/usr/local/lib/libhts.a
/usr/local/include/htslib
```

#### [sra-tools](https://github.com/ncbi/sra-tools)
Ubuntu LTS version: 2.11.3<br/>
GitHub version: 3.0.0

These are useful tools for using data in the GenBank Sequence Read Archives.
They provide precompiled binaries on their GitHub page. We'll use these since they're a bit newer than what's in the Ubuntu repositories.

```
sudo su -
cd /usr/local/src
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-ubuntu64.tar.gz
tar xvzf sratoolkit.3.0.0-ubuntu64.tar.gz 
cd sratoolkit.3.0.0-ubuntu64
# check the docs
less README.md
less README-blastn
less README-vdb-config

# link in the binaries - only take the parents, not the version-named ones
cd bin
for i in `ls -1 | grep -v '\.3$' | grep -v '3.0.0$'`; do \
  ln -s /usr/local/src/sratoolkit.3.0.0-ubuntu64/bin/$i /usr/local/bin; 
done
# this needs to be configured once (both root and ubuntu), even if accepting the defaults
vdb-config -i
# type 'x' to just exit with default configuration
# now go back to being the ubuntu user, probably exit from the su shell
vdb-config -i
# type 'x' to exit with default configuration
# this should leave a file in /home/ubuntu/.ncbi/user-settings.mkfg
```

#### [SLC Closet](#SLC-Closet)
These are a lot of old utility scripts I've accumulated over the years.
This doesn't mean they're all useful or nice though!

A few useful biology-related utilities:
* `slc-blastn` (and blastp, blastx, etc. versions) - this can work directly with the `.tgz` files output by `SLC-wgs.pl`
* `blast-mlst.pl` - will pick up the SRST2 MLST databases if you follow the [instructions below](#SRST2)
* `blast-vf.pl` - will pick up the SRST2 gene databases if you follow the [instructions below](#SRST2)
* `SLC-wgs.pl` - clean wrapper for `velvet`, `SPAdes`, and `SKESA` assembly, will do a `prokka` annotation and put everything into a `.tgz` file
* `fasta-length.pl`

Some non-biology related tools I use all the time:
* `mean`
* `ncols`
* `randomize.pl`
* `switchcols`
```
sudo su -
cd /usr/local/src
git clone https://github.com/swainechen/closet
cd closet/bin
for i in *; do ln -s /usr/local/src/closet/bin/$i /usr/local/bin; done
ln -s /usr/local/src/closet/lib/slchen.pm /usr/local/lib/site_perl
ln -s /usr/local/src/closet/lib/Slctree.pm /usr/local/lib/site_perl

# some convenient links
cd /usr/local/bin
for i in blastp blastx tblastn tblastx; do ln -s slc-blastn slc-$i; done
```

### Read processing
* [bowtie](#bowtie)
* [bwa](#bwa)
* [deML](#deML)
* [fastp](#fastp)
* [lacer](#lacer)
* [minimap2](#minimap2)
* [Porechop](#Porechop)
* [poretools](#poretools)
* [seqmagick](#seqmagick)
* [seqtk](#seqtk)
* [STAR](#STAR)
* [Trimmomatic](#Trimmomatic)

#### [bowtie2](https://github.com/BenLangmead/bowtie2)
Ubuntu LTS version: 2.4.4<br/>
GitHub version: 2.4.5

This is one of the well known short read mappers.
[SRST2](#SRST2) uses this but has some version requirements (2.2.9).
The Ubuntu focal LTS repositories include bowtie2 at version 2.4.4.
The latest (April 2022) online at GitHub is version 2.4.5.

We'll install the latest release from GitHub as the default. Installing the older 2.2.9 will be dealt with in the [SRST2](#SRST2) section.

Install from the Ubuntu repositories:
```
sudo apt install bowtie2

# check where it is and the version - this should be /usr/bin/bowtie2
which bowtie2
bowtie2 --version
```

Install from the GitHub release (**Recommended**):
```
sudo su -
cd /usr/local/src
wget https://github.com/BenLangmead/bowtie2/releases/download/v2.4.5/bowtie2-2.4.5-linux-x86_64.zip
unzip bowtie2-2.4.5-linux-x86_64.zip
cd bowtie2-2.4.5-linux-x86_64

# always check the docs
less README.md
for i in bowtie*; do ln -s /usr/local/src/bowtie2-2.4.5-linux-x86_64/$i /usr/local/bin; done

# check where it is and the version - this should be /usr/local/bin/bowtie2
which bowtie2
bowtie2 --version
```

#### [bwa](https://github.com/lh3/bwa)
This isn't being updated so much so I prefer to use the Ubuntu repositories, which have the same main version as the current release (as of May 2021) on the GitHub repository (0.7.17).

From the Ubuntu repositories (**Recommended**):
```
sudo apt install bwa

# check where it is and the version
which bwa
bwa
```

From the GitHub release tarball:
```
sudo su -
cd /usr/local/src
wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
tar xvzf bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make
ln -s /usr/local/src/bwa/bwa /usr/local/bin

# check where it is and the version
which bwa
bwa
```

From the GitHub repository:
```
sudo su -
cd /usr/local/src
git clone https://github.com/lh3/bwa.git
cd bwa
make
ln -s /usr/local/src/bwa/bwa /usr/local/bin

# check where it is and the version
which bwa
bwa
```

#### [deML](https://github.com/grenaud/deML)
GitHub release version: 1.1.3

This is a maximum likelihood demultiplexer that is useful for when you have designed a custom multiplexing strategy. This happens a lot with Tn-seq experiments, for example.

There seem to be some updates since the last release, and the 1.1.3 release tarball doesn't seem to compile, seems to be some issues with changing paths. So we'll install from a clone of the current GitHub repository.

On newer systems, this requires an older compiler (gcc-9 and g++-9 for Ubuntu). I won't go into the details here, but the short version is to install gcc-9 and g++-9 with apt, change the `/usr/bin/gcc` and `/usr/bin/g++` links temporarily to point to those, then do the compile. Remember to revert the `/usr/bin` links back afterwards to use the up-to-date compiler for other things.
```
sudo su -
cd /usr/local/src
git clone https://github.com/grenaud/deML.git
cd deML
# check the docs
less README.md

# build and install
make
ln -s /usr/local/src/deML/src/deML /usr/local/bin

# check the binary runs ok
deML
```

#### [fastp](https://github.com/OpenGene/fastp)
Ubuntu LTS version: 0.20.1<br/>
GitHub version: 0.23.2

This is a utility for preprocessing fastq files. We'll install the GitHub version. (The Ubuntu version of course is straightforward with an `apt install fastp`).


```
sudo su -
cd /usr/local/src
wget https://github.com/OpenGene/fastp/archive/refs/tags/v0.23.2.tar.gz
tar xvzf v0.23.2.tar.gz
cd fastp-0.23.2/
# check the docs
less README.md

# Recent versions seem to need deflate and isal libraries
apt install libdeflate-dev libisal-dev
make
ln -s /usr/local/src/fastp-0.23.2/fastp /usr/local/bin

# check the binary runs ok
fastp
```

#### [lacer](https://github.com/swainechen/lacer)
GitHub version: 0.426

This is a base quality score recalibrator (the "Q" in FASTQ files).
It does not require knowledge of common SNPs and therefore is the only program that can generally recalibrate base quality scores on any organism (programs like GATK's BaseRecalibrator were designed primarily for human sequencing data, for example, and therefore only work well on data from a limited set of organisms).
This means that base quality recalibration, leading to more accurate trimming and SNP calling, can now be used on any organism (specifically including bacteria).

The main program is just a perl script. There's a fast and stripped-down `lacepr` program that does the basics of what GATK's PrintReads does. It's faster, is a bit more flexible to use to manage read groups, and also can recalibrate fastq files directly (so you can map "enough" to generate a recalibration, then recalibrate the fastq file, then do a full mapping with recalibrated data).

There are some Perl dependencies - if you've done everything in the [initial system setup](system.md), all of these should be available already. Of course, all of these are already set up on the CHENLAB-PUBLIC AMI.
```
sudo su -
cd /usr/local/src
git clone https://github.com/swainechen/lacer.git
cd lacer
# check the docs
less README.md

ln -s /usr/local/src/lacer/lacer.pl /usr/local/bin
```

To get lacepr to compile and install, we need an older version of samtools (up to 1.9). There were changes starting in samtools/htslib version 1.10 that haven't been incorporated into lacepr yet.
```
sudo su -
cd /usr/local/src/lacer/lacepr
# compile and install lacepr
# as is, this works for samtools 1.6 - 1.9
# check the docs first
less INSTALL

# get and compile samtools 1.9
wget https://sourceforge.net/projects/samtools/files/samtools/1.9/samtools-1.9.tar.bz2/download -O samtools-1.9.tar.bz2
tar xvjf samtools-1.9.tar.bz2
cd samtools-1.9
# this compile also require switching to an older gcc-9, as for deML above
make

# for the actual lacepr compile, you can switch back to gcc-11
cd /usr/local/src/lacer/lacepr
SAMTOOLS=/usr/local/src/lacer/lacepr/samtools-1.9
HTSLIB=/usr/local/src/lacer/lacepr/samtools-1.9/htslib-1.9
gcc -I$SAMTOOLS -I$HTSLIB lacepr.c -L$SAMTOOLS -L$HTSLIB -lbam -l:libhts.a -lz -lpthread -lm -llzma -lbz2 -ldeflate -lcurl -lcrypto -o lacepr
ln -s /usr/local/src/lacer/lacepr/lacepr /usr/local/bin
```

#### [minimap2](https://github.com/lh3/minimap2)
Ubuntu LTS version: 2.24<br/>
GitHub version: 2.24

The Ubuntu Focal (20.04) LTS repositories have version 2.17.
The current release version (May 2021) in GitHub is 2.24.
As this is under active development, we'll have to keep up with new releases (similar instructions to below).

From the Ubuntu repositories:
```
sudo apt install minimap2

# check where it is and the version
which minimap2
minimap2
```

From the GitHub release tarball (**Recommended**):
```
sudo su -
cd /usr/local/src
wget https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2
tar xvjf minimap2-2.24_x64-linux.tar.bz2
ln -s /usr/local/src/minimap2-2.24_x64-linux/minimap2 /usr/local/bin
ln -s /usr/local/src/minimap2-2.24_x64-linux/minimap2.1 /usr/local/man/man1/

# check where it is and the version
which minimap2
minimap2
```

#### [Porechop](https://github.com/rrwick/Porechop)
Ubuntu LTS version: 2.4<br/>
GitHub version: 2.4

This is a tool for adapter trimming for Oxford Nanopore sequencing data.
It seems to work, but officially is no longer supported.
Since it hasn't been updated in a while, the Ubuntu version is up-to-date with what's on GitHub, so use that version.

Install from the Ubuntu Focal LTS repositories (**Recommended**):
```
sudo apt install porechop

# check it runs ok
porechop --version
```

Install from GitHub:
```
sudo su -
cd /usr/local/src
wget https://github.com/rrwick/Porechop/archive/refs/tags/v0.2.4.tar.gz
tar xvzf v0.2.4.tar.gz
cd Porechop-0.2.4
# check the docs
less README.md

# the install script defaults to /usr/local already
python3 setup.py install

# check it runs ok
porechop --version
```

#### [poretools](https://github.com/arq5x/poretools)
These are utilities to work with Oxford Nanopore sequencing data.
This hasn't been updated in a while, so the Ubuntu repositories have the same version as the latest releast on GitHub.

Install from the Ubuntu Focal LTS repositories:
```
sudo apt install poretools

# check it runs ok
poretools --version
```

#### [seqmagick](https://fhcrc.github.io/seqmagick/)
Ubuntu LTS version: 0.8.4<br/>
Pip version: 0.8.4

This is a useful utility for converting between sequence formats.
This can be installed with pip (pip3 for python3), but the same version is available from the Ubuntu LTS repositories, so we'll use that.
```
sudo apt install seqmagick
seqmagick -V
```

#### [seqtk](https://github.com/lh3/seqtk)
Ubuntu LTS version: 1.3<br/>
GitHub release version: 1.3

This is another useful utility for mangling sequence files.
This seems to not be updated so frequently as well. Therefore, the Ubuntu LTS repositories have the latest version that's available as a release on the GitHub site (version 1.3). So I'd recomment using the Ubuntu packaged version.

From the Ubuntu repositories (**Recommended**):
```
sudo apt install seqtk

# check where it is and the version
which seqtk
seqtk
```

From the GitHub release tarball:
```
sudo su -
cd /usr/local/src
wget https://github.com/lh3/seqtk/archive/refs/tags/v1.3.tar.gz
tar xvzf v1.3.tar.gz
cd seqtk-1.3
# check the docs
less README.md

# make and install
make
ln -s /usr/local/src/seqtk-1.3/seqtk /usr/local/bin

# check where it is and the version
which seqtk
seqtk
```

#### [STAR](https://github.com/alexdobin/STAR)
This is a fast aligner that is splicing-aware, which is good for mapping host-pathogen combined datasets.

We'll install the latest GitHub release version.

```
sudo su -
cd /usr/local/src
wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.10a.tar.gz
tar xvzf 2.7.10a.tar.gz
cd STAR-2.7.10a
# check the docs
less README.md

# link the binary over
cd bin
ln -s /usr/local/src/STAR-2.7.10a/bin/Linux_x86_64_static/STAR* /usr/local/bin

# test it
STAR --help
```

#### [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
Ubuntu focal LTS version: 0.39<br/>
Online version: 0.39

Since the latest version is already in the Ubuntu repositories, we'll use that.
This provides bash wrappers to run `trimmomatic.jar`.
```
sudo apt install trimmomatic

# check where it is and the version
TrimmomaticPE
TrimmomaticSE

# can be useful to see what files were included - for ex. what adapter sequence files are included and where they are stored
dpkg -L trimmomatic
```

#### [BBMap](https://sourceforge.net/projects/bbmap/)
Ubuntu focal LTS version: 38.95<br/>
SourceForge version: 38.96

This is another mapper (BBMap) that has a read trimmer as well (BBDuk), along with a few other utilities.
Since the SourceForge version is newer, we'll prefer that.

From the Ubuntu repositories:
```
sudo apt install bbmap

# check where it is and the version
bbmap.sh
bbduk.sh
```

From the SourceForge release tarball:
```
sudo su -
cd /usr/local/src
# the URL again is generated
# you may need to get a new URL by clicking from https://sourceforge.net/projects/bbmap/
wget 'https://downloads.sourceforge.net/project/bbmap/BBMap_38.96.tar.gz?ts=gAAAAABiSjQCnvRNlvYnDls__gSGWZQaEiWukvpmwTVKm-aJxhDNAq7uYEjTdw-vo6EOemmH_Oy5bUL9V-DKBJesqFCn9Jp3EA%3D%3D&r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fbbmap%2Ffiles%2Flatest%2Fdownload' -O BBMap_38.96.tar.gz
tar xvzf BBMap_38.96.tar.gz
cd bbmap
# check the docs
less README.md
less docs/readme.txt

# we follow what the Ubuntu maintainers did and only link the scripts for bbduk, bbmap, bbnorm, bloomfilter, dedupe, and reformat
for i in bbduk bbmap bbnorm bloomfilter dedupe reformat; do ln -s /usr/local/src/bbmap/$i.sh /usr/local/bin; done
```

### Assembly
* [A5 assembler](#A5-assembler)
* [canu](#canu)
* [Flye](#Flye)
* [miniasm](#miniasm)
* [raven](#raven)
* [redbean](#redbean)
* [SGA](#SGA)
* [skesa](#skesa)
* [SPAdes](#SPAdes)
* [velvet](#velvet)
* [Trycycler](#Trycycler)
* [Unicycler](#Unicycler)
* [OPERA-LG](#OPERA-LG)
* [FinIS](#FinIS)
* [Contiguity](#Contiguity)

#### [A5 assembler](https://sourceforge.net/p/ngopt/wiki/A5PipelineREADME/)
This is an assembler designed for Illumina short read sequencing, stitching together a few tools.
The requirements are covered in other sections of this guide:
* [bwa](#bwa)
* [samtools](#samtools)
* [SGA](#SGA)
* [bowtie](#bowtie2)
* [Trimmomatic](#Trimmomatic)

```
sudo su -
cd /usr/local/src
# this project is on SourceForge, so you may need to regenerate the direct link below by clicking starting at https://sourceforge.net/projects/ngopt/files/latest/download
wget 'https://downloads.sourceforge.net/project/ngopt/a5_miseq_linux_20160825.tar.gz?ts=gAAAAABgnAz6WvH-pOx65oLmdQNesYbL38ifqyQIc0Cd2se1eQlL1LsJq-Et2BhXxvzP9I2t6RF9I3FDBJ4UxY_oLsq9FbXdUg%3D%3D&r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fngopt%2Ffiles%2Flatest%2Fdownload' -O a5_miseq_linux_20160825.tar.gz
tar xvzf a5_miseq_linux_20160825.tar.gz
# check the docs
less README.txt

# this figures out the real directory and will find the binaries it needs when run, so we just need to link the main perl script to /usr/local/bin
ln -s /usr/local/src/a5_miseq_linux_20160825/bin/a5_pipeline.pl /usr/local/bin

# test it
cd /usr/local/src/a5_miseq_linux_20160825
./test.a5.sh
```

#### [canu](https://github.com/marbl/canu)
Ubuntu LTS version: 2.0<br/>
GitHub version: 2.2

This is an assembler designed for noisy long-read sequencing (i.e. Oxford Nanopore and PacBio).
We'll install the GitHub version since it's newer than the one in the Ubuntu repositories.

Install from Ubuntu repositories:
```
sudo apt install canu
```

Install from GitHub (**Recommended**):
```
sudo su -
cd /usr/local/src
wget https://github.com/marbl/canu/releases/download/v2.2/canu-2.2.tar.xz
tar xvJf canu-2.2.tar.xz
cd canu-2.2
# check the docs
less README.md

# build and link the binaries
cd /usr/local/src/canu-2.2/src
make
# this will find the full path despite the symlink, where it will then find the other binaries it needs, so we only need the canu binary linked into /usr/local/bin
ln -s /usr/local/src/canu-2.2/build/bin/canu /usr/local/bin

# test it
canu
```

#### [Flye](https://github.com/fenderglass/Flye)
This is another assembler for long read sequencing (i.e. Oxford Nanopore and PacBio).

```
sudo su -
cd /usr/local/src
wget https://github.com/fenderglass/Flye/archive/refs/tags/2.9.tar.gz
tar xvzf 2.9.tar.gz
cd Flye-2.9
# check the docs
less README.md

# build and install into /usr/local - default for the setup script
python3 setup.py install

# test it
flye
python3 /usr/local/src/Flye-2.9/flye/tests/test_toy.py
```

#### [miniasm](https://github.com/lh3/miniasm)
Ubuntu LTS version: 0.3<br/>
GitHub version: 0.3

This is an assembler designed for noisy long reads (i.e. Oxford Nanopore and PacBio) with speed in mind.
This hasn't been updated in a while so we'll go with the Ubuntu repository, as it has the most current version.

Install from Ubuntu repositories (**Recommended**):
```
sudo apt install miniasm
```

Install from GitHub:
```
sudo su -
cd /usr/local/src
wget https://github.com/lh3/miniasm/archive/refs/tags/v0.3.tar.gz
tar xvzf v0.3.tar.gz
cd miniasm-0.3
# check the docs
less README.md

# make and link the binary
make
ln -s /usr/local/src/miniasm-0.3/miniasm /usr/local/bin
ln -s /usr/local/src/miniasm-0.3/minidot /usr/local/bin
```

#### [raven](https://github.com/lbcb-sci/raven)
This is another assembler for long reads (Oxford Nanopore and PacBio).

```
sudo su -
cd /usr/local/src
wget https://github.com/lbcb-sci/raven/archive/refs/tags/1.8.1.tar.gz
tar xvzf 1.8.1.tar.gz
cd raven-1.8.1
# check the docs
less README.md

# build and link the binaries
cmake -S ./ -B./build -DRAVEN_BUILD_EXE=1 -DCMAKE_BUILD_TYPE=Release
cmake --build build
ln -s /usr/local/src/raven-1.8.1/build/bin/raven /usr/local/bin

# test it
raven
```

#### [redbean](https://github.com/ruanjue/wtdbg2)
This is another assembler for long reads (Oxford Nanopore and PacBio).

Installing the latest release (2.5) from GitHub - this is just a binary package:
```
sudo su -
cd /usr/local/src
wget https://github.com/ruanjue/wtdbg2/releases/download/v2.5/wtdbg-2.5_x64_linux.tgz
tar xvzf wtdbg-2.5_x64_linux.tgz
cd wtdbg-2.5_x64_linux
# no docs to see here!
# link the binaries and test
for i in *; do ln -s /usr/local/src/wtdbg-2.5_x64_linux/$i /usr/local/bin; done
```

#### [SGA](https://github.com/jts/sga)
Ubuntu LTS version: 0.10.15<br/>
GitHub version: 0.10.15

Since this hasn't been updated in a while, and the versions are the same, we'll go with the Ubuntu version.

```
sudo apt install sga
```

#### [skesa](https://github.com/ncbi/SKESA)
Ubuntu LTS version: 2.3.0<br/>
GitHub version: 2.4.0

This is a set of assemblers (SKESA and SAUTE) with some companion programs. SKESA in particular was designed for microbial genomes.

We'll install the latest GitHub release version (tagged as "Update for NGS v.2.11.0" as of May, 2021 - but this is v2.4.0 of SKESA (the 2.11.0 is for the NGS bit)).
We'll take the source and compile.
```
sudo su -
cd /usr/local/src
wget https://github.com/ncbi/SKESA/archive/refs/tags/skesa.2.4.0_saute.1.3.0_2.tar.gz
tar xvzf skesa.2.4.0_saute.1.3.0_2.tar.gz
cd SKESA-skesa.2.4.0_saute.1.3.0_2
# check the docs
less README.md

# build - there are some assumptions for the full make where it tries to
# build ngs-sdk, even though we already have it installed, leading to a compile
# error. For now use the nongs version
make -f Makefile.nongs

# link the binaries
for i in skesa saute saute_prot gfa_connector kmercounter; do ln -s /usr/local/src/SKESA-skesa.2.4.0_saute.1.3.0_2/$i /usr/local/bin; done

# test that it works, at least for showing the help screen
skesa
saute
```

#### [SPAdes](https://cab.spbu.ru/software/spades/)
Ubuntu LTS version: 3.13.1<br/>
Online version: 3.15.4

We'll install the latest version from https://cab.spbu.ru/software/spades/. Note that since we're on a newer kernel / system / libc, we'll have to compile from source.
```
sudo su -
cd /usr/local/src
wget https://cab.spbu.ru/files/release3.15.4/SPAdes-3.15.4.tar.gz
tar xvzf SPAdes-3.15.4.tar.gz
cd SPAdes-3.15.4
# check the docs
less README.md

# compile and install
PREFIX=/usr/local ./spades_compile.sh

# check it works ok
spades.py --test
```

#### [velvet](https://github.com/dzerbino/velvet/tree/master)
This is a popular assembler.
It hasn't been updated in a while, so the Ubuntu repositories have the latest release version (1.2.10).
However, the compile parameters are important and not great in the Ubuntu packaged version (they use `MAXKMERLENGTH=31` up to `MAXKMERLENGTH=127`).
You can install the `velvet-long` package to have it compiled with the `LONGSEQUENCES` option.
However, all the binaries now have suffixes like `_127` and `_long`.
So we'll compile this from the GitHub release tarball.

From the Ubuntu repositories:
```
sudo apt install velvet velvet-long

# check where it is and the version
which velvetg
ls `which velvetg`*
velvetg
```

From the GitHub repository (**Recommended**):
```
sudo su -
cd /usr/local/src
wget https://github.com/dzerbino/velvet/archive/refs/tags/v1.2.10.tar.gz
tar xvzf v1.2.10.tar.gz
cd velvet-1.2.10
# check the docs
less README.txt
# the info on compile options is in the pdf manual, online at https://github.com/dzerbino/velvet/blob/master/Manual.pdf
# we are not specifying BIGASSEMBLY, which is generally not needed for bacterial assemblies...
make 'MAXKMERLENGTH=255' 'LONGSEQUENCES=1'
for i in velvetg velveth; do ln -s /usr/local/src/velvet-1.2.10/$i /usr/local/bin; done
# also pull in some contrib utilities, one of which needs a perl library
ln -s /usr/local/src/velvet-1.2.10/contrib/VelvetOptimiser-2.2.4/VelvetOptimiser.pl /usr/local/bin
ln -s /usr/local/src/velvet-1.2.10/contrib/shuffleSequences_fasta/shuffleSeq* /usr/local/bin
ln -s /usr/local/src/velvet-1.2.10/contrib/VelvetOptimiser-2.2.4/VelvetOpt/ /usr/local/lib/site_perl/
```
The last line for the VelvetOpt perl libraries goes to `/usr/local/lib/site_perl`.
This should work on a stock Ubuntu install.
You can check the options for this by running `perl -e 'print join ("\n", @INC), "\n";'` - this is the default search path for perl modules, and it just needs to be in one of those places (though something under `/usr/local/` is a good idea).

#### [Trycycler](https://github.com/rrwick/Trycycler/wiki)
This is a method to produce a consensus assembly. It was also designed to help with assembly of long read sequencing (Oxford Nanopore and PacBio).

Looking at the installation methods, I prefer to take the latest release to build and install. We'll use `pip3` for the installation itself.
```
sudo su -
cd /usr/local/src
wget https://github.com/rrwick/Trycycler/archive/refs/tags/v0.5.3.tar.gz
tar xvzf v0.5.3.tar.gz
cd Trycycler-0.5.3
# check the docs
less README.md
less requirements.txt

# install (defaults to /usr/local/bin for the binary)
pip3 install /usr/local/src/Trycycler-0.5.3/

# test it
trycycler --help
```

#### [Unicycler](https://github.com/rrwick/Unicycler)
Ubuntu LTS version: 0.4.8<br/>
GitHub version: 0.5.0

We'll use the latest GitHub version.

Ubuntu version:
```
sudo apt install unicycler

# test it out
unicycler
```

GitHub version (**Recommended**):
```
sudo su -
cd /usr/local/src
wget https://github.com/rrwick/Unicycler/archive/refs/tags/v0.5.0.tar.gz
tar xvzf v0.5.0.tar.gz
cd Unicycler-0.5.0
# check the docs
less README.md

# the install automatically goes to /usr/local/bin
python3 setup.py install

# test it out
unicycler
unicycler --help_all
```

#### [OPERA-LG](https://sourceforge.net/p/operasf/wiki/The%20OPERA%20wiki/)
This is a scaffolder that can be run after you do your primary assembly (with SOAP, velvet, etc.).

```
sudo su -
cd /usr/local/src
# this is on SourceForge, again you might have to get your own direct link by clicking from https://sourceforge.net/projects/operasf/files/OPERA-LG%20version%202.0.6/OPERA-LG_v2.0.6.tar.gz/download
wget 'https://downloads.sourceforge.net/project/operasf/OPERA-LG%20version%202.0.6/OPERA-LG_v2.0.6.tar.gz?ts=gAAAAABgm8W1f35V0V9XSqvQRztfQyjbS0e2qaEtdS660mhxxJ4GBG8Io1tgWoM2kOAKkS0HfeDnkQFCbnEGzx33NcJGIIJaUQ%3D%3D&r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Foperasf%2Ffiles%2FOPERA-LG%2520version%25202.0.6%2FOPERA-LG_v2.0.6.tar.gz%2Fdownload' -O OPERA-LG_v2.0.6.tar.gz
tar xvzf OPERA-LG_v2.0.6.tar.gz
cd OPERA-LG_v2.0.6
# check the docs
less README

# there are already precompiled binaries, so just link them in
cd bin
for i in *; do ln -s /usr/local/src/OPERA-LG_v2.0.6/bin/$i /usr/local/bin; done

# since we're doing this not great thing of running around as root, we need to fix permissions
chmod 755 /usr/local/src/OPERA-LG_v2.0.6
cd /usr/local/src/OPERA-LG_v2.0.6
chmod -R a+r *
find . -mindepth 1 -type d | xargs chmod 755
find . -mindepth 1 -executable -type f | xargs chmod 755

# this also has a bug in the preprocess script which is due to samtools version differences:
sed -i -e 's/-\\@ 20//' /usr/local/src/OPERA-LG_v2.0.6/bin/preprocess_reads.pl

# this can now be run as (obviously change all arguments in <>):
preprocess_reads.pl --contig <assembled contigs> --illumina-read1 <R1 fastq.gz> --illumina-read2 <R2 fastq.gz> --out preprocess.bam --map-tool bwa --samtools-dir /usr/local/src/samtools-0.1.18/
OPERA-LG <assembled contigs> preprocess.bam <output_dir> /usr/local/src/samtools-0.1.18/
```
Note that the wiki says that this depends on samtools <=0.1.19 and blasr <=1.3.1.
Instructions for keeping an old version of samtools (0.1.18) can be found in the [SRST2](#SRST2) section, and the location of the binaries should be `/usr/local/src/samtools-0.1.18`.

#### [FinIS](https://sourceforge.net/p/finis/wiki/FinIS%20wiki/)
FinIS is an assembly finisher, which generally is used after the [OPERA-LG](#OPERA-LG) scaffolder.
This is on SourceForge at v0.3, and a binary is included.
It requires [MOSEK](https://www.mosek.com/) 6 (this is currently at 9).

This really only works at this point with a velvet assembly, and it needs the `LastGraph` file from a velvet run (specify `-clean no` for `velvetg`).
```
sudo su -
cd /usr/local/src
# FinIS is on SourceForge, you may need to generate your own direct link
# by clicking starting from https://sourceforge.net/projects/finis/files/v0.3/
wget 'https://downloads.sourceforge.net/project/finis/v0.3/v0.3.tar.gz?ts=gAAAAABgl49QRj_nAX5gccQ9FeGwccNWKR_QnNWBbNTwV1HGd2cpWueqjuolVNhg8C27rHgE0kQdAd0IeYzGRkmTyTvowuvqsQ%3D%3D&r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Ffinis%2Ffiles%2Fv0.3%2Fv0.3.tar.gz%2Fdownload' -O v0.3.tar.gz
tar xvzf v0.3.tar.gz
mv v0.3 FinIS-0.3
cd FinIS-0.3
# link the binary
ln -s /usr/local/src/FinIS-0.3/finis /usr/local/bin

# download MOSEK 6 - https://www.mosek.com/downloads/6/
cd /usr/local/src
wget https://download.mosek.com/stable/6/mosektoolslinux64x86.tar.gz
tar xvzf mosektoolslinux64x86.tar.gz
# this has some binaries and libraries
# technically only the libraries are needed though for FinIS
# binaries
for i in lmgrd lmutil mampl mosek MOSEKLM moseksi mskdgopt mskexpopt mskscopt msktestlic; do ln -s /usr/local/src/mosek/6/tools/platform/linux64x86/bin/$i /usr/local/bin/; done
# libraries
for i in libiomp5.so libmosek64.so libmosek64.so.6.0 libmosekglb64.so.6.0 libmosekglbnoomp64.so.6.0 libmosekjava6_0.so libmoseknoomp64.so libmoseknoomp64.so.6.0 libmosekxx6_0.so libscopt.so; do ln -s /usr/local/src/mosek/6/tools/platform/linux64x86/bin/$i /usr/local/lib/; done
# then LD_LIBRARY_PATH needs to get set
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
echo "# for MOSEK 6" >> /home/ubuntu/.bashrc
echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib" >> /home/ubuntu/.bashrc

# we can verify this with test data that comes with the FinIS source
# there are two test datasets provided
finis test_dataset/velvet/conf.config

# the soap config file needs to drop the num_threads and mosek_runtime lines
# these seem to have changed as well - comment these two out and the test
# should run fine - note this requires more memory (~4GB) than
# a t3a.small instance has (2GB) however
finis test_dataset/soap/conf.config
```

### Post-processing, variant calling
* [BLASR](#BLASR)
* [breseq](#breseq)
* [GATK](#GATK)
* [GraphMap2](#GraphMap2)
* [lofreq](#lofreq)
* [mash](#mash)
* [medaka](#medaka)
* [nanopolish](#nanopolish)
* [pilon](#pilon)
* [racon](#racon)
* [snippy](#snippy)

#### [BLASR]
Ubuntu LTS version: 5.3.5

There used to be a Github repository for this (https://github.com/PacificBiosciences/blasr) but now it gives a 404 error. We can still install from the Ubuntu repositories.

```
sudo apt install blasr

# test it
blasr --help
```

#### [breseq](https://barricklab.org/twiki/bin/view/Lab/ToolsBacterialGenomeResequencing)
This is a pipeline for finding mutations in a short read sequencing data set relative to a reference sequence.

```
sudo su -
cd /usr/local/src
wget https://github.com/barricklab/breseq/releases/download/v0.37.0/breseq-0.37.0-Linux-x86_64.tar.gz
tar xvzf breseq-0.37.0-Linux-x86_64.tar.gz
cd breseq-0.37.0-Linux-x86_64
# check the docs
less README
less documentation/installation.html

# update owners and permissions because we're running around as root
chown -R root:root /usr/local/src/breseq-0.37.0-Linux-x86_64
cd /usr/local/src/breseq-0.37.0-Linux-x86_64
chmod -R a+r *
find . -mindepth 1 -type d | xargs chmod 755
find . -mindepth 1 -executable -type f | xargs chmod 755

# link over the binary and the share files
cd /usr/local/src/breseq-0.37.0-Linux-x86_64/bin
for i in *; do ln -s /usr/local/src/breseq-0.37.0-Linux-x86_64/bin/$i /usr/local/bin; done
ln -s /usr/local/src/breseq-0.37.0-Linux-x86_64/share/breseq /usr/local/share/

```

#### [GATK](https://gatk.broadinstitute.org/hc/en-us)
```
sudo su -
cd /usr/local/src
wget https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip
unzip gatk-4.2.6.1.zip
cd gatk-4.2.6.1
# check the docs
less README.md

# the main gatk invocation script will figure out the original directory, so we can just link this in to /usr/local/bin
ln -s /usr/local/src/gatk-4.2.6.1/gatk /usr/local/bin

# test it
gatk --help
```

#### [GraphMap2](https://github.com/lbcb-sci/graphmap2)
GitHub release version: 0.6.4
This is a mapper that was designed for Oxford Nanopore and PacBio reads.
This pulls in some other modules with git, so it's easier to clone the repository as recommended.

```
sudo su -
cd /usr/local/src
git clone https://github.com/lbcb-sci/graphmap2 
cd graphmap2
# read the docs
less README.md
less INSTALL.md

# compile
make modules
make

# link the binary
ln -s /usr/local/src/graphmap2/bin/Linux-x64/graphmap2 /usr/local/bin
```

#### [lofreq](https://csb5.github.io/lofreq/)
This is a very good variant caller with a strong theoretical basis (uses all quality information in a model-based algorithm). We'll install from the release tarball.
```
sudo su -
cd /usr/local/src
wget https://github.com/CSB5/lofreq/raw/master/dist/lofreq_star-2.1.5_linux-x86-64.tgz
tar xvzf lofreq_star-2.1.5_linux-x86-64.tgz
cd lofreq_star-2.1.5_linux-x86-64/bin
for i in *; do ln -s /usr/local/src/lofreq_star-2.1.5_linux-x86-64/bin/$i /usr/local/bin; done
```

#### [mash](https://mash.readthedocs.io/en/latest/)
Ubuntu LTS version: 2.3<br/>
GitHub release version: 2.3

Mash is a tool to estimate sequence distance using MinHash. Since the Ubuntu LTS version is the same, we'll use that.

Install from Ubuntu repositories (**Recommended**):
```
sudo apt install mash
```

Install from Github:
```
sudo su -
cd /usr/local/src
wget https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar
tar xvf mash-Linux64-v2.3.tar
cd mash-Linux64-v2.3
# check the docs online, there is just one binary
ln -s /usr/local/src/mash-Linux64-v2.3/mash /usr/local/bin
```

#### [medaka](https://github.com/nanoporetech/medaka)
** As of this writing, medaka seems to require an update to work with Python 3.0, which is what is distributed with Ubuntu 22.04 **

This is a sequence polisher for Oxford Nanopore data.
There are a few ways to install this.
The recommended build uses a virtual environment in python.
Of the choices, though, for a single user machine, I would prefer using `pip` over this.

```
# the following will automatically pull in quite a few other dependencies
sudo pip3 install medaka

# check where it went
which medaka
ls -lrt /usr/local/bin

# test it
medaka
medaka_consensus
medaka_variant
medaka_haploid_variant
```

#### [nanopolish](https://github.com/jts/nanopolish)
Ubuntu LTS version: 0.13.3<br/>
GitHub version: 0.14.0

This is a signal-level polisher for Oxford Nanopore data.
We'll use the GitHub version since it's newer.
This seems to be intended to be built from a cloned git repository, so that the source for dependencies is automatically pulled in.

Ubuntu version:
```
sudo apt install nanopolish
```

GitHub version (**Recommended**):
```
sudo su -
cd /usr/local/src
git clone --recursive https://github.com/jts/nanopolish.git
cd nanopolish
# check the docs
less README.md

# build and link the binary
make
ln -s /usr/local/src/nanopolish/nanopolish /usr/local/bin

# test it out
nanopolish
```

#### [pilon](https://github.com/broadinstitute/pilon/wiki)
Ubuntu LTS version: 1.24<br/>
GitHub version: 1.24

This is a sequence polisher that is popular for refining assemblies with short read sequencing data.
This comes as just a `jar` file.
There are differing opinions as to where to put `jar` files to keep things organized.
As the ones related to bioinformatics are usually like "programs", I typically put them in /usr/local/bin with other executables.
Another option would be to write a simple shell script that calls `java` with all the right options, including the full path to the `jar` file, but that will hide some of the details and require some translation from online documentation.
The `jar` files are usually specified with a full path anyway so it doesn't really matter that much, so long as you know where it is.

From the Ubuntu repositories (**Recommended**):
```
sudo apt install pilon
```
N.B. This Ubuntu version puts the jar file in `/usr/share/java/pilon.jar`.
There's also a script at `/usr/bin/pilon`:
```
#! /bin/sh
set -e

# export JAVA_HOME=/usr/lib/jvm/default-java
export JAVA_CMD=java

# Include the wrappers utility script
. /usr/lib/java-wrappers/java-wrappers.sh

# For memory setting see https://github.com/rrwick/Unicycler/issues/63
run_java -Xms128M -Xmx16384m -jar /usr/share/java/pilon.jar "$@"
```

From the website:
```
sudo su -
cd /usr/local/src
wget https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar
# the instructions you'll see online refer to pilon.jar, so we'll link it as such
ln -s /usr/local/src/pilon-1.24.jar /usr/local/bin/pilon.jar

# test it out
java -jar /usr/local/bin/pilon.jar
```
Unicycler is looking for a binary called `pilon`, so we can modify the Ubuntu package's script:
```
sudo su -
# you can copy and paste what's below, I've escaped the necessary characters

#### BEGIN copy-paste for pilon script ####
cat << EOF > /usr/local/bin/pilon
#! /bin/sh
set -e

# export JAVA_HOME=/usr/lib/jvm/default-java
export JAVA_CMD=java

# Include the wrappers utility script
. /usr/lib/java-wrappers/java-wrappers.sh

# For memory setting see https://github.com/rrwick/Unicycler/issues/63
run_java -Xms128M -Xmx16384m -jar /usr/local/bin/pilon.jar "\$@"
EOF
#### END copy-paste for pilon script ####

# make it executable
chmod +x /usr/local/bin/pilon

# test it
pilon
```
However, note that calling the `pilon` script asks for 16GB of heap space, which you may have to modify depending on your machine's RAM.
That said, this script was for Unicycler and that seems to need the space as per the [referenced link](https://github.com/rrwick/Unicycler/issues/63).

#### [racon](https://github.com/lbcb-sci/racon)
Ubuntu LTS version: 1.5.0<br/>
GitHub version: 1.5.0

This is a polisher / consensus module for uncorrected long reads (i.e. Oxford Nanopore and PacBio).
We'll use the GitHub version since it's newer than the one in the Ubuntu repositories.

Install from Ubuntu repositories (**Recommended**):
```
sudo apt install racon
```

Install from GitHub:
```
sudo su -
cd /usr/local/src
wget https://github.com/lbcb-sci/racon/archive/refs/tags/1.5.0.tar.gz
tar xvzf 1.5.0.tar.gz
cd racon-1.5.0
# check the docs
less README.md

# build and link the binaries
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
ln -s /usr/local/src/racon-1.5.0/build/bin/racon /usr/local/bin

# test it
racon
```

#### [snippy](https://github.com/tseemann/snippy)
This is a tool to call SNPs based on a reference sequence.
It can also output a core SNP alignment, which can then be used to make a phylogenetic tree.

```
sudo su -
cd /usr/local/src
wget https://github.com/tseemann/snippy/archive/refs/tags/v4.6.0.tar.gz
tar xvzf v4.6.0.tar.gz
cd snippy-4.6.0
# check the docs
less README.md

# link in the binaries
cd /usr/local/src/snippy-4.6.0/bin
for i in *; do ln -s /usr/local/src/snippy-4.6.0/bin/$i /usr/local/bin; done

# quick check
snippy --check
```

### Annotation and classification
* [abricate](#abricate)
* [EzClermont](#EzClermont)
* [GBS-SBG](#GBS-SBG)
* [Kleborate](#Kleborate)
* [Kraken 2](#Kraken-2)
* [prokka](#prokka)
* [referenceseeker](#referenceseeker)
* [Roary](#Roary)
* [SeqSero](#SeqSero)
* [SeqSero2](#SeqSero2)
* [SRST2](#SRST2)

#### [abricate](https://github.com/tseemann/abricate)
This is a tool to predict resistances from assemblies.

```
sudo su -

# this requires any2fasta as a dependency
cd /usr/local/src
wget https://github.com/tseemann/any2fasta/archive/refs/tags/v0.4.2.tar.gz
tar xvzf v0.4.2.tar.gz
cd any2fasta-0.4.2
# check the docs
less README.md

# link the binary
ln -s /usr/local/src/any2fasta-0.4.2/any2fasta /usr/local/bin

# now for abricate
cd /usr/local/src
wget https://github.com/tseemann/abricate/archive/refs/tags/v1.0.0.tar.gz
tar xvzf v1.0.0.tar.gz
cd abricate-1.0.0
# check the docs
less README.md

# link the binary, set it up
ln -s /usr/local/src/abricate-1.0.0/bin/abricate /usr/local/bin
abricate --check
abricate --setupdb
```

#### [EzClermont](https://github.com/nickp60/EzClermont)
This is a tool to predict phylotypes for E. coli using assemblies.

```
sudo pip3 install ezclermont
```

#### [GBS-SBG](https://github.com/swainechen/GBS-SBG)
This is a tool to predict serotypes for Group B Streptococcus (S. agalactiae) using short reads or assemblies.

```
sudo su -
cd /usr/local/src
git clone https://github.com/swainechen/GBS-SBG
cd GBS-SBG
# check the docs
less README.md

# link the binary
ln -s /usr/local/src/GBS-SBG/GBS-SBG.pl /usr/local/bin

GBS-SBG.pl -help
```

#### [Kleborate](https://github.com/katholt/Kleborate/wiki)
This is a tool that integrates multiple analyses for Klebsiella.

This is under active development, so we'll install from the latest GitHub repository.

```
sudo su -
cd /usr/local/src
git clone --recursive https://github.com/katholt/Kleborate.git
cd Kleborate/kaptive
git pull https://github.com/katholt/Kaptive master
cd ..
# this will put the main kleborate in /usr/local/bin by default
python3 setup.py install

# test it
kleborate -h
```

#### [Kraken 2](https://ccb.jhu.edu/software/kraken2/)
Ubuntu LTS version: 2.1.2<br/>
GitHub Release version: 2.1.2

This is a popular tool to do species / taxonomic classification of short reads using k-mers.

This takes quite a bit of space due to the database.
It also requires quite a bit of RAM, at least 32GB for the standard library.
Note some databases are now available as part of the [Repository of Open Data on AWS](https://aws.amazon.com/opendata/open-data-sponsorship-program/): https://registry.opendata.aws/jhu-indexes/.
We'll install the Standard 8GB database here, which has archaea, bacteria, viruses, plasmids, human, and UniVec_Core - refer to the [documentation](https://benlangmead.github.io/aws-indexes/k2) for more details or other databases.
We'll put these in /usr/local/lib/Kraken2

The Ubuntu repositories are up-to-date with the latest GitHub release, so we'll prefer to install from the repositories.

Install from Ubuntu repositories (**Recommended**):
```
sudo apt install kraken2
```

Install from Github:
```
sudo su -
cd /usr/local/src
wget https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.2.tar.gz
tar xvzf v2.1.2.tar.gz
cd kraken2-2.1.2
# check the docs
less README.md
less docs/MANUAL.markdown

# compile then link over the binaries
./install_kraken2.sh /usr/local/src/kraken2-2.1.2/bin
for i in kraken2 kraken2-build kraken2-inspect; do ln -s /usr/local/src/kraken2-2.1.2/bin/$i /usr/local/bin; done
```

Kraken needs libraries to run - we'll put these in `/usr/local/lib/Kraken2`, and they will need to be specified when you run, as noted below:
```
mkdir /usr/local/lib/Kraken2
cd /usr/local/lib/Kraken2
# install the updated standard 8GB library
mkdir /usr/local/lib/Kraken2/k2_standard_8gb_20210517
cd /usr/local/lib/Kraken2/k2_standard_8gb_20210517
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_8gb_20210517.tar.gz
tar xvzf k2_standard_8gb_20210517.tar.gz

# then this can be run as (again watch out for RAM - you probably realistically should have 16GB total)
kraken2 --db /usr/local/lib/Kraken2/k2_standard_8gb_20210517 fq_1.gz fq2.gz
```

#### [prokka](https://github.com/tseemann/prokka)
Ubuntu LTS version: 1.14.6<br/>
GitHub Release version: 1.14.5

This is a popular genome annotation tool.
This hasn't been updated in about 1-2 years according to GitHub, but somehow the Ubuntu repositories have a newer one! We'll install the Ubuntu version.

```
sudo apt install prokka

# this has another issue with licensing with Debian/Ubuntu. Check /usr/share/doc/prokka/README.Debian
# note the default directory for these is ~/.local/lib/prokka/db
# so run this as the ubuntu user, not as root
sh /usr/share/doc/prokka/get-additional-data
prokka --list
```
For reference, these libraries get stored in `/var/lib/prokka/db`
This throws an error because it assumes a numeric version number for BioPerl, but the version is ok and this still runs.

#### [referenceseeker](https://github.com/oschwengers/referenceseeker)
This program uses a kmer-based approach to find a good/close reference sequence for mapping, variant calling, etc.

We'll install the latest Github release.

```
sudo su -
cd /usr/local/src
wget https://github.com/oschwengers/referenceseeker/archive/refs/tags/v1.8.0.tar.gz
tar xvzf v1.8.0.tar.gz
cd referenceseeker-1.8.0
# check the docs
less README.md

# we need to install xopen
pip3 install xopen

# link the binaries over
cd bin
for i in *; do ln -s /usr/local/src/referenceseeker-1.8.0/bin/$i /usr/local/bin; done

# test it
referenceseeker --help
referenceseeker /usr/local/src/referenceseeker-1.8.0/test/db /usr/local/src/referenceseeker-1.8.0/test/data/Salmonella_enterica_CFSAN000189.fasta

# this requires some quite large databases. Here's an example for archaea (the smallest one
# I recommend putting them in a "standard" location like /usr/local/lib/referenceseeker
# The AMI we've provided has all the RefSeq databases installed already
cd /usr/local/lib
mkdir referenceseeker
cd referenceseeker
wget https://zenodo.org/record/4415843/files/archaea-refseq.tar.gz
tar xvzf archaea-refseq.tar.gz
```

The CHENLAB-PUBLIC AMI has all the RefSeq databases downloaded and installed already from [https://github.com/oschwengers/referenceseeker#databases](https://github.com/oschwengers/referenceseeker#databases).
These can be used by specifying one of the following as the db location (first argument) for `referenceseeker`:
* /usr/local/lib/referenceseeker/archaea-refseq
* /usr/local/lib/referenceseeker/bacteria-refseq
* /usr/local/lib/referenceseeker/fungi-refseq
* /usr/local/lib/referenceseeker/protozoa-refseq
* /usr/local/lib/referenceseeker/viral-refseq

#### [Roary](https://github.com/sanger-pathogens/Roary)
Ubuntu LTS version: 3.13.0<br/>
GitHub Release version: 3.13.0

This is a popular pan genome prediction tool.
This hasn't been updated in about 18 months, and the Ubuntu and GitHub versions are the same. So we'll install from the Ubuntu repositories.

(Note previously for Ubuntu LTS 20.04, the Ubuntu packages seemed to depend strictly on R version 3, which is why we previously preferred the GitHub release; but this seems to be resolved in 22.04.)

Install from Ubuntu repositories (**Recommended*):
```
sudo apt install roary

# check it out
roary -a
```

Install from GitHub:
```
sudo su -
cd /usr/local/src
wget https://github.com/sanger-pathogens/Roary/archive/refs/tags/v3.13.0.tar.gz
tar xvzf v3.13.0.tar.gz
cd Roary-3.13.0
# check the docs
less README.md

# install some dependencies
apt install mcl
cpanm  Array::Utils Bio::Perl Exception::Class File::Basename File::Copy File::Find::Rule File::Grep File::Path File::Slurper File::Spec File::Temp File::Which FindBin Getopt::Long Graph Graph::Writer::Dot List::Util Log::Log4perl Moose Moose::Role Text::CSV PerlIO::utf8_strict Devel::OverloadInfo Digest::MD5::File

# link the binaries and libraries
cd /usr/local/src/Roary-3.13.0/bin
for i in *; do ln -s /usr/local/src/Roary-3.13.0/bin/$i /usr/local/bin; done
ln -s /usr/local/src/Roary-3.13.0/lib/* /usr/local/lib/site_perl

# check it out
roary -a
```
This complains about cd-hit being not the right version - it asks for min 4.6, but we have 4.8.1 installed.

#### [SeqSero](https://github.com/denglab/SeqSero)
This still works well, but has been superseded by [SeqSero2](#SeqSero2).

This a serotype predictor for Salmonella.
This has some syntax that is now deprecated, and it also will probably soon need a port to Python 3.
However, the [patch here](https://github.com/denglab/SeqSero/issues/6) should allow it to run for a while still.

```
sudo su -
cd /usr/local/src
wget https://github.com/denglab/SeqSero/archive/refs/tags/v1.0.1.tar.gz
tar xvzf v1.0.1.tar.gz
cd SeqSero-1.0.1
# check the docs
less README.md
# patch the file libs/Otype_determine_analysis.py as in https://github.com/denglab/SeqSero/issues/6
wget https://github.com/denglab/SeqSero/files/6603725/SeqSero.patch.txt
patch -p1 < SeqSero.patch.txt
# as of Python 3.10 it also needs this touch up
sed -i -e 's/print /print(/' -e '/print/s/$/)/' SeqSero.py

# link the binary
ln -s /usr/local/src/SeqSero-1.0.1/SeqSero.py /usr/local/bin
```
Note that this gives a lot of `FASTA-Reader: Title ends with at least 20 valid nucleotide characters` errors from some of the blast steps due to inclusion of primer sequences in the fasta headers from isPcr.

#### [SeqSero2](https://github.com/denglab/SeqSero2)
This is an updated Salmonella serotype predictor.

```
sudo su -
cd /usr/local/src
wget https://github.com/denglab/SeqSero2/archive/refs/tags/v1.2.1.tar.gz
tar xvzf v1.2.1.tar.gz
cd SeqSero2-1.2.1
# check the docs
less README.md
# install - we're doing this as root as discussed above
python3 -m pip install .

# quick test, at least for the help screen
SeqSero2_package.py
```

#### [SRST2](https://github.com/katholt/srst2)
This is a popular short read analysis program, good for calling MLSTs, resistances, and serotypes directly from short reads.
The README on GitHub recommends installing from the git repository directly.
This requires Python 2, most obviously due to changes in the `print` syntax.
`pip` is not available for Python 2 in Ubuntu Focal, so we have to install this a bit manually.
We will also have to install the Python 2 version of scipy, which isn't in the Ubuntu repositories any more.
This also requires older samtools and bowtie2 versions.
Then there's some customization with environment variables and setting up the reference databases.
```
sudo su -
cd /usr/local/src
wget https://bootstrap.pypa.io/pip/2.7/get-pip.py
apt install python2 python2-dev
python2 get-pip.py
pip install scipy
git clone https://github.com/katholt/srst2
pip install srst2/

# install samtools 0.1.18
sudo su -
cd /usr/local/src
wget https://github.com/samtools/samtools/archive/refs/tags/0.1.18.tar.gz
tar xvzf 0.1.18.tar.gz
cd samtools-0.1.18
# check the docs
less INSTALL
make
# don't link the binaries to /usr/local/bin - this is only going to be used for SRST2

# install bowtie2 2.2.9
sudo su -
cd /usr/local/src
# this is the binaries, so no compilation needed
wget https://github.com/BenLangmead/bowtie2/releases/download/v2.2.9/bowtie2-2.2.9-linux-x86_64.zip
unzip bowtie2-2.2.9-linux-x86_64.zip
cd bowtie2-2.2.9
# check the docs
less MANUAL
# don't link the binaries to /usr/local/bin - this is only going to be used for SRST2

# environment variables
export SRST2_SAMTOOLS=/usr/local/src/samtools-0.1.18/samtools
export SRST2_BOWTIE2=/usr/local/src/bowtie2-2.2.9/bowtie2
export SRST2_BOWTIE2_BUILD=/usr/local/src/bowtie2-2.2.9/bowtie2-build

# environment variables for the ubuntu user
echo "# For SRST2" >> /home/ubuntu/.bashrc
echo "export SRST2_SAMTOOLS=/usr/local/src/samtools-0.1.18/samtools" >> /home/ubuntu/.bashrc
echo "export SRST2_BOWTIE2=/usr/local/src/bowtie2-2.2.9/bowtie2" >> /home/ubuntu/.bashrc
echo "export SRST2_BOWTIE2_BUILD=/usr/local/src/bowtie2-2.2.9/bowtie2-build" >> /home/ubuntu/.bashrc
```

SRST2 requires reference libraries.
Some are included in the source distribution.
To keep with following "standard" locations, I put these into `/usr/local/lib/SRST2`.
I also have a few utility scripts to help manage these. These are included with the [SLC Closet](#SLC-Closet) scripts.
The data files released with SRST2 we keep in a version directory.
The MLST and VFDB sequences we have to update manually, and they get their own directories.
Finally, individual organisms have their own directory, which helps to customize.
The overall structure looks like:
```
/usr/local/lib/SRST2/0.2.0
                    /MLST
                    /VFDB
                    /Ecoli
                    /Sagalactiae
                    /<other species>
```

Now, let's set it all up - remember to use the older versions of Bowtie2 and Samtools as per the SRST2 docs:
```
sudo su -
SRST2LIB=/usr/local/lib/SRST2
VERSION=0.2.0
BTBUILD=/usr/local/src/bowtie2-2.2.9/bowtie2-build
SAMTOOLS=/usr/local/src/samtools-0.1.18/samtools

# Setting up SRST2 reference libraries
mkdir -p $SRST2LIB/$VERSION
ln -s /usr/local/src/srst2/data/* $SRST2LIB/$VERSION

# MLST - use the basic SRST2 scripts - just some examples
mkdir -p $SRST2LIB/MLST
cd $SRST2LIB/MLST
# E. coli
getmlst.py --species "Escherichia coli#1" && mv profiles_csv ecoli.txt
$BTBUILD "Escherichia_coli#1.fasta" "Escherichia_coli#1.fasta"
$SAMTOOLS faidx "Escherichia_coli#1.fasta"
# S. agalactiae
getmlst.py --species "Streptococcus agalactiae" && mv profiles_csv sagalactiae.txt
$BTBUILD Streptococcus_agalactiae.fasta Streptococcus_agalactiae.fasta
$SAMTOOLS faidx Streptococcus_agalactiae.fasta
# add others as appropriate

# VFDB - follow the SRST2 directions
mkdir -p $SRST2LIB/VFDB
cd $SRST2LIB/VFDB
wget http://www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz
gunzip VFDB_setB_nt.fas.gz
# we need an old version of biopython for python 2.7, which the SRST2 scripts require
pip install biopython==1.76
python2 /usr/local/src/srst2/database_clustering/VFDBgenus.py --infile VFDB_setB_nt.fas
for i in *.fsa; do 
  cd-hit -i $i -o ${i%.*}_cdhit90 -c 0.9 > ${i%.*}_cdhit90.stdout &&
  python2 /usr/local/src/srst2/database_clustering/VFDB_cdhit_to_csv.py --cluster_file ${i%.*}_cdhit90.clstr --infile $i --outfile ${i%.*}_cdhit90.csv &&
  python2 /usr/local/src/srst2/database_clustering/csv_to_gene_db.py -t ${i%.*}_cdhit90.csv -o ${i%.*}_VF_clustered.fasta -s 5 &&
  $BTBUILD ${i%.*}_VF_clustered.fasta ${i%.*}_VF_clustered.fasta &&
  $SAMTOOLS faidx ${i%.*}_VF_clustered.fasta
done
```
For individual organisms, SRST2 can be used to do lots of things - serotyping, for example.
This requires some customization of the databases though.
I have an update script that I use, this is put into `/usr/local/lib/SRST2/update.sh`:
```
# Call this from the lower directory!
# i.e.:  cd /usr/local/lib/SRST2/Ecoli && bash ../update.sh
#
SRST2=/usr/local/lib/SRST2
VERSION=0.2.0
BOWTIEBUILD=/usr/local/src/bowtie2-2.2.9/bowtie2-build
SAMTOOLS=/usr/local/src/samtools-0.1.18/samtools
BASE=${PWD##*/}
DATE=`date +%F`
if [ $BASE == "SRST2" ]; then
  echo "Call this from the organism directory"
  echo "For example from within /usr/local/lib/SRST2/Ecoli"
  exit
else
  # this script comes from the SLC Closet respository
  combine-SRST2-fasta.pl -srst2 $SRST2 -version $VERSION > $BASE-combined-$DATE.fasta
  $BOWTIEBUILD $BASE-combined-$DATE.fasta $BASE-combined-$DATE.fasta
  $SAMTOOLS faidx $BASE-combined-$DATE.fasta
fi
```
This script will look for a file called `combine.db` in the organism directory, merge the fasta files, and put a timestamp on the combined fasta file.
All that's needed is to generate the `combine.db` file, and additional fasta files (formatted for SRST2 already) can then be added to further customize.
An example for E. coli follows.

Get the [fimH allele database](https://bitbucket.org/genomicepidemiology/fimtyper_db):
```
wget https://bitbucket.org/genomicepidemiology/fimtyper_db/raw/e11ab5129f6235ae5b40210b04bbe2f4091ba793/fimH.fsa
# fix the fasta headers
perl -i -ne 'chomp; if (/^>/) { $i++; s/>/>1__fimHType__/; s/$/__$i/; } print "$_\n";' fimH.fsa
```

File `/usr/local/lib/SRST2/Ecoli/combine.db`
```
$SRST2/$VERSION/ARGannot_r3.fasta
$SRST2/$VERSION/EcOH.fasta
$SRST2/$VERSION/Plasmid18Replicons.fasta
$SRST2/$VERSION/PlasmidFinder.fasta
$SRST2/VFDB/Citrobacter_VF_clustered.fasta
$SRST2/VFDB/Edwardsiella_VF_clustered.fasta
$SRST2/VFDB/Enterobacter_VF_clustered.fasta
$SRST2/VFDB/Erwinia_VF_clustered.fasta
$SRST2/VFDB/Escherichia_VF_clustered.fasta
$SRST2/VFDB/Klebsiella_VF_clustered.fasta
$SRST2/VFDB/Proteus_VF_clustered.fasta
$SRST2/VFDB/Salmonella_VF_clustered.fasta
$SRST2/VFDB/Serratia_VF_clustered.fasta
$SRST2/VFDB/Shigella_VF_clustered.fasta
$SRST2/VFDB/Yersinia_VF_clustered.fasta
fimH.fsa
```

Now create the database:
```
sudo su -
cd /usr/local/lib/SRST2/Ecoli
bash ../update.sh
```

Now this can be run to call resistances, E. coli seroptypes, plasmid predictions, virulence factor predictions, and fimH types (check the date of the file):
```
srst2 --input_pe strainA_1.fastq.gz strainA_2.fastq.gz --output strainA_typing --log --gene_db /usr/local/lib/SRST2/Ecoli/Ecoli-combined-2021-05-24.fasta
# output goes to:
#   strainA_typing__genes__Ecoli-combined-2021-05-24__results.txt
#   strainA_typing__fullgenes__Ecoli-combined-2021-05-24__results.txt
```

### Visualization
* [bandage](#bandage)
* [BRIG](#BRIG)
* [Circos](#Circos)
* [EasyFig](#EasyFig)
* [SeqFindr](#SeqFindr)
* [slcview](#slcview)

#### [Bandage](https://rrwick.github.io/Bandage/)
Ubuntu LTS version: 0.9.0<br/>
GitHub Release version: 0.9.0

This is a program that allows visualization of assembly graphs. This is a graphical application so you'll have to sort out forwarding X packets if you're on a remote connection (like to an AWS EC2 instance).

If you're on Linux and have a window manager, this can be as trivial as adding `-X` or `-Y` to your ssh command. This may be more complicated if you're on Windows or MacOS, though there are quite a few solutions out there.

From the Ubuntu repositories (**Recommended**):
```
sudo apt install bandage
```

From GitHub:
```
sudo su -
cd /usr/local/src
wget https://github.com/rrwick/Bandage/releases/download/v0.9.0/Bandage_Ubuntu-x86-64_v0.9.0_AppImage.zip
mkdir Bandage-0.9.0
cd Bandage-0.9.0
unzip ../Bandage_Ubuntu-x86-64_v0.9.0_AppImage.zip
ln -s /usr/local/src/Bandage-0.9.0/Bandage_Ubuntu-x86-64_v0.9.0.AppImage /usr/local/bin/Bandage
```

#### [BRIG](https://sourceforge.net/projects/brig/)
Ubuntu LTS version: 0.95<br/>
GitHub Release version: 0.95

BRIG stands for BLAST Ring Image Generator, and is useful for generating the typical bacterial genome comparison figures you see in many papers.

This hasn't been updated in a while, so the version in the Ubuntu respositories is the latest. We'll use that.

This also has a graphical interface, so you'll have to sort out forwarding X packets if you're on an AWS instance.

```
sudo apt install brig

```

#### [Circos](http://circos.ca/)
Ubuntu LTS version: 0.69-9<br/>
Website Release version: 0.69-9

This is another program for generating circular images and is quite general and widely used.

The Ubuntu repositories have the same version as the latest release so we'll use that.

```
sudo apt install circos
```

#### [EasyFig](https://mjsull.github.io/Easyfig/)
This is a program for drawing figures of genetic comparisons between two sequences.

This has a graphical interface so you'll have to sort out forwarding X packets if you're on an AWS instance.

```
sudo su -
cd /usr/local/src
wget https://github.com/mjsull/Easyfig/releases/download/2.2.2/Easyfig_2.2.2_linux.tar.gz
tar xvzf Easyfig_2.2.2_linux.tar.gz
ln -s /usr/local/src/Easyfig_2.2.2_linux/Easyfig /usr/local/bin
```

On the 20.04 LTS version of Ubuntu, this gives a lot of fontconfig warnings, but it seems to run ok (see https://askubuntu.com/questions/1098809/applications-load-slowly-on-ubuntu-18-10-seems-related-to-fontconfig)

#### [SeqFindr](http://mscook.github.io/SeqFindR/)
This is a tool to draw presence-absence plots, such as for virulence genes.

```
sudo su -
# we generally have python3, so force python2 for the install
pip2 install SeqFindr
# this gives some errors on trying to compile matplotlib due to C++ standards, but it seems to install using some legacy method and imports fine into python2 after this

# check it, at least for the help screen
SeqFindr -h
```

#### [slcview](https://github.com/swainechen/slcview)
This is a program to generate clustergrams, the traditional heatmaps coming from .cdt, .gtr, and .atr files from Mike Eisen's Cluster and Gavin Sherlock's XCluster. There is also an slctree.pl program in my [Closet repository](https://github.com/swainechen/closet) that builds on this - slctree allows scriptable and customizable tree generation, though this has generally been superseded by many good packages in R and other languages.

```
sudo su -
cd /usr/local/src
git clone https://github.com/swainechen/slcview.git
cd slcview
# check the docs
less README

for i in bon.pl slcview.pl; do ln -s /usr/local/src/slcview/$i /usr/local/bin; done
ln -s /usr/local/src/slcview/Slcview.pm /usr/local/lib/site_perl

# test it, at least for the help screen
slcview.pl
```
