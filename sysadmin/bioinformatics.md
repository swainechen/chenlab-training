# Bioinformatics software
Bioinformatics is moving very fast. It's hard to keep up with version changes, and because of some details of academic funding and training, there is a lot of variation in organization and design.

From a setup point of view, there are two categories of software out there. Some software that can be installed in a relatively "standard" unix-like manner, while others require some custom updates or custom configuration (sometimes for compilation, nonstandard paths, or for needing special libraries / data to be downloaded).
* [Software that follows a standard installation](#standard-installs)
* [Software requiring some customization](#customized-installs)

Generally the strategy is to keep all the downloads and sources in /usr/local/src and then link into /usr/local/bin or other appropriate directories in /usr/local.

### Another reminder about using `sudo` and the root account:
The commands below are written to do most things as root - this is not the best practice. The [standard recommendation](https://tldp.org/HOWTO/Software-Building-HOWTO-3.html) is to only use sudo or the root account when it's absolutely needed. There are a couple ways to do this:
- download into /home/ubuntu/src, build there, then sudo install
- make a `staff` or `admin` or `src` group, add the `ubuntu` user to that group, `chown -R root:staff /usr/local/src` (or even the whole `/usr/local` tree), then do all the builds as ubuntu +/- the install as root

Here, though, we're doing the simple version to just run around as root. This takes away a layer of potential problems and if you're on a fresh cloud machine, it's not a huge deal to mess something up since there really shouldn't be other users you're concerned about and you can just make a new instance.


Note that the version numbers are included in the commands below - you'll have to update those to whatever the current version is at the time.

## Overview
### Standard installs
#### General
* [BLAST+](#BLAST)
* [HYPHY](#HYPHY)
* [Kingfisher](#Kingfisher)
* [meme](#meme)
* [samtools](#samtools-including-htslib)
* [sra-tools](#sra-tools)

#### Read processing
* [bowtie](#bowtie)
* [bwa](#bwa)
* [deML](#deML)
* [fastp](#fastp)
* [lacer](#lacer)
* [minimap](#minimap)
* [Porechop](#Porechop)
* [poretools](#poretools)
* [seqmagick](#seqmagick)
* [seqtk](#seqtk)
* [Trimmomatic](#Trimmomatic)

#### Assembly
* [a5 assembler](#a5-assembler)
* [canu](#canu)
* [miniasm](#miniasm)
* [racon](#racon)
* [sga](#sga)
* [skesa](#skesa)
* [SPAdes](#SPAdes)
* [velvet](#velvet)
* [Trycycler](#Trycycler)
* [Unicycler](#Unicycler)
* [OPERA](#OPERA)
* [GapCloser](#GapCloser)
* [Contiguity](#Contiguity)

#### Post-processing, variant calling
* [GATK](#GATK)
* [graphmap](#graphmap)
* [lofreq](#lofreq)
* [pilon](#pilon)
* [nanopolish](#nanopolish)

#### Annotation and classification
* [abricate](#abricate)
* [Kraken](#Kraken)
* [prokka](#prokka)
* [SeqSero](#SeqSero)

#### Visualization
* [bandage](#bandage)
* [BRIG](#BRIG)
* [Circos](#Circos)
* [EasyFig](#EasyFig)
* [SeqFindr](#SeqFindr)
* [slcview](#slcview)

### Customized installs
* [ASCP](#ASCP)
* [FinIS](#FinIS)
* [genome-tools](#genome-tools)
* [SRST2](#SRST2)

## Detailed instructions

### General
* [BLAST+](#BLAST)
* [HYPHY](#HYPHY)
* [Kingfisher](#Kingfisher)
* [meme](#meme)
* [samtools](#samtools-including-htslib)
* [sra-tools](#sra-tools)

#### [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
The Ubuntu repositories have this, but it's at version 2.9.0-2. The current version right now is at 2.11.0-1. Generally for most users, BLAST has been pretty stable, but here's how to update it if you need the latest version (you'll have to repeat these for each update as well after you do this manual install).

First check on what we already have and where it is
```
# check the current version
blastn -h
# this should be 2.9.0+ if it's from the Ubuntu repositories
# this will likely be 2.11.0+ or higher if you manually updated it as below

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
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.11.0+-x64-linux.tar.gz
tar xvzf ncbi-blast-2.11.0+-x64-linux.tar.gz
cd ncbi-blast-2.11.0+
# always good to look at the README, though not much in this one
less README
cd /usr/local/src/ncbi-blast-2.11.0+/bin
for i in *; do ln -s /usr/local/src/ncbi-blast-2.11.0+/bin/$i /usr/local/bin; done
```

You can verify that the default version is what you expect with the commands in the first block above. You may also want to remove the older Ubuntu version to avoid any potential confusion: `sudo apt purge ncbi-blast+ ncbi-blast+-legacy blast2`.

#### [HYPHY](https://github.com/veg/hyphy)
I haven't had great luck with conda, and this setup is intended to be used for a single person or to drive a specific pipeline. Therefore, we'll compile HYPHY from source using the latest release tarball.
```
sudo su -
cd /usr/local/src
wget https://github.com/veg/hyphy/archive/refs/tags/2.5.31.tar.gz
tar xvzf 2.5.31.tar.gz
cd hyphy-2.5.31
# check the instructions
less README.md

# then do the installation. We have openMPI installed from the initial system setup
cmake .
make HYPHYMPI
# this defaults to install to /usr/local
make install
```

#### [meme](http://meme-suite.org/meme/)
Current version is 5.3.3.
```
sudo su -
cd /usr/local/src
wget https://meme-suite.org/meme/meme-software/5.3.3/meme-5.3.3.tar.gz
tar xvzf meme-5.3.3.tar.gz
cd meme-5.3.3

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
echo 'export PATH=$PATH:/usr/local/libexec/meme-5.3.3' >> /home/ubuntu/.bashrc
```

#### [Kingfisher](https://github.com/wwood/kingfisher-download)
This is a convenient tool for downloading public data sets, such as from GenBank or ENA.
Again, I haven't had a lot of good experience with conda.
Fortunately, this only requires one additional library to install (extern).
```
sudo su -
pip3 install extern
cd /usr/local/src
git clone https://github.com/wwood/kingfisher-download
cd kingfisher-download
# check the docs
less README.md
ln -s /usr/local/src/kingfisher-download/bin/kingfisher /usr/local/bin
```
If you've installed [ASCP](#ASCP) as described in this guide, you can use the `-m ena-ascp` method as well as the more standard methods.

#### [samtools](http://www.htslib.org/) (including htslib)

The HTSlib / samtools suite provides core tools and libraries that are also used by many other bioinformatics software tools. Again there is a version in the Ubuntu repositories (1.10), but the current version as of this writing is 1.12. This may or may not matter for what you're doing. If it does, here's how to update it.

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
wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2
tar xvjf samtools-1.12.tar.bz2
cd samtools-1.12
# always check the basic docs
less README
less INSTALL
# do the main configuration and installation
./configure --enable-configure-htslib
make
# the default prefix is to install to /usr/local already
make install

# htslib is configured - go make and install
cd /usr/local/src/samtools-1.12/htslib-1.12
make
make install
```

Then bcftools:
```
sudo su -
cd /usr/local/src
wget https://github.com/samtools/bcftools/releases/download/1.12/bcftools-1.12.tar.bz2
tar xvjf bcftools-1.12.tar.bz2
cd bcftools-1.12
# always check the basic docs
less README
less INSTALL
# do the main configuration and installation
./configure
make
# the default prefix is to install to /usr/local already
make install
```

Again, you can go back above to verify what version you're now using by default (you may need to do a `hash -r` to have bash update program locations). Also, you can remove the version from the Ubuntu repositories if you want to avoid potential confusion: `sudo apt purge samtools libhts3 bcftools`.

Note a couple useful pieces of information for the hts libraries:
```
/usr/local/lib/libhts.a
/usr/local/include/htslib
```

#### [sra-tools](https://github.com/ncbi/sra-tools)
These are useful tools for using data in the GenBank Sequence Read Archives.
They provide precompiled binaries on their GitHub page.
```
sudo su -
cd /usr/local/src
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.0/sratoolkit.2.11.0-ubuntu64.tar.gz
tar xvzf sratoolkit.2.11.0-ubuntu64.tar.gz 
cd sratoolkit.2.11.0-ubuntu64
# check the docs
less README.md
less README-blastn
less README-vdb-config

# link in the binaries - only take the parents, not the version-named ones
for i in `ls -1 | grep -v '\.2$' | grep -v '2.11.0$'`; do \
  ln -s /usr/local/src/sratoolkit.2.11.0-ubuntu64/bin/$i /usr/local/bin; 
done
# this needs to be configured once (both root and ubuntu), even if accepting the defaults
vdb-config -i
# type 'x' to just exit with default configuration
# now go back to being the ubuntu user, probably exit from the su shell
vdb-config -i
# type 'x' to exit with default configuration
# this should leave a file in /home/ubuntu/.ncbi/user-settings.mkfg
```

### Read processing
* [bowtie](#bowtie)
* [bwa](#bwa)
* [deML](#deML)
* [fastp](#fastp)
* [lacer](#lacer)
* [minimap](#minimap)
* [Porechop](#Porechop)
* [poretools](#poretools)
* [seqmagick](#seqmagick)
* [seqtk](#seqtk)
* [Trimmomatic](#Trimmomatic)

#### [bowtie2](https://github.com/BenLangmead/bowtie2)
This is one of the well known short read mappers.
[SRST2](#SRST2) uses this but has some version requirements (2.2.9).
The Ubuntu focal LTS repositories include bowtie2 at version 2.3.5.1.
The latest (May 2021) online at GitHub is version 2.4.2.

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
wget https://github.com/BenLangmead/bowtie2/releases/download/v2.4.2/bowtie2-2.4.2-linux-x86_64.zip
unzip bowtie2-2.4.2-linux-x86_64.zip
cd bowtie2-2.4.2-linux-x86_64

# always check the docs
less README.md
for i in bowtie*; do ln -s /usr/local/src/bowtie2-2.4.2-linux-x86_64/$i /usr/local/bin; done

# check where it is and the version - this should be /usr/local/bin/bowtie2
which bowtie2
bowtie2 --version
```

#### [bwa](https://github.com/lh3/bwa)
This isn't being updated so much though so I prefer to use the Ubuntu repositories, which have the same main version as the current release (as of May 2021) on the GitHub repository (0.7.17).

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
```
sudo su -
cd /usr/local/src
git clone https://github.com/grenaud/deML.git
cd deML
# check the docs
less README.md

make
ln -s /usr/local/src/deML/src/deML /usr/local/bin

# check the binary runs ok
deML
```

#### [fastp](https://github.com/OpenGene/fastp)
Ubuntu LTS version: 0.20.0
GitHub version: 0.20.1

We'll install the GitHub version. (The Ubuntu version of course is straightforward with an `apt install fastp`).

```
sudo su -
cd /usr/local/src
wget https://github.com/OpenGene/fastp/archive/refs/tags/v0.20.1.tar.gz
tar xvzf v0.20.1.tar.gz
cd fastp-0.20.1/
# check the docs
less README.md

make
ln -s /usr/local/src/fastp-0.20.1/fastp /usr/local/bin

# check the binary runs ok
fastp
```

#### [lacer](https://github.com/swainechen/lacer)
GitHub version: 0.424

This is a base quality score recalibrator (the "Q" in FASTQ files).
It does not require knowledge of common SNPs and therefore is the only program that can generally recalibrate base quality scores on any organism (programs like GATK's BaseRecalibrator were designed primarily for human sequencing data, for example, and therefore only work well on data from a limited set of organisms).
This means that base quality recalibration, leading to more accurate trimming and SNP calling, can now be used on any organism (specifically including bacteria).

The main program is just a perl script. There's a fast and stripped-down `lacepr` program that does the basics of what GATK's PrintReads does.

There are some Perl dependencies - if you've done everything in the [initial system setup](system.md), all of these should be available already.
```
sudo su -
cd /usr/local/src
git clone https://github.com/swainechen/lacer.git
cd lacer
# check the docs
less README.md

ln -s /usr/local/src/lacer/lacer.pl /usr/local/bin
```

To get lacepr to compile and install, we need an older version of samtools (up to 1.9). There were many changes starting in samtools/htslib version 1.10 that haven't been incorporated into lacepr yet.
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
make
cd /usr/local/src/lacer/lacepr
SAMTOOLS=/usr/local/src/lacer/lacepr/samtools-1.9
HTSLIB=/usr/local/src/lacer/lacepr/samtools-1.9/htslib-1.9
gcc -I$SAMTOOLS -I$HTSLIB lacepr.c -L$SAMTOOLS -L$HTSLIB -lbam -l:libhts.a -lz -lpthread -lm -llzma -lbz2 -o lacepr
ln -s /usr/local/src/lacer/lacepr/lacepr /usr/local/bin
```

#### [minimap](https://github.com/lh3/minimap2)
The Ubuntu Focal (20.04) LTS repositories have version 2.17.
The current release version (May 2021) in GitHub is 2.18.
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
wget https://github.com/lh3/minimap2/releases/download/v2.18/minimap2-2.18_x64-linux.tar.bz2
tar xvjf minimap2-2.18_x64-linux.tar.bz2
ln -s /usr/local/src/minimap2-2.18_x64-linux/minimap2 /usr/local/bin
ln -s /usr/local/src/minimap2-2.18_x64-linux/minimap2.1 /usr/local/man/man1/

# check where it is and the version
which minimap2
minimap2
```

#### [Porechop](https://github.com/rrwick/Porechop)
This is a tool for adapter trimming for Oxford Nanopore sequencing data.
It seems to work, but officially is no longer supported.
Since it hasn't been updated in a while, the Ubuntu version is up-to-date with what's on GitHub, so use that version.

Install from the Ubuntu Focal LTS repositories:
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
This is a useful utility for converting between sequence formats.
This can be installed with pip (pip3 for python3).
Installing as root as for other software, this will also be available for the ubuntu user (and others).
```
sudo pip3 install seqmagick
seqmagick -V
```

#### [seqtk](https://github.com/lh3/seqtk)
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

#### [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
Ubuntu focal LTS version: 0.39
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
Ubuntu focal LTS version: 38.79
SourceForge version: 38.90

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
wget 'https://downloads.sourceforge.net/project/bbmap/BBMap_38.90.tar.gz?ts=gAAAAABglrnyrwR2Q1MnHROncBIoOt3unY2XNT7spAMQTg8KopS73pB07gOriP7rcnn1yqmDMctWFxcJPSV9x0-K8AUaV6hXJw%3D%3D&r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fbbmap%2Ffiles%2Flatest%2Fdownload' -O BBMap_38.90.tar.gz
tar xvzf BBMap_38.90.tar.gz
cd bbmap
# check the docs
less README.md
less docs/readme.txt

# we follow what the Ubuntu maintainers did and only link the scripts for bbduk, bbmap, bbnorm, bloomfilter, dedupe, and reformat
for i in bbduk bbmap bbnorm bloomfilter dedupe reformat; do ln -s /usr/local/src/bbmap/$i.sh /usr/local/bin; done
```

### Assembly
* [a5 assembler](#a5-assembler)
* [canu](#canu)
* [miniasm](#miniasm)
* [racon](#racon)
* [sga](#sga)
* [skesa](#skesa)
* [SPAdes](#SPAdes)
* [velvet](#velvet)
* [Trycycler](#Trycycler)
* [Unicycler](#Unicycler)
* [OPERA](#OPERA)
* [GapCloser](#GapCloser)
* [Contiguity](#Contiguity)

#### a5 assembler
#### canu
#### miniasm
#### racon
#### sga

#### [skesa](https://github.com/ncbi/SKESA)
Ubuntu LTS version: 2.3.0
GitHub version: 2.4.0

This is a set of assemblers (SKESA and SAUTE) with some companion programs. SKESA in particular was designed for microbial genomes.

We'll install the latest GitHub release version (tagged as "Update for NGS v.2.11.0" as of May, 2021 - but this is v2.4.0 of SKESA (the 2.11.0 is for the NGS bit)).
We'll take the source and compile.
```
sudo su -
cd /usr/local/src
wget https://github.com/ncbi/SKESA/archive/refs/tags/skesa.2.4.0_saute.1.3.0_1.tar.gz
tar xvzf skesa.2.4.0_saute.1.3.0_1.tar.gz
cd SKESA-skesa.2.4.0_saute.1.3.0_1
# check the docs
less README.md

# build - there are some assumptions for the full make where it tries to
# build ngs-sdk, even though we already have it installed, leading to a compile
# error. For now use the nongs version
make -f Makefile.nongs

# link the binaries
for i in skesa saute saute_prot gfa_connector kmercounter; do ln -s /usr/local/src/SKESA-skesa.2.4.0_saute.1.3.0_1/$i /usr/local/bin; done

# test that it works
skesa
saute
```

#### [SPAdes](https://cab.spbu.ru/software/spades/)
Ubuntu LTS version: 3.13.1
Online version: 3.15.2
However, Unicycler needs a version no later than 3.13.0. We'll install the latest online version for regular use, then do a second install with the older version for Unicycler.

The latest version from https://cab.spbu.ru/software/spades/ as default (i.e. in `/usr/local/bin`)
```
sudo su -
cd /usr/local/src
wget https://cab.spbu.ru/files/release3.15.2/SPAdes-3.15.2-Linux.tar.gz
tar xvzf SPAdes-3.15.2-Linux.tar.gz
cd SPAdes-3.15.2-Linux
# check the docs
less share/spades/README.md

# link in the binaries
cd bin
for i in *; do ln -s /usr/local/src/SPAdes-3.15.2-Linux/bin/$i /usr/local/bin; done

# check it works ok
```

The older 3.13.0 version for Unicycler - this can be downloaded at the [SPAdes GitHub repository](https://github.com/ablab/spades)
We will keep this in /usr/local/src and will need to specify the full path in order to use it:
```
sudo su -
cd /usr/local/src
wget https://github.com/ablab/spades/releases/download/v3.13.0/SPAdes-3.13.0-Linux.tar.gz
tar xvzf SPAdes-3.13.0-Linux.tar.gz
cd SPAdes-3.13.0-Linux

# this was a binary distribution, so it should be done - just check the version
/usr/local/src/SPAdes-3.13.0-Linux/bin/spades.py
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
ln -s /usr/local/src/velvet/contrib/VelvetOptimiser-2.2.4/VelvetOptimiser.pl /usr/local/bin
ln -s /usr/local/src/velvet/contrib/shuffleSequences_fasta/shuffleSeq* /usr/local/bin
ln -s /usr/local/src/velvet/contrib/VelvetOptimiser-2.2.4/VelvetOpt/ /usr/local/lib/site_perl/
```
The last line for the VelvetOpt perl libraries goes to `/usr/local/lib/site_perl`.
This should work on a stock Ubuntu install.
You can check the options for this by running `perl -e 'print join ("\n", @INC), "\n";'` - this is the default search path for perl modules, and it just needs to be in one of those places (though something under `/usr/local/` is a good idea).

#### [Trycycler](https://github.com/rrwick/Trycycler/wiki)

#### [Unicycler](https://github.com/rrwick/Unicycler)
Ubuntu LTS version: 0.4.8
GitHub version: 0.4.9

We'll use the latest GitHub version.

#### [OPERA](https://github.com/CSB5/OPERA-MS)
#### GapCloser
#### Contiguity

### Post-processing, variant calling
* [GATK](#GATK)
* [graphmap](#graphmap)
* [lofreq](#lofreq)
* [pilon](#pilon)
* [nanopolish](#nanopolish)

#### GATK
#### graphmap

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

#### pilon
#### nanopolish

### Annotation and classification
* [abricate](#abricate)
* [Kraken](#Kraken)
* [prokka](#prokka)
* [SeqSero](#SeqSero)

#### abricate
#### Kraken
#### prokka
#### SeqSero

### Visualization
* [bandage](#bandage)
* [BRIG](#BRIG)
* [Circos](#Circos)
* [EasyFig](#EasyFig)
* [SeqFindr](#SeqFindr)
* [slcview](#slcview)

#### [Bandage](https://rrwick.github.io/Bandage/)
#### [BRIG]
#### [Circos](http://circos.ca/)
#### [EasyFig]
#### [SeqFindr]
#### [slcview](https://github.com/swainechen/slcview)

### Customized installs
* [ASCP](#ASCP)
* [FinIS](#FinIS)
* [genome-tools](#genome-tools)
* [SRST2](#SRST2)

#### [ASCP](http://downloads.asperasoft.com/connect2/)
This is useful for fast downloads, such as from ENA or Genbank. This seems to be best (and easiest) to install as the user (ubuntu).
```
# make sure you're in a user (ubuntu) shell and not root
cd /home/ubuntu
wget https://d3gcli72yxqn2z.cloudfront.net/connect_latest/v4/bin/ibm-aspera-connect-3.11.2.63-linux-g2.12-64.tar.gz
tar xvzf ibm-aspera-connect-3.11.2.63-linux-g2.12-64.tar.gz
./ibm-aspera-connect-3.11.2.63-linux-g2.12-64.sh
# update your path to include this (for ex. in .bashrc)
echo '# path for ascp' >> /home/ubuntu/.bashrc
echo 'export PATH=$PATH:/home/ubuntu/.aspera/connect/bin' >> /home/ubuntu/.bashrc
```
Note that, if needed, the standard key required is at `/home/ubuntu/.aspera/connect/etc/asperaweb_id_dsa.openssh`.

#### [FinIS](https://sourceforge.net/p/finis/wiki/FinIS%20wiki/)
FinIS is an assembly finisher, which generally is used after the [OPERA](#OPERA) scaffolder.
This is on SourceForge at v0.3.
This is somewhat complex to get running, as it hasn't been updated in a while.
It seems to suffer from some changes in C++ conventions since 2014, and also was last tested on [MOSEK](https://www.mosek.com/) 6 (this is currently at 9).
I've provided a compiled binary that will run on the Ubuntu images at AWS. Note this still needs MOSEK 6:
[FinIS binary]()
```
sudo su -
cd /usr/local/bin
wget #link_to_FinIS

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
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/bin
echo "# for MOSEK 6" >> /home/ubuntu/.bashrc
echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib" >> /home/ubuntu/.bashrc

# we can verify this with test data that comes with the FinIS source
# again this is on SourceForge, you may need to generate your own direct link
# by clicking starting from https://sourceforge.net/projects/finis/files/v0.3/
cd /usr/local/src
wget 'https://downloads.sourceforge.net/project/finis/v0.3/v0.3.tar.gz?ts=gAAAAABgl49QRj_nAX5gccQ9FeGwccNWKR_QnNWBbNTwV1HGd2cpWueqjuolVNhg8C27rHgE0kQdAd0IeYzGRkmTyTvowuvqsQ%3D%3D&r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Ffinis%2Ffiles%2Fv0.3%2Fv0.3.tar.gz%2Fdownload' -O v0.3.tar.gz
tar xvzf v0.3.tar.gz
mv v0.3 FinIS-0.3
cd FinIS-0.3

# there are two test datasets provided
FinIS test_dataset/velvet/conf.config

# the soap config file needs to drop the num_threads and mosek_runtime lines
# these seem to have changed as well - comment these two out and the test
# should run fine - note this requires more memory (~4GB) than
# a t3a.small instance has (2GB) however
FinIS test_dataset/soap/conf.config
```

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
python2 get-pip.py
pip install scipy
git clone https://github.com/katholt/srst2
pip install srst2/

# install samtools 0.1.18
sudo su -
cd /usr/local/src
wget https://github.com/samtools/samtools/archive/refs/tags/0.1.18.tar.gz
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
I also have a few utility scripts to help manage these.
```
# Setting up SRST2 reference libraries
```
