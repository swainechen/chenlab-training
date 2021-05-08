# Bioinformatics software
Bioinformatics is moving very fast. It's hard to keep up with version changes, and because of some details of academic funding and training, there is a lot of variation in organization and design.

I've split this section into two parts - software that can be installed a relatively "standard" unix-like manner, and those that require some custom updates or custom configuration.
* [Software that follows a standard installation](#standard-installs)
* [Software requiring some customization](#customized-installs)

Generally the strategy is to keep all the downloads and sources in /usr/local/src and then link into /usr/local/bin or other appropriate directories in /usr/local.

Note that the version numbers are included in the commands below - you'll have to update those to whatever the current version is at the time.

## Standard installs
### General
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

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

* HYPHY

* meme

* [samtools](http://www.htslib.org/) (including htslib)

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

### Read processing
* bowtie
* [bwa](https://github.com/lh3/bwa)

This isn't being updated so much though so I prefer to use the Ubuntu repositories, which have the same main version as the current release (as of May 2021) on the github repository (0.7.17).

From the Ubuntu repositories (*Recommended*):
```
sudo apt install bwa

# check where it is and the version
which bwa
bwa
```

From the github release tarball:
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

From the github repository:
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

* deML
* fastp
* [minimap](https://github.com/lh3/minimap2)

The Ubuntu Focal (20.04) LTS repositories have version 2.17.
The current release version (May 2021) in github is 2.18.
As this is under active development, we'll have to keep up with new releases (similar instructions to below).

From the Ubuntu repositories:
```
sudo apt install minimap2

# check where it is and the version
which minimap2
minimap2
```

From the github release tarball (*Recommended*):
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

* porechop
* poretools
* seqmagick
* seqtk
* Trimmomatic

### Assembly
* a5 assembler
* canu
* miniasm
* racon
* sga
* SPAdes
* velvet
* unicycler
* OPERA
* GapCloser
* Contiguity

### Post-processing, variant calling
* GATK
* graphmap
* lofreq
* pilon
* nanopolish

### Annotation and classification
* abricate
* Kraken
* prokka
* SeqSero

### Visualization
* slcview
* BRIG
* EasyFig
* SeqFindr

## Customized installs
* [ASCP](from http://downloads.asperasoft.com/connect2/)

This is useful for fast downloads, such as from ENA or Genbank. This seems to be easiest to install as the user (ubuntu).
```
# make sure you're in a user (ubuntu) shell and not root
cd /home/ubuntu
wget https://d3gcli72yxqn2z.cloudfront.net/connect_latest/v4/bin/ibm-aspera-connect-3.11.2.63-linux-g2.12-64.tar.gz
tar xvzf ibm-aspera-connect-3.11.2.63-linux-g2.12-64.tar.gz
./ibm-aspera-connect-3.11.2.63-linux-g2.12-64.sh
# update your path to include this (for ex. in .bashrc)
echo 'export PATH=$PATH:/home/ubuntu/.aspera/connect/bin' >> .bashrc
```
Note that, if needed, the standard key required is at `/home/ubuntu/.aspera/connect/etc/asperaweb_id_dsa.openssh`.

* FinIS

* genome-tools

* SRST2
