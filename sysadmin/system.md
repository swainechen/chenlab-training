# Initial "New Car" Setup and General Sysadmin
Starting from a bare system image means we have to do some initial updates and get all the software installed. General system stuff will be done here, as well as a few preparatory steps that will be needed for installing some more specific bioinformatics tools later.

First, initial updates (this is the standard updating I do on machines, you can stick this into cron, but I like to make sure updates aren't happening during some long-running process. Userdata for AWS instances is also a good place to do an initial update on startup):
```
sudo su -
apt update
apt upgrade -y
```
Everything after this will be run as root - we're doing sysadmin! You can also add `sudo` in front of all these commands below if you're staying as the default `ubuntu` user.

Add some repositories we'll need later:
```
apt install software-properties-common dirmngr
apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
apt update
```

Install regular ubuntu packages that we'll need for later. This took about 30 min on a t3a.small AWS instance.
```
apt install -y automake awscli cmake cython evince gnuplot-nox imagemagick \
  libbio-samtools-perl libcairo2-dev libcurl4-openssl-dev \
  libdatetime-format-dateparse-perl libdatetime-format-dbi-perl \
  libfile-type-perl libfile-which-perl libhdf5-dev libhtml-template-perl \
  libimage-magick-perl libjemalloc-dev libjson-perl liblzma-dev \
  libmariadbclient-dev libmemory-usage-perl libmodule-build-perl \
  libopenmpi-dev libsparsehash-dev libssl-dev libterm-progressbar-perl \
  libtext-csv-perl libv8-dev libxml-compile-perl libxml-compile-wsdl11-perl \
  libxml2-dev libxslt1-dev mlocate mysql-client openjdk-8-jdk parallel pdl \
  prodigal python python-numpy snakemake zlib1g-dev

# some initial software
apt install -y cd-hit clonalframeml fastdnaml fastqc fasttree \
  gubbins kraken njplot ncbi-entrez-direct ncbi-tools-bin paml soapdenovo2 \
  vcftools mauve-aligner mummer
```

Generally we would try to install system packages and relatively stable software from the Ubuntu repositories - this helps with security updates and keeping things organized on the system. One quirk at the time of this writing (May 2021) is that the Ubuntu awscli package is at 1.17.14-1, and there ends up being a library issue (see [here](https://github.com/boto/boto3/issues/2596)) - so we'll use pip3 to do the install (currently at version 1.19.69). (Note that version 2 of the awscli suite seems to require manual install (i.e. therefore requiring manual updates...).)
```
pip3 install awscli
```
