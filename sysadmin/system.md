# Initial "New Car" Setup and General Sysadmin
```
sudo su -
apt update
apt upgrade -y
```

Add some repositories we'll need later
```
apt install software-properties-common dirmngr
apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
apt update
```

Install regular ubuntu packages that we'll need for later. This took about 30 min on a t3a.small AWS instance.
```
apt install -y automake cmake cython evince gnuplot-nox imagemagick libbio-samtools-perl libcairo2-dev libcurl4-openssl-dev libdatetime-format-dateparse-perl libdatetime-format-dbi-perl libfile-type-perl libfile-which-perl libhdf5-dev libhtml-template-perl libimage-magick-perl libjemalloc-dev libjson-perl liblzma-dev libmariadbclient-dev libmemory-usage-perl libmodule-build-perl libopenmpi-dev libsparsehash-dev libssl-dev libterm-progressbar-perl libtext-csv-perl libv8-dev libxml2-dev libxslt1-dev mysql-client openjdk-8-jdk parallel pdl prodigal python python-numpy snakemake zlib1g-dev

# some initial software
apt install -y bcftools cd-hit clonalframeml fastdnaml fastqc fasttree gubbins kraken njplot ncbi-entrez-direct ncbi-tools-bin paml soapdenovo2 vcftools mauve-aligner mummer
```
