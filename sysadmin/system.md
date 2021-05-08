# Initial "New Car" Setup and General Sysadmin
Starting from a bare system image means we have to do some initial updates and get all the software installed. General system stuff will be done here, as well as a few preparatory steps that will be needed for installing some more specific bioinformatics tools later.

First, initial updates (this is the standard updating I do on machines, you can stick this into cron, but I like to make sure updates aren't happening during some long-running process. Userdata for AWS instances is also a good place to do an initial update on startup):
```
sudo su -
apt update
apt upgrade -y
```
Everything after this will be run as root - we're doing sysadmin! You can also add `sudo` in front of all these commands below if you're staying as the default `ubuntu` user.

### Use of `sudo` and the root account
I should note that running all this stuff as root isn't the best practice. The [standard recommendation](https://tldp.org/HOWTO/Software-Building-HOWTO-3.html) is to only use sudo or the root account when it's absolutely needed. There are a couple ways to do this:
- download into /home/ubuntu/src, build there, then sudo install
- make a `staff` or `admin` or `src` group, add the `ubuntu` user to that group, `chown -R root:staff /usr/local/src` (or even the whole `/usr/local` tree), then do all the builds as ubuntu +/- the install as root

Here, though, we're doing the simple version to just run around as root. This takes away a layer of potential problems and if you're on a fresh cloud machine, it's not a huge deal to mess something up since there really shouldn't be other users you're concerned about and you can just make a new instance.

## Software Repositories and PPAs
Add some repositories we'll need later:
```
apt install software-properties-common dirmngr
apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
apt update
```

## Standard Ubuntu Packages (from the repositories)
Install regular ubuntu packages that we'll need for later. This took about 30 min on a t3a.small AWS instance.
```
apt install -y automake awscli cmake cpanminus cython evince gnuplot-nox \
  imagemagick libbio-samtools-perl libcairo2-dev libcurl4-openssl-dev \
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

## Perl 

Many modules are available in the Ubuntu repositories. These may be older versions, but generally for our purposes these will be "stable enough" for us, and the auto-updating and security fixes are important.

For those modules that aren't available already, it's useful to set up CPAN to manage the installs. [CPAN Instructions](http://www.cpan.org/modules/INSTALL.html)

Here we'll install just two modules that help with some of the other software later (i.e. meme, lacer): `Math::CDF` `PDL::Parallel::threads`

First, `Math::CDF`:
```
sudo su -
cpanm Math::CDF

# check it works - the following command should give no error if things are ok
perl -e 'use Math::CDF'
```

The next module, `PDL::Parallel::threads`, is required for `lacer`. This has a little bug and requires some tweaking, documented [here](https://github.com/run4flat/PDL-Parallel-threads/issues/1).
```
sudo su -
cnapm PDL::Parallel::threads
# will see an error. It should point you to a build.log, go to that directory
# the directory below is what I had - the last part will be different for you
cd /root/.cpanm/work/1620464035.165032/
# edit Build.PL as per https://github.com/run4flat/PDL-Parallel-threads/issues/1
# this was on line 8 for me

perl Build.PL
./Build
./Build test
./Build install

# test to see if it works - no error and no output for the following if ok
perl -e 'use PDL::Parallel::threads'
```

## AWS Utilities

Generally we would try to install system packages and relatively stable software from the Ubuntu repositories - this helps with security updates and keeping things organized on the system. One quirk at the time of this writing (May 2021) is that the Ubuntu awscli package is at 1.17.14-1, and there ends up being a library issue (see [here](https://github.com/boto/boto3/issues/2596)) - so we'll use pip3 to do the install (currently at version 1.19.69). (Note that version 2 of the awscli suite seems to require manual install (i.e. therefore requiring manual updates...).)
```
pip3 install awscli
```
