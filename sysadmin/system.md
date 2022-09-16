# Initial "New Car" Setup and General Sysadmin
Starting from a bare system image means we have to do some initial updates and get all the software installed. General system stuff will be done here, as well as a few preparatory steps that will be needed for installing some more specific bioinformatics tools later.

## Sections:
* [First steps](#First-steps)
* [Software Repositories and PPAs](#Software-Repositories-and-PPAs)
* [Standard Ubuntu Packages](#Standard-Ubuntu-Packages)
* [Perl](#Perl)
* [AWS Utilities](#AWS-Utilities)
* [User environment](#User-environment)
  - [Bash](#Bash)
    - [Command Prompt](#Command-Prompt)
    - [Bash options](#Bash-options)
    - [Path](#Path)
    - [Aliases](#Aliases)
    - [Editors](#Editors)
    - [Screen](#Screen)
    - [Non-interactive shells](#Non-interactive-shells)
  - [SSH Keys](#SSH-Keys)
  - [Extra data storage](#Extra-data-storage)

## First steps
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
# For R
sudo su -
apt install --no-install-recommends software-properties-common dirmngr
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
apt update

# For Docker
sudo mkdir -p /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg
echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
apt update
```

## Standard Ubuntu Packages
Install regular Ubuntu packages from the base repositories that we'll need for later. This took about 30 min on a t3a.small AWS instance.
```
####################
# for Ubuntu 20.04 #
####################
apt install -y automake cmake cpanminus cpanoutdated cython evince fig2dev \
  gnuplot-nox imagemagick libbio-samtools-perl libboost-all-dev \
  libcairo2-dev libcurl4-openssl-dev libdatetime-format-dateparse-perl \
  libdatetime-format-dbi-perl libdb-dev libfile-type-perl libfile-which-perl \
  libhdf5-dev libhtml-template-perl libimage-magick-perl libjemalloc-dev \
  libjson-perl liblzma-dev libmariadbclient-dev libmemory-usage-perl \
  libmodule-build-perl libopenmpi-dev libsparsehash-dev libssl-dev \
  libsys-meminfo-perl libterm-progressbar-perl libtext-csv-perl libv8-dev \
  libxml-compile-perl libxml-compile-wsdl11-perl libxml2-dev libxslt1-dev \
  mlocate mysql-client openjdk-8-jdk parallel pdl prodigal python \
  python-numpy python3-pip snakemake swig xfig zlib1g-dev

# some initial software
apt install -y cd-hit clonalframeml fastdnaml fastqc fasttree \
  gubbins njplot ncbi-entrez-direct ncbi-tools-bin paml snp-sites \
  soapdenovo2 vcftools mauve-aligner

####################
# for Ubuntu 22.04 #
####################
apt install -y automake clang cmake cpanminus cpanoutdated docker-ce \
  docker-ce-cli containerd.io docker-compose-plugin evince fig2dev \
  gnuplot-nox imagemagick libbio-samtools-perl libboost-all-dev libbz2-dev \
  libcairo2-dev libcurl4-openssl-dev libdatetime-format-dateparse-perl \
  libdatetime-format-dbi-perl libdb-dev libfile-map-perl libfile-type-perl \
  libfile-which-perl libgeos-dev libhdf5-dev libhtml-template-perl \
  libimage-magick-perl libjemalloc-dev libjson-perl liblapack-dev \
  liblzma-dev libmariadb-dev libmemory-usage-perl libmodule-build-perl \
  libncurses-dev libopenmpi-dev libsparsehash-dev libssl-dev \
  libsys-info-perl libsys-meminfo-perl libterm-progressbar-perl \
  libtext-csv-perl libv8-dev libxml-compile-perl \
  libxml-compile-wsdl11-perl libxml2-dev libxslt1-dev meson mlocate \
  mysql-client openjdk-8-jdk parallel pdl prodigal python-is-python3 \
  python3-pip snakemake swig xfig zlib1g-dev

# some initial software
apt install -y cd-hit clonalframeml fastdnaml fastqc fasttree \
  gubbins njplot ncbi-entrez-direct ncbi-tools-bin paml snp-sites \
  soapdenovo2 vcftools mauve-aligner

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

The next module, `PDL::Parallel::threads`, is required for `lacer`.

_N.B.: Until recently (Feb 2022) this module had a little bug that required some tweaking, documented [here](https://github.com/run4flat/PDL-Parallel-threads/issues/1). As of version 0.04, however, it compiles and installs cleanly with just the commands below._
```
sudo su -
cpanm PDL::Parallel::threads

# test to see if it works - no error and no output for the following if ok
perl -e 'use PDL::Parallel::threads'
```

## AWS Utilities

Generally we would try to install system packages and relatively stable software from the Ubuntu repositories - this helps with security updates and keeping things organized on the system. One quirk at the time of this writing (May 2021) is that the Ubuntu awscli package is at 1.17.14-1, and there ends up being a library issue (see [here](https://github.com/boto/boto3/issues/2596)) - so we'll use pip3 to do the install (currently at version 1.19.69). (Note that version 2 of the awscli suite seems to require manual install (i.e. therefore requiring manual updates...).)
```
# for version 1
pip3 install awscli

# for version 2 (installed on CHENLAB-PUBLIC v2.1x AMIs)
sudo su -
cd /usr/local/src
wget https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip
unzip awscli-exe-linux-x86_64.zip
./aws/install
```

## User environment

Many of us spend a lot of time on the command line. It's worth it for your wrists and for your productivity to optimize your environment. Here are a few things I have set up in my shell.

### Bash
#### Command Prompt
A good command prompt can be very helpful.
I picked this one up from an old slashdot.org post long ago and have only made minor tweaks.
You can play with the [ANSI color codes](https://gist.github.com/fnky/458719343aabd01cfb17a3a4f7296797) (here's a [howto](https://www.lihaoyi.com/post/BuildyourownCommandLinewithANSIescapecodes.html)) and wrap this in an `case \`hostname -s` in` statement to give you a different color on different machines if you're using a few regularly!

This just can go into your `.bashrc` (there's an if statement at lines 59-63 in the default Ubuntu `.bashrc`, you'll probably want to modify that if statement as well). 
```
PS1="\[\033[00m\]\!\\[\033[0m\]|\`if [[ \$? = "0" ]]; then echo "\\[\\033[32m\\]"; else echo "\\[\\033[31m\\]"; fi\`\u@\h:\[\033[36m\]\`if [[ `pwd|wc -c|tr -d " "` > 18 ]]; then echo "\\W"; else echo "\\w"; fi\`\$\[\033[0m\] "; echo -ne "\033]0;`hostname -s`:`pwd`\007"
```
This prompt gives you command number, username, host, and current path.
The color of the user@host part is green if the last command had an exit code of 0 (success) and red if not.

One thing to note is that you should respect the prompt for a non-interactive shell. See [below](#Non-interactive-shells) - make sure this is only set for interactive shells. Otherwise, you can run into problems like [this](https://superuser.com/questions/395356/scp-doesnt-work-but-ssh-does) with `scp` and `rsync`.

#### Bash options
The standard one to set is `noclobber` for safety!
There are quite a few others that you should probably explore - [tldp link](https://tldp.org/LDP/abs/html/options.html)

```
echo "# for safety" >> /home/ubuntu/.bashrc
echo "set -o noclobber" >> /home/ubuntu/.bashrc
```

#### Path
Many software packages outside the standard repositories (i.e. like most biology software) would try to avoid problems with unknown potential configurations by just putting all files in one directory, then telling you to add something to your path.

This is fine from a practical point of view at first, but with a lot of software can get unwieldy to manage and slow to search (in extreme cases).
All the setup in this guide aims to keep the path pretty minimal and standard.
Where it seems more sensible to add something to the `PATH`, this is included in the instructions in this guide, and those are done at the end of `/home/ubuntu/.bashrc`.

The default `PATH` is typically:
```
echo $PATH
# /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin
```

Some of the software in this guide requires some additional paths due to difficulties in getting them "nicely" installed into a standard tree.

Two other things you sometimes see on a path are:
* `/home/ubuntu/bin` (or equivalent, may be set as `~/bin`)
* `.` (i.e. the current directory)

On a multiuser system, I often would include my own personal bin directory.
This setup guide is assuming you're on the cloud and basically have your own machine, and root permissions - in which case, I prefer to just manage /usr/local/bin.
The inclusion of `.`, the current directory, can be convenient when you're starting out, but leads to some [potential](http://www.faqs.org/faqs/unix-faq/faq/part2/section-13.html) [attacks](https://superuser.com/questions/156582/why-is-not-in-the-path-by-default).
Most of these again are more relevant in a multiuser setting, but it's good to have good habits, and you may find yourself in such a setting at some point anyway.
I don't include `.` in my path, and you just get used to adding the `./` whenever you need it explicitly.

#### Aliases
These can save you a lot of time. One I picked up that I use constantly is `m`, which is an alias for `less`. Lots of people use common `ls` aliases like `l`, `ll`, `lrt`, or `lrtail`.

Here's a list of the common ones I have - these go into `/home/ubuntu/.bash_aliases` (if you use this filename, the default Ubuntu `.bashrc` will load them automatically):
```
eval `dircolors -b`
alias md=mkdir
alias rd=rmdir

alias ll='ls --show-control-chars --color=auto -F -l'
alias la='ls --show-control-chars -A'
alias l='ls --show-control-chars --color=auto -CF'
alias lrtail='ls --show-control-chars --color=auto -F -ltr | tail'
alias dir='ls --show-control-chars --color=auto --format=vertical'
alias vdir='ls --show-control-chars --color=auto --format=long'
alias m='less'
alias lrt='ls --show-control-chars --color=auto -F -ltr'
alias ..='cd ..'
alias ...='cd ..;cd ..'
alias pu='pushd'
alias po='popd'
alias pd='pushd'

alias rm='rm -i'
alias mv='mv -i'
alias cp='cp -i'

alias log='history > `logfilename.pl`'
```
I also usually make this file for the root user too.

I have a simple `logfilename.pl` script to set filenames - this can go into `/usr/local/bin`, make sure it's executable (this is also installed as part of the [SLC Closet](bioinformatics.md#SLC-Closet) scripts):
```
#!/usr/bin/perl
#
# log command history
#
use warnings;
use strict;

my $date = `date +%F`;
chomp $date;

my $i = 1;
while (-f "$date"."_$i.log") {
  $i++;
}

print "$date", "_$i.log";
```

With these changes, you can save a command log by just typing `log` before actually logging off.

#### Editors
The editor is one of the most commonly used programs when working on a command line.
This is relatively easy to set and you should do it - if you're administering your own machine, I recommend setting this up for both your default user (ubuntu) and the root account.
Commands below are for the ubuntu user only.
```
echo "# Editor" >> /home/ubuntu/.bashrc
# Pick one of the following three lines - or modify accordingly if you use another one
# note that the stock Ubuntu image doesn't have emacs installed
# that can be installed with "sudo apt install emacs"
echo "export EDITOR=/usr/bin/vim"
echo "export EDITOR=/usr/bin/emacs"
echo "export EDITOR=/usr/bin/nano"
```
Then also be sure to set your own preferences, whether that's in `/home/ubuntu/.vimrc`, `/home/ubuntu/.emacs`, `/home/ubuntu/.nanorc`, or the appropriate other configuration file for your editor of choice.

#### Screen
This is a hugely useful tool for working on a remote machine.
It solves lots of problems with connection timeouts on long-running processes and provides many other features.
I highly recommend you learn to use this.
Basically, as soon as I start on a new machine, I run `screen` first, and then every time I log in it's a simple `screen -dr` and I'm back to where I left off.

I refer you to the [screen](https://www.gnu.org/software/screen/) [documentation](https://www.gnu.org/software/screen/manual/) and [these](https://opensource.com/article/17/3/introduction-gnu-screen) [useful](https://linuxize.com/post/how-to-use-linux-screen/) [primers](https://tpaschalis.github.io/gnu-screen-primer/#:~:text=GNU%20Screen%20is%20a%20command,processes%20independently%20from%20each%20other.).

#### Non-interactive shells
This can be a big "gotcha" when you're trying to automate tasks. You set everything up, test it all while logged in, make an AMI, and then launch it with some userdata to run your scripts / pipeline. And it doesn't work.

One of the subtleties is that (by default) a different environment is set when you are in an interactive vs. noninteractive shell. For more information, here are some [links](https://tldp.org/LDP/abs/html/intandnonint.html).

Again, on a machine that you fully control, and that you're going to replicate, I find it easier to just remove this layer of complexity. I only do this for the ubuntu user, as some system / maintenance scripts are designed to run non-interactively.

On the default Ubuntu 20.04 (Focal) AMI at AWS (which this whole site is based on), there are just a few lines at the top of the default `.bashrc` that detect whether you're on an interactive shell and shorcuts out if you are not - these are lines 5-9 in the version I'm looking at:
```
# If not running interactively, don't do anything
case $- in
    *i*) ;;
      *) return;;
esac
```
You can just comment these out so that there is no difference between interactive and non-interactive shells.

(Brief note, sometimes you'll see this on machines where you can log in and run something successfully, but if you do an `ssh user@host -c 'command'` then that will fail. One of the first things to check is the environment, and specifically your `PATH` or `LD_LIBRARY_PATH` or other variables.)

### SSH keys
The basics are well laid out [in](https://www.digitalocean.com/community/tutorials/how-to-set-up-ssh-keys--2) [other](https://upcloud.com/community/tutorials/use-ssh-keys-authentication/) [websites](https://www.ssh.com/academy/ssh/public-key-authentication).
Here are the common tasks, for a more full explanation check out the links in the previous sentence.
#### Create a new public/private key pair:
```
ssh-keygen
# for more options
ssh-keygen -f <output_keyfile>
# i.e. ssh-keygen -f mykeyfile
# this will generate two files:
# mykeyfile - this is the private key
# mykeyfile.put - this is the public key
chmod 0400 <output_keyfile>
```
#### Add a new public key to enable login
```
cat <output_keyfile.pub> >> /home/ubuntu/.ssh/authorized_keys
```

### Extra data storage
Data just keeps growing. You will likely need to add more storage at some point.
The AWS console has it's own instructions on how to create an EBS volume and attach it to a running instance.
Once it's attached to your machine, the common things are:
#### Figure out the device name
```
lsblk
# Look under "NAME" for things like sda, sda1, sdb, sdc, nvme0n1, nvme0n1p1
# remember that those are under /dev, i.e. /dev/sda, /dev/nvme0n1p1, etc.
# sda is canonically the whole disk (in physical terms), sda1 the partition
# similarly, nvme0n1 is like a disk, and nvmen1p1 a partition, though this
# loses some meaning when everything is virtual
```

#### Make a filesystem
Note this only applies to a new EBS volume. Don't do this to a volume restored from a snapshot (unless you really know what you're doing!).

Also, traditionally you would have to also partition your drive.
This is less relevant with the virtualization on the cloud, but you still can do it.
I won't cover it here, but the difference between devices like `/dev/sda` and `/dev/sda1` is that the former is the whole drive device, and the latter is a partition (and you could have `/dev/sda2`, `/dev/sda3`, etc.).
The newer machines on AWS now generally use devices named like `/dev/nvme0n1` (whole drive) and `/dev/nvme0n1p1` (partition).
```
# check if there's already a file system
# note that <device name> would be something like /dev/nvme0n1p1
sudo file -s <device name>
# if you see "<device name>: data" then there's no file system yet
# BE VERY CAREFUL YOU HAVE THE RIGHT DEVICE NAME HERE!!
# Otherwise you may blow away your root drive or (worse) another data drive!
sudo mkfs.ext4 <device name>
sudo file -s <device name>
# now you should see filesystem info
```

#### Mount a filesystem
The traditional place to put these things is in `/mnt`. It's good to plan ahead and use some sustainable / sensible names, so I call things `/mnt/volume1`, `/mnt/volume2`, etc.
```
ls -l /mnt
# if needed, make a directory for a mount point
sudo mkdir /mnt/volume1
sudo mount <device name> /mnt/volume1
# this would typically be "sudo mount /dev/nvme0n1p1 /mnt/volume1
# if this is the first time, you may need to give permission to the ubuntu user
sudo chown ubuntu:ubuntu /mnt/volume1
ls -l /mnt/volume1
```

#### Adding swap
In general the idea is you should spin up an instance with enough RAM for what you want to do, so the AWS instances by default have no swap space. However, sometimes you just need a little extra headroom for a temporary burst in memory usage (sometimes forking short processes can eat a lot of RAM, for example). Here's how to add some swap space.
```
# Many of these need root privileges
sudo su -
# You can change the file location and name (the of= option)
# You can also change the size - here it's 4 million blocks of 1024 bytes = 4GB
dd if=/dev/zero of=/tmp/swap0 bs=1024 count=4M
chmod 0400 /tmp/swap0
mkswap /tmp/swap0
swapon /tmp/swap0
```
