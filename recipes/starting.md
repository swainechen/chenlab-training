# Getting started

Everything on this website is based on the public AMI developed by the [Chen Lab](https://swainechen.github.io) at the National Univesrity of Singapore and the Genome Institute of Singapore.
You can use this on any AWS account.
The AMI is currently built in the `ap-southeast-1` region.
You can search for `CHENLAB-PUBLIC` in "Community AMIs" when launching an instance on AWS, click on the results in "Community AMIs", and pick the most recent version.

I recommend using a machine with 16 or 32GB of RAM if you're doing any assembly work (for example, an `r5a.large` or `r5a.xlarge` instance type on AWS).

A good thing to do (if you just have your own workstation - if you're automating, be careful about automatic updates for version control and potential new bugs) is to update things regularly:
```
# System updates
sudo apt update
sudo apt upgrade
sudo updatedb

# Perl
cpan-outdated -p | sudo cpanm

# local software, depending on what you have installed, some examples below
sudo su -
cd /usr/local/src/closet
git pull
cd /usr/local/src/kingfisher
git pull
# etc. etc.
```

Update locally installed R packages - unless you've configured user access to the system directories, probably easiest to run `R` as root, then:
```
# taken from https://www.r-bloggers.com/2014/11/update-all-user-installed-r-packages-again/
install.packages( 
    lib  = lib <- .libPaths()[1],
    pkgs = as.data.frame(installed.packages(lib), stringsAsFactors=FALSE)$Package,
    type = 'source'
)

# another alternative is rvcheck, which will also take care of BioConductor and devtools installs
install.packges("rvcheck")
rvcheck::update_all(check_R = TRUE, which = c("CRAN", "BioC", "github"))
```

Python packages seem a bit trickier, especially since we need some older versions and Python2.7 versions for specific bioinformatics software.
You should probably update these only when you have a good reason to do so.
