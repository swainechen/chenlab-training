# Setting up R and Bioconductor
R is one of the big engines for bioinformatics research these days. This will install the latest stable version of R along with the latest Bioconductor (as of May, 2022, this is R 4.2.0 with Bioconductor 3.15).

The default Ubuntu repositories for the LTS releases are not designed to keep up with the latest versions of software like R. So we'll follow the [advice of CRAN](https://cran.r-project.org/bin/linux/ubuntu/) and use their repositories.

```
sudo su -
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
apt install --no-install-recommends r-base
# then we'll run R as root
R
```

It's good to check where things are going to go first. Then do a basic initial update. This takes about 5 min on a t3a.small instance:
```
.libPaths()
update.packages(ask=F)
```

Install some regular R packages (this takes about 30 min):
```
install.packages(c("apcluster", "d3r", "devtools", "fields", "fpc",
  "genoPlotR", "gplots", "hexbin", "PopGenome", "r2d3", "randomcoloR",
  "RMySQL", "Rtsne", "rvcheck", "scatterD3", "Seurat", "shiny",
  "shinydashboard", "svglite", "tidyverse", "tsne"))
install_github("prabhakarlab/RCAv2")
```

Install Bioconductor (following the [standard instructions](https://www.bioconductor.org/install/))
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.15")
```

Then install some Bioconductor packages:
```
BiocManager::install(c("Biostrings", "ComplexHeatmap", "DESeq2", "edgeR",
  "ggtree", "IRanges", "SNPRelate", "SRAdb", "Sushi", "bambu", "ggbio"))
```
