## README

This folder contains scripts used for analysis of single-cell RNA-seq data in our recent definitive endoderm screen. (For atacTFAP analysis, there is a separate repo.)

You can find the data at [GEO: GSE127202](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127202). You will need a couple of key R packages: `Seurat` version 2 (from the Satija lab) and `thymusatlastools2` available [here](https://github.com/maehrlab/thymusatlastools2). We also use [`freezr`](https://github.com/ekernf01/freezr) to save code and session info and to track processed data for use in downstream scripts. To install the correct versions, try the installation script provided in `scripts/install.R`.

`main.R` is the master script for the paper draft. To reproduce the paper, 

- in `main.R`, edit the variable `proj_dir` to point to this repo, 
- in `clustering/setup.Rmd`, edit the lines that read in the data. 
- in `clustering/integrate_gAmp.Rmd`, edit the lines that read in the data. 

Then run `main.R` line by line. 

We're passionate about reproducibility, but we understand that exact reproduction across environments is very difficult. If you have trouble, please file an issue on [github](https://github.com/maehrlab/de_screen_analysis).
