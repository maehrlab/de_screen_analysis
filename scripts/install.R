setwd("PATH/TO/THIS/FOLDER")
# You need renv first if you don't have it yet
if(!requireNamespace("renv", quietly = T)) install.packages("renv")
# Isolate this library from the system library
renv::activate()
# Get some fundamental tools
install.packages("remotes")
install.packages("versions")
# It works better to use new ggplot2 with old R code than old ggplot2 with new R language versions.
install.packages("ggplot2")
# Install the necessary packages as they were towards the end of this project
library(versions)
versions::install.dates(
  pkgs = c("Seurat", "devtools", "tidyverse", "quantreg"),
  dates = "2019-01-01"
)
# Install packages from this project
remotes::install_github("maehrlab/thymusatlastools2")
remotes::install_github("ekernf01/freezr")
