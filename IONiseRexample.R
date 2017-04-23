#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Script File: IONiseRexample.R
# Date of creation: 5 Feb 2017
# Date of last modification: 23 Apr 2017
# Author: Seraya Maouche <seraya.maouche@iscb.org>
# Short Description: This script provides functions to QC mobile sequecing data from 
#                     Oxford Nanopore MinION
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Install IONiseR Bioconductor package
# http://bioconductor.org/packages/release/bioc/html/IONiseR.html
source("https://bioconductor.org/biocLite.R")
biocLite("IONiseR")

# Load IONiseR and required packages
library(tibble)
library(hwriter)
library(rhdf5)
library(ShortRead)
library(IONiseR)
library(ggplot2)
library(gridExtra)


# Read data
dataDir <- "M:/Professional/MyCode/MobileSequencing/data"
fast5files <- list.files(dataDir, pattern = '.fast5$')
setwd(dataDir)
example.summary <- readFast5Summary(fast5files)
str(example.summary)
