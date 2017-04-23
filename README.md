
# Quality Assessment of mobile sequecing data from the Oxford Nanopore MinION, a portable and real time device for DNA and RNA sequencing.
By Seraya Maouche |
Latest update : 23 April 2017

The [IONiseR Bioconductor package](http://bioconductor.org/packages/release/bioc/html/IONiseR.html) provides functions for quality assessment (QC) of mobile sequencing experiments.
Data used in this example was initialy published by [Ashton et al (2015)](http://www.nature.com/nbt/journal/v33/n3/full/nbt.3103.html) in Nature Biotechnology. This Salmonella Typhi dataset is available as part of the Bioconductor data package [minionSummaryData](https://bioconductor.org/packages/release/data/experiment/html/minionSummaryData.html).


1- Install IONiseR and required packages
```R
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

# Install the minionSummaryData dataset
biocLite("minionSummaryData")
library(minionSummaryData)

```

2- Read summary data example
```R
dataDir <- "M:/Professional/MyCode/MobileSequencing/data"
fast5files <- list.files(dataDir, pattern = '.fast5$')
setwd(dataDir)

# ****** Extract data
#Extract readInfo slot
# The function readInfo() accesses the readInfo slot stored in an 
# object derived from the Fast5Summary class.
if(require(minionSummaryData) ) {
  data(s.typhi.rep2, package = 'minionSummaryData')
  readInfo( s.typhi.rep2 )
}

# Extract rawData slot
# Accesses the rawData slot stored in an object derived from 
# the Fast5Summary class.
data(s.typhi.rep2, package = 'minionSummaryData')
rawData( s.typhi.rep2 )

# Accesses the fastq slot stored in an object 
# derived from the Fast5Summary class.
data(s.typhi.rep2, package = 'minionSummaryData')
fastq( s.typhi.rep2 )

```

3 - Visualization of metrics about sequecing experiments
```R
# Plot the accumulation of reads over the duration of the experiment.
data(s.typhi.rep2, package = 'minionSummaryData')
plotReadAccumulation( s.typhi.rep2 )

``` 

<p align="center">
  <img src="" width=""/>
</p>


### Links

* [IONiseR Package on Bioconductor](http://bioconductor.org/packages/release/bioc/html/IONiseR.html)
* [The minionSummaryData package](https://bioconductor.org/packages/release/data/experiment/html/minionSummaryData.html)

    
