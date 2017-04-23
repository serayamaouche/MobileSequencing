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

# Install the minionSummaryData dataset
biocLite("minionSummaryData")
library(minionSummaryData)

# Read summary data example
dataDir <- "M:/Professional/MyCode/MobileSequencing/data"
fast5files <- list.files(dataDir, pattern = '.fast5$')
setwd(dataDir)

# The function readFast5Summary() reads one or more fast5 files and 
# collects summary information about them.
example.summary <- readFast5Summary(fast5files)
str(example.summary)

# Load data from the minionSummaryData package
library(minionSummaryData)
data(s.typhi.rep1)
s.typhi.rep1

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

# Extract 2D reads
if( require(minionSummaryData) ) {
  data(s.typhi.rep2, package = 'minionSummaryData')
  fastq2D( s.typhi.rep2 )
}

# Extract complement reads
data(s.typhi.rep2, package = 'minionSummaryData')
fastqComplement( s.typhi.rep2 )

# Extract baseCalled slot
# Accesses the baseCalled slot stored in an object derived from the Fast5Summary
# class using the function baseCalled().
call <- baseCalled(s.typhi.rep1)
call
baseCalled(s.typhi.rep1[1:2])

# Plot the accumulation of reads over the duration of the experiment.
data(s.typhi.rep2, package = 'minionSummaryData')
plotReadAccumulation( s.typhi.rep2 )

# Plot the proportion of template, complement and 2D reads found a
# dataset using plotReadCategoryCount()
# and Visualise the mean base quality of each read using the
# plotReadCategoryQuals() function
p1 <- plotReadCategoryCounts(s.typhi.rep1)
p2 <- plotReadCategoryQuals(s.typhi.rep1)
grid.arrange(p1, p2, ncol = 2)

# Plot the number of active channels for each minute of run time
p2 <- plotActiveChannels(s.typhi.rep1)
grid.arrange(p1, ncol = 1)

# View changes in signal against run time using plotReadTypeProduction()
data(s.typhi.rep2, package = 'minionSummaryData')
p1 <- plotReadTypeProduction(s.typhi.rep1)
p2 <- plotReadTypeProduction(s.typhi.rep2)
grid.arrange(p1, p2, ncol = 2)

# Plot the mean rate at which events occur
data(s.typhi.rep2, package = 'minionSummaryData')
p1 <- plotEventRate(s.typhi.rep1)
p2 <- plotEventRate(s.typhi.rep2)
grid.arrange(p1, p2, ncol = 2)

# channelActivityPlot() can be used to visualise a specified metric over all channels over time.
require(dplyr)
data(s.typhi.rep3, package = 'minionSummaryData')
## we will plot the median raw signal for each read on z-axis
z_scale = select(rawData(s.typhi.rep2), id, median_signal)
channelActivityPlot( s.typhi.rep2, zScale = z_scale )


# Create layout plot of flowcell
p1 <- layoutPlot(s.typhi.rep1, attribute = "nreads")
p2 <- layoutPlot(s.typhi.rep1, attribute = "kb")
grid.arrange(p1, p2, ncol = 2)

# heatmap example
library(dplyr)
read_count_2D <- baseCalled(s.typhi.rep1) %>% ## start with base called reads
  filter(strand == 'template') %>% ## keep template so we don't count things twice
  left_join(readInfo(s.typhi.rep1), by = 'id') %>% ## channel stored in @readInfo slot, match by id column
  group_by(channel) %>% ## group according to channel
  summarise(d2_count = length(which(full_2D == TRUE)), ## count those with full 2D status
  d2_prop = length(which(full_2D == TRUE)) / n()) ## divide by total count of reads from channel

p1 <- channelHeatmap(read_count_2D, zValue = 'd2_count')
p2 <- channelHeatmap(read_count_2D, zValue = 'd2_prop')
grid.arrange(p1, p2, ncol = 2)

#library(dplyr)
if( require(minionSummaryData) ) {
  data(s.typhi.rep2, package = 'minionSummaryData')
  ## calculate and plot the mean number of events recorded by each channel
  avgEvents <- left_join(readInfo(s.typhi.rep2), rawData(s.typhi.rep2), by = 'id') %>%
    group_by(channel) %>%
    summarise(mean_nevents = mean(num_events))
  channelHeatmap(avgEvents, zValue = 'mean_nevents')
}

# Display correlation between pentemer proportions in two time windows
plotKmerFrequencyCorrelation( s.typhi.rep1, only2D = FALSE )

# Write FASTQ-formatted file using the ShortRead package
# library(ShortRead)
# The fastq object is accessed using fastq() function
writeFastq( fastq( s.typhi.rep1 ), file = tempfile() )
writeFastq( fastq2D( s.typhi.rep1 ), file = tempfile() 


