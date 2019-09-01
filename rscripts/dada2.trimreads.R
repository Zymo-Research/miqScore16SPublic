# Title     : zrdocker-dada2:builderrormodels.R
# Objective : Takes in parameters, runs dada2 error model building functionality
# Created by: mweinstein
# Created on: 3/1/2019

library(getopt)
library(dada2)
optionSpecs = matrix(c(
    'R1', 'f', 1, "character", "Paired-end 1 reads fastq file",
    'R2', 'r', 1, "character", "Paired-end 2 reads fastq file",
    'R1_output', 'x', 1, "character", "Output file for filtered paired-end 1 reads",
    'R2_output', 'y', 1, "character", "Output file for filtered paired-end 2 reads",
    'R1_primer_length', 'a', 1, "numeric", "Length of primer for paired-end 1",
    'R2_primer_length', 'b', 1, "numeric", "Length of primer for paired-end 2",
    'R1_truncate_length', 'c', 1, "numeric", "Base position to truncate paired-end 1 reads",
    'R2_truncate_length', 'd', 1, "numeric", "Base position to truncate paired-end 2 reads",
    'truncateQuality', 'q', 1, "numeric", "Truncate read after a base of this quality score",
    'R1_maxEE', 'm', 1, "numeric", "Paired-end 1 maximum expected error",
    'R2_maxEE', 'n', 1, "numeric", "Paired-end 2 maximum expected error"
    ), byrow=TRUE, ncol=5
)

print("Checking trim R script options")
options = getopt(optionSpecs)

print("Starting dada2 trimmer function")
filterAndTrim(options$R1, options$R1_output, options$R2, options$R2_output, trimLeft=c(options$R1_primer_length, options$R2_primer_length), truncLen=c(options$R1_truncate_length, options$R2_truncate_length), maxEE=c(options$R1_maxEE, options$R2_maxEE), truncQ=options$truncateQuality, rm.phix=TRUE, compress=TRUE, verbose=TRUE, multithread=FALSE)
print("Completed dada2 trimmer function")
