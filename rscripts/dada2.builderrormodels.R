# Title     : zrdocker-dada2:builderrormodels.R
# Objective : Takes in parameters, runs dada2 error model building functionality
# Created by: mweinstein
# Created on: 3/1/2019

library(getopt)
library(dada2)
optionSpecs = matrix(c(
    "input", 'i', 1, "character", "Filtered fastq from previous step",
    "output", 'o', 1, "character", "Error model data output file"
    ), byrow=TRUE, ncol=5
)

print("Checking error model building R script options")
options = getopt(optionSpecs)

print("Starting dada2 error model build function")
filteredReads <- readLines(options$input)
error_model <- learnErrors(filteredReads, multithread=TRUE)
print("Completed error model build, printing results")
saveRDS(error_model, options$output)
print("Completed dada2 error model build function")


