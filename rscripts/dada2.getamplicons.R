# # Title     : zrdocker-dada2:getamplicons.R
# # Objective : Takes in parameters, infers amplicon variants using filtered reads and an error model from previous functions
# # Created by: mweinstein
# # Created on: 3/1/2019

library(getopt)
library(dada2)

default_min_overlap <- 10
default_max_mismatch <- 2

optionSpecs = matrix(c(
    'R1', 'f', 1, "character", "Trimmed ead 1 fastq file",
    'R2', 'r', 1, "character", "Trimmed read 2 fastq file",
    'R1_error_model', 'x', 1, 'character', "Read 1 error model from previous step",
    'R2_error_model', 'y', 1, 'character', "Read 2 error model from previous step",
    'output_file', 'o', 1, 'character', "Output file path for sequence table in CSV",
    'output_taxa', 't', 1, 'character', "Output for taxa CSV",
    'chimera_free_output', 'k', 1, 'character', "Chimera-free sequence table output CSV",
    'output_sequence_table', 's', 1, 'character', "Output file for chimera-free sequence table in native format",
    'rdp_database', 'd', 1, 'character', "Microbial sequence reference database",
    'min_overlap', 'l', 2, 'integer', "Minimum overlap between paired ends",
    "max_mismatch", 'm', 2, 'integer', "Maximum mismatch between sequences that can be grouped"
    ), byrow=TRUE, ncol=5
)

print("Checking amplicon analysis R script options")
options = getopt(optionSpecs)

if (is.null(options$min_overlap)){
    min_overlap <- default_min_overlap
} else {
    min_overlap <- options$min_overlap
}

if (is.null(options$max_mismatch)){
    max_mismatch <- default_max_mismatch
} else {
    max_mismatch <- options$max_mismatch
}    

print("Starting dada2 amplicon analysis functions")

print("Reading in error models")
R1_error_model <- readRDS(file=options$R1_error_model)
R2_error_model <- readRDS(file=options$R2_error_model)
print("Error models loaded")
print("Beginning read dereplication")
R1_dereplicated <- derepFastq(options$R1)
R2_dereplicated <- derepFastq(options$R2)
print("Dereplication completed")
print("Starting divisive amplicon denoising algorithm")
R1_denoise <- dada(R1_dereplicated, err=R1_error_model, multithread=TRUE)
R2_denoise <- dada(R2_dereplicated, err=R2_error_model, multithread=TRUE)
print("Denoising completed")
print("Merging read pairs")
merged_reads <- mergePairs(R1_denoise, R1_dereplicated, R2_denoise, R2_dereplicated, minOverlap = min_overlap, maxMismatch = max_mismatch)
print("Done merging reads")
print("Making sequence table")
sequence_table <- makeSequenceTable(merged_reads)
print("Made sequence table")
print("Finding chimeras")
sequence_table_chimera_free <- removeBimeraDenovo(sequence_table, method="consensus", multithread=TRUE, verbose=TRUE)
print("Chimera removal completed")
print("Assigning taxa")
taxa <- assignTaxonomy(sequence_table_chimera_free, options$rdp_database, multithread=TRUE)
print("Assigned taxa")
print("Saving results")
saveRDS(sequence_table_chimera_free, options$output_sequence_table)
write.csv(sequence_table, options$output_file)
write.csv(sequence_table_chimera_free, options$chimera_free_output)
write.csv(taxa, options$output_taxa)
print("Results saved")
message("Completed dada2 amplicon analysis functions")