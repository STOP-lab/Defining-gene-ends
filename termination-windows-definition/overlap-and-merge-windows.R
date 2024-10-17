#!/usr/bin/env Rscript

library(GenomicRanges)
library(rtracklayer)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: overlap_and_merge.R <input_file1.bed> <input_file2.bed> <output_file.bed>")
}

input_file1 <- args[1]
input_file2 <- args[2]
output_file <- args[3]

# Import bed files
file1 <- import.bed(input_file1)
file2 <- import.bed(input_file2)

# Find overlaps
overlaps <- findOverlaps(file1, file2)

# Extract the overlapping ranges
overlapping_ranges <- c(file1[queryHits(overlaps)], file2[subjectHits(overlaps)])

# Merge the overlapping ranges
merged_overlaps <- reduce(overlapping_ranges)

# Export 
export.bed(merged_overlaps, output_file)





