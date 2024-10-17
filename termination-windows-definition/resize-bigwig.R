#!/usr/bin/env Rscript

# Load required libraries
library(rtracklayer)

# Get the path to the bedGraph file as an argument
args <- commandArgs(trailingOnly = TRUE)
bigwig_file <- args[1]

# Load the bedGraph file
bw <- import.bw(bigwig_file)

# Resize bedGraph file
bw_resized <- resize(bw, width = 150, fix = "center")

# Get coverage file
cov <- coverage(bw_resized, weight = "score")
grcov <- GRanges(cov) 

# Create a filename for the output file
output_file <- paste0(getwd(), "/", basename(bigwig_file), "_resizedTo150.bw")

# Export the modified bedGraph file to the current directory
export(grcov, con = output_file, format = "bigWig")

# Print a message indicating the location of the exported file
cat("Modified bedGraph file exported to:", getwd(), "/", output_file, "\n")
