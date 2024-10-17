#!/usr/bin/env Rscript

# usage /multiply_bigwig_by_given_number.R path-to-bw-file.bw scaling-factor


# Load required libraries
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)

# Get the path to the BigWig file as an argument
args <- commandArgs(trailingOnly = TRUE)
bigwig_file <- args[1]
scaling_factor <- as.numeric(args[2])

# Load the BigWig file
bw <- import.bw(bigwig_file)

# Multiply the scores by given scaling factor
bw$score <- scaling_factor * bw$score

# Create a filename for the output file
output_file <- paste0(basename(bigwig_file), "_multiplied.bw")

# Export the modified BigWig file to the current directory
export.bw(bw, output_file)

# Print a message indicating the location of the exported file
cat("Multiplied BigWig file exported to:", getwd(), "/", output_file, "\n")
