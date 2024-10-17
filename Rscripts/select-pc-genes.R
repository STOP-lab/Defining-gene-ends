# Script to extract from Gencode annotation genes that:
# - are protein-coding
# - have the closest downstream neighbor localized +6kb from their 3' end 
# - don't overlap another annotated transcript on the same strand

library(plyranges)
library(rtracklayer)
library(dplyr)

# Upload annotation file
anno <- readGFF("/dysk2/groupFolders/magdak/annotation/gencode.v40.basic.annotation.gtf.gz")

# Extract all genes from annotation
anno_g <- anno[which(anno$type == "gene"),]
rm(anno)

# Make GRanges from GTF file
annotated_genes <- makeGRangesFromDataFrame(anno_g, keep.extra.columns = TRUE)
# Extract genes that have another annotated gene on the same strand
genes_overlapping <- annotated_genes[countOverlaps(annotated_genes, annotated_genes) > 1L]

# Extract downstream genes and remove NAs
ds <- precede(annotated_genes)
RM <- which(is.na(ds))
ds <- ds[-RM]
downstream <- annotated_genes[ds]
annotated_genes <- annotated_genes[-RM]

# Extract genes that have downstream neighbors localized at least 6 kb downstream of the annotated gene end
dd <- distance(annotated_genes, downstream)
g6kb <- annotated_genes[which(dd >= 6000)]

# Extract genes that don't overlap another annotated transcript at the same strand
g6kb_no <- g6kb[countOverlaps(g6kb, genes_overlapping) <= 1L ] 

# Extract protein-coding genes
pc_genes <- g6kb_no[g6kb_no$gene_type == "protein_coding",]
pc_genes_f <- pc_genes[strand(pc_genes) == "+"]
pc_genes_r <- pc_genes[strand(pc_genes) == "-"]
pc_genes_r <- 

# Export as .txt / .bed
export.bed(pc_genes_f, "/dysk2/groupFolders/magdak/t4-review/pc_genes_F.bed")
export.bed(pc_genes_r, "/dysk2/groupFolders/magdak/t4-review/pc_genes_R.bed")

