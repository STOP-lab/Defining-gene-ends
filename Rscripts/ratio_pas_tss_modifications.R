library(rtracklayer)
library(BRGenomics)
library(plyranges)
library(tidyr)
library(ggplot2)
library(dplyr)

file_paths <- c(
  "/dysk2/groupFolders/magdak/t4-review/nascent-RNA/GSM2357382_total_mnetseq_F.bw_hg38.bw",
  "/dysk2/groupFolders/magdak/t4-review/nascent-RNA/GSM2357382_total_mnetseq_R.bw_hg38.bw",
  "/dysk2/groupFolders/magdak/t4-review/nascent-RNA/pol2mods/GSM2165964_s5ph-mnetseq_17098_02_F.bw_hg38.bw",
  "/dysk2/groupFolders/magdak/t4-review/nascent-RNA/pol2mods/GSM2165964_s5ph_mnetseq_17098_02_R.bw_hg38.bw",
  "/dysk2/groupFolders/magdak/t4-review/nascent-RNA/pol2mods/GSM2165966_s2ph_mnetseq_2080_04_F.bw_hg38.bw",
  "/dysk2/groupFolders/magdak/t4-review/nascent-RNA/pol2mods/GSM2165966_s2ph_mnetseq_2080_04_R.bw_hg38.bw",
  "/dysk2/groupFolders/magdak/t4-review/nascent-RNA/pol2mods/GSM2165968_y1ph_mnetseq_2080_02_F.bw_hg38.bw",
  "/dysk2/groupFolders/magdak/t4-review/nascent-RNA/pol2mods/GSM2165968_y1ph_mnetseq_2080_02_R.bw_hg38.bw",
  "/dysk2/groupFolders/magdak/t4-review/nascent-RNA/pol2mods/GSM2165958_s7ph_17098_05_F_hg38.bw",
  "/dysk2/groupFolders/magdak/t4-review/nascent-RNA/pol2mods/GSM2165958_s7ph_17098_05_R_hg38.bw",
  "/dysk2/groupFolders/magdak/t4-review/nascent-RNA/pol2mods/GSM2165961_t4ph_1904_11_F.bw_hg38.bw",
  "/dysk2/groupFolders/magdak/t4-review/nascent-RNA/pol2mods/GSM2165961_t4ph_1904_11_R.bw_hg38.bw"
)

# Import BigWig files
bigwig_list <- lapply(file_paths, import.bw)

# Iterate over the list of imported BigWig files
for (i in seq_along(bigwig_list)) {
  if (i %% 2 == 0) {  # Check if index is even
    bigwig_list[[i]] <- transform(bigwig_list[[i]], score = -1 * score)  # Multiply scores by -1 for even indexed files
    strand(bigwig_list[[i]]) <- "-"  # Set strand as "-"
  } else {
    strand(bigwig_list[[i]]) <- "+"  # Set strand as "+" for odd indexed files
  }
}

df <- read.table("/dysk2/groupFolders/magdak/analysis/subsets/hg38_genes6kb_basicannotation_v40.txt", sep = "\t", header = TRUE)

genes <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)

total <- c(bigwig_list[[1]], bigwig_list[[2]])
s5 <- c(bigwig_list[[3]], bigwig_list[[4]])
s2 <- c(bigwig_list[[5]], bigwig_list[[6]])
y1 <- c(bigwig_list[[7]], bigwig_list[[8]])
s7 <- c(bigwig_list[[9]], bigwig_list[[10]])
t4 <- c(bigwig_list[[11]], bigwig_list[[12]])

tss <- resize(genes, width = 200, fix = "start", ignore.strand = FALSE)
pas <- resize(genes, width = 200, fix ="end", ignore.strand = FALSE)


tss_total <- getCountsByRegions(total, tss)
pas_total <- getCountsByRegions(total, pas)
tss_s5 <- getCountsByRegions(s5, tss)
pas_s5 <- getCountsByRegions(s5, pas)
tss_s2 <- getCountsByRegions(s2, tss)
pas_s2 <- getCountsByRegions(s2, pas)
tss_y1 <- getCountsByRegions(y1, tss)
pas_y1 <- getCountsByRegions(y1, pas)
tss_s7 <- getCountsByRegions(s7, tss)
pas_s7 <- getCountsByRegions(s7, pas)
tss_t4 <- getCountsByRegions(t4, tss)
pas_t4 <- getCountsByRegions(t4, pas)

total_ratio <- pas_total / tss_total
s5_ratio <- pas_s5 / tss_s5
s2_ratio <- pas_s2 / tss_s2
y1_ratio <- pas_y1 / tss_y1
s7_ratio <- pas_s7 / tss_s7
t4_ratio <- pas_t4 / tss_t4


genes$total <- total_ratio
genes$s5 <- s5_ratio
genes$s2 <- s2_ratio
genes$y1 <- y1_ratio
genes$s7 <- s7_ratio
genes$t4 <- t4_ratio

genes_df <- as.data.frame(genes)
genes_short <- genes_df %>% select(-c(seqnames, start, end, strand, width, pointseq_avg))

genes_mods <- pivot_longer(genes_short, cols = c("total", "s5", "s2", "y1", "s7", "t4"), names_to = "modification", values_to = "ratio")
genes_mods$modification <- as.factor(genes_mods$modification)
genes_mods$ratio_log10 <- log(genes_mods$ratio)
genes_mods$ratio_log2 <- log2(genes_mods$ratio)

order <- c("total", "s5", "s2", "y1", "s7", "t4")
p <- genes_mods %>% ggplot(aes(x = factor(modification, level = order), y = ratio_log2, fill = modification))

pp <- p + 
geom_boxplot(width = 0.5, outlier.alpha = 0.2) +
theme_minimal()

ggsave("ratio_mods_pas_to_tss_withs7.pdf", plot = pp, device = "pdf", path = "/dysk2/groupFolders/magdak/t4-review/figs", width = 22, height = 14, units = "cm")

genes_mods %>% group_by(modification) %>% summarize(median = median(ratio_log2, na.rm = TRUE), n = n())

  modification median     n
  <fct>         <dbl> <int>
1 s2            -1.37 10387
2 s5            -2.62 10387
3 s7            -2.95 10387
4 t4             2.58 10387
5 total         -5.81 10387
6 y1            -5.23 10387



genes_mods %>% group_by(modification) %>% summarize(median = median(ratio, na.rm = TRUE), n = n())


  modification median     n
  <fct>         <dbl> <int>
1 s2           0.388  10387
2 s5           0.163  10387
3 s7           0.130  10387
4 t4           6.00   10387
5 total        0.0179 10387
6 y1           0.0266 10387















