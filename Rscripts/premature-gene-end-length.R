### Extract gene-end, premature and other termination windows based on non-overlapping, 6kb genes

library(rtracklayer)
library(ggplot2)
library(dplyr)

p1 <- "/dysk2/groupFolders/magdak/t4-review/termination-windows-kinga-schlackow_F_hg38.bed"
p2 <- "/dysk2/groupFolders/magdak/t4-review/termination-windows-kinga-schlackow_R_hg38.bed"
p3 <- "/dysk2/groupFolders/magdak/analysis/subsets/hg38_genes6kb_basicannotation_v40.txt"

b1 <- import.bed(p1)
b2 <- import.bed(p2)
strand(b1) <- "+"
strand(b2) <- "-"

windows <- c(b1, b2)

df <- read.table(p3, header = TRUE, sep = "\t")
genes <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)

# Extract premature termination windows - subset windows within gene bodies
premature <- subsetByOverlaps(windows, genes, type = "within")
prem_F <- premature[strand(premature) == "+"]
prem_R <- premature[strand(premature) == "-"]

# Export bed files
export.bed(prem_F, "/dysk2/groupFolders/magdak/t4-review/windows-subsets-proper/premature_F.bed")
export.bed(prem_R, "/dysk2/groupFolders/magdak/t4-review/windows-subsets-proper/premature_R.bed")

# Function to expand range
expandRange = function(x, upstream=0, downstream=5000) {
  strand_is_minus = strand(x) == "-"
  on_plus = which(!strand_is_minus)
  on_minus = which(strand_is_minus)
  start(x)[on_plus] = start(x)[on_plus] - upstream
  start(x)[on_minus] = start(x)[on_minus] - downstream
  end(x)[on_plus] = end(x)[on_plus] + downstream
  end(x)[on_minus] = end(x)[on_minus] + upstream
  x
}

# Expand ranges 6kb downstream to extract gene-end terminators
genes_exp <- expandRange(genes, 0, 6000)

# Resize expanded genes to 6kb regions corresponding to artificial termination regions +6kb after annotated PAS
tregs <- resize(genes_exp, 6000, fix = "end")

# Extract gene-end termination windows as overlapping +6kb regions after annotated PAS
normal <- subsetByOverlaps(windows, tregs, type = "any")
normal_F <- normal[strand(normal) == "+"]
normal_R <- normal[strand(normal) == "-"]

# Export bed files
export.bed(normal_F, "/dysk2/groupFolders/magdak/t4-review/windows-subsets-proper/normal_F.bed")
export.bed(normal_R, "/dysk2/groupFolders/magdak/t4-review/windows-subsets-proper/normal_R.bed")

# Define other termination windows as windows that don't overlap with the extended genes 
other <- subsetByOverlaps(windows, genes_exp, invert = TRUE, type = "any")
other_F <- other[strand(other) == "+"]
other_R <- other[strand(other) == "-"]

# Export bed files
export.bed(other_F, "/dysk2/groupFolders/magdak/t4-review/windows-subsets-proper/other_F.bed")
export.bed(other_R, "/dysk2/groupFolders/magdak/t4-review/windows-subsets-proper/other_R.bed")


# Converr GRanges object to data frame to plot length of the windows
normal_df <- as.data.frame(normal)
normal_df$type = "normal"
premature_df <- as.data.frame(premature)
premature_df$type <- "premature"
other_df <- as.data.frame(other)
other_df$type <- "other"

# Merge dataframes to plot all termination windows together
df <- rbind(normal_df, premature_df, other_df)

# Plot
p <- df %>% group_by(type) %>% ggplot(aes(x = width))

p + geom_histogram(binwidth = 200) +
theme_minimal() + 
xlim(0, 15000) +
facet_grid(factor(type, levels = c("normal", "premature", "other")) ~ .)

df %>% group_by(type) %>% summarise(median = median(width), n = n())

#
#  type      median     n
#  <chr>      <dbl> <int>
#1 normal      2900  3566
#2 other       2106  8546
#3 premature   1559  1376

# Save the plot
ggsave("localization_termination_windows_normal_prem_other.pdf", device = "pdf", units = "cm", width = 23, height = 23, path = "/dysk2/groupFolders/magdak/t4-review/figs/")

# Plot only gene-end and premature termination windows
p <- df %>% filter(type != "other") %>% group_by(type) %>% ggplot(aes(x = width))

p + geom_histogram(binwidth = 200) +
theme_minimal() + 
xlim(0, 15000) +
facet_grid(type ~ .)

ggsave("localization_termination_windows_normal_prem.pdf", device = "pdf", units = "cm", width = 23, height = 17, path = "/dysk2/groupFolders/magdak/t4-review/figs/")

