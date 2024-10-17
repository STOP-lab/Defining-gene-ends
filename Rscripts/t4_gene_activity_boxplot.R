library(rtracklayer)
library(BRGenomics)
library(plyranges)
library(fabricatr)

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

path1 <- "/dysk2/groupFolders/magdak/t4-review/nascent-RNA/GSM2357382_total_mnetseq_F.bw_hg38.bw"
path2 <- "/dysk2/groupFolders/magdak/t4-review/nascent-RNA/GSM2357382_total_mnetseq_R.bw_hg38.bw"
path3 <- "/dysk2/groupFolders/magdak/t4-review/nascent-RNA/pol2mods/GSM2165961_t4ph_1904_11_F.bw_hg38.bw"
path4 <- "/dysk2/groupFolders/magdak/t4-review/nascent-RNA/pol2mods/GSM2165961_t4ph_1904_11_R.bw_hg38.bw"


# Import mNETseq total files and merge Fwd and Rev
total_F <- import.bw(path1)
strand(total_F) <- "+"
total_R <- import.bw(path2)
strand(total_R) <- "-"
total_R$score <- -1*total_R$score
total <- c(total_F, total_R)

# Import mNETseq T4ph file and merge Fwd and Rev
T4_F <- import.bw(path3)
strand(T4_F) <- "+"
T4_R <- import.bw(path4)
strand(T4_R) <- "-"
T4_R$score <- -1*T4_R$score
T4_mnetseq <- c(T4_F, T4_R)


# Load selected genes (see - "select-pc-genes.R" ) 
anno_g <- read.table("/dysk2/groupFolders/magdak/analysis/subsets/hg38_genes6kb_basicannotation_v40.txt", sep = "\t", header = TRUE)
annotated_genes <- makeGRangesFromDataFrame(anno_g, keep.extra.columns = TRUE)

# Expand genes and calculate their length
ext <- expandRange(annotated_genes, 0, 5000)
l <- width(ext)

# Sum of total & T4ph mnet-seq signal within extended genes
mnetseq_total <- getCountsByRegions(total, ext) / l
T4 <- getCountsByRegions(T4_mnetseq, ext) / l

# Make data frame
data_pc <- data.frame(chrom = seqnames(annotated_genes), pos = ranges(annotated_genes), strand = strand(annotated_genes), name = annotated_genes$gene_name, pol_total = mnetseq_total, T4ph = T4)

# Divide genes into 5 groups based on Pol2 total mNETseq signal
data_pc$activity <- split_quantile(data_pc$pol_total, type = 5)
sum(data_pc$activity == "5")


# Boxplot gene activity vs T4ph
box <- data_pc %>% ggplot(aes(x = activity, y = T4ph, fill = activity)) +
  geom_boxplot(width = 0.5, outlier.alpha = 0.2) +
  ylim(0,1) +
  labs(x = "Gene activity",
       y = "Average Pol II T4ph mNET-seq signal (RPM)",
       title = "Pol II T4ph signal vs gene activity") +
  scale_x_discrete(labels = c("No activity", "Very low", "Low", "Medium", "High")) +
  scale_fill_brewer(palette = "Blues") +
  guides(fill = "none")  +
  #scale_fill_viridis_d() +
  theme_minimal() +
  theme(legend.position = "bottom", 
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 20),
        #plot.subtitle = element_text(size = 18),
        axis.text.y  = element_text(size = 14),
        axis.text.x = element_text(size = 18),
        strip.text = element_text(size = 18),
        legend.title = element_blank())
 
# Save the plot  
ggsave("t4_vs_activity_boxplot.pdf", plot = box, device = "pdf", path = "/dysk2/groupFolders/magdak/t4-review/figs", width = 22, height = 14, units = "cm")

