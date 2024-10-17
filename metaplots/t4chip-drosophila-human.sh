### Convert T4ph ChIP-seq from https://pubmed.ncbi.nlm.nih.gov/28782042/ to bedGraph and get Drosophila genes

R

library(rtracklayer)
library(AnnotationHub)


p1 <- "/dysk2/groupFolders/magdak/t4-review/drosophila-chip/Thr4P.BG3.mock.ave.log.enrichment.sgr"

t<-read.table(p1, sep = "\t", header = TRUE)
colnames(t) <- c("seqnames", "end", "score")
df <- data.frame(seqnames = t$seqnames, start = t$end - 1, end = t$end, score = t$score)
df$seqnames<-sub("^","chr",df$seqnames)
gr <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)

export.bedGraph(gr, "/dysk2/groupFolders/magdak/t4-review/drosophila-chip/Thr4P.mock.ave.log.enrichmet.bedGraph")


ah <- AnnotationHub()
dm3 <- ah[["AH52252"]]
genes <- genes(dm3)

export.bed(genes, "/dysk2/groupFolders/magdak/t4-review/drosophila-chip/drosophila_dm3_from_annohub.bed")

q()


### Convert bedGraph to bigWig
bedGraphToBigWig /dysk2/groupFolders/magdak/t4-review/drosophila-chip/Thr4P.mock.ave.log.enrichmet.bedGraph http://hgdownload.soe.ucsc.edu/goldenPath/dm3/bigZips/dm3.chrom.sizes /dysk2/groupFolders/magdak/t4-review/drosophila-chip/Thr4P.mock.ave.log.enrichmet.bw 



### Compute matrix and plot profile for T4ph ChIP-seq in Drosophila

bigwig1="/dysk2/groupFolders/magdak/t4-review/drosophila-chip/Thr4P.mock.ave.log.enrichmet.bw"
region1="/dysk2/groupFolders/magdak/t4-review/drosophila-chip/drosophila_dm3_from_annohub.bed"
output="/dysk2/groupFolders/magdak/t4-review/matrix/drosophila-t4ph_annohub.gz"
computeMatrix scale-regions -S ${bigwig1} \
                              -R ${region1} \
                              --beforeRegionStartLength 5000 \
                              --regionBodyLength 10000 \
                              --afterRegionStartLength 5000 \
                              --skipZeros -o ${output} \
                              --startLabel TSS \
                              --endLabel PAS -p max/2 

matrix="/dysk2/groupFolders/magdak/t4-review/matrix/drosophila-t4ph_annohub.gz"
output="/dysk2/groupFolders/magdak/t4-review/figs/drosophila-t4ph_annohub.pdf"
plotProfile -m ${matrix} -out ${output} --plotTitle "T4ph ChIP-seq signal" --startLabel TSS --endLabel PAS --regionsLabel "genes" --yMin 0 



### Compute matrix and plot profile for T4ph ChIP-seq in human (lifted to hg38 using CrossMap.py)

CrossMap.py bigwig /dysk2/groupFolders/magdak/liftchain/hg19ToHg38.over.chain.gz GSM920945_t4ph_chip_rep1_hg19.bw GSM920945_t4ph_chip_rep1_hg38

bigwig1="/dysk2/groupFolders/magdak/t4-review/chip/GSM920945_t4ph_chip_rep1_hg38.bw"
region1="/dysk2/groupFolders/magdak/analysis/subsets/hg38_genes6kb_basicannotation_v40.bed"
output="/dysk2/groupFolders/magdak/t4-review/matrix/t4ph-chip-hintermair.gz"
computeMatrix scale-regions -S ${bigwig1} \
                              -R ${region1} \
                              --beforeRegionStartLength 5000 \
                              --regionBodyLength 10000 \
                              --afterRegionStartLength 5000 \
                              --skipZeros -o ${output} \
                              --startLabel TSS \
                              --endLabel PAS -p max/2 


matrix="/dysk2/groupFolders/magdak/t4-review/matrix/t4ph-chip-hintermair.gz"
output="/dysk2/groupFolders/magdak/t4-review/figs/t4ph-chip-hintermair.pdf"
plotProfile -m ${matrix} -out ${output} --plotTitle "T4ph ChIP-seq signal" --startLabel TSS --endLabel PAS 



