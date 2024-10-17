### 1. Create average files from 3 biological replicates from GSE123105 and 2 biological replicates from GSE81662 
wiggletools mean GSM2165961_WTCHG_1904_11_F.bw GSM2357387_WTCHG_31850_03_F.bw > siLUC_FS_avg.bw
wiggletools mean GSM3495849_T4ph_mNETseq_siLUC_rep1_Fwd.bw GSM3495850_T4ph_mNETseq_siLUC_rep2_Fwd.bw GSM3495851_T4ph_mNETseq_siLUC_rep3_Fwd.bw > siLUC_Favg.bw
# Repeat the same procedure for bigwigs from the reverse strand


### 2. Multiply reverse files by -1, and files from GSE81662 by 0.01 (to have the same normalization)
# Use the R script "multiply_bigwig_by_given_number_R":
# usage multiply_bigwig_by_given_number.R path-to-bw-file.bw scaling-factor
Rscript multiply_bigwig_by_given_number.R siLUC_FS_avg.bw 0.01
Rscript multiply_bigwig_by_given_number.R siLUC_RS_avg.bw -0.01
Rscript multiply_bigwig_by_given_number.R siLUC_Ravg.bw -1


### 3. Resize to 150bp 
# Use the R script "resize_bigwig_150.R" for 4 files (two papers, forward and reverse strand)
for file in *avg*; do Rscript resize_bigwig_150.R ${file} & done


### 4. Lift the resized files to hg38 and convert to bedGraphs
# Download the chain file : https://pythonhosted.org/CrossMap/
for file in *resized*; do CrossMap.py bigwig path/to/lift/chain/file/hg19ToHg38.over.chain.gz ${file} ${file%.*}_hg38 & done 
for file in *hg38*; do bigWigToBedGraph ${file} ${file%.*}.bedGraph & done


### 5. Find the enrichment windows
for file in *bedGraph; do macs2 bdgbroadcall -i ${file} -o /magdak/t4-review/t4ph-mnetseq/${file%.*}.bed & done

### 4. Find overlapping ranges and merge
# Use the R script "overlap_and_merge.R" : Usage: overlap_and_merge.R <input_file1.bed> <input_file2.bed> <output_file.bed>
Rscript overlap_and_merge.R siLUC_FS_avg.bed siLUC_Favg.bed overlap_merge_F.bed
