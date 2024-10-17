# 1. Plot nascent RNA files over genes (lifted to hg38 by CrossMap.py)

bigwig1="/dysk2/groupFolders/magdak/t4-review/nascent-RNA/scaled/GSM4826609_F_POINTseq_HeLa_Untr_Rep2.bw"
bigwig2="/dysk2/groupFolders/magdak/t4-review/nascent-RNA/scaled/GSM7119959_L_EGFP_rep1_tt_corr_ff_noJncReads_plus.bw"
bigwig3="/dysk2/groupFolders/magdak/t4-review/nascent-RNA/scaled/GSM2692352_PRO-Seq-DMSO-FW.bw.hg38.bw"
bigwig4="/dysk2/groupFolders/magdak/t4-review/nascent-RNA/scaled/GSM2357382_total_mnetseq_F.bw_hg38.bw"
bigwig5="/dysk2/groupFolders/magdak/t4-review/nascent-RNA/scaled/multiplied_to_genes/GSM2165961_t4ph_1904_11_F.bw_hg38.bw"
region1="/dysk2/groupFolders/magdak/analysis/subsets/hg38_genes6kb_basicannotation_v40_F.bed"
output="/dysk2/groupFolders/magdak/t4-review/matrix/unscaled_4techniques_genes_F.gz"
computeMatrix scale-regions -S ${bigwig1} ${bigwig2} ${bigwig3} ${bigwig4} ${bigwig5} \
                              -R ${region1} \
                              --beforeRegionStartLength 5000 \
                              --regionBodyLength 10000 \
                              --afterRegionStartLength 5000 \
                              --missingDataAsZero -o ${output} \
                              --startLabel TSS \
                              --endLabel PAS -p max/2 


matrix="/dysk2/groupFolders/magdak/t4-review/matrix/unscaled_4techniques_genes_F.gz"
plotData="/dysk2/groupFolders/magdak/t4-review/figs/plotData/unscaled_4techniques_genes_F.tab"
output="/dysk2/groupFolders/magdak/t4-review/figs/unscaled_4techniques_genes_F.pdf"
plotProfile -m ${matrix} -out ${output} --outFileNameData ${plotData} --plotTitle " " --startLabel TSS --endLabel PAS --regionsLabel "genes " --perGroup --samplesLabel  "POINT-seq" "TT-seq" "PRO-seq" "Total Pol II mNET-seq" "T4ph mNET-seq"


# 2. Open .tab file and find the max value for each bigwig (did it manually in the spreadsheet)

# 3. Calculate the scaling factor = 1/max
#					max			1/max
#POINT-seq			38.1162148	0.02623555
#TT-seq				19.8261267	0.0504385
#PRO-seq 			15.1559241	0.0659808
#Total Pol mNET-seq	1.76232479	0.56743229
#T4ph mNET-seq		0.54521299	1.8341456

# 4. Use multiply_bigwig_by_given_number.R Rscript to get scaled bigwigs, eg:

Rscript multiply_bigwig_by_given_number.R GSM4826609_F_POINTseq_HeLa_Untr_Rep2.bw 0.02623555

# 5. Once again calculate matrix on multiplied bigwigs and plot

bigwig1="/dysk2/groupFolders/magdak/t4-review/nascent-RNA/scaled/multiplied_to_genes/GSM4826609_F_POINTseq_HeLa_Untr_Rep2.bw_multiplied.bw"
bigwig2="/dysk2/groupFolders/magdak/t4-review/nascent-RNA/scaled/multiplied_to_genes/GSM7119959_L_EGFP_rep1_tt_corr_ff_noJncReads_plus.bw_multiplied.bw"
bigwig3="/dysk2/groupFolders/magdak/t4-review/nascent-RNA/scaled/multiplied_to_genes/GSM2692352_PRO-Seq-DMSO-FW.bw.hg38.bw_multiplied.bw"
bigwig4="/dysk2/groupFolders/magdak/t4-review/nascent-RNA/scaled/multiplied_to_genes/GSM2357382_total_mnetseq_F.bw_hg38.bw_multiplied.bw"
bigwig5="/dysk2/groupFolders/magdak/t4-review/nascent-RNA/scaled/multiplied_to_genes/GSM2165961_t4ph_1904_11_F.bw_hg38.bw_multiplied.bw"
region1="/dysk2/groupFolders/magdak/analysis/subsets/hg38_genes6kb_basicannotation_v40_F.bed"
output="/dysk2/groupFolders/magdak/t4-review/matrix/scaledto1_4tech_t4_genes_F.gz"
computeMatrix scale-regions -S ${bigwig1} ${bigwig2} ${bigwig3} ${bigwig4} ${bigwig5} \
                              -R ${region1} \
                              --beforeRegionStartLength 5000 \
                              --regionBodyLength 10000 \
                              --afterRegionStartLength 5000 \
                              --missingDataAsZero -o ${output} \
                              --startLabel TSS \
                              --endLabel PAS -p max/2 


matrix="/dysk2/groupFolders/magdak/t4-review/matrix/scaledto1_4tech_t4_genes_F.gz"
output="/dysk2/groupFolders/magdak/t4-review/figs/scaledto1_4techniques_t4_genes_F.pdf"
outdata="/dysk2/groupFolders/magdak/t4-review/figs/plotData/scaledto1_4techniques_t4_genes_F.tab"
plotProfile -m ${matrix} -out ${output} --outFileNameData ${outdata} --plotTitle " " --startLabel TSS --endLabel PAS --regionsLabel "Genes" --perGroup --samplesLabel  "POINT-seq" "TT-seq" "PRO-seq" "Total Pol II mNET-seq" "T4ph mNET-seq"

# 6. The same workflow was used to plot profiles over gene-end and premature termination windows. To plot Pol II modifications instead of plotProfile, plotHeatmap function was used. 

