
# Defining gene ends: RNA polymerase II CTD threonine 4 phosphorylation marks transcription termination regions genome-wide

This repository contains code used for the analyses in the paper 'Defining gene ends: RNA polymerase II CTD threonine 4 phosphorylation marks transcription termination regions genome-wide' published in Nucleic Acids Research on xxx.

## Required tools
 - [CrossMap.py](https://crossmap.readthedocs.io/en/latest/)
 - [wiggletools](https://github.com/Ensembl/WiggleTools)
 - [macs2](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html)
 - [deeptools](https://deeptools.readthedocs.io/en/develop/)
 - bedGraphToBigWig and bigWigToBedGraph from [kentUtils](https://github.com/ENCODE-DCC/kentUtils)


## Termination windows definition

First, we selected genes that:
- are protein-coding,
- do not have any annotated trasncript on the same strand,
- their 3' end is isolated from a downstream neighbor at least 6kb.   
(see [`select-pc-genes.R`](https://github.com/STOP-lab/T4ph-review/blob/main/Rscripts/select-pc-genes.R))

To define termination windows follow [`definition-of-windows.sh`](https://github.com/STOP-lab/T4ph-review/blob/main/workflows/windows-definition.sh).   
All needed Rscripts are [here](https://github.com/STOP-lab/T4ph-review/tree/main/Rscripts).     

We defined premature and gene-end termination windows as described [here](https://github.com/STOP-lab/T4ph-review/blob/main/Rscripts/premature-gene-end-length.R).

## Metaplots

All metaplots were plotted using [deeptools](https://deeptools.readthedocs.io/en/develop/).     
To plot T4ph ChIP-seq metaplots follow [`metaplots-t4ph-chip.sh`](https://github.com/STOP-lab/T4ph-review/blob/main/workflows/metaplots-t4ph-chip.sh).    
To plot nascent RNA signal over genes/termination windows follow [`metaplots-nascentRNA.sh`](https://github.com/STOP-lab/T4ph-review/blob/).    

## Other analyses 

To calculate PAS/TSS ratio of Pol II modifications follow [`ratio_pas_tss_modifications.R`](https://github.com/STOP-lab/T4ph-review/blob/main/Rscripts/ratio_pas_tss_modifications.R).    
To plot average T4ph mNET-seq signal based on genes activity follow [`t4_gene_activity_boxplot.R`](https://github.com/STOP-lab/T4ph-review/blob/main/Rscripts/t4_gene_activity_boxplot.R).    
To plot termination windows length follow [`premature-gene-end-length.R`](https://github.com/STOP-lab/T4ph-review/blob/main/Rscripts/premature-gene-end-length.R).      



