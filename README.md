
# Defining gene ends: RNA polymerase II CTD threonine 4 phosphorylation marks transcription termination regions genome-wide

This repository contains code used for the analyses in the paper 'Defining gene ends: RNA polymerase II CTD threonine 4 phosphorylation marks transcription termination regions genome-wide' published in Nucleic Acids Research on xxx.

## Required tools
 - [CrossMap.py](https://crossmap.readthedocs.io/en/latest/)
 - [wiggletools](https://github.com/Ensembl/WiggleTools)
 - [macs2](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html)
 - [deeptools](https://deeptools.readthedocs.io/en/develop/)


## Termination windows definition

First, we selected genes that:
- are protein-coding,
- do not have any annotated trasncript on the same strand,
- their 3' end is isolated from a downstream neighbor at least 6kb.
(see [`select-pc-genes.R`](https://github.com/STOP-lab/T4ph-review/blob/main/Rscripts/select-pc-genes.R))

To define termination windows follow [`definition-of-windows.sh`](https://github.com/STOP-lab/T4ph-review/blob/main/workflows/windows-definition.sh).
All needed Rscripts are [here](https://github.com/STOP-lab/T4ph-review/tree/main/Rscripts).


