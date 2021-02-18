# Calculating codon coverage and visualizing coverage distribution of four drug resistant genes of _Plasmodium falciparum_ for pooled and individual samples (Haiti)

The four genes which are studies here are:
- **PfDHFR:** dihydrofolate reductase
- **PfMDR1:** multidrug resistance gene
- **PfDHPS:** dihydropteroate synthase
- **PfCRT:** chloroquine resistance transporter

## Prerequistes

1. NeST (Next-generation Sequence-analysis Toolkit) output:
	1.1. Install NeST from https://github.com/CDCgov/NeST.
	1.2. Follow the guidelines to get _output_FM_SR_DD_RG.bam_ files for each gene for each sample type. This file is located in _NeST output folder > sample > alignments > output_FM_SR_DD_RG.bam_.
2. Samtools depth:
	1. http://www.htslib.org/doc/samtools-depth.html
	2. This is used to get nucleotide depth.
