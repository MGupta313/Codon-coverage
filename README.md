# Calculating codon coverage and visualizing coverage distribution of four drug resistant genes of _Plasmodium falciparum_ for pooled and individual samples (Haiti)

The four genes which are studies here are:
- **PfDHFR:** dihydrofolate reductase
- **PfMDR1:** multidrug resistance gene
- **PfDHPS:** dihydropteroate synthase
- **PfCRT:** chloroquine resistance transporter

## Prerequistes

1. NeST (Next-generation Sequence-analysis Toolkit) output:
	1. Install NeST from https://github.com/CDCgov/NeST.
	2. Follow the guidelines to get _output_FM_SR_DD_RG.bam_ files for each gene for each sample type. This file is located in _NeST output folder > sample > alignments > output_FM_SR_DD_RG.bam_.
2. Samtools depth:
	1. [Samtools depth documentation](http://www.htslib.org/doc/samtools-depth.html)
	2. This is used to get nucleotide depth.
3. Reference mdr.bed file:
	1. This contains information about 6 genes of _Plasmodium falciparum_
	2. This is used to get the start and end positions of the exonic regions of the four genes.
4. Geneious Prime:
	1. Geneious Prime is a powerful bioinformatics software solution packed with fundamental molecular biology and sequence analysis tools. [Geneious Prime](https://www.geneious.com/)
	2. This is used to map sample sequences to reference gene sequence (CDS region) to get rid of the intronic regions, if any.
	3. Sample sequences mapped to reference gene sequence are always required to get start and end positions in cases of PfCRT and PfDHPS. They are sometimes required for PfDHFR and PfMDR1 as well.
	

