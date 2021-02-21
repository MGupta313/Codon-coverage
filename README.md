# Calculating codon coverage and visualizing coverage distribution of four drug resistant genes of _Plasmodium falciparum_ for pooled and individual samples (Haiti)

The four genes which are studies here are:
- **PfDHFR:** dihydrofolate reductase
- **PfMDR1:** multidrug resistance gene
- **PfDHPS:** dihydropteroate synthase
- **PfCRT:** chloroquine resistance transporter

## Prerequistes

1. NeST (Next-generation Sequence-analysis Toolkit) output:
	1. Install NeST from https://github.com/CDCgov/NeST.
	2. Follow the guidelines to get _output_FM_SR_DD_RG.bam_ files for each gene for each sample type. (_NeST output folder > sample > alignments > output_FM_SR_DD_RG.bam_)
2. Samtools depth:
	1. [Samtools depth documentation](http://www.htslib.org/doc/samtools-depth.html)
	2. This is used to get nucleotide depth.
3. Reference mdr.bed file:
	1. This contains information about 6 genes of _Plasmodium falciparum_
	2. This is used to get the start and end positions of the exonic regions of the four genes.
4. Geneious Prime:
	1. [Geneious Prime](https://www.geneious.com/) is a powerful bioinformatics software solution packed with fundamental molecular biology and sequence analysis tools.
	2. This is used to map sample sequences to reference gene sequence (CDS region) to get rid of the intronic regions, if any.
	3. Sample sequences mapped to reference gene sequence are always required to get start and end positions in cases of PfCRT and PfDHPS. They are sometimes required for PfDHFR and PfMDR1 as well.
	4. This will generate a single document containing all the sequences assembled to the selected gene.

## Codon coverage calculation:

`python calc_cov_ind.py <sample_type> <gene_name> <nest_output_path> <bed_file_path> <geneious_document_path>`

User inputs required:
- sample_type: Individual or Pooled
- gene_name: PfDHFR or PfMDR1 or PfDHPS or PfCRT
- nest_output_path: /data/Haiti_ind_samples/
- bed_file_path: /data/mdr.bed
- geneious_document_path: /data/individual/PfDHFR/91\ documents\ from\ PfDHFR\ Assembled.txt

**NOTE:** Use the exact options as mentioned in sample_type and gene_name because the input is case sensitive.

Individual samples:
`python calc_cov_ind.py Individual PfDHFR <nest_output_path> <bed_file_path> <geneious_document_path>`

