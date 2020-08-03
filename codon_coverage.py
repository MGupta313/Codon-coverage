"""
Calculating codon coverage:

Step1: Running NeST on raw reads of 10 pooled samples

## python3 nest.py -i ~/Desktop/pooled_sample_fq -a ref/pfalciparum/adapters.fa -r ref/pfalciparum/mdr.fa -o local/Haiti_pooled_samples -b ref/pfalciparum/mdr.bed -m bowtie2 --varofint ref/pfalciparum/Reportable_SNPs.csv

Step2: Run Samtools depth on the bam file generated from NeST to get nucleotide depth

## samtools depth <bamfile> <output.txt>

Step3: Using mdr.bed file for getting the exonic and intronic regions plus codon start and stop positions for reference for PfDHFR and PfMDR1 genes

Step4: Using the output.txt file to extract exonic region for PfDHFR gene. Check if the start/stop positions match bed file. If not, go through geneious assembled file as reference to get exonic region for the gene.

Step5: Calculate codon coverage by adding nucleotide depth in groups of 3. Remove the codon positions where codon coverage is less than 5.

Step6: Calculate average codon coverage.

Step7: Using jupyter notebook, categorize the codon positions in groups of 30 for PfDHFR so that the plot is not too crowded. Generate a boxenplot to see the distribution of codon coverage.

Step8: Save the figure in TIF format.
"""

import os
import subprocess
import re

def cal_nucl_depth():
	directory = r'/Users/mansi/Desktop/Haiti_pooled_bam/' #input directory
	for filename in os.listdir(directory):
		if filename != ".DS_Store":
			filepath = os.path.join(directory, filename)
			bamfile_path = filepath+"/alignments/output_FM_SR_DD_RG.bam"
			bamoutfile = filename[10:13]+"_depth.txt"
			return bamoutfile

cal_nucl_depth()


			#print("running samtools depth on", bamfile)
			#f=open(bamoutfile, "w")
			#subprocess.call(['samtools', 'depth', bamfile], stdout = f)














