import os
import subprocess
import re


directory = r'/Users/mansi/Desktop/Codon_coverage_data/Haiti_pooled_bam/' #input directory
for filename in os.listdir(directory):
	if filename != ".DS_Store":
		#print(filename)
		sample_num = filename[10:13]
		filepath = os.path.join(directory, filename)
		bamfile_path = filepath+"/alignments/output_FM_SR_DD_RG.bam"
		new_path = "/Users/mansi/Desktop/Haiti_bam_files/"+sample_num+".bam"
		subprocess.call(['cp', bamfile_path, new_path])