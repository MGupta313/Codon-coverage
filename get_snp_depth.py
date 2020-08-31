"""
import os
import re

directory = r'/Users/mansi/Desktop/Codon_coverage_data (PfDHFR)/Coverage_results'
for filename in os.listdir(directory):
	if filename != ".DS_Store":
		if re.search(pattern = "avg", string = filename):
			with open(filename) as f:
				lines = f.readlines()
				print(lines)
"""

import re
import subprocess

avg_filename = "snp_depth_010_PfMDR1.csv"
file = '/Users/mansi/Desktop/Codon_coverage_data (PfMDR1)/Coverage_results/avg_010_PfMDR1_result.txt'
with open(file) as f:
	lines = f.readlines()
	for l in lines:
		new = l.strip()
		new_lines = new.split('\t')
		if new_lines[0] == '86':
			print(new_lines[0], new_lines[1])
			with open(avg_filename, "a+") as fo:
				fo.write("Codon_pos"+"\t"+ "Codon_coverage"+"\n")
				fo.write(new_lines[0]+"\t"+new_lines[1]+"\n")
		if new_lines[0] == "184":
			print(new_lines[0], new_lines[1])
			with open(avg_filename, "a+") as fo:
				fo.write(new_lines[0]+"\t"+new_lines[1]+"\n")
		if new_lines[0] == "1034":
			print(new_lines[0], new_lines[1])
			with open(avg_filename, "a+") as fo:
				fo.write(new_lines[0]+"\t"+new_lines[1]+"\n")
		if new_lines[0] == "1042":
			print(new_lines[0], new_lines[1])
			with open(avg_filename, "a+") as fo:
				fo.write(new_lines[0]+"\t"+new_lines[1]+"\n")
		if new_lines[0] == "1246":
			print(new_lines[0], new_lines[1])
			with open(avg_filename, "a+") as fo:
				fo.write(new_lines[0]+"\t"+new_lines[1]+"\n")
			

subprocess.call(['mv', avg_filename, "/Users/mansi/Desktop/Codon_coverage_data (PfMDR1)/Coverage_results/"])
	