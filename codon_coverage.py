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
	directory = r'/Users/mansi/Desktop/Codon_coverage_data/Haiti_pooled_bam/' #input directory
	for filename in os.listdir(directory):
		if filename != ".DS_Store":
			filepath = os.path.join(directory, filename)
			bamfile_path = filepath+"/alignments/output_FM_SR_DD_RG.bam"
			depth_output = filename[10:13]+"_depth.txt"
			
			print("running samtools depth on", filename)
			f=open(depth_output, "w")
			subprocess.call(['samtools', 'depth', bamfile_path], stdout = f)
			subprocess.call(['mv', depth_output, '/Users/mansi/Desktop/Codon_coverage_data/Samtools_output/'])

			print("calculated nucleotide depth for", filename, "\n")
cal_nucl_depth()

def process_bed_file(input_bed_file, gene):
	with open(input_bed_file) as f:
		lines = f.readlines()
		final_lines = []
		for l in lines:
			line = l.strip()
			line = line.split("\t")
			final_lines.append(line)
		#print(final_lines)

		list1 = []
		list2 = []
		list3 = []
		list4 = []
		bed_file = {}
		for i in range(len(final_lines)):
			list1.append(final_lines[i][11].strip('"').rstrip(","))
			list3.append(final_lines[i][10].strip('"').rstrip(","))
		for j in range(len(list1)):
			list2.append(list1[j].split(",")) #start sites
			list4.append(list3[j].split(",")) #exon size

		#print(list2)
		#print(list4)

		for k in range(len(final_lines)):
			bed_file[final_lines[k][0]] = []
		 
		for m in range(len(list2)):
			for n in range(len(list2[m])):
				list5 = []
				list5.append(list2[m][n])
				list5.append(list4[m][n])
				fin_tup = tuple(list5)
				bed_file[final_lines[m][0]].append(fin_tup)

		for key in list(bed_file.keys()):
			if key != gene:
				del(bed_file[key])

		print("Processed bed file for", gene, "\n")
		return bed_file

process_bed_file("/Users/mansi/Desktop/CDC/mdr.bed", "PfDHFR")


def process_geneious_file(input_geneious_file):

	subprocess.call(['mkdir', '/Users/mansi/Desktop/Codon_coverage_data/Geneious_assembled/'])

	with open(input_geneious_file) as f:
		lines_list = f.readlines()
		converted_list = []
		final_lines = []
		for element in lines_list:
			converted_list.append(element.strip())
		#print(converted_list)
		while("" in converted_list) : 
	   		converted_list.remove("")
	   	for element in converted_list:
	   		new_lines = element.replace(" ", '')
	   		final_lines.append(new_lines)
	#print(final_lines)

	new_list = []
    for i in range(0, len(final_lines)):
    	if final_lines[i][0].isdigit():
    		new_list.append(final_lines[i])
    	else:
    		new_list.append(re.sub(pattern = '\d', repl= '', string= final_lines[i])) 
    
    filename = ""
	for element in new_list:
		if element[0].isdigit():
			filename = element.split("assembled")[0]+"_assembled.txt"
			with open(filename, 'w') as fo:
				fo.write(element+"\n")
		else:
			#print("else", filename)
			with open(filename, 'a') as fo:
				fo.write(element)
		subprocess.call(['mv', filename, "/Users/mansi/Desktop/Codon_coverage_data/Geneious_assembled"])

	print("Processed geneious file for", gene, "\n")

process_geneious_file("/Users/mansi/Desktop/CDC/10 documents from PFDHFRpooled.txt", "PfDHFR")


def process_depth_file(input_depth_file):
	with open(input_depth_file) as f2:
		lines2 = f2.readlines()
		final_lines2 = []
		for l in lines2:
			line2 = l.rstrip("\n")
			line2 = line2.split("\t")
			final_lines2.append(line2)
		#print(final_lines2)
		
		file_txt = {}
		for i in range(len(final_lines2)):
			file_txt[final_lines2[i][0]] = []
			
		for i in range(len(final_lines2)):
			list7 = []
			list7.append(final_lines2[i][1])
			list7.append(final_lines2[i][2])
			final_tup = tuple(list7)
			file_txt[final_lines2[i][0]].append(final_tup)
	return file_txt

process_depth_file("/Users/mansi/Desktop/Codon_coverage_data/Samtools_output/001_depth.txt")

def get_exon(depth_file, bed_ref_file, geneious_ref_file):
	final_dict = {}
	for key in bed_ref_file:
		if key in depth_file:
			temp = depth_file[key]
			temp_list = []
			for i in range(len(bed_ref_file[key])):
				start_from_bed_file = int(bed_ref_file[key][i][0])
				if start_from_bed_file != 1:
					end_from_bed_file = start_from_bed_file + int(bed_ref_file[key][i][1])
				else:
					end_from_bed_file = int(bed_ref_file[key][i][1])
				#add if the start codon matches the mdr start pos, 
				#if not open assembled file and get positions from there
				if int(temp[0][0]) == start_from_bed_file:
					print("Start position from bed file matches start position of sample seq")
					print("start_from_bed_file", start_from_bed_file)
					print("end_from_bed_file", end_from_bed_file)
					for j in range(len(temp)):
						if int(temp[j][0]) >= start_from_bed_file and int(temp[j][0]) <= end_from_bed_file:
							temp_list.append(int(temp[j][1]))
					#if they don't match it means the sample seq is truncated
					#check for start/stop using geneious assembled seq
				else:
					print("start position does not match bed file start position")
					print("going through geneious assembled file now")
					with open("/Users/mansi/Desktop/CDC/geneious_output/001-PfDHFR-assembled.txt") as f3:
						print("opening geneious assembled file")
						lines3 = f3.readlines()
						for i in range(0, len(lines3[1]), 3):
							if re.search(pattern="\w{3}", string = lines3[1][i:i+3]):
								start_pos = i+1
								break
						print("start_pos", start_pos)

						for i in range(0, len(lines3[1]), 3):
							if re.search(pattern="\w{3}", string = lines3[1][i:i+3]):
								#print(lines3[1][i:i+3], i)
								end_pos=i+1
						print("end_pos", end_pos)
						for j in range(len(temp)):
							if int(temp[j][0]) >= start_pos and int(temp[j][0]) <= end_pos:
								temp_list.append(int(temp[j][1]))

			final_dict[key] = temp_list
	#print(final_dict)

get_exon(file_txt, bed_file, )




