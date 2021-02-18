"""
Haiti Codon Coverage Distribution:

Step 0: Running NeST on raw reads of 100 individual samples

## python3 nest.py -i ~/data/individual_sample_fq -a ref/pfalciparum/adapters.fa -r ref/pfalciparum/mdr.fa -o local/Haiti_individual_samples
					-b ref/pfalciparum/mdr.bed -m bowtie2 --varofint ref/pfalciparum/Reportable_SNPs.csv

Step 1: Run Samtools depth on the bam file generated from NeST to get nucleotide depth

## samtools depth <bamfile> <output.txt>

Step 2: Using mdr.bed file for getting the exonic and intronic regions plus codon start and end positions for reference for PfDHFR, PfMDR1, 
PfDHPS and PfCRT genes.

Step 3: Using the output.txt file to extract exonic region for PfDHFR gene. Check if the start/end positions match bed file. 

Step 4: If start/end do not match, then go through geneious assembled file as reference to get exonic region for the gene.

Step 5: Calculate codon coverage by adding nucleotide depth in groups of 3. Remove the codon positions where codon coverage is less than 5. 
Calculate average codon coverage, which is what will be used for further studies.

Step 6: Using jupyter notebook, categorize the codon positions in groups of 30 for PfDHFR so that the plot is not too crowded. Generate plots to see the distribution of codon coverage.

###################################################################
###################### Usage of this script: ######################

python calc_cov_ind.py <sample_type> <gene_name> <nest_output_path> <bed_file_path> <geneious_document_path>

Options/Examples:
 - sample_type: Individual or Pooled
 - gene_name: PfDHFR or PfMDR1 or PfDHPS or PfCRT
 - nest_output_path: /data/Haiti_ind_samples/
 - bed_file_path: /data/mdr.bed
 - geneious_document_path: /data/Individual/PfDHFR/91\ documents\ from\ PfDHFR\ Assembled.txt
"""

import subprocess
import sys
import os
import re


"""
Calculating nucleotide depth using Samtools depth.

Running on output_FM_SR_DD_RG.bam files for each sample generated from NeST.
Some samples did not generate any files to calculate nucleotide depth so those were removed.
For all other files, depth was calculated
"""

def cal_nucl_depth(directory): 
	# directory = r'/Users/mansi/Desktop/haiti-ind/Haiti_ind/' #input directory which has all files generated from NeST
	for filename in os.listdir(directory):
		if filename != ".DS_Store":
			# print(filename)
			if sample_type == "Individual":
				sample_num = filename[9:15].split("P")[0] # getting the sample number from the filename
			if sample_type == "Pooled":
				sample_num = filename[10:15].split("P")[0]
			# print(sample_num)
			
			filepath = os.path.join(directory, filename)
			sample_dir = os.listdir(filepath)
			if len(sample_dir) == 0: # checking for empty files, samples with no output
				print("Empty directory {}".format(filename))
				with open("empty_files.txt", "a") as fo: # recording the empty files/samples in a text file
					fo.write(filename+"\n")
			else:
				bamfile_path = filepath+"/alignments/output_FM_SR_DD_RG.bam" # non-empty files used for depth calculation
				depth_output = sample_num+"_"+gene_name+"_depth.txt"
				print("Running samtools depth on:  {}".format(filename))
				f=open(depth_output, "w")
				subprocess.call(['samtools', 'depth', bamfile_path], stdout = f)
				subprocess.call(['mv', depth_output, depth_dir]) # moving the files to a new folder
				print("Calculated nucleotide depth for {} \n".format(sample_num))
			
"""
Processing reference bed file to obtain start and stop sites for genes.

Just tranforming data from bed file to usable lists.
Finally converting into a dictionary with informtion for all genes.
"""
def process_bed_file(input_bed_file, gene):
	with open(input_bed_file) as f: # processing the bed file
		lines = f.readlines()
		final_lines = []
		for l in lines:
			line = l.strip() # removing any new line characters from either end
			line = line.split("\t") # removing the tab characters
			final_lines.append(line) # has all the rows as list of lists with one list for each gene
		#print(final_lines)

		list1 = [] # single list of exon start positions for all genes 
		list3 = [] # single list of exon sizes for all genes

		list2 = [] # list of lists containing start sites for each gene
		list4 = [] # list of lists containing exon sizes for each gene
		

		for i in range(len(final_lines)):
			list1.append(final_lines[i][11].strip('"').rstrip(",")) # index 11 = exon start positions
			list3.append(final_lines[i][10].strip('"').rstrip(",")) # index 10 = exon sizes
		for j in range(len(list1)):
			list2.append(list1[j].split(",")) #start sites from list1 directly
			list4.append(list3[j].split(",")) #exon sizes from list3 directly

		# bed_file dictionary --> gene: [exon start position, exon size]
		bed_file = {}
		for k in range(len(final_lines)):
			bed_file[final_lines[k][0]] = [] # adding genes as keys
		 
		for m in range(len(list2)):
			for n in range(len(list2[m])):
				list5 = []
				list5.append(list2[m][n]) # appending start sites
				list5.append(list4[m][n]) # appending exon sizes
				fin_tup = tuple(list5)    # creating a tuple to make them a pair
				bed_file[final_lines[m][0]].append(fin_tup) # appending it to the dictionary corresponding to each gene

		for key in list(bed_file.keys()):
			if key != gene: # keeping information for just the gene in use eg. PfDHFR
				del(bed_file[key])

		#final dictionary --> {'PfDHFR': [('1', '1827')]}

		#print("Processed bed file for", gene, "\n")
		print("Processed bed file for {} \n".format(gene))
	return bed_file

"""
Processing files generated using Geneious.

When samples are mapped to reference gene sequence using Geneious,
the final output comes as a single file containing the assembled
sequences for all the samples.

That single document is being processed below to separate into the
respective sequence files.
"""

def process_geneious_file(input_geneious_file, gene):

	with open(input_geneious_file) as f:
		lines_list = f.readlines()
		converted_list = []
		final_lines = []
		for element in lines_list:
			converted_list.append(element.strip()) # removing new line characters
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
		file_list = []
		for element in new_list:
			if element[0].isdigit():
				if sample_type == "Individual":
					filename = element.split("assembled")[0][9:15].split("P")[0]+"_"+gene_name+"_assembled.txt"
				if sample_type == "Pooled":
					filename = element.split("assembled")[0][10:15].split("P")[0]+"_"+gene_name+"_assembled.txt"
				file_list.append(filename)
				#print("if -", filename)
				with open(filename, 'w') as fo:
					fo.write(element+"\n")
			else:
				#print("else", filename)
				with open(filename, 'a') as fo:
					fo.write(element)

		for item in file_list:
			subprocess.call(['mv', item, geneious_dir])

"""
Processing the depth files generated from samtools depth above.

Storing start and end positions based on bed file, sample depth files and geneious files, as required.
Storing depth values between start and end positions.
"""

def process_depth_file_dhfr_mdr1(input_depth_file, ref_bed_file, input_geneious_filepath):
	#print("Processing", input_depth_file+"\n")
	print("Processing {} \n".format(input_depth_file))
	with open(input_depth_file) as f2:
		lines2 = f2.readlines()
		final_lines2 = []
		for l in lines2:
			line2 = l.rstrip("\n")
			line2 = line2.split("\t")
			final_lines2.append(line2) # a list of lists with 'gene', 'nucleotide position', 'nucleotide depth' as values
		#print(final_lines2)
		
		# file_txt dictionary --> {'PfMDR1': [], 'PfCRT': [], 'PfDHFR': [], 'PfDHPS': []}
		file_txt = {}
		for i in range(len(final_lines2)):
			file_txt[final_lines2[i][0]] = [] # creating a dictionary with genes as keys
			
		for i in range(len(final_lines2)):
			list7 = []
			list7.append(final_lines2[i][1])
			list7.append(final_lines2[i][2])
			final_tup = tuple(list7)
			file_txt[final_lines2[i][0]].append(final_tup) # Appending position and depth pairs for each gene
		# file_txt --> {'PfMDR1': [(pos, depth), (pos, depth), (pos, depth)..]} for all genes
	
	"""
	Using the information from both nucleotide depth file and bed file dictionary.
	Checking if the start and end sites from bed file match the positions from sample sequences.
	
	For some genes, the exonic region starts from position 1 (eg. PfDHFR) so the exon size directly becomes
	the end site whereas for others (eg. PfCRT), exonic region starts from some other position so in that case,
	to get the the end position, exon size was added to exon start site.
	"""
	final_dict = {}
	for key in ref_bed_file:
		if key in file_txt: # taking the information from nucleotide depth file just for the gene being run
			temp = file_txt[key]
			temp_list = []
			for i in range(len(ref_bed_file[key])):
				start_from_bed_file = int(ref_bed_file[key][i][0])
				if start_from_bed_file != 1: # if exon starts from position 1 or not
					end_from_bed_file = start_from_bed_file + int(ref_bed_file[key][i][1]) # end site = exon size + exon start position
				else:
					end_from_bed_file = int(ref_bed_file[key][i][1])
				"""
				If the first nucleotide position for the given gene in the nucleotide depth file matches
				the bed file start position, then append all the depth values, for start to end positions,
				to a temp list. These start and end psoitions are stored and used for further calculations.

				If the first nucleotide position for the given gene in the nucleotide depth file does not
				match the bed file start position, then use the Geneious assembled file.

				Geneious assembled file is processed first and then the start and end positions from that
				sequence file are stored and used for further calculations.
				"""
				#add if the start codon matches the mdr start pos, 
				#if not open assembled file and get positions from there
				if int(temp[0][0]) == start_from_bed_file:
					print("Start position from bed file matches start position of sample seq")
					#print("start_from_bed_file", str(start_from_bed_file))
					print("start_from_bed_file: {}".format(start_from_bed_file))
					#print("end_from_bed_file", str(end_from_bed_file) + "\n")
					print("end_from_bed_file: {} \n".format(end_from_bed_file))
					for j in range(len(temp)):
						if int(temp[j][0]) >= start_from_bed_file and int(temp[j][0]) <= end_from_bed_file:
							temp_list.append(int(temp[j][1]))
					#if they don't match it means the sample seq is truncated
					#check for start/stop using geneious assembled seq
				else:
					print("start position does not match bed file start position"+"\n")
					print("going through geneious assembled file now")
					with open(input_geneious_filepath) as f3:
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

			final_dict[key] = temp_list # list of depth values for nucleotide positions between start and end positions.
	
	return final_dict

"""
Processing the depth files generated from samtools depth above.

Storing start and end positions based on bed file, sample depth files and geneious files, as required.
Storing depth values between start and end positions.
"""
def process_depth_file_dhps_crt(input_depth_file, ref_bed_file, input_geneious_filepath):
	#print("Processing", input_depth_file+"\n")
	print("Processing {} \n".format(input_depth_file))
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
	"""
	Using the information from both nucleotide depth file and bed file dictionary.
	Checking if the start and end sites from bed file match the positions from sample sequences.

	Both PfDHPS and PfCRT contain introns, therefore, they directly use geneious assembled files in which samples
	have been mapped to the CDS region of the respective genes so that only the exonic region remains.
	"""
	final_dict = {}
	for key in ref_bed_file:
		if key in file_txt:
			temp = file_txt[key]
			temp_list = []
			print("start position does not match bed file start position"+"\n")
			print("going through geneious assembled file now")
			with open(input_geneious_filepath) as f3:
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
	return final_dict

"""
Calculating codon coverage from nucleotide depth selected from above.

1 codon = 3 nucleotides
codon coverage = nucleotide depth 1+ nucleotide depth 2+ nucleotide depth 3
Coverage cut off = 5
"""
def cal_codon_cov(input_dict):
	codon_cov_dict={}
	for key in input_dict:
		sum1=0
		temp_list2=[]
		for i in range(len(input_dict[key])):
			sum1 = sum1 + int(input_dict[key][i])
			if (i+1)%3 == 0: # codon consideration: 3 nucleotides make 1 codon
				if sum1 > 5: # threshold set = 5
					temp_list2.append(sum1)
					sum1=0
		codon_cov_dict[key] = temp_list2
		#avg_codon_cov[key] = avg_list
	print("Calculated codon coverage"+"\n")
	return codon_cov_dict

"""
Actual codon coverage = average of 3 nucleotide depths.

For any further calculation/visualization, average coverage was used.
"""

def avg_codon_cov(input_dict):
	#codon_cov_dict={}
	avg_codon_cov = {}
	for key in input_dict:
		sum1=0
		avg = 0
		temp_list2=[]
		avg_list = []
		for i in range(len(input_dict[key])):
			sum1 = sum1 + int(input_dict[key][i])
			if (i+1)%3 == 0:
				if sum1 > 5: # threshold set = 5
					avg = sum1//3
					temp_list2.append(sum1)
					avg_list.append(avg)
					sum1=0
		#codon_cov_dict[key] = temp_list2
		avg_codon_cov[key] = avg_list

	print("Calculated average codon coverage"+"\n")
	return avg_codon_cov

"""
Creating a tab separated file with two columns - Codon position and Codon coverage

Creating 2 types of files:
  1.) Coverage results
  2.) Avergae coverage results (*for any further study, this is to be used)
"""

def create_result_file(codon_cov, avg_cov, sample, gene):
	print("Creating codon coverage and average coverage results files"+"\n")
	
	filename = sample+"_"+gene+"_result.txt"
	with open(filename, "w") as fo:
		fo.write("Codon_pos"+"\t"+ "Codon_coverage"+"\n")
		for key in codon_cov:
			for i in range(len(codon_cov[key])):
				fo.write(str(i+1)+"\t"+str(codon_cov[key][i])+"\n")

	subprocess.call(['mv', filename, coverage_dir])
	
	avg_filename = "avg_"+sample+"_"+gene+"_result.txt"
	with open(avg_filename, "w") as fo:
		fo.write("Codon_pos"+"\t"+ "Codon_coverage"+"\n")
		for key in avg_cov:
			for i in range(len(avg_cov[key])):
				fo.write(str(i+1)+"\t"+str(avg_cov[key][i])+"\n")

	subprocess.call(['mv', avg_filename, avg_coverage_dir])


"""
Creating a file with just the coverage information separated by comma. No information about position.

Creating category files in 2 types:
  1.) Coverage results
  2.) Avergae coverage results (*for any further study, this is to be used)
"""

def create_category(codon_cov, avg_cov, sample, gene):
	print("Creating category files for codon coverage and average coverage"+"\n")
	
	filename = sample+"_"+gene+"_cat_cov.txt"
	
	with open(filename, "w") as fo:
		intermed_list = []
		for key in codon_cov:
			for i in range(len(codon_cov[key])):
				fo.write(str(codon_cov[key][i])+",")

	subprocess.call(['mv', filename, category_dir])
	
	avg_filename = "avg_"+sample+"_"+gene+"_cat_cov.txt"

	with open(avg_filename, "w") as fo:
		intermed_list = []
		for key in avg_cov:
			for i in range(len(avg_cov[key])):
				fo.write(str(avg_cov[key][i])+",")

	subprocess.call(['mv', avg_filename, avg_category_dir])

#########################################################################
#Main run code#
#Prequistes:
		# 1.) bam files from NeST
		# 2.) bed file in the working directory
		# 3.) assembled sequences from geneious
#########################################################################

sample_type = sys.argv[1]
gene_name = sys.argv[2]
nest_output_path = sys.argv[3]
bed_file_path = sys.argv[4]
geneious_document_path = sys.argv[5]


depth_dir = 'Haiti_coverage/'+sample_type+'/'+gene_name+'/Samtools_output/'
geneious_dir = 'Haiti_coverage/'+sample_type+'/'+gene_name+'/Geneious_assembled/'
coverage_dir = 'Haiti_coverage/'+sample_type+'/'+gene_name+'/Coverage_results/'
category_dir = 'Haiti_coverage/'+sample_type+'/'+gene_name+'/Category_results/'
avg_coverage_dir = 'Haiti_coverage/'+sample_type+'/'+gene_name+'/Coverage_results/avg_cat_cov/'
avg_category_dir = 'Haiti_coverage/'+sample_type+'/'+gene_name+'/Category_results/avg_cov/'

subprocess.call(['mkdir', 'Haiti_coverage'])
subprocess.call(['mkdir', 'Haiti_coverage/'+sample_type])
subprocess.call(['mkdir', 'Haiti_coverage/'+sample_type+'/'+gene_name])
subprocess.call(['mkdir', depth_dir])
subprocess.call(['mkdir', geneious_dir])
subprocess.call(['mkdir', coverage_dir])
subprocess.call(['mkdir', category_dir])
subprocess.call(['mkdir', avg_coverage_dir])
subprocess.call(['mkdir', avg_category_dir])

#Step 1: Calculate nucleotide depth
cal_nucl_depth(nest_output_path)

#Step 2: Process bed file
bed_file = process_bed_file(bed_file_path, gene_name)

#Step 3: Process geneious document
process_geneious_file(geneious_document_path, gene_name)

#Step 4: Process nucleotide depth file and get start/end positions
"""
If gene is either PfDHFR or MDR1 then bed file is used first then if required geneious files

If gene is either PfDHPS or PfCRT then bed file is not used since these two genes contain introns. 
They directly use geneious files which are mapped against the CDS region.
"""

if gene_name == "PfDHFR" or "PfMDR1":	
	for filename in os.listdir(depth_dir):
		for file in os.listdir(geneious_dir):
			if filename[0:6].split("_")[0] == file[0:6].split("_")[0]:
				sample_id = filename[0:6].split("_")[0]
				# Processing each sample with its respective geneious file
				depth_filepath = os.path.join(depth_dir, filename)
				geneious_filepath = os.path.join(geneious_dir, file)
				depth_dict = process_depth_file_dhfr_mdr1(depth_filepath, bed_file, geneious_filepath)
				# recording if there was information (reads) for a sample
				if len(depth_dict) == 0:
					print("No output from Samtools depth for ", filename, "\n")
					with open("No_samtools_output.txt", "a") as fo:
						fo.write(depth_filepath+"\n")
				else:
					if len(depth_dict[gene_name]) == 0:
						print("Not enough reads for Geneious assembly", filename, "\n")
						with open("No_geneious.txt", "a") as fo:
							fo.write(depth_filepath+"\n")

	#Step 5: Calculating codon coverage		
					else: # if there is some information (reads) for a sample then calculating codon coverage
						codon_cov_dict_final = cal_codon_cov(depth_dict)
						avg_codon_cov_final = avg_codon_cov(depth_dict)
						create_result_file(codon_cov_dict_final, avg_codon_cov_final, sample_id, gene_name)
						create_category(codon_cov_dict_final, avg_codon_cov_final, sample_id, gene_name)
						print("Complete all calculations for", filename, "\n")

if gene_name == "PfDHPS" or "PfCRT":
	for filename in os.listdir(depth_dir):
		for file in os.listdir(geneious_dir):
			if filename[0:6].split("_")[0] == file[0:6].split("_")[0]:
				sample_id = filename[0:6].split("_")[0]
				# Processing each sample with its respective geneious file
				depth_filepath = os.path.join(depth_dir, filename)
				geneious_filepath = os.path.join(geneious_dir, file)
				depth_dict = process_depth_file_dhfr_mdr1(depth_filepath, bed_file, geneious_filepath)
				# recording if there was information (reads) for a sample
				if len(depth_dict) == 0:
					print("No output from Samtools depth for ", filename, "\n")
					with open("No_samtools_output.txt", "a") as fo:
						fo.write(depth_filepath+"\n")
				else:
					if len(depth_dict[gene_name]) == 0:
						print("Not enough reads for Geneious assembly", filename, "\n")
						with open("No_geneious.txt", "a") as fo:
							fo.write(depth_filepath+"\n")

	#Step 5: Calculating codon coverage		
					else: # if there is some information (reads) for a sample then calculating codon coverage
						codon_cov_dict_final = cal_codon_cov(depth_dict)
						avg_codon_cov_final = avg_codon_cov(depth_dict)
						create_result_file(codon_cov_dict_final, avg_codon_cov_final, sample_id, gene_name)
						create_category(codon_cov_dict_final, avg_codon_cov_final, sample_id, gene_name)
						print("Complete all calculations for", filename, "\n")









