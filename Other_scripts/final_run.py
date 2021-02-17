import os
import subprocess
import re

#print("Opening mdr.bed file for processing")
with open("/Users/mansi/Desktop/CDC/mdr.bed") as f:
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
		if key != "PfDHFR":
			del(bed_file[key])

	#print(bed_file)


#with open("/Users/mansi/Desktop/CDC/scripts/001-PfDHFR-NeST.txt") as f2:
with open("/Users/mansi/Desktop/CDC/samtools_output/001.bam.txt") as f2:
	print("001")
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
#print(file_txt)

final_dict = {}
for key in bed_file:
	if key in file_txt:
		temp = file_txt[key]
		temp_list = []
		for i in range(len(bed_file[key])):
			start_from_bed_file = int(bed_file[key][i][0])
			if start_from_bed_file != 1:
				end_from_bed_file = start_from_bed_file + int(bed_file[key][i][1])
			else:
				end_from_bed_file = int(bed_file[key][i][1])
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

def cal_codon_cov(input_dict):
	codon_cov_dict={}
	#avg_codon_cov = {}
	for key in input_dict:
		sum1=0
		#avg = 0
		temp_list2=[]
		#avg_list = []
		for i in range(len(input_dict[key])):
			sum1 = sum1 + int(input_dict[key][i])
			if (i+1)%3 == 0:
				#avg = sum1//3
				temp_list2.append(sum1)
				#avg_list.append(avg)
				sum1=0
		codon_cov_dict[key] = temp_list2
		#avg_codon_cov[key] = avg_list
	return codon_cov_dict

codon_cov_dict_final = cal_codon_cov(final_dict)
#print(codon_cov_dict_final)

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
				avg = sum1//3
				temp_list2.append(sum1)
				avg_list.append(avg)
				sum1=0
		#codon_cov_dict[key] = temp_list2
		avg_codon_cov[key] = avg_list
	return avg_codon_cov

avg_codon_cov_final = avg_codon_cov(final_dict)

with open("/Users/mansi/Desktop/CDC/test/001-PfDHFR-category.txt", "w") as fo:
	intermed_list = []
	for key in codon_cov_dict_final:
		for i in range(len(codon_cov_dict_final[key])):
			fo.write(str(codon_cov_dict_final[key][i])+",")

with open("/Users/mansi/Desktop/CDC/test/avg_001-PfDHFR-category.txt", "w") as fo:
	intermed_list = []
	for key in avg_codon_cov_final:
		for i in range(len(avg_codon_cov_final[key])):
			fo.write(str(avg_codon_cov_final[key][i])+",")


#codon coevrage for each codon position

with open("/Users/mansi/Desktop/CDC/test/001-PfDHFR-result.txt", "w") as fo:
	fo.write("Codon_pos"+"\t"+ "Codon_coverage"+"\n")
	for key in codon_cov_dict_final:
		for i in range(len(codon_cov_dict_final[key])):
			fo.write(str(i+1)+"\t"+str(codon_cov_dict_final[key][i])+"\n")

with open("/Users/mansi/Desktop/CDC/test/avg_001-PfDHFR-result.txt", "w") as fo:
	fo.write("Codon_pos"+"\t"+ "Codon_coverage"+"\n")
	for key in avg_codon_cov_final:
		for i in range(len(avg_codon_cov_final[key])):
			fo.write(str(i+1)+"\t"+str(avg_codon_cov_final[key][i])+"\n")

print("All done")


"""
def creating_categories(input_dict2):
	temp_dict = {}
	for key in input_dict2:
		sum2=0
		temp_list3 = []
		for i in range(len(input_dict2[key])):
			sum2= sum2+int(input_dict2[key][i])
			if (i+1)%20 == 0:
				temp_list3.append(sum2)
				sum2=0
		temp_dict[key] = temp_list3
	return temp_dict

#category_dict = creating_categories(codon_cov_dict_final)

#codon coverage added for each 20 codon positions - to categorize them
with open("001-PfDHFR-category.txt", "w") as fo:
	fo.write("codon_category"+"\t"+ "codon_coverage"+"\n")
	for key in category_dict:
		for i in range(len(category_dict[key])):
			fo.write("category "+str(i+1)+"\t"+str(category_dict[key][i])+"\n")

"""






