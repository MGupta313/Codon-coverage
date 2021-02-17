import re 
import os
import subprocess


directory = r'/Users/mansi/Desktop/Codon_coverage_data/Geneious_output/' #input directory
for filename in os.listdir(directory):
	sample_num = filename[10:13]
	filepath = os.path.join(directory, filename)
	print(filename)
	with open(filepath) as f:
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


	def remove(input_list): 
		#pattern = '[0-9]'
		new_list = []
		for i in range(0, len(input_list)):
			if input_list[i][0].isdigit():
				new_list.append(input_list[i])
			else:
				new_list.append(re.sub(pattern = '\d', repl= '', string= input_list[i])) 
		return new_list

	temp_list = remove(final_lines)
	#print(temp_list)

	def create_new_file(num_list):
		for element in temp_list:
			if element[0].isdigit():
				new_filepath = "/Users/mansi/Desktop/Codon_coverage_data/Geneious_assembled/"+filename
				with open(new_filepath, 'w') as fo:
					fo.write(filename+"\n")
			else:
				with open(new_filepath, 'a') as fo:
					fo.write(element)

	create_new_file(temp_list)









