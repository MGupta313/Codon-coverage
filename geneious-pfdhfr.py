import re 

with open("10 documents from PFDHFRpooled.txt") as f:
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
	filename = ""
	for element in temp_list:
		if element[0].isdigit():
			filename = element.split("assembled")[0]+"-assembled.txt"
			#print("if -", filename)
			with open(filename, 'w') as fo:
				fo.write(element+"\n")
		else:
			#print("else", filename)
			with open(filename, 'a') as fo:
				fo.write(element)

create_new_file(temp_list)

with open('001-PfDHFR-assembled.txt') as sf:
	sf_lines = sf.readlines()
	#print(sf_lines[1])
	#find first ATG codon after ?
	#first eleminate ?
	#from there start +3 sequences
	#then find first ATG
	#then find first TAA,TAG,TGA
	for i in range(len(sf_lines[1])):
		if re.search(pattern = "[^?]", string = sf_lines[1][i]):
			temp_pos = i
			break
	for j in range(temp_pos, len(sf_lines[1]), 3):
		if re.search(pattern = "ATG", string = sf_lines[1][j:j+3]):
			start_pos = j+1
			break
	#print("start_pos", start_pos, sf_lines[1][start_pos-1])
	for j in range(temp_pos, len(sf_lines[1])):
		if re.search(pattern = "[?]", string = sf_lines[1][j]):
			end_pos = j
			break
	#print("end_pos", end_pos, sf_lines[1][end_pos-1])

