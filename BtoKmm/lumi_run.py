import os
import csv
import ast
import json

with open('file_run.txt', 'r') as file:
	doc = file.readlines()


query = 'das_client --query="lumi file='

selectedLumi= {705, 513, 580, 869, 645, 678, 1257, 1033, 713, 1163, 1164, 1129, 1041, 568, 1597}

file_lumi=[]
k=0
for i,line in enumerate(doc):
	if '317392' in line:
		print('\n\n\n', i+1, '/', len(doc))
		query_ = query+line.split(',')[0].strip()
		query_+='" --json'
		stream = os.popen(query_)
		output = stream.read()
		#print(query_, '\n')
		#print(output)
		#x = ast.literal_eval(output.strip())
		#file_lumi.append(line+','+output.strip())
		x = json.loads(output.strip())
		for entrada in x:
			lumisection = entrada['lumi'][0]['lumi_section_num']
			for ele in lumisection:
				if ele in selectedLumi:
					print(query_,'\n\t->' ,ele)
					print('Te encontre! ----> ', line.split(',')[0].strip())
				#else:
					#print(ele, ' No esta en las lumi sections elegidas')
	#if i==1000: break


#print(file_lumi)
#print(k) 
	#print(i+1, '/', len(doc))
	#query_ = query+line.strip()
	#query_+='"'
	#stream = os.popen(query_)
	#output = stream.read()
	#if '317392' in output:
	#	print('Found ya!\n', line.strip())
	#file_run.append([line.strip(), output.strip().replace('"', '')])

#with open('file_run.txt', 'w+') as csv_file:
#	writer = csv.writer(csv_file, delimiter=',')
#	for l in file_run:
#		writer.writerow(l)
	
