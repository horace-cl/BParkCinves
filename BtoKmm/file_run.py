import os
import csv

with open('files.txt', 'r') as file:
	doc = file.readlines()


query = 'das_client --query="run file='

file_run=[]
for i,line in enumerate(doc):
	#if i==5: break
	query_ = query+line.strip()
	query_+='"'
	stream = os.popen(query_)
	output = stream.read()
	if '317392' in output:
		print('Found ya!\n', line.strip())
	file_run.append([line.strip(), output.strip().replace('"', '')])

with open('file_run.txt', 'w+') as csv_file:
	writer = csv.writer(csv_file, delimiter=',')
	for l in file_run:
		writer.writerow(l)
	
