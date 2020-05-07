with open('crab.log', 'r') as file:
	doc = file.readlines()

print(len(doc))

lumis = []
events = []
# for l in doc:
# 	if 'BasicSingleVertexState' in l: 
# 		print(len(l),'\n\n' ,l[:10],'\n\n')
# 		spli = l.split(':')
# 		for i,var in enumerate(spli):
# 			if 'lumi' in var:
# 				lumis.append(int(spli[i+1].strip().split(' ')[0]))
# 				events.append(int(spli[i+2].strip().split(' ')[0]))
# 				print('->\t', spli[i+1].strip().split(' ')[0])
# 				print('----->\t',spli[i-10:i+10], '\n')
# 		print('\n\n\n\n\n')


# print(len(lumis), set(lumis),  len(set(lumis)))
# print(len(events), set(events),  len(set(events)))
strings = []
for l in doc:
	if 'BasicSingleVertexState' in l:
		separado = l.split(' ')
		for i,var in enumerate(separado):
			if 'run:' in var:
				string = separado[i+1]+':'+separado[i+3]+':'+separado[i+5]
				strings.append(string)
				print(separado[i-1:i+7])

for k in set(strings):
	print(k)