import sys

i=0
for line in sys.stdin:
	if line[0:6]=="# Gene":
		i+=1
		continue

	try:
		bg,end,sign=line.replace("\n","").split("\t")
		print(str(i)+"\t"+sign+"\t"+bg+"\t"+end)
	except: 
		pass	

