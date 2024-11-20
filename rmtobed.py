# Script to convert repeatmasker .out format to .bed format
# Usage python rmtobed.py <repeatmasker.out>

import sys
with open(sys.argv[1],'r') as f:
	counter=0
	for line in f:
		counter+=1
		if counter<4:
			continue

		line=tuple(line.split())
		if line[8]=="C":
			print(line[4]+"\t"+line[5]+"\t"+line[6]+"\t"+line[9]+","+line[10]+","+line[12]+","+line[13])
		else:
			print(line[4]+"\t"+line[5]+"\t"+line[6]+"\t"+line[9]+","+line[10]+","+line[11]+","+line[12])

f.close()