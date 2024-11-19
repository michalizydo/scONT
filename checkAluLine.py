# Simple check for counting SINE/Alu and LINE/L1-containing variants
# Usage python checkaluline.py <repeatmasker_output.out>

def main():
	import re
	count=0
	insa=[]
	insl=[]
	import sys
	rmfile=sys.argv[1]
	with open(rmfile,'r') as f:
		for line in f:

			count+=1
			if count<4:
				continue
			line=re.sub(' +', ' ', line).lstrip(" ")

			if "LINE/L1" in line:
				insl.append(tuple(line.split(" "))[4])
			elif "SINE/Alu" in line:
				insa.append(tuple(line.split(" "))[4])
	insl=set(insl)
	insa=set(insa)
	print(rmfile+" Alu: "+str(len(insa))+" LINE/L1: "+str(len(insl)))
	
if __name__=="__main__":
    main()
