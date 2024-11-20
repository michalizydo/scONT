# Counts breadth of coverage using mosdepth .bed output. 
# Usage: cat <mosd.output.bed> | python countcoveragebreadth.py  

def main():
       	chromosomes=[]
       	lastchrom=''
       	import sys
       	for line in sys.stdin:
       		liner=tuple(line.split("\t"))  
       		if liner[0]==lastchrom:
       			pass
       		else:
       			lastchrom=liner[0]
       			chromosomes.append([lastchrom])
       		chromosomes[-1].append(tuple([int(liner[2])-int(liner[1]),int(liner[3])]))    
       	for i in chromosomes:
       		#print(i[0])
       		totalsum=0
       		covsum=0
       		for j in i[1:]:
       			totalsum+=j[0]
       			covsum+=(j[1]*j[0])
       		print(i[0],end=' ')
       		print(float(covsum)/float(totalsum))
if __name__=="__main__":
	main()


