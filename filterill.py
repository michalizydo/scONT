# Script to filter deletions genotyped in PicoPLEX samples.
import sys
#var=[]
for i in sys.stdin:
	if i.startswith("#") or i.endswith("0\n") or i.endswith(".\n"):
		continue
	else:		
		print(i,end='')
#print(len(set(var)))
