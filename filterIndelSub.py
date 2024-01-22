import sys
for i in sys.stdin:
	if i.startswith("#"):
		print(i,end='')
	else:
		j=tuple(i.split("\t"))
		if len(j[3])==1 and len(j[4])==1:
			print(i,end='')

